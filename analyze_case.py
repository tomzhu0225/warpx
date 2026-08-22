#!/usr/bin/env python3
"""Analyze and plot one WarpX radiation analytic verification case."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import yt
from scipy.constants import Boltzmann, c, elementary_charge, physical_constants
from scipy.integrate import solve_ivp

RADIATION_CONSTANT = 4.0 * physical_constants["Stefan-Boltzmann constant"][0] / c


def wrapped_title(text: str, width: int = 58) -> str:
    return "\n".join(textwrap.wrap(text, width=width))


def load_field(plotfile: Path, name: str) -> np.ndarray:
    """Load one level-zero plotfile field without assuming a yt field type."""
    yt.funcs.mylog.setLevel("ERROR")
    ds = yt.load(str(plotfile))
    key = next((field for field in ds.field_list if field[1] == name), None)
    if key is None:
        raise KeyError(f"field {name!r} not found in {plotfile}: {ds.field_list}")
    grid = ds.covering_grid(
        level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
    )
    return np.squeeze(np.asarray(grid[key]))


def plotfiles(raw: Path) -> tuple[Path, Path]:
    files = sorted(path for path in raw.glob("plt*") if (path / "Header").exists())
    if len(files) < 2:
        raise RuntimeError(f"expected initial and final plotfiles under {raw}")
    return files[0], files[-1]


def cell_geometry(config: dict) -> tuple[list[np.ndarray], np.ndarray]:
    geometry = config["geometry"]
    lo = np.asarray(config["domain_lo_m"], dtype=float)
    hi = np.asarray(config["domain_hi_m"], dtype=float)
    cells = np.asarray(config["cells"], dtype=int)
    widths = (hi - lo) / cells
    centers = [
        lo[axis] + (np.arange(cells[axis]) + 0.5) * widths[axis]
        for axis in range(cells.size)
    ]
    if geometry == "1d":
        volumes = np.full(cells[0], widths[0])
    elif geometry == "2d":
        volumes = np.full(tuple(cells), widths[0] * widths[1])
    elif geometry == "rcyl":
        edges = lo[0] + np.arange(cells[0] + 1) * widths[0]
        volumes = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    else:
        raise ValueError(f"unknown geometry {geometry}")
    return centers, volumes


def reshape_field(values: np.ndarray, volumes: np.ndarray) -> np.ndarray:
    if values.shape == volumes.shape:
        return values
    if values.size != volumes.size:
        raise ValueError(f"field shape {values.shape} does not match {volumes.shape}")
    return values.reshape(volumes.shape)


def save_json(path: Path, values: dict) -> None:
    path.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n")


def analyze_beer_lambert(case: Path, config: dict) -> dict:
    raw = case / "results" / "raw"
    plots = case / "plots"
    results = case / "results"
    plots.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    initial_plotfile, final_plotfile = plotfiles(raw)
    deposition = load_field(final_plotfile, "radiation_material_energy")
    centers, volumes = cell_geometry(config)
    deposition = reshape_field(deposition, volumes)

    particle_energy = np.atleast_2d(np.loadtxt(raw / "particle_energy.txt"))
    # ParticleEnergy columns 2 and 3 are the weighted total and the sole
    # photon-species energy; the trailing columns are unweighted means.
    initial_energy = float(particle_energy[0, 2])
    final_energy = float(particle_energy[-1, 2])
    alpha = float(config["alpha_m_inv"])
    path_length = float(config["path_length_m"])
    expected_final = initial_energy * np.exp(-alpha * path_length)
    deposited = float(np.sum(deposition))
    expected_deposited = initial_energy - expected_final
    transmission_error = abs(final_energy - expected_final) / expected_final
    deposition_error = abs(deposited - expected_deposited) / expected_deposited
    ledger_error = abs(final_energy + deposited - initial_energy) / initial_energy

    geometry = config["geometry"]
    start = np.asarray(config["start_m"], dtype=float)
    direction = np.asarray(config["direction"], dtype=float)
    direction /= np.linalg.norm(direction)

    if geometry in ("1d", "rcyl"):
        coordinate = centers[0]
        mask = (coordinate >= start[0]) & (
            coordinate <= start[0] + path_length + (coordinate[1] - coordinate[0])
        )
        order_coordinate = coordinate[mask]
        order_deposition = deposition[mask]
        projection = order_coordinate - start[0]
        xlabel = "z [m]" if geometry == "1d" else "r [m]"

        fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
        ax.axvspan(start[0], start[0] + path_length, color="#d8e6f3", alpha=0.8)
        ax.axvline(start[0], color="#c0392b", lw=2, label="packet start")
        ax.annotate(
            "propagation",
            xy=(start[0] + 0.75 * path_length, 0.55),
            xytext=(start[0] + 0.2 * path_length, 0.55),
            arrowprops={"arrowstyle": "->", "lw": 2},
        )
        ax.set(xlabel=xlabel, ylabel="normalized packet energy", ylim=(0, 1.05))
        ax.set_title(wrapped_title(config["title"] + " -- initial"), fontsize=12)
        ax.legend(loc="upper right")
        fig.savefig(plots / "initial.png", dpi=180)
        plt.close(fig)

        density = deposition / volumes
        fig, ax = plt.subplots(figsize=(7.0, 4.0), constrained_layout=True)
        ax.plot(coordinate, density, color="#1f77b4", lw=1.8)
        ax.axvline(
            start[0] + path_length, color="#c0392b", ls="--", label="packet final"
        )
        ax.set(xlabel=xlabel, ylabel="deposited energy density [J m$^{-3}$]")
        ax.set_title(
            wrapped_title(config["title"] + " -- final deposition"), fontsize=12
        )
        ax.legend()
        fig.savefig(plots / "final.png", dpi=180)
        plt.close(fig)
    else:
        x, z = centers
        xgrid, zgrid = np.meshgrid(x, z, indexing="ij")
        projection_all = (xgrid - start[0]) * direction[0] + (
            zgrid - start[1]
        ) * direction[1]
        active = deposition > max(float(np.max(deposition)) * 1.0e-14, 0.0)
        order = np.argsort(projection_all[active])
        projection = projection_all[active][order]
        order_deposition = deposition[active][order]

        fig, ax = plt.subplots(figsize=(6.0, 5.2), constrained_layout=True)
        ax.set_xlim(config["domain_lo_m"][0], config["domain_hi_m"][0])
        ax.set_ylim(config["domain_lo_m"][1], config["domain_hi_m"][1])
        ax.plot(start[0], start[1], "o", color="#c0392b", label="packet start")
        end = start + path_length * direction
        ax.annotate(
            "",
            xy=end,
            xytext=start,
            arrowprops={"arrowstyle": "->", "lw": 2.5, "color": "#1f77b4"},
        )
        ax.set(xlabel="x [m]", ylabel="z [m]", aspect="equal")
        ax.set_title(wrapped_title(config["title"] + " -- initial"), fontsize=12)
        ax.legend()
        fig.savefig(plots / "initial.png", dpi=180)
        plt.close(fig)

        density = deposition / volumes
        fig, ax = plt.subplots(figsize=(6.1, 5.2), constrained_layout=True)
        image = ax.pcolormesh(x, z, density.T, shading="auto", cmap="magma")
        ax.plot([start[0], end[0]], [start[1], end[1]], color="cyan", ls="--", lw=1.2)
        ax.plot(end[0], end[1], "o", color="cyan", ms=4)
        ax.set(xlabel="x [m]", ylabel="z [m]", aspect="equal")
        ax.set_title(
            wrapped_title(config["title"] + " -- final deposition"), fontsize=12
        )
        fig.colorbar(image, ax=ax, label="deposited energy density [J m$^{-3}$]")
        fig.savefig(plots / "final.png", dpi=180)
        plt.close(fig)

    order_mask = (projection >= -1.0e-14) & (projection <= path_length * 1.01)
    projection = projection[order_mask]
    order_deposition = order_deposition[order_mask]
    order = np.argsort(projection)
    projection = projection[order]
    order_deposition = order_deposition[order]
    numerical_transmission = 1.0 - np.cumsum(order_deposition) / initial_energy
    analytic_transmission = np.exp(-alpha * np.maximum(projection, 0.0))
    dense_s = np.linspace(0.0, path_length, 400)

    np.savetxt(
        results / "profile.csv",
        np.column_stack((projection, numerical_transmission, analytic_transmission)),
        delimiter=",",
        header="path_m,numerical_transmission,analytic_transmission",
        comments="",
    )
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    ax.plot(dense_s, np.exp(-alpha * dense_s), "k-", lw=2, label="Beer-Lambert")
    ax.step(
        projection,
        numerical_transmission,
        where="post",
        color="#d35400",
        lw=1.4,
        label="WarpX cumulative deposition",
    )
    ax.plot(path_length, final_energy / initial_energy, "o", color="#1f77b4", ms=6)
    ax.set(
        xlabel="path length [m]",
        ylabel="$E/E_0$",
        ylim=(0.0, 1.03),
        title=wrapped_title(config["title"] + " -- analytic comparison"),
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(plots / "comparison.png", dpi=180)
    plt.close(fig)

    summary = {
        "case": case.name,
        "family": config["family"],
        "geometry": geometry,
        "initial_photon_energy_J": initial_energy,
        "simulated_final_photon_energy_J": final_energy,
        "analytic_final_photon_energy_J": expected_final,
        "simulated_deposited_energy_J": deposited,
        "analytic_deposited_energy_J": expected_deposited,
        "relative_transmission_error": transmission_error,
        "relative_deposition_error": deposition_error,
        "relative_energy_ledger_error": ledger_error,
        "passed": bool(
            transmission_error < 2.0e-11
            and deposition_error < 2.0e-11
            and ledger_error < 2.0e-11
        ),
    }
    return summary


def analyze_exchange(case: Path, config: dict) -> dict:
    raw = case / "results" / "raw"
    plots = case / "plots"
    results = case / "results"
    plots.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    initial_plotfile, final_plotfile = plotfiles(raw)
    centers, volumes = cell_geometry(config)
    final_radiation = reshape_field(
        load_field(final_plotfile, "radiation_diffusion_energy"), volumes
    )
    final_material = reshape_field(
        load_field(final_plotfile, "radiation_material_energy"), volumes
    )
    initial_te = reshape_field(load_field(initial_plotfile, "Te"), volumes)
    final_te = reshape_field(load_field(final_plotfile, "Te"), volumes)

    active = np.ones(volumes.shape, dtype=bool)
    if "active_region_m" in config:
        active_lo, active_hi = config["active_region_m"]
        active = (centers[0] >= active_lo) & (centers[0] < active_hi)

    radiation_table = np.atleast_2d(np.loadtxt(raw / "radiation_energy.txt"))
    time = radiation_table[:, 1]
    simulated_energy = radiation_table[:, 4]
    cumulative_material = radiation_table[:, 6]
    total_volume = float(np.sum(volumes[active]))
    simulated_density = simulated_energy / total_volume

    density = float(config["electron_density_m3"])
    gamma = float(config["gamma"])
    alpha = float(config["planck_alpha_m_inv"])
    initial_temperature = (
        float(config["initial_electron_temperature_eV"]) * elementary_charge / Boltzmann
    )
    heat_capacity = density * Boltzmann / (gamma - 1.0)
    total_density = heat_capacity * initial_temperature

    def rhs(_time: float, radiation_density: np.ndarray) -> np.ndarray:
        temperature = np.maximum(
            (total_density - radiation_density[0]) / heat_capacity, 0.0
        )
        return np.array(
            [c * alpha * (RADIATION_CONSTANT * temperature**4 - radiation_density[0])]
        )

    reference = solve_ivp(
        rhs,
        (float(time[0]), float(time[-1])),
        np.array([0.0]),
        t_eval=time,
        rtol=2.0e-12,
        atol=1.0e-12,
        method="DOP853",
    )
    if not reference.success:
        raise RuntimeError(reference.message)
    analytic_density = reference.y[0]
    analytic_temperature = (total_density - analytic_density) / heat_capacity
    inferred_sim_temperature = (total_density - simulated_density) / heat_capacity

    final_density_field = final_radiation / volumes
    final_reference_density = float(analytic_density[-1])
    final_reference_temperature = float(analytic_temperature[-1])
    volume_weighted_te = float(
        np.sum(final_te[active] * volumes[active]) / total_volume
    )
    trajectory_error = float(
        np.max(np.abs(simulated_density - analytic_density)) / np.max(analytic_density)
    )
    final_density_error = (
        abs(simulated_density[-1] - final_reference_density) / final_reference_density
    )
    final_temperature_error = (
        abs(volume_weighted_te - final_reference_temperature)
        / final_reference_temperature
    )
    uniformity_error = float(
        np.std(final_density_field[active]) / np.mean(final_density_field[active])
    )
    conservation_error = (
        abs(simulated_energy[-1] + cumulative_material[-1]) / simulated_energy[-1]
    )
    current_material_exchange = float(radiation_table[-1, 5])
    plotfile_current_exchange_error = abs(
        float(np.sum(final_material)) - current_material_exchange
    ) / max(abs(current_material_exchange), np.finfo(float).tiny)

    simulated_te_eV = inferred_sim_temperature * Boltzmann / elementary_charge
    analytic_te_eV = analytic_temperature * Boltzmann / elementary_charge
    simulated_tr_eV = (
        (np.maximum(simulated_density, 0.0) / RADIATION_CONSTANT) ** 0.25
        * Boltzmann
        / elementary_charge
    )
    analytic_tr_eV = (
        (np.maximum(analytic_density, 0.0) / RADIATION_CONSTANT) ** 0.25
        * Boltzmann
        / elementary_charge
    )
    np.savetxt(
        results / "history.csv",
        np.column_stack(
            (
                time,
                simulated_density,
                analytic_density,
                simulated_te_eV,
                analytic_te_eV,
                simulated_tr_eV,
                analytic_tr_eV,
            )
        ),
        delimiter=",",
        header=(
            "time_s,simulated_Er_J_m3,analytic_Er_J_m3,"
            "simulated_Te_eV,analytic_Te_eV,simulated_Tr_eV,analytic_Tr_eV"
        ),
        comments="",
    )

    geometry = config["geometry"]
    initial_tr = np.zeros_like(initial_te)
    final_tr = (
        (np.maximum(final_density_field, 0.0) / RADIATION_CONSTANT) ** 0.25
        * Boltzmann
        / elementary_charge
    )
    initial_te_eV_field = initial_te * Boltzmann / elementary_charge
    final_te_eV_field = final_te * Boltzmann / elementary_charge

    if geometry in ("1d", "rcyl"):
        coordinate = centers[0]
        xlabel = "z [m]" if geometry == "1d" else "r [m]"
        plotted = active if geometry == "rcyl" else np.ones_like(active)
        for path, te_field, tr_field, label in (
            (plots / "initial.png", initial_te_eV_field, initial_tr, "initial"),
            (plots / "final.png", final_te_eV_field, final_tr, "final"),
        ):
            fig, ax = plt.subplots(figsize=(7.0, 4.2), constrained_layout=True)
            ax.plot(coordinate[plotted], te_field[plotted], lw=2, label="$T_e$")
            ax.plot(coordinate[plotted], tr_field[plotted], lw=2, label="$T_r$")
            if "active_region_m" in config:
                active_lo, active_hi = config["active_region_m"]
                ax.axvspan(
                    active_lo,
                    active_hi,
                    color="#d8e6f3",
                    alpha=0.35,
                    label="active plasma annulus",
                )
                ax.set_xlim(config["domain_lo_m"][0], config["domain_hi_m"][0])
            ax.set(xlabel=xlabel, ylabel="temperature [eV]")
            ax.set_title(wrapped_title(config["title"] + f" -- {label}"), fontsize=12)
            ax.grid(alpha=0.25)
            ax.legend()
            fig.savefig(path, dpi=180)
            plt.close(fig)
    else:
        x, z = centers
        maximum_temperature = float(config["initial_electron_temperature_eV"])
        for path, te_field, tr_field, label in (
            (plots / "initial.png", initial_te_eV_field, initial_tr, "initial"),
            (plots / "final.png", final_te_eV_field, final_tr, "final"),
        ):
            fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), constrained_layout=True)
            for ax, field, title in zip(axes, (te_field, tr_field), ("$T_e$", "$T_r$")):
                image = ax.pcolormesh(
                    x,
                    z,
                    field.T,
                    shading="auto",
                    cmap="inferno",
                    vmin=0.0,
                    vmax=maximum_temperature,
                )
                ax.set(xlabel="x [m]", ylabel="z [m]", title=title, aspect="equal")
                fig.colorbar(image, ax=ax, label="temperature [eV]")
            fig.suptitle(wrapped_title(config["title"] + f" -- {label}"), fontsize=12)
            fig.savefig(path, dpi=180)
            plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    time_ps = time * 1.0e12
    axes[0].plot(time_ps, analytic_te_eV, "k-", lw=2, label="ODE $T_e$")
    axes[0].plot(time_ps, simulated_te_eV, "o", ms=3, markevery=3, label="WarpX")
    axes[0].set(xlabel="time [ps]", ylabel="$T_e$ [eV]")
    axes[1].plot(time_ps, analytic_tr_eV, "k-", lw=2, label="ODE $T_r$")
    axes[1].plot(time_ps, simulated_tr_eV, "o", ms=3, markevery=3, label="WarpX")
    axes[1].set(xlabel="time [ps]", ylabel="$T_r$ [eV]")
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend()
    fig.suptitle(
        wrapped_title(config["title"] + " -- analytic exchange history", 80),
        fontsize=13,
    )
    fig.savefig(plots / "comparison.png", dpi=180)
    plt.close(fig)

    summary = {
        "case": case.name,
        "family": config["family"],
        "geometry": geometry,
        "total_volume_m3_per_suppressed_length": total_volume,
        "simulated_final_radiation_energy_J": float(simulated_energy[-1]),
        "analytic_final_radiation_energy_J": final_reference_density * total_volume,
        "simulated_final_electron_temperature_eV": volume_weighted_te
        * Boltzmann
        / elementary_charge,
        "analytic_final_electron_temperature_eV": final_reference_temperature
        * Boltzmann
        / elementary_charge,
        "relative_history_Linf_error": trajectory_error,
        "relative_final_radiation_error": final_density_error,
        "relative_final_temperature_error": final_temperature_error,
        "relative_spatial_nonuniformity": uniformity_error,
        "relative_reduced_diagnostic_conservation_error": conservation_error,
        "relative_plotfile_current_exchange_error": plotfile_current_exchange_error,
        "passed": bool(
            trajectory_error < 3.0e-2
            and final_density_error < 2.0e-2
            and final_temperature_error < 2.0e-2
            and uniformity_error < 3.0e-2
            and conservation_error < 2.0e-10
            and plotfile_current_exchange_error < 2.0e-10
        ),
    }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    case = args.case.resolve()
    config = json.loads((case / "case.json").read_text())
    if config["family"] == "beer_lambert":
        summary = analyze_beer_lambert(case, config)
    elif config["family"] == "matter_radiation_exchange":
        summary = analyze_exchange(case, config)
    else:
        raise ValueError(config["family"])
    save_json(case / "results" / "summary.json", summary)
    print(json.dumps(summary, indent=2, sort_keys=True))
    if not summary["passed"]:
        raise SystemExit(f"analytic acceptance failed for {case.name}")


if __name__ == "__main__":
    main()
