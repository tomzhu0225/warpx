#!/usr/bin/env python3
"""Run all or one of the radiation analytic verification cases."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CASES = (
    "1d_searchlight_beer_lambert",
    "1d_matter_radiation_exchange",
    "2d_searchlight_beer_lambert",
    "2d_matter_radiation_exchange",
    "rcyl_radial_beer_lambert",
    "rcyl_matter_radiation_exchange",
)
DEFAULT_EXECUTABLES = {
    "1d": Path(
        "/home/tomzhu0225/src/warpx-radiation-transport/"
        "build-boundary-pr-1d/bin/warpx.1d.NOMPI.OMP.DP.PDP.EB"
    ),
    "2d": Path(
        "/home/tomzhu0225/src/warpx-radiation-transport/"
        "build-four-force-public-2d/bin/warpx.2d.NOMPI.OMP.DP.PDP.EB"
    ),
    "rcyl": Path(
        "/home/tomzhu0225/src/warpx-radiation-transport/"
        "build-boundary-pr-rcyl/bin/warpx.rcylinder.NOMPI.OMP.DP.PDP.EB"
    ),
}
ENV_KEYS = {"1d": "WARPX_1D", "2d": "WARPX_2D", "rcyl": "WARPX_RCYL"}


def run_case(name: str, python: Path) -> dict:
    case = ROOT / name
    config = json.loads((case / "case.json").read_text())
    geometry = config["geometry"]
    executable = Path(os.environ.get(ENV_KEYS[geometry], DEFAULT_EXECUTABLES[geometry]))
    if not executable.is_file():
        raise FileNotFoundError(f"missing {geometry} WarpX executable: {executable}")

    raw = case / "results" / "raw"
    if raw.exists():
        shutil.rmtree(raw)
    raw.mkdir(parents=True)
    (case / "plots").mkdir(parents=True, exist_ok=True)
    log = case / "results" / "run.log"
    environment = os.environ.copy()
    environment.setdefault("OMP_NUM_THREADS", "1")
    print(f"[{name}] {executable.name}", flush=True)
    with log.open("w") as stream:
        completed = subprocess.run(
            [str(executable), "inputs"],
            cwd=case,
            env=environment,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        tail = "\n".join(log.read_text(errors="replace").splitlines()[-40:])
        raise RuntimeError(f"WarpX failed for {name}:\n{tail}")
    subprocess.run(
        [str(python), str(ROOT / "analyze_case.py"), str(case)],
        cwd=case,
        env=environment,
        check=True,
    )
    return json.loads((case / "results" / "summary.json").read_text())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", choices=CASES, action="append")
    parser.add_argument(
        "--python",
        type=Path,
        default=Path(os.environ.get("PYTHON", sys.executable)),
        help="Python containing numpy, scipy, matplotlib and yt",
    )
    args = parser.parse_args()
    selected = args.case or list(CASES)
    # Do not resolve a virtual-environment interpreter symlink: executing the
    # resolved system binary would discard the venv's site-packages.
    summaries = [run_case(name, args.python) for name in selected]
    passed = all(summary["passed"] for summary in summaries)
    (ROOT / "summary.json").write_text(
        json.dumps({"passed": passed, "cases": summaries}, indent=2) + "\n"
    )
    print(f"completed {len(summaries)} cases: passed={passed}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
