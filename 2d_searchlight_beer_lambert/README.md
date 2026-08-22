# Cartesian 2-D Beer-Lambert searchlight

A photon packet crosses the Cartesian mesh at 45 degrees while a uniform
3 m^-1 extinction coefficient attenuates it over 0.7 m. The test simultaneously
checks oblique packet streaming, localized deposition and the exact
`exp(-alpha L)` transmission law.

Run `./run.sh`. `plots/final.png` visualizes the diagonal deposition track and
`plots/comparison.png` overlays the cumulative WarpX transmission on the
analytic curve. Transmission is exact at the reported precision and the total
energy-ledger error is 2.0e-16.
