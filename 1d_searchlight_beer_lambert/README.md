# Planar 1-D Beer-Lambert searchlight

A single photon packet propagates 0.7 m through a uniform extinction
coefficient of 3 m^-1. The exact transmission is `exp(-2.1)`. This combines
the standard searchlight streaming test with an exact absorbing-medium
solution.

Run `./run.sh`. The numerical and analytic transmission profiles are in
`results/profile.csv`; `plots/` contains the initial geometry, final deposition
and analytic comparison. The committed double-precision result has zero
transmission, deposition and energy-ledger error at the reported precision.
