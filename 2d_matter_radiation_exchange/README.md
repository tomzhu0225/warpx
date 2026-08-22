# Cartesian 2-D radiation-matter exchange

This is the uniform Turner-Stone-style grey equilibration benchmark on a
32-by-32 Cartesian mesh. It verifies that the material-radiation source and
suppressed-length area accounting are independent of transverse position.

Run `./run.sh`. The final radiation error is 6.45e-5, the full-history error is
4.53e-3, spatial nonuniformity is 4.2e-15 and the reduced-diagnostic energy
ledger closes to 1.2e-14. The initial and final plots share a physical 0–10 eV
color scale so roundoff is not visually exaggerated.
