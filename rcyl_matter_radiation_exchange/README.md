# RCYL radiation-matter exchange

The uniform grey equilibration benchmark is applied to a stationary plasma
annulus from 2.5 to 7.5 mm. Keeping the active material away from the axis and
particle boundary isolates cylindrical volume accounting from coordinate and
boundary artifacts.

Run `./run.sh`. At 256 radial cells, the final radiation/history error is
0.455%, the final electron-temperature error is 0.128% and the conservative
energy ledger closes to 4.0e-16. A 64/128/256-cell study produced 1.84%, 0.916%
and 0.455% radiation errors, showing first-order convergence.
