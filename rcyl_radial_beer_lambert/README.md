# RCYL radial Beer-Lambert searchlight

A radially propagating photon packet is attenuated over 0.7 m by a uniform
3 m^-1 extinction coefficient. The integrated solution remains
`E(r)=E0 exp(-alpha r)` while the deposited energy density includes exact
cylindrical annular volumes.

Run `./run.sh`. The comparison plot verifies the integrated exponential law;
the final plot intentionally shows the radial dilution of deposited energy
density. Transmission, deposition and total energy close exactly at the
reported precision.
