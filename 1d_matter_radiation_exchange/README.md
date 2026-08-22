# Planar 1-D radiation-matter exchange

A stationary, uniform 10 eV electron plasma emits into an initially empty grey
radiation field. With transport disabled, the exact reference is the coupled
Turner-Stone-style material-radiation ODE described in the parent README.

Run `./run.sh`. `results/history.csv` contains both WarpX and high-accuracy ODE
histories for radiation energy density, electron temperature and radiation
temperature. The final radiation error is 6.45e-5 and the conservative energy
ledger closes to 2.0e-16 relative error.
