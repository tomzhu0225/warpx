# QDSMC electron-entropy remap investigation

These figures support the WarpX Discussion about cumulative numerical
diffusion in the QDSMC electron-entropy gather/deposit cycle.

- `Te_foot_spacetime_propagation.png` compares the radial and temporal
  electron-temperature structure of the experimental RCYL Hybrid-PIC run
  with the implicit kinetic reference.
- `E_series_foot_mechanism.png` compares the full QDSMC control with a
  no-marker-advection control, a run with half as many remaps, and a radial
  marker-volume control at 40 ns.

The underlying RCYL case uses experimental geometry support on
`feature/hybrid-implicit-mag-diffusion`. The proposed gather/deposit mechanism
is geometry-independent; the Discussion asks whether it is intended QDSMC
behavior and requests guidance on a minimal supported-geometry regression.
