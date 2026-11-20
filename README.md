# YAeHMOP: Yet Another extended Hueckel Molecular Orbital Package

This software was original written by Greg Landrum during his time in the Hoffmann Group at Cornell University. Modifications to support COHP analysis and f-orbitals were provided by Wingfield Glassey.

The license has been updated to the standard BSD 2-clause license for the github version of the code.


## Fork additions (ac-dias/yaehmop, branch adias-dev)

- Spin-polarized runs with one parameter table per channel:
  `bind --parm-up up.dat --parm-down down.dat input` (outputs `input_up.*`, `input_down.*`).
- Spin-orbit coupling (on-site lambda L.S): `bind --soc soc.dat input params.dat`
  (2n spinor bands, spinor H(k)/S(k) dumps, `.SOC` file with H_SOC).

Build with LAPACK (`cmake ../tightbind -DUSE_BLAS_LAPACK=ON`). Details, file formats,
limitations and tests: [docs/SPIN_SOC.md](docs/SPIN_SOC.md); example in `examples/spin_soc/`.
