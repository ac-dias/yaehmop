# YAeHMOP: Yet Another extended Hueckel Molecular Orbital Package

This software was original written by Greg Landrum during his time in the Hoffmann Group at Cornell University. Modifications to support COHP analysis and f-orbitals were provided by Wingfield Glassey.

The license has been updated to the standard BSD 2-clause license for the github version of the code.


## Fork additions (ac-dias/yaehmop, branch adias-dev)

This fork adds spin polarization and spin-orbit coupling to the `tightbind` program
(`bind`). Runs without the new options give the same output as the original program.

**Developer documentation: [tex/advances_in_yaehmop.pdf](tex/advances_in_yaehmop.pdf)**
("Advances in YAeHMOP"). It explains every new and modified file, how each change works,
what changed in the inputs and outputs, how the changes were tested, and two worked
examples. It was written to support a pull request to the original repository. The LaTeX
source is in `tex/` (`make` builds it).

New features:

- **Spin-polarized runs** with one parameter table per spin channel:
  `bind --parm-up up.dat --parm-down down.dat input`. The two channels run as two
  independent calculations, and every output gets `_up` / `_down` before the input
  extension. Optional: `--nelec-up/--nelec-down` (electrons per channel) and
  `--const-up/--const-down` (Wolfsberg-Helmholz K per channel).
- **Spin-orbit coupling** (on-site lambda L.S): `bind --soc soc.dat input params.dat`.
  Gives 2n spinor bands, the Fermi level and total DOS, spinor H(k)/S(k) in the `.HAM`/`.OV`
  dumps, and a `.SOC` file with H_SOC.
- **`THE CONST DIAG <K>`**: new input keyword for the K of the diagonal H_ii(k) term (an
  orbital with its own images). The original code had 1.75 hardwired here, and the default
  is still 1.75.
- Fix: the `.HAM` dump is truncated when it is opened (it used to be appended to).

Build with LAPACK (`cmake ../tightbind -DUSE_BLAS_LAPACK=ON`); SOC needs it.

| where | what |
|-------|------|
| [tex/advances_in_yaehmop.pdf](tex/advances_in_yaehmop.pdf) | developer documentation of all the changes |
| [docs/SPIN_SOC.md](docs/SPIN_SOC.md) | user guide: options, file formats, limitations |
| [new_tests/](new_tests/) | band structures of a MoS2 monolayer with SOC and of bcc Fe with one table per spin channel: inputs, outputs, figures, README |
| [examples/spin_soc/](examples/spin_soc/) | small diamond example |
| [tests/spin_soc/run_tests.py](tests/spin_soc/run_tests.py) | automated tests (`python3 tests/spin_soc/run_tests.py --bind build/bind`) |

### License

The fork keeps the original BSD 2-clause license ([license.txt](license.txt)). The new
code (`tightbind/soc.c` and the changes to the other files) is under the same license.
`make install` now also installs `license.txt` (in `share/yaehmop/`), as the license
requires for binary redistributions. The parameter sets used in `new_tests/` and
`tests/spin_soc/data/` come from
[yaehmop-bandfit](https://github.com/tromer-unb/yaehmop-bandfit) (MIT license). When you
publish results, cite the YAeHMOP authors (G. A. Landrum and W. V. Glassey) and this fork.
