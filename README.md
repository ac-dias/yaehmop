# YAeHMOP: Yet Another extended Hueckel Molecular Orbital Package

This software was original written by Greg Landrum during his time in the Hoffmann Group at Cornell University. Modifications to support COHP analysis and f-orbitals were provided by Wingfield Glassey.

The license has been updated to the standard BSD 2-clause license for the github version of the code.


## Fork additions (ac-dias/yaehmop)

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
  `--common-fermi` gives one Fermi level for both channels. Every spin-polarized run
  writes a summary `<input>.spin` (Fermi level, N_up, N_down, moments per atom).
- **Spin-orbit coupling** (on-site lambda L.S): `bind --soc soc.dat input params.dat`.
  Gives 2n spinor bands, spinor H(k)/S(k) in the `.HAM`/`.OV` dumps, a `.SOC` file with
  H_SOC, and every option of the program with SOC: projected DOS (also resolved in spin,
  keyword `Spin Resolved DOS`), COOP/COHP, charges, populations, modified Mulliken, charge
  iteration, zeta variation, Walsh, molecules (electrostatic energy, average properties),
  FMO/FCO, `MO print` for viewkel, sparse dumps and netCDF, plus spin moments.
- **SOC on top of two spin channels**: `bind --soc soc.dat --parm-up up.dat --parm-down
  down.dat input`, H = diag(H_up, H_down) + H_SOC. The tables may have different Hii and
  zetas (S = diag(S_up, S_down) when the zetas differ); charge iteration and zeta variation
  change both tables alike.
- **Fixes of old defects** found on the way: the `vary` list of `Charge Iteration` crashed
  (uninitialized pointer); `my_realloc` read past the old block; `calc_occupations` read
  outside the eigenvalue array; the Walsh item of the reduced charge matrix used a wrong
  index; netCDF energies and coefficients were written from the first k point only; the OP
  matrix was not allocated when only the modified Mulliken analysis was asked for.
- **viewkel**: a file named `*_up.*` is shown together with its `*_down.*` partner
  (bands, DOS, COOP).
- **`THE CONST DIAG <K>`**: new input keyword for the K of the diagonal H_ii(k) term (an
  orbital with its own images). The original code had 1.75 hardwired here, and the default
  is still 1.75.
- **Double-zeta p orbitals**: a p shell may now have two Slater exponents, as d and f
  shells already could. In the parameter table, set `Nzeta = 2` on the p line and give
  `exp1 exp2 coeff1 coeff2`; the contraction is normalized as for d/f. A contracted
  inner function plus a diffuse tail gives, for example, interlayer coupling across a van
  der Waals gap that a single exponent cannot reach. Single-zeta tables give exactly the
  original results. The overlap integrals (bands, DOS, COOP, H/S dumps) use the contracted
  function, and so do the post-processing routines: the molecular electrostatic energy
  (exact potential of the contracted density), the zeta variation (both exponents
  shifted, contraction renormalized), the `.MO` file (three extra fields on a double-zeta
  p line) with viewkel's MO surfaces and contours, and the netCDF output. The Muller
  iteration uses its own single-zeta p, as it already did for d, and says so in the
  output. The netCDF writer also got four old fixes (wrong parameter index, int indices
  passed as `size_t*`, file closed after the memory was freed, swapped `sprintf`
  arguments in the Walsh-step file name). Double-zeta s is not
  supported. Adapted from the research patch of yaehmop-bandfit
  (`experimental/carbon_2p_double_zeta.patch`).
- **Structure converter** `interface/kel_crystal_conversion.py`: reads the crystal structure of VASP
  (POSCAR/CONTCAR), FHI-aims (geometry.in), SIESTA (.fdf, .XV, .STRUCT_OUT), DFTB+ (.gen, dftb_in.hsd)
  or Quantum ESPRESSO (pw.x input) and writes a periodic YAeHMOP input (`Geometry Crystallographic`,
  `Lattice`, `Crystal Spec`; optional `THE CONST`, `Electrons`, `Band`, `K points`). The periodic
  dimension is detected from the vacuum and slabs wrapped across the vacuum are made contiguous.
- Fixes: the `.HAM` dump is truncated when it is opened (it used to be appended to); the
  projected DOS uses the right k-point weights (the original curves change slightly).

Build with LAPACK (`cmake ../tightbind -DUSE_BLAS_LAPACK=ON`); SOC needs it.

| where | what |
|-------|------|
| [tex/advances_in_yaehmop.pdf](tex/advances_in_yaehmop.pdf) | developer documentation of all the changes |
| [docs/SPIN_SOC.md](docs/SPIN_SOC.md) | user guide: options, file formats, limitations |
| [new_tests/](new_tests/) | band structures of a MoS2 monolayer with SOC and of bcc Fe with one table per spin channel: inputs, outputs, figures, README |
| [interface/](interface/) | `yaehmop2wtb.py` (lean) and `yaehmop2wtb-full.py`: convert the exported H and S (spinless, SOC or spin polarized) to the [WanTiBEXOS](https://github.com/ac-dias/wantibexos) tight-binding format, with the basis file and, with `--rmatrix`, the position matrix of the Slater basis (`BSE_CENTER_FILE`); the full version also writes the `wtb.x` input and k path. MoS2 example validated with `wtb.x` |
| [interface/kel_crystal_conversion.py](interface/kel_crystal_conversion.py) | crystal structure of VASP, FHI-aims, SIESTA, DFTB+ or Quantum ESPRESSO → periodic YAeHMOP input (`kel_crystal_conversion.py --format VASP --in POSCAR`; `--help` for all options); single file, standard library only |
| [examples/spin_soc/](examples/spin_soc/) | small diamond example |
| [tests/spin_soc/run_tests.py](tests/spin_soc/run_tests.py) | automated tests (`python3 tests/spin_soc/run_tests.py --bind build/bind`) |

### License

The fork keeps the original BSD 2-clause license ([license.txt](license.txt)). The new
code (`tightbind/soc.c` and the changes to the other files) is under the same license.
`make install` now also installs `license.txt` (in `share/yaehmop/`), as the license
requires for binary redistributions. The parameter sets used in `new_tests/` and
`tests/spin_soc/data/` come from yaehmop-bandfit, our parameter-fitting tool (its
repository is still private; MIT license). When you
publish results, cite the YAeHMOP authors (G. A. Landrum and W. V. Glassey) and this fork.
