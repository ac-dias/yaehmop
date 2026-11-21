# new_tests: band structures with SOC (MoS2) and spin polarization (Fe)

These two examples use the parametrizations fitted with yaehmop-bandfit (our
parameter-fitting tool; its repository is still private) and the fork's new options
(fork `ac-dias/yaehmop`, see `docs/SPIN_SOC.md`).

| folder      | system                | parametrization                                                          | new options                       |
|-------------|-----------------------|--------------------------------------------------------------------------|-----------------------------------|
| `MoS2_SOC/` | MoS2 monolayer        | EHT fit to FHI-aims HSE06+SOC bands + on-site lambda (Mo d, S p)         | `--soc`                           |
| `Fe_spin/`  | bcc Fe (ferromagnet)  | one EHT fit per spin channel to VASP PBE (ISPIN=2) bands                 | `--parm-up/--parm-down`, `--const-*`, `--nelec-*` |

Each folder has the YAeHMOP input (`*.bind`) and the parameter files. Its `out/`
subfolder holds the outputs of the runs below, a copy of the input, the plotting script
`plot_bands.py` and the figure.

## Requirements

Build the fork with LAPACK, or load the module installed on our machine:

```
module load yaehmop-dev        # provides the command `yaehmop` (bash already has a builtin called `bind`)
```

To build it yourself instead:

```
mkdir build && cd build
cmake ../tightbind -DUSE_BLAS_LAPACK=ON && make
# then use ./bind in place of `yaehmop`
```

The figures need Python 3 with numpy and matplotlib.

## MoS2 monolayer with spin-orbit coupling (`MoS2_SOC/`)

Files:
- `MoS2_bands.bind`: structure (3 atoms; atoms 4 and 5 only define the lattice vectors a
  and b), THE CONST 1.7221, band path Γ-M-K-Γ (40 points per line), 18 electrons,
  `Average Properties` on a 12×12×1 Γ-centred mesh (Fermi level and DOS), and `Dump Hamil`
  / `Dump Overlap` (H(k), S(k) on that mesh).
- `MoS2_params.dat`: EHT parameters (Mo s, p, d double zeta; S s, p).
- `MoS2_soc.dat`: lambda (eV), Mo d 0.063352 and S p 0.197577.

Run:

```
cd MoS2_SOC/out
cp ../MoS2_bands.bind .
yaehmop --soc ../MoS2_soc.dat MoS2_bands.bind ../MoS2_params.dat     # with SOC

mkdir -p nosoc && cd nosoc                                           # same model, no SOC
cp ../../MoS2_bands.bind .
yaehmop MoS2_bands.bind ../../MoS2_params.dat
cd ..

python3 plot_bands.py               # -> bands_soc.png (options: --emin, --emax)
```

Outputs in `out/` (SOC run):
- `MoS2_bands.bind.band`: 34 spinor bands (17 orbitals × 2) along the path.
- `MoS2_bands.bind.HAM` and `.OV`: spinor H(k) and S(k), 34×34, on the 144 mesh points.
  These are binary files: int nk, int n, then per k an n×n double array. Element (i,j)
  with i<j is `m[i][j] - i*m[j][i]`. The basis is orbitals 1..17 spin up, then 18..34
  spin down.
- `MoS2_bands.bind.SOC`: on-site H_SOC in text form.
- `MoS2_bands.bind.out`: lambdas used per atom, `#SOC_STATES: 34`, total DOS and
  `#Fermi_Energy` (the top of the valence band for an insulator).
- `MoS2_bands.bind.status`: progress log.
- `nosoc/`: the same files from the run without SOC (17 bands, 17×17 matrices).

Expected results (they match the yaehmop-bandfit fit):

| quantity                            | without SOC | with SOC |
|-------------------------------------|-------------|----------|
| band gap (eV)                       | 2.3130      | 2.2681   |
| valence band splitting at K (eV)    | 0           | 0.1650   |
| Fermi level / VBM (eV)              | −9.0366     | −9.0108  |

## bcc Fe, spin polarized (`Fe_spin/`)

Files:
- `Fe_bands.bind`: primitive bcc cell (1 atom; atoms 2-4 only define the lattice
  vectors), band path Γ-H-N-Γ-P-H (the path of the DFT reference, 40 points per line),
  `Average Properties` on a 10×10×10 mesh, and `Dump Hamil` / `Dump Overlap`.
- `Fe_up.dat`, `Fe_down.dat`: EHT parameters of each spin channel (Fe 4s, 4p, 3d double
  zeta). Both the Hii and the zetas differ between the channels.

Each channel was fitted with its own Wolfsberg-Helmholz constant K, and the DFT run has
5.119 spin-up and 2.881 spin-down electrons (magnetic moment 2.238 μB). These values go on
the command line:

```
cd Fe_spin/out
cp ../Fe_bands.bind .
yaehmop --parm-up ../Fe_up.dat --parm-down ../Fe_down.dat \
        --const-up 1.752806708980162 --const-down 1.7490786283756976 \
        --nelec-up 5.119 --nelec-down 2.881 \
        Fe_bands.bind

python3 plot_bands.py --emin -8 --emax 6     # -> bands_spin.png
```

The two channels run as two independent calculations. Every output gets `_up` or `_down`
before the input extension:
- `Fe_bands_up.bind.band`, `Fe_bands_down.bind.band`: 9 bands per channel.
- `Fe_bands_up.bind.HAM` / `.OV` and the `_down` pair: H(k), S(k) of each channel, 9×9,
  1000 mesh points. The two `.OV` files differ, because the zetas differ.
- `Fe_bands_up.bind.out`, `Fe_bands_down.bind.out`: header with the channel
  (`#SPIN_CHANNEL`), its parameter file, `#ELECTRONS_THIS_CHANNEL`, `#THE_CONST`, the
  DOS and the `#Fermi_Energy` of the channel.
- `*.status`: progress logs.

Expected Fermi levels: −14.332 eV (up) and −15.055 eV (down). Each channel has its own
Fermi level because its electron count is fixed. The fit put the DFT Fermi level at
−14.259 eV (up) and −14.895 eV (down). The small differences come from the 10×10×10 mesh
without smearing. In the figure each channel is drawn relative to its own Fermi level.

## Notes

- License: the parametrizations come from yaehmop-bandfit (MIT license); the program is
  under the BSD 2-clause license of YAeHMOP (`license.txt` at the top of the repository).

- The diagonal term of H(k), which couples an orbital with its own images in the other
  cells, uses K = 1.75 by default, as in the original program. The new keyword
  `THE CONST DIAG <K>` changes it. Both fits were made with the default, so these inputs
  don't use the keyword.
- The plotting script works in any `out/` folder. It finds the kind of run from the file
  names: `*_up.*band` + `*_down.*band` means spin polarized; `<name>.band` with a
  `<name>.SOC` file means SOC.
