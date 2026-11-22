# interface: YAeHMOP H and S → WanTiBEXOS, and DFT structures → YAeHMOP

**Structure converter.** `kel_crystal_conversion.py` writes a periodic YAeHMOP input from the crystal
structure of another code (one file, Python standard library only):

```bash
python3 kel_crystal_conversion.py --format VASP   --in POSCAR             # -> POSCAR.bind
python3 kel_crystal_conversion.py --format AIMS   --in geometry.in --out MoS2.bind
python3 kel_crystal_conversion.py --format SIESTA --in MoS2.fdf            # also .XV, .STRUCT_OUT
python3 kel_crystal_conversion.py --format DFTB   --in geo_end.gen         # or dftb_in.hsd
python3 kel_crystal_conversion.py --format QE     --in scf.in --band "G 0 0 0, M 0.5 0 0, K 1/3 1/3 0, G 0 0 0" \
        --band-points 40 --kmesh 11 11 1 --average-properties --const 1.75
python3 kel_crystal_conversion.py --help                                   # every option
```

The periodic dimension comes from the vacuum (an empty slab of 6 Å or more across a lattice direction;
`--dimension` forces it); non-periodic directions are moved to the end, and atoms wrapped across the
vacuum are brought back next to their slab, since YAeHMOP has no images along those directions.
Checked on a 4L TiTe2 slab (with one Te wrapped across the vacuum), anatase TiO2 and rutile TiS2: the
YAeHMOP bands of the inputs converted from VASP (direct and Cartesian), FHI-aims (atom and atom_frac),
DFTB+ (.gen, inline hsd), QE (crystal and angstrom) and SIESTA (.fdf Ang and Fractional, .XV,
.STRUCT_OUT) agree with each other to the 1e-5 eV printed in the .band file.

---

# YAeHMOP H and S → WanTiBEXOS

`yaehmop2wtb.py` reads the Hamiltonian and overlap written by YAeHMOP (this fork, `yaehmop-dev`)
and writes them in the non-orthogonal tight-binding format of
[WanTiBEXOS](https://github.com/ac-dias/wantibexos), the same layout as `siesta2wtb.py`,
read with `DFT= S`. It handles the three kinds of run of the fork:

| YAeHMOP run | outputs used | WanTiBEXOS files |
|-------------|--------------|------------------|
| spinless (`bind RUN.bind params.dat`) | `RUN.bind.HAM`, `.OV`, `.out` | `huckel-NP.dat`, `basis_set-NP` |
| SOC (`bind --soc soc.dat RUN.bind params.dat`) | same, with 2n×2n spinor matrices | `huckel-SOC.dat`, `basis_set-SOC` |
| spin polarized (`bind --parm-up ... --parm-down ... RUN.bind`) | `RUN_up.bind.*`, `RUN_down.bind.*` | `huckel-SP.dat`, `basis_set-SP` |

It also writes `input-<TYPE>.dat` and `kpoints.dat` (the band path of `RUN.bind`), ready for
`wtb.x < input-<TYPE>.dat`, and a log `yaehmop2wtb-<TYPE>.log`.

Requirements: Python 3 with numpy and scipy.

## 1. Prepare the YAeHMOP run

The input must contain `Dump Hamil`, `Dump Overlap` and a `K points` block with a
**full Γ-centred mesh** `N1×N2×N3`. YAeHMOP writes H(k) and S(k) at those points, and the
script turns them into H(R), S(R) by an inverse Fourier transform:

    H(R) = (1/N) Σ_k exp(−2πi k·R) H(k),    |R_i| ≤ r_i

YAeHMOP sums the lattice up to `r` cells in each direction (second line of the `Lattice`
section, e.g. `5 5 0`). The transform is therefore **exact** when

- `N_i ≥ 2 r_i + 1` along each periodic direction (11×11×1 for `5 5 0`, 11×11×11 for
  `5 5 5`). With a smaller mesh the script stops; `--allow-aliasing` forces an
  approximate conversion.
- the k points are written with **12 or more decimals** (e.g. `0.090909090909091`).
  YAeHMOP uses them exactly as written, and with 6-8 decimals the rebuilt H(k) is only good
  to ~10⁻⁵-10⁻⁶ eV. The script warns when the points are rounded.

A mesh block can be generated with a few lines of Python:

```python
N = 11
pts = [(i/N, j/N, 0.0) for i in range(N) for j in range(N)]       # 2D; add k/N for 3D
print("K points"); print(len(pts))
for a, b, c in pts: print(f"{a:.15f} {b:.15f} {c:.15f} 1.0")
```

Add a `Band` block too, if you want the script to write the k path for WanTiBEXOS.

## 2. Convert

```
python3 yaehmop2wtb.py path/to/RUN.bind -o wtb_folder [--zero auto|fermi|midgap|<eV>]
```

`RUN.bind` is the input file of the run, with the outputs next to it. The kind of run is
detected from the files: `.SOC` states in the `.out` mean SOC, and `RUN_up`/`RUN_down` files
mean spin polarized. The orbitals of each element come from the `.out` (Extended Hückel
Parameters; double-zeta p shells, listed as `3P DZ:`, included). The structure, the lattice and the k points come from `RUN.bind`
(crystallographic or Cartesian geometry; `*` atoms with parameters given in the input are
not supported).

The script checks itself: it rebuilds H(k) and S(k) on the mesh from H(R), S(R) (error
~10⁻¹² eV when the conditions above hold). It also checks that H(R) is real (spinless and
spin polarized) and that S(R) is real. The log lists these numbers.

**Energy zero.** `wtb.x` diagonalizes (H(k) − E_F·I, S(k)), which is not a rigid shift
when S ≠ I. The file therefore holds H(R) − E0·S(R) with Fermi level 0, so `wtb.x` gives
exactly the YAeHMOP eigenvalues minus E0. These are the choices for E0:

- `auto` (default): the midgap of the mesh eigenvalues if there is a gap, otherwise the
  Fermi level of the `.out` (`Average Properties` needed);
- `fermi`, `midgap`, or a value in eV.

In a spin-polarized run each channel uses its own zero, so that both channels have E_F at
0, as when they were fitted to a DFT calculation.

**File layout** (as `siesta2wtb.py`): `NP`/`SOC`/`SP`, scissors (0), Fermi level (0),
three lattice vectors [Å], basis size, number of R vectors, a comment line, then one line
`Rx Ry Rz [Å]  i  j  ReH  ImH  S` per element. The basis file has
`index species x y z l m spin`. With SOC and SP the orbitals 1..n are spin up (spin 1)
and n+1..2n spin down (spin −1). In SP the up-down blocks are zero.

## 3. Run WanTiBEXOS

```
module load wtb-teste
cd wtb_folder && mkdir -p out
wtb.x < input-SOC.dat          # bands in out/bands.dat (energies relative to the VBM)
```

## Example: MoS2 monolayer (`example_MoS2/`)

The parametrization is the yaehmop-bandfit fit to FHI-aims HSE06+SOC bands (the same one as
`new_tests/MoS2_SOC`). The input `MoS2.bind` has the band path Γ-M-K-Γ (40 points per line)
and the 11×11×1 mesh with 15 decimals. `run_example.sh` does everything:

```
cd example_MoS2
./run_example.sh      # yaehmop-dev (with and without SOC) -> yaehmop2wtb.py -> wtb.x -> compare_bands.py
```

| folder / file | content |
|---------------|---------|
| `MoS2.bind`, `MoS2_params.dat`, `MoS2_soc.dat` | YAeHMOP input, EHT parameters, SOC constants |
| `yaehmop_soc/`, `yaehmop_nosoc/` | YAeHMOP outputs (`.band`, `.HAM`, `.OV`, `.out`, `.SOC`) |
| `wtb_soc/`, `wtb_nosoc/` | converted files (`huckel-*.dat`, `basis_set-*`), `input-*.dat`, `kpoints.dat`, conversion log, `wtb.x` outputs in `out/` |
| `compare_bands.py` | comparison YAeHMOP × WanTiBEXOS, point by point on the same k points |
| `validation.txt`, `bands_yaehmop_vs_wtb.png` | results |

Results (energies relative to the VBM, bands within VBM−6 … VBM+4 eV):

| | without SOC | with SOC |
|---|---|---|
| converted file diagonalized in double precision vs YAeHMOP | 5.4×10⁻⁷ eV | 5.4×10⁻⁷ eV |
| WanTiBEXOS (`wtb.x`) vs YAeHMOP | 4.9×10⁻⁵ eV | 6.6×10⁻⁵ eV |
| band gap, YAeHMOP / WanTiBEXOS | 2.312981 / 2.312991 eV | 2.268108 / 2.268106 eV |
| SOC splitting at K, YAeHMOP / WanTiBEXOS | — | 0.165011 / 0.165014 eV |

![bands](example_MoS2/bands_yaehmop_vs_wtb.png)

The remaining differences come from `wtb.x`, which works in single precision. The converted
file itself reproduces YAeHMOP to 5×10⁻⁷ eV; this is the 8-digit precision of the YAeHMOP
`.band` file. Over all bands the differences reach 10⁻³ eV for `wtb.x`, in states 100 eV or
more above the gap.

**Spin-polarized check (bcc Fe, not included as files).** The Fe fit of `new_tests/Fe_spin`
was run on an 11×11×11 mesh and converted to `huckel-SP.dat`. The rebuilt H(k) error was
10⁻¹² eV. `wtb.x` read the SP file and gave the bands of both channels: the double-precision
bands of the file agree with YAeHMOP to 5×10⁻⁷ eV and the `wtb.x` bands to 10⁻³ eV in the window
E_F−8 … E_F+6 eV. The larger `wtb.x` error than in MoS2 comes from single precision with the
nearly singular S(k) of Fe. Far above E_F (4p states at ~5700 eV) single precision fails
completely, which does not affect the bands near E_F.

## Runs with the newer options of yaehmop-dev

- `bind --soc soc.dat --parm-up up.dat --parm-down down.dat RUN.bind` (SOC on top of two
  spin channels) writes one 2n x 2n spinor dump and is converted as a SOC run
  (`huckel-SOC.dat`).
- `bind --common-fermi --parm-up ... --parm-down ... RUN.bind`: the script reads the common
  Fermi level from `RUN.bind.spin` and uses it as the energy zero of both channels, so the
  spin-up and spin-down bands keep their relative position. Without `--common-fermi`
  (fixed moment) each channel keeps its own Fermi level as zero, as before.
