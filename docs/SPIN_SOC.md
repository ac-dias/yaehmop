# Spin-polarized runs and spin-orbit coupling (adias-dev)

Two additions to `bind`, made in the fork `ac-dias/yaehmop` (branch `adias-dev`).
Runs without the new options give the same output as the original program.

```
bind [options] <inputfile> [paramfile]

  --parm-up FILE --parm-down FILE   spin-polarized run, one parameter table per channel
  --nelec-up X   --nelec-down X     electrons in each channel (default: half of the
                                    electrons of the input in each one)
  --soc FILE                        spin-orbit coupling (on-site lambda L.S)
  -h, --help                        list the options
```

## Spin-polarized runs

```
bind --parm-up up.dat --parm-down down.dat MoS2.bind
```

Each spin channel is an independent extended Hückel calculation that uses its own
parameter table (same format as `eht_parms.dat`), so the channels can have different
Hii values, and different zetas too. The two channels run in two separate processes
at the same time.

**Outputs.** The suffix goes before the extension of the input name. Every output gets it:

| spinless              | spin up                  | spin down                  |
|-----------------------|--------------------------|----------------------------|
| `MoS2.bind.out`       | `MoS2_up.bind.out`       | `MoS2_down.bind.out`       |
| `MoS2.bind.band`      | `MoS2_up.bind.band`      | `MoS2_down.bind.band`      |
| `MoS2.bind.HAM`, `.OV`| `MoS2_up.bind.HAM`, `.OV`| `MoS2_down.bind.HAM`, `.OV`|
| `MoS2.bind.status`    | `MoS2_up.bind.status`    | `MoS2_down.bind.status`    |

(An input without an extension, `diamond`, gives `diamond_up.out`, `diamond_up.band`, ...)

**Occupations.** Each channel holds one electron per state. By default each channel gets
half of the electrons of the input (`Electrons` or `Charge` keyword). To fix the moment,
use `--nelec-up` and `--nelec-down`. Each channel has its own Fermi level
(`#Fermi_Energy` in its `.out`). A Fermi level shared by both channels is not
implemented yet.

**Charges.** Net charges are given per channel: half of the atom's valence electrons
minus the channel population. Add up and down to get the total. Charge iteration is
refused in spin-polarized runs.

**H and S.** With `Dump Hamil` / `Dump Overlap` in the input, each channel writes its own
`.HAM`/`.OV` (n×n, same binary format as without spin). If both tables use the same zetas,
the two `.OV` files are identical.

The `.out` header of each channel records the channel (`#SPIN_CHANNEL: up`), its
parameter file and its electron count (`#ELECTRONS_THIS_CHANNEL`).

## Spin-orbit coupling

```
bind --soc soc.dat MoS2.bind params.dat
```

`soc.dat` has one line per element and shell, with lambda in eV. Lines starting with `#`
or `;` are comments. This is the same format as the `soc_parameters.txt` written by
yaehmop-bandfit:

```
# element shell lambda_eV
Mo  d  0.063352
S   p  0.197577
```

Model (the same as `parameter_fitter/soc.py` in yaehmop-bandfit). The coupling acts only
between orbitals of the same shell on the same atom:

    H_SOC = sum over atoms A and shells l >= 1 of  lambda_{A,l} L.S
    H_spinor(k) = I2 (x) H(k) + H_SOC,      S_spinor(k) = I2 (x) S(k)

The spinor basis has 2n states: the n orbitals with spin up (YAeHMOP order), then the same
n orbitals with spin down. A p shell splits into j=3/2 (+lambda/2) and j=1/2 (-lambda). A d
shell splits into j=5/2 (+lambda) and j=3/2 (-3lambda/2). Each spinor state holds one
electron.

**Outputs with SOC**

- `.band`: 2n bands per k point. The header says `2n orbitals in the unit cell.`, so
  existing readers work unchanged.
- `.HAM`, `.OV` (with `Dump Hamil` / `Dump Overlap` and a `K points` block): spinor
  H(k) and S(k), 2n×2n, in the usual binary format. That format is: int nk, int dim, then
  per k a dim×dim double array, where element (i,j) with i<j is `m[i][j] - i*m[j][i]`.
  The `read_matrix_dump` reader of yaehmop-bandfit reads these files directly.
- `.SOC`: text copy of the on-site H_SOC (nonzero elements `i j Re Im`). With it, the
  spinor H(R) can be rebuilt from a spinless run: `I2 (x) H(R) + delta_{R,0} H_SOC`.
- `.out`: the lambdas used for each atom, `#SOC_STATES: 2n`, the Fermi level and the
  total DOS (with `Average Properties`), and the energies per k point (with
  `print energy levels`).

**Not implemented with SOC yet.** The program stops with a message naming the option:
projected DOS, COOP/COHP, FMO/FCO, populations and charges, charge iteration, Walsh
diagrams, MO/wave-function printing, printing of H/S in the `.out`, molecules, THIN mode,
and SOC together with spin polarization. SOC needs the build with LAPACK
(`-DUSE_BLAS_LAPACK=ON`).

## Tests done (2026-11-20)

The tests are automated in `tests/spin_soc/run_tests.py` and need only numpy and scipy:

```
python3 tests/spin_soc/run_tests.py --bind build/bind [--ref-bind <bind of the original program>]
```

The SOC reference inside the script is an independent numpy version of the lambda L.S
model. The MoS2 inputs (`tests/spin_soc/data/`) come from the yaehmop-bandfit validation
fit `aims/MoS2_monolayer_hse06_soc` (EHT parameters, lambdas and structure).

- Without the new options: `.out` and `.band` of the examples diamond, XeF4, Te2Br and
  CoNb4Si are byte-identical to the original program.
- Spin polarization with the same table in both channels: each channel's `.band` and
  `.HAM` equal the spinless run. Energies, occupations and net charges are exactly half.
  With a different Co d level in the down channel (CoNb4Si), only that channel changes.
- SOC with lambda = 0: each spinless band appears twice, and the spinor `.HAM`/`.OV` are
  exactly I2 (x) H(k), I2 (x) S(k). The Fermi level is unchanged.
- SOC on MoS2 (yaehmop-bandfit fit `aims/MoS2_monolayer_hse06_soc`): the spinor H(k)
  and S(k) from the C code are identical (difference 0.0) to `soc.py` applied to the
  spinless dumps. Band differences are up to 2e-5 eV in the valence and conduction
  window; this is the same size as the original program's own difference to scipy, and
  comes from LAPACK and the 8-digit `.band` format. The VB splitting at K is 0.1650 eV
  and the gap is 2.2681 eV, the same values as the bandfit `fit_summary.json`.

## Example

`examples/spin_soc/`: diamond with a carbon table per channel (`C_up.dat`,
`C_down.dat`) and `soc.dat`:

```
bind diamond.bind C_up.dat                              # spinless
bind --soc soc.dat diamond.bind C_up.dat                # SOC
bind --parm-up C_up.dat --parm-down C_down.dat diamond.bind   # spin polarized
```
