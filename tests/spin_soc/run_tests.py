#!/usr/bin/env python3
"""Tests of the spin-polarized (--parm-up/--parm-down) and SOC (--soc) runs.

    python3 tests/spin_soc/run_tests.py --bind build/bind [--ref-bind build-ref/bind]

Needs numpy and scipy. Every run happens in a temporary directory.
--ref-bind (a build of the original, unmodified program) adds a regression
test: without the new options the outputs must be byte-identical.

The SOC reference is an independent numpy version of the on-site lambda L.S
model (the same one as parameter_fitter/soc.py in yaehmop-bandfit), applied to
the spinless H(k), S(k) dumps of the program.
"""
import argparse, filecmp, shutil, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
import scipy.linalg as sl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = HERE / "data"
EXAMPLES = ROOT / "examples"
TABLE = ROOT / "tightbind" / "eht_parms.dat"
FAILS = []


def check(name, ok, detail=""):
    print(f"[{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def run(bind, args, cwd, expect_ok=True):
    r = subprocess.run([str(bind)] + [str(a) for a in args], cwd=cwd,
                       capture_output=True, text=True)
    if expect_ok and r.returncode != 0:
        raise RuntimeError(f"bind {args} failed:\n{r.stdout}\n{r.stderr}")
    return r


def read_band(path):
    text = Path(path).read_text()
    nst = int([l for l in text.splitlines() if "orbitals in the unit cell" in l][0].split()[0])
    body = text.split("; Begin band data.")[1].split("#END")[0].strip().splitlines()
    kp, ev, cur = [], [], []
    for l in body:
        if l.startswith("; K point"):
            if cur:
                ev.append(cur)
            cur = []
            kp.append([float(x) for x in l.split(":")[1].split()])
        else:
            cur.append(float(l))
    ev.append(cur)
    ev = np.array(ev)
    assert ev.shape[1] == nst
    return np.array(kp), ev


def read_dump(path):
    """.HAM/.OV: int nk, int n, per k an n*n double array; (i,j), i<j, is m[i,j] - i m[j,i]."""
    raw = Path(path).read_bytes()
    nk, n = np.frombuffer(raw[:8], dtype=np.int32)
    data = np.frombuffer(raw[8:8 + 8 * nk * n * n], dtype=np.float64).reshape(nk, n, n)
    iu = np.triu_indices(n, 1)
    out = np.zeros((nk, n, n), complex)
    for q in range(nk):
        m = data[q]
        H = np.diag(np.diag(m)).astype(complex)
        H[iu] = m[iu] - 1j * m.T[iu]
        H[(iu[1], iu[0])] = np.conj(H[iu])
        out[q] = H
    return out


def grep(path, key):
    return [l for l in Path(path).read_text().splitlines() if key in l]


def fermi(path):
    return float(grep(path, "#Fermi_Energy")[-1].split()[-1])


def net_charges(path):
    lines = Path(path).read_text().splitlines()
    i = lines.index("# AVERAGE NET CHARGES")
    out = []
    for l in lines[i + 1:]:
        p = l.split()
        if len(p) != 3 or not p[0].isdigit():
            break
        out.append(float(p[2]))
    return np.array(out)


# ---------------- independent SOC model (numpy) ----------------
REAL = {"px": (1, "c"), "py": (1, "s"), "pz": (0, "0"),
        "dx2-y2": (2, "c"), "dz2": (0, "0"), "dxy": (2, "s"), "dxz": (1, "c"), "dyz": (1, "s"),
        "fz3": (0, "0"), "fxz2": (1, "c"), "fyz2": (1, "s"), "fxyz": (2, "s"),
        "fz(x2-y2)": (2, "c"), "fx(x2-3y2)": (3, "c"), "fy(3x2-y2)": (3, "s")}
SHELL = {"s": ["s"], "p": ["px", "py", "pz"], "d": ["dx2-y2", "dz2", "dxy", "dxz", "dyz"],
         "f": ["fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)", "fx(x2-3y2)", "fy(3x2-y2)"]}


def angular_momentum(l, names):
    U = np.zeros((len(names), 2 * l + 1), complex)
    for a, nm in enumerate(names):
        m, kind = REAL[nm]
        if kind == "0":
            U[a, l] = 1
        elif kind == "c":
            U[a, l + m] = (-1) ** m / np.sqrt(2); U[a, l - m] = 1 / np.sqrt(2)
        else:
            U[a, l - m] = 1j / np.sqrt(2); U[a, l + m] = -1j * (-1) ** m / np.sqrt(2)
    m = np.arange(-l, l + 1)
    Lp = np.zeros((2 * l + 1, 2 * l + 1), complex)
    for i, mm in enumerate(m[:-1]):
        Lp[i + 1, i] = np.sqrt(l * (l + 1) - mm * (mm + 1))
    Ls = ((Lp + Lp.conj().T) / 2, (Lp - Lp.conj().T) / 2j, np.diag(m).astype(complex))
    return [U.conj() @ L @ U.T for L in Ls]


def soc_reference(symbols, shells_of, lambdas):
    """H_SOC (2n x 2n, spin up first); shells_of[el] = shells in YAeHMOP order."""
    blocks, n = [], 0
    for at, el in enumerate(symbols):
        for sh in shells_of[el]:
            blocks.append((el, sh, n)); n += len(SHELL[sh])
    H = np.zeros((2 * n, 2 * n), complex)
    for el, sh, s0 in blocks:
        lam = lambdas.get((el.upper(), sh), 0.0)
        if sh == "s" or lam == 0.0:
            continue
        Lx, Ly, Lz = angular_momentum("spdf".index(sh), SHELL[sh])
        u = np.arange(s0, s0 + len(SHELL[sh])); d = u + n
        H[np.ix_(u, u)] += lam / 2 * Lz;  H[np.ix_(d, d)] -= lam / 2 * Lz
        H[np.ix_(u, d)] += lam / 2 * (Lx - 1j * Ly); H[np.ix_(d, u)] += lam / 2 * (Lx + 1j * Ly)
    return H


def shells_from_table(path):
    out = {}
    for l in Path(path).read_text().splitlines():
        p = l.split()
        if len(p) >= 6 and not l.startswith(";") and p[5] in "spdf":
            out.setdefault(p[0].upper(), []).append(p[5])
    return {el: [s for s in "spdf" if s in v] for el, v in out.items()}


# ---------------- tests ----------------
def diamond_input(dst):
    """examples/spin_soc/diamond.bind (C from the table, H/S dumps, energy levels)."""
    shutil.copy(EXAMPLES / "spin_soc" / "diamond.bind", dst / "diamond.bind")
    for f in ("C_up.dat", "C_down.dat", "soc.dat"):
        shutil.copy(EXAMPLES / "spin_soc" / f, dst / f)


def test_regression(bind, ref_bind, tmp):
    for name in ("diamond", "XeF4", "Te2Br", "CoNb4Si"):
        out = {}
        for tag, exe in (("ref", ref_bind), ("new", bind)):
            d = tmp / f"reg_{tag}"; d.mkdir(exist_ok=True)
            shutil.copy(EXAMPLES / name, d / name)
            run(exe, [name, TABLE], d)
            out[tag] = d
        for ext in (".out", ".band"):
            a, b = out["ref"] / (name + ext), out["new"] / (name + ext)
            if a.exists():
                check(f"regression {name}{ext} identical to the original program",
                      filecmp.cmp(a, b, shallow=False))


def test_spin_same_tables(bind, tmp):
    d = tmp / "sp_same"; d.mkdir()
    diamond_input(d)
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "diamond.bind", ref)
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--parm-up", TABLE, "--parm-down", TABLE, "diamond.bind"], d)
    for ch in ("up", "down"):
        for ext in (".band", ".HAM", ".OV"):
            check(f"spin, same tables: {ch} {ext} = spinless",
                  filecmp.cmp(d / f"diamond_{ch}.bind{ext}", ref / f"diamond.bind{ext}", shallow=False))
        check(f"spin, same tables: {ch} Fermi = spinless",
              fermi(d / f"diamond_{ch}.bind.out") == fermi(ref / "diamond.bind.out"))
        e = float(grep(d / f"diamond_{ch}.bind.out", "#ELECTRONS_THIS_CHANNEL")[0].split()[-1])
        check(f"spin: {ch} channel has half of the electrons", e == 4.0, f"{e}")


def test_spin_different_tables(bind, tmp):
    d = tmp / "sp_diff"; d.mkdir()
    up, dn = d / "up.dat", d / "dn.dat"
    shutil.copy(TABLE, up)
    txt = TABLE.read_text()
    old = "CO      27      9      2     3     d  -13.180"
    assert old in txt
    dn.write_text(txt.replace(old, "CO      27      9      2     3     d  -12.180"))
    src = (EXAMPLES / "CoNb4Si").read_text()
    src = src.replace("\nprint\n", "\nprint\nnet charges\n").replace("\ncharge \n", "\nDump Hamil\nDump Overlap\ncharge \n")
    (d / "CoNb4Si.bind").write_text(src)
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "CoNb4Si.bind", ref)
    run(bind, ["CoNb4Si.bind", up], ref)
    run(bind, ["--parm-up", up, "--parm-down", dn, "CoNb4Si.bind"], d)
    check("spin, different tables: up .band = spinless (up table)",
          filecmp.cmp(d / "CoNb4Si_up.bind.band", ref / "CoNb4Si.bind.band", shallow=False))
    check("spin, different tables: up .HAM = spinless (up table)",
          filecmp.cmp(d / "CoNb4Si_up.bind.HAM", ref / "CoNb4Si.bind.HAM", shallow=False))
    check("spin, different tables: down .HAM differs",
          not filecmp.cmp(d / "CoNb4Si_down.bind.HAM", ref / "CoNb4Si.bind.HAM", shallow=False))
    check("spin, different tables: down .OV = spinless (same zetas)",
          filecmp.cmp(d / "CoNb4Si_down.bind.OV", ref / "CoNb4Si.bind.OV", shallow=False))
    q_up, q_ref = net_charges(d / "CoNb4Si_up.bind.out"), net_charges(ref / "CoNb4Si.bind.out")
    check("spin: net charges of the up channel = half of the spinless ones",
          np.allclose(2 * q_up, q_ref, atol=2e-6), f"max {np.abs(2*q_up-q_ref).max():.1e}")


def test_spin_const(bind, tmp):
    d = tmp / "sp_const"; d.mkdir()
    diamond_input(d)
    ref = d / "ref"; ref.mkdir()
    txt = (d / "diamond.bind").read_text()
    (ref / "diamond.bind").write_text(txt.replace("\nprint\n", "\nTHE CONST 1.60\nprint\n", 1))
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--parm-up", TABLE, "--parm-down", TABLE, "--const-up", "1.60",
               "diamond.bind"], d)
    check("spin: --const-up 1.60 = spinless run with THE CONST 1.60",
          filecmp.cmp(d / "diamond_up.bind.band", ref / "diamond.bind.band", shallow=False))
    check("spin: down channel keeps the K of the input",
          not filecmp.cmp(d / "diamond_down.bind.band", ref / "diamond.bind.band", shallow=False))


def test_const_diag(bind, tmp):
    d = tmp / "kdiag"; d.mkdir()
    diamond_input(d)
    txt = (d / "diamond.bind").read_text()
    (d / "k.bind").write_text(txt.replace("\nprint\n", "\nTHE CONST 1.60\nTHE CONST DIAG 1.60\nprint\n", 1))
    (d / "k175.bind").write_text(txt.replace("\nprint\n", "\nTHE CONST DIAG 1.75\nprint\n", 1))
    run(bind, ["k.bind", TABLE], d)
    run(bind, ["k175.bind", TABLE], d)
    run(bind, ["diamond.bind", TABLE], d)
    H, S = read_dump(d / "k.bind.HAM"), read_dump(d / "k.bind.OV")
    hii = np.array([-21.4, -11.4, -11.4, -11.4] * 2)      # C s, p (2 atoms)
    K = 1.60
    diff = max(np.abs(np.diag(h).real - hii * (1 - K + K * np.diag(s).real)).max() for h, s in zip(H, S))
    check("THE CONST DIAG 1.60: H_ii(k) = Hii (1 - K + K S_ii(k)) with K = 1.60", diff < 1e-9, f"{diff:.1e}")
    check("THE CONST DIAG 1.75 = default (original program)",
          filecmp.cmp(d / "k175.bind.band", d / "diamond.bind.band", shallow=False))


def test_soc_lambda0(bind, tmp):
    d = tmp / "soc0"; d.mkdir()
    diamond_input(d)
    (d / "soc0.dat").write_text("C p 0.0\n")
    ref = d / "ref"; ref.mkdir()
    shutil.copy(d / "diamond.bind", ref)
    run(bind, ["diamond.bind", TABLE], ref)
    run(bind, ["--soc", "soc0.dat", "diamond.bind", TABLE], d)
    _, e0 = read_band(ref / "diamond.bind.band")
    _, e1 = read_band(d / "diamond.bind.band")
    diff = np.abs(np.sort(np.repeat(e0, 2, axis=1), 1) - e1).max()
    check("SOC lambda=0: 2n bands = spinless bands twice", diff == 0.0, f"{diff:.1e}")
    I2 = np.eye(2)
    for ext in (".HAM", ".OV"):
        A, B = read_dump(ref / f"diamond.bind{ext}"), read_dump(d / f"diamond.bind{ext}")
        diff = max(np.abs(np.kron(I2, a) - b).max() for a, b in zip(A, B))
        check(f"SOC lambda=0: spinor {ext} = I2 x spinless", diff == 0.0, f"{diff:.1e}")
    check("SOC lambda=0: Fermi level unchanged",
          fermi(d / "diamond.bind.out") == fermi(ref / "diamond.bind.out"))


def test_soc_mos2(bind, tmp):
    d = tmp / "mos2"; d.mkdir()
    params, socf = DATA / "MoS2_params.dat", DATA / "MoS2_soc.dat"
    shutil.copy(DATA / "MoS2_path.bind", d / "path.bind")
    run(bind, ["path.bind", params], d)
    _, e0 = read_band(d / "path.bind.band")
    run(bind, ["--soc", socf, "path.bind", params], d)
    kp, esoc = read_band(d / "path.bind.band")
    # the same k points as a K points list: spinless and spinor dumps
    head = (d / "path.bind").read_text().split("\nBand\n")[0]
    klist = "\n".join(f"{k[0]:.6f} {k[1]:.6f} {k[2]:.6f} 1.0" for k in kp)
    (d / "kl.bind").write_text(head + f"\n\nDump Hamil\nDump Overlap\n\nK points\n{len(kp)}\n{klist}\n")
    run(bind, ["kl.bind", params], d)
    H, S = read_dump(d / "kl.bind.HAM"), read_dump(d / "kl.bind.OV")
    run(bind, ["--soc", socf, "kl.bind", params], d)
    Hc, Sc = read_dump(d / "kl.bind.HAM"), read_dump(d / "kl.bind.OV")
    lam = {}
    for l in socf.read_text().splitlines():
        if l.strip() and not l.startswith("#"):
            el, sh, v = l.split(); lam[(el.upper(), sh)] = float(v)
    Hsoc = soc_reference(["MO", "S", "S"], shells_from_table(params), lam)
    I2 = np.eye(2)
    dH = max(np.abs(hc - (np.kron(I2, h) + Hsoc)).max() for hc, h in zip(Hc, H))
    dS = max(np.abs(sc - np.kron(I2, s)).max() for sc, s in zip(Sc, S))
    check("SOC MoS2: spinor H(k) = I2 x H(k) + H_SOC (numpy model)", dH < 1e-10, f"{dH:.1e}")
    check("SOC MoS2: spinor S(k) = I2 x S(k)", dS < 1e-12, f"{dS:.1e}")
    epy = np.array([sl.eigh(np.kron(I2, h) + Hsoc, np.kron(I2, s), eigvals_only=True)
                    for h, s in zip(H, S)])
    nocc = 18                       # Mo 6 + 2 x S 6 electrons, one per spinor state
    vbm = esoc[:, :nocc].max()
    win = (esoc > vbm - 6) & (esoc < vbm + 4)
    dE = np.abs(esoc - epy)[win].max()
    check("SOC MoS2: bands (VBM-6..VBM+4 eV) = numpy model", dE < 1e-4, f"{dE:.1e} eV")
    iK = np.argmin(np.linalg.norm(kp - np.array([1 / 3, 1 / 3, 0]), axis=1))
    split = esoc[iK, nocc - 1] - esoc[iK, nocc - 2]
    gap0 = e0[:, nocc // 2:].min() - e0[:, :nocc // 2].max()
    gap = esoc[:, nocc:].min() - esoc[:, :nocc].max()
    check("SOC MoS2: VB splitting at K = 0.165 eV (yaehmop-bandfit fit)", abs(split - 0.1650) < 1e-3, f"{split:.4f}")
    check("SOC MoS2: gap 2.313 -> 2.268 eV (yaehmop-bandfit fit)",
          abs(gap0 - 2.3130) < 1e-3 and abs(gap - 2.2681) < 1e-3, f"{gap0:.4f} -> {gap:.4f}")
    check("SOC MoS2: .SOC file written", (d / "kl.bind.SOC").exists())


def test_errors(bind, tmp):
    d = tmp / "err"; d.mkdir()
    diamond_input(d)
    shutil.copy(EXAMPLES / "XeF4", d / "XeF4")
    r = run(bind, ["--parm-up", "C_up.dat", "diamond.bind"], d, expect_ok=False)
    check("error: --parm-up without --parm-down", r.returncode != 0 and "both" in r.stderr)
    r = run(bind, ["--soc", "soc.dat", "--parm-up", "C_up.dat", "--parm-down", "C_down.dat",
                   "diamond.bind"], d, expect_ok=False)
    check("error: SOC + spin polarization refused", r.returncode != 0 and "not implemented" in r.stderr)
    run(bind, ["--soc", "soc.dat", "XeF4", TABLE], d, expect_ok=False)
    check("error: SOC in a molecule refused", any("FATAL" in l for l in grep(d / "XeF4.status", "FATAL")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bind", required=True, help="bind built from this branch")
    ap.add_argument("--ref-bind", help="bind of the original program (regression test)")
    ap.add_argument("--keep", action="store_true", help="keep the temporary directory")
    a = ap.parse_args()
    bind = Path(a.bind).resolve()
    tmp = Path(tempfile.mkdtemp(prefix="yaehmop_spin_soc_"))
    try:
        if a.ref_bind:
            test_regression(bind, Path(a.ref_bind).resolve(), tmp)
        test_spin_same_tables(bind, tmp)
        test_spin_different_tables(bind, tmp)
        test_spin_const(bind, tmp)
        test_const_diag(bind, tmp)
        test_soc_lambda0(bind, tmp)
        test_soc_mos2(bind, tmp)
        test_errors(bind, tmp)
    finally:
        if a.keep:
            print("outputs in", tmp)
        else:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{len(FAILS)} failure(s)" if FAILS else "\nall tests passed")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
