#!/usr/bin/env python3
"""YAeHMOP H and S -> WanTiBEXOS tight-binding file and basis file (lean version).

    python3 yaehmop2wtb.py RUN.bind [-o OUTDIR] [--efermi EV] [--rmatrix [--params PARMS [PARMS_DOWN]]]

Reads only the results of a YAeHMOP run made with 'Dump Hamil', 'Dump Overlap'
and a full Gamma-centred 'K points' mesh (outputs next to RUN.bind) and writes

  huckel-<TYPE>.dat    H(R) - E_F S(R) and S(R) in the non-orthogonal WanTiBEXOS
                       layout (as siesta2wtb.py, DFT= S), Fermi level 0
  basis_set-<TYPE>     index species x y z l m spin
  huckel-<TYPE>_r.dat  with --rmatrix: the position matrix <m,0|r|n,R> of the
                       Slater basis (Angstrom), in the layout of siesta2wtb.py
                       --rmatrix (Wannier90 seedname_r.dat; BSE_CENTER_FILE)

TYPE is NP (spinless), SOC (bind --soc, 2n x 2n spinors) or SP (RUN_up.bind.*
and RUN_down.bind.*), detected from the files. yaehmop2wtb-full.py does the
same plus the energy-zero choices (midgap), the input and k-path files for
wtb.x and a log.

H(R) = (1/N) sum_k exp(-2 pi i k.R) H(k) for |R_i| <= r_i (the overlap radius
of the Lattice section), exact when N_i >= 2 r_i + 1 and the k points are
written with 12 or more decimals. The energy zero E_F is the Fermi level of the
.out ('Average Properties'), the common one of RUN.bind.spin for a spin
polarized run with bind --common-fermi, or --efermi.

--rmatrix: X_ij(R) = <phi_i,0| r |phi_j,R>, absolute coordinates in the frame
of basis_set-*, spin diagonal. The orbitals are YAeHMOP's normalized Slater
functions, c1 r^(n-1) e^(-z1 r) [+ c2 r^(n-1) e^(-z2 r)] times real spherical
harmonics, integrated numerically (bipolar coordinates for two centres,
spherical for one) for the atom pairs with a non-zero overlap in YAeHMOP. The
same quadrature gives S, compared with YAeHMOP's (~1e-7 with --params). The
.out prints the exponents with 4 decimals (S then agrees to ~5e-5); --params
gives the parameter file of the run for the full precision (SP: the up and the
down file, since the two channels have their own orbitals, and the same for
SOC on two tables with different zetas, where it is required; not with charge
iteration, where the .out has the final exponents). ~10-40 s for a small cell.

Needs numpy.
"""
import argparse, itertools, re, sys
from math import factorial
from pathlib import Path
import numpy as np

BOHR = 1 / 1.889644746             # Angstrom: the AUI of YAeHMOP's overlaps (R_overlap_mat.c),
                                   # not its BOHR = 0.52918 (bind.h)
ORBITALS = {"s": ["s"], "p": ["px", "py", "pz"],
            "d": ["dx2-y2", "dz2", "dxy", "dxz", "dyz"],
            "f": ["fz3", "fxz2", "fyz2", "fxyz", "fz(x2-y2)", "fx(x2-3y2)", "fy(3x2-y2)"]}
ORBITAL_LM = {"s": (0, 0), "px": (1, 1), "py": (1, -1), "pz": (1, 0),
              "dx2-y2": (2, 2), "dz2": (2, 0), "dxy": (2, -2), "dxz": (2, 1), "dyz": (2, -1),
              "fz3": (3, 0), "fxz2": (3, 1), "fyz2": (3, -1), "fxyz": (3, -2),
              "fz(x2-y2)": (3, 2), "fx(x2-3y2)": (3, 3), "fy(3x2-y2)": (3, -3)}


# ------------------------------------------------------------------ YAeHMOP files
def read_bind(path):
    """Atoms, lattice vectors [A], periodic dimension, overlap radius, k points."""
    lines = [l.split(";")[0].strip() for l in Path(path).read_text().splitlines()]
    lines = [l for l in lines if l]
    up = [l.upper() for l in lines]
    info = dict(kpoints=None, crystal=None)
    i = 0
    while i < len(lines):
        u = up[i]
        if u.startswith("GEOMETRY"):
            info["frac"] = "CRYSTAL" in u
            n = int(lines[i + 1].split()[0])
            info["atoms"] = [(p[1], np.array(p[2:5], float))
                             for p in (l.split() for l in lines[i + 2:i + 2 + n])]
            i += 2 + n
        elif u == "LATTICE":
            dim = int(lines[i + 1].split()[0])
            info.update(dim=dim, overlaps=[int(x) for x in lines[i + 2].split()[:3]],
                        ends=[tuple(int(x) for x in lines[i + 3 + j].split()[:2]) for j in range(dim)])
            i += 3 + dim
        elif u.startswith("CRYSTAL SPEC"):
            info["crystal"] = [[float(x) for x in lines[i + j].split()[:3]] for j in (1, 2)]
            i += 3
        elif u.startswith("K POINTS") or u == "KPOINTS":
            n = int(lines[i + 1].split()[0])
            info["kpoints"] = np.array([l.split()[:3] for l in lines[i + 2:i + 2 + n]], float)
            i += 2 + n
        else:
            i += 1
    if "atoms" not in info or "dim" not in info:
        raise SystemExit(f"{path}: needs Geometry and Lattice sections")
    if info["kpoints"] is None:
        raise SystemExit(f"{path}: no 'K points' block (the H/S dumps are written on those k points)")
    atoms, dim = info["atoms"], info["dim"]
    if info["frac"]:
        (a, b, c), angles = info["crystal"]
        al, be, ga = np.radians(angles)
        cx, cy = c * np.cos(be), c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
        cell = np.array([[a, 0, 0], [b * np.cos(ga), b * np.sin(ga), 0],
                         [cx, cy, np.sqrt(max(c * c - cx * cx - cy * cy, 0.0))]])
        atoms = [(s, f @ cell) for s, f in atoms]
    else:
        cell = np.zeros((3, 3))
        for j, (b, e) in enumerate(info["ends"]):
            cell[j] = atoms[e - 1][1] - atoms[b - 1][1]
        pos = np.array([p for _, p in atoms[:len(atoms) - dim]])
        for j in range(dim, 3):              # vacuum directions: orthogonal, 15 A of vacuum
            if j == 1:
                t = np.array([0.0, 0.0, 1.0])
                if abs(cell[0] @ t) > 0.9 * np.linalg.norm(cell[0]):
                    t = np.array([1.0, 0.0, 0.0])
                v = np.cross(cell[0], t)
            else:
                v = np.cross(cell[0], cell[1])
            v /= np.linalg.norm(v)
            cell[j] = v * max(20.0, (np.ptp(pos @ v) if len(pos) > 1 else 0.0) + 15.0)
    real = [(s, p) for s, p in atoms[:len(atoms) - dim] if not s.startswith("&")]
    if any(s.startswith("*") for s, _ in real):
        raise SystemExit("atoms with parameters given in the input ('*') are not supported")
    info.update(real=real, cell=cell)
    return info


def read_out(path):
    """Shells of each element (or atom, with charge iteration), Fermi level, SOC flag."""
    text = Path(path).read_text()
    shells, cur = {}, None
    for l in text.splitlines():
        m = re.match(r"ATOM:\s+(\S+)\s+Atomic number:", l)
        if m:
            cur = m.group(1).upper()
            shells.setdefault(cur, {})
            continue
        m = re.match(r"\s+(\d+)([SPDF])(?: DZ)?:\s+(.*)", l)
        if cur and m:
            v = [float(x) for x in m.group(3).split()]
            sh = m.group(2).lower()
            # s, p: Hii z1;  p DZ, d, f: Hii c1 z1 c2 z2  (coefficients already normalized)
            shells[cur][sh] = dict(n=int(m.group(1)), stos=[(1.0, v[1])] if len(v) == 2 else
                                   [(c, z) for c, z in ((v[1], v[2]), (v[3], v[4])) if c != 0.0])
            continue
        if cur and not l.startswith("\t"):
            cur = None
    fermi = re.findall(r"#Fermi_Energy:\s*(\S+)", text)
    return dict(shells=shells, fermi=float(fermi[-1]) if fermi else None,
                soc=bool(re.search(r"#SOC_STATES:", text)),
                # SOC on two parameter tables with different zetas: S = diag(S_up, S_down)
                soc_two_zetas=bool(re.search(r"#SOC_SPIN_ZETAS:\s*different", text)))


def read_params(path):
    """Slater functions from a YAeHMOP parameter file, normalized as bind does."""
    shells = {}
    for l in Path(path).read_text().splitlines():
        p = l.split(";")[0].split()
        if len(p) < 11 or p[0].upper() == "END":
            continue
        sym, nzeta, nq, ang = p[0].upper(), int(p[3]), int(p[4]), p[5].lower()
        hii, z1, z2, c1, c2 = (float(x) for x in p[6:11])
        if hii == 0.0:
            continue
        if nzeta == 2:
            s12 = (4 * z1 * z2 / (z1 + z2) ** 2) ** (nq + 0.5)
            norm = np.sqrt(c1 * c1 + c2 * c2 + 2 * s12 * c1 * c2)
            stos = [(c1 / norm, z1), (c2 / norm, z2)]
        else:
            stos = [(1.0 if ang in "sp" else c1, z1)]
        shells.setdefault(sym, {}).setdefault(ang, dict(n=nq, stos=stos))
    return shells


def read_dump(path):
    """H(k) or S(k): int nk, int n, then per k an n x n doubles array in which
    element (i,j), i<j, is m[i][j] - i m[j][i]."""
    raw = Path(path).read_bytes()
    nk, n = (int(x) for x in np.frombuffer(raw[:8], dtype=np.int32))
    m = np.frombuffer(raw[8:8 + 8 * nk * n * n], dtype=np.float64)
    if m.size != nk * n * n:
        raise SystemExit(f"{path}: shorter than its header says ({nk} k points, n = {n})")
    m = m.reshape(nk, n, n)
    iu = np.triu_indices(n, 1)
    out = np.zeros((nk, n, n), complex)
    out[:, np.arange(n), np.arange(n)] = m[:, np.arange(n), np.arange(n)]
    out[:, iu[0], iu[1]] = m[:, iu[0], iu[1]] - 1j * m[:, iu[1], iu[0]]
    out[:, iu[1], iu[0]] = np.conj(out[:, iu[0], iu[1]])
    return out


# ------------------------------------------------------------------ H(R), S(R)
def mesh_size(kpts, dim, ov):
    N = []
    for d in range(3):
        if d >= dim:
            N.append(1)
            continue
        vals = np.unique(np.round(np.mod(kpts[:, d] + 1e-9, 1.0), 6))
        if not np.allclose(vals, np.arange(len(vals)) / len(vals), atol=2e-6):
            raise SystemExit(f"the k points are not a uniform Gamma-centred mesh along direction {d + 1}")
        N.append(len(vals))
    if len(kpts) != np.prod(N):
        raise SystemExit(f"{len(kpts)} k points, but a {N[0]}x{N[1]}x{N[2]} mesh has {np.prod(N)}")
    if any(N[d] < 2 * ov[d] + 1 for d in range(dim)):
        raise SystemExit(f"the mesh {N} is too small for the overlap radius {ov}: "
                         f"use N_i >= 2 r_i + 1 (or yaehmop2wtb-full.py --allow-aliasing)")
    dev = max(np.abs(kpts[:, d] - np.round(kpts[:, d] * N[d]) / N[d]).max() for d in range(dim))
    if dev > 1e-9:
        print(f"WARNING: k points rounded by up to {dev:.1e}: write them with 12 or more decimals "
              "for an exact H(R)", file=sys.stderr)
    return N


def real_space(Hk, Sk, kpts, dim, ov, e0, cutoff):
    """[(R, H(R) - e0 S(R), S(R))] for the R vectors with some |H| or |S| > cutoff,
    and the largest error of H(k), S(k) rebuilt on the mesh."""
    rng = [range(-ov[d], ov[d] + 1) if d < dim else (0,) for d in range(3)]
    cells = []
    Hr, Sr = np.zeros_like(Hk), np.zeros_like(Sk)
    for R in itertools.product(*rng):
        R = np.array(R)
        ph = np.exp(-2j * np.pi * (kpts @ R))
        h = np.tensordot(ph, Hk, 1) / len(kpts)
        s = np.tensordot(ph, Sk, 1) / len(kpts)
        Hr += np.conj(ph)[:, None, None] * h
        Sr += np.conj(ph)[:, None, None] * s
        if max(np.abs(h).max(), np.abs(s).max()) > cutoff:
            cells.append((R, h - e0 * s, s.real))
    if max(np.abs(s.imag).max() for _, _, s in cells) > 1e-6:
        raise SystemExit("S(R) is not real")
    return cells, np.abs(Hr - Hk).max(), np.abs(Sr - Sk).max()


# ------------------------------------------------------------------ position matrix
def _ylm(orb, v, r):
    """Real spherical harmonics (the orbital names of YAeHMOP) at the points v (|v| = r)."""
    x, y, z = (v[:, i] / r for i in range(3))
    pi = np.pi
    return {"s": 0 * x + np.sqrt(1 / (4 * pi)),
            "px": np.sqrt(3 / (4 * pi)) * x, "py": np.sqrt(3 / (4 * pi)) * y,
            "pz": np.sqrt(3 / (4 * pi)) * z,
            "dx2-y2": np.sqrt(15 / (16 * pi)) * (x * x - y * y),
            "dz2": np.sqrt(5 / (16 * pi)) * (3 * z * z - 1),
            "dxy": np.sqrt(15 / (4 * pi)) * x * y, "dxz": np.sqrt(15 / (4 * pi)) * x * z,
            "dyz": np.sqrt(15 / (4 * pi)) * y * z,
            "fz3": np.sqrt(7 / (16 * pi)) * z * (5 * z * z - 3),
            "fxz2": np.sqrt(21 / (32 * pi)) * x * (5 * z * z - 1),
            "fyz2": np.sqrt(21 / (32 * pi)) * y * (5 * z * z - 1),
            "fxyz": np.sqrt(105 / (4 * pi)) * x * y * z,
            "fz(x2-y2)": np.sqrt(105 / (16 * pi)) * z * (x * x - y * y),
            "fx(x2-3y2)": np.sqrt(35 / (32 * pi)) * x * (x * x - 3 * y * y),
            "fy(3x2-y2)": np.sqrt(35 / (32 * pi)) * y * (3 * x * x - y * y)}[orb]


class Atom:
    """Slater orbitals of one atom (Angstrom units) and their extent."""

    def __init__(self, shells, tol=1e-7):
        self.orbs, self.radial, rc = [], [], 0.0
        for sh in "spdf":
            if sh not in shells:
                continue
            n = shells[sh]["n"]
            # (c, N, z) per Slater function, z in 1/A, N = (2z)^(n+1/2) / sqrt((2n)!)
            stos = [(c, (2 * z / BOHR) ** (n + 0.5) / np.sqrt(factorial(2 * n)), z / BOHR)
                    for c, z in shells[sh]["stos"]]
            for orb in ORBITALS[sh]:
                self.orbs.append(orb)
                self.radial.append((n, stos))
            r = np.linspace(0.01, 60.0, 6000)
            f = np.abs(sum(c * N * r ** (n - 1) * np.exp(-z * r) for c, N, z in stos)) * r
            rc = max(rc, r[np.nonzero(f > tol)[0][-1]] if np.any(f > tol) else 1.0)
        self.rc = rc

    def values(self, v):
        r = np.maximum(np.linalg.norm(v, axis=1), 1e-12)
        out = []
        for orb, (n, stos) in zip(self.orbs, self.radial):
            rad = sum(c * N * r ** (n - 1) * np.exp(-z * r) for c, N, z in stos)
            out.append(rad * _ylm(orb, v, r))
        return np.array(out)


def _gauss(lo, hi, h, n, cuts=()):
    """Gauss-Legendre on [lo, hi] in pieces split at the cuts (kinks of the
    integrand), graded near r = 0 (0.05, 0.1, 0.2, ... A, for the compact
    components of d and f shells) and of at most h elsewhere."""
    e = np.concatenate([0.05 * 2.0 ** np.arange(int(np.log2(h / 0.05)) + 1),
                        np.arange(h, hi + h, h), np.asarray(cuts, float)])
    e = np.unique(np.concatenate([[lo, hi], e[(e > lo + 1e-10) & (e < hi - 1e-10)]]))
    x, w = np.polynomial.legendre.leggauss(n)
    half = 0.5 * (e[1:] - e[:-1])[:, None]
    return (half * x + 0.5 * (e[1:] + e[:-1])[:, None]).ravel(), (half * w).ravel()


def _quadrature(A, B, Ra, Rb, h, n, nphi=12):
    """Points and weights over the region where both atoms' orbitals are non-zero:
    spherical about A for A = B, bipolar (r_a, r_b, phi) about the A-B axis
    otherwise (volume element r_a r_b / d)."""
    phi = 2 * np.pi * np.arange(nphi) / nphi
    d = np.linalg.norm(B - A)
    if d < 1e-8:
        r, wr = _gauss(0.0, min(Ra, Rb), h, n)
        ct, wt = np.polynomial.legendre.leggauss(10)
        st = np.sqrt(1 - ct ** 2)
        u = np.stack(np.broadcast_arrays(st[:, None] * np.cos(phi), st[:, None] * np.sin(phi),
                                         ct[:, None] + 0 * phi), -1).reshape(-1, 3)
        w = (wt[:, None] + 0 * phi).ravel() * 2 * np.pi / nphi
        return (A + r[:, None, None] * u[None]).reshape(-1, 3), (wr[:, None] * r[:, None] ** 2 * w).ravel()
    ez = (B - A) / d
    t = np.array([1.0, 0, 0]) if abs(ez[0]) < 0.9 else np.array([0, 1.0, 0])
    ex = t - ez * (t @ ez)
    ex /= np.linalg.norm(ex)
    circle = np.cos(phi)[:, None] * ex + np.sin(phi)[:, None] * np.cross(ez, ex)
    P, W = [], []
    # the r_b range [|r_a - d|, min(r_a + d, Rb)] has kinks at r_a = d and r_a = Rb - d
    ra, wa = _gauss(0.0, Ra, h, n, (d, Rb - d))
    for x, wx in zip(ra, wa):
        lo, hi = abs(x - d), min(x + d, Rb)
        if hi - lo < 1e-10:
            continue
        rb, wb = _gauss(lo, hi, h, n)
        z = (x * x - rb * rb + d * d) / (2 * d)
        rho = np.sqrt(np.clip(x * x - z * z, 0, None))
        P.append((A + z[:, None, None] * ez + rho[:, None, None] * circle).reshape(-1, 3))
        W.append(np.repeat(wx * x * rb * wb / d * 2 * np.pi / nphi, nphi))
    return np.concatenate(P), np.concatenate(W)


def position_matrix(real, sites, cell, Rs, Sref, h=2.0, n=10):
    """X[iR, i, j] = <phi_i,0| r |phi_j,R> and the overlap Sq of the same quadrature,
    for the atom pairs whose block of YAeHMOP's S(R) is non-zero."""
    first = np.cumsum([0] + [len(a.orbs) for a in sites])
    no = first[-1]
    X, Sq = np.zeros((len(Rs), no, no, 3)), np.zeros((len(Rs), no, no))
    for iR, R in enumerate(Rs):
        for a, (_, A) in enumerate(real):
            ia = slice(first[a], first[a + 1])
            for b, (_, B0) in enumerate(real):
                ib = slice(first[b], first[b + 1])
                if not np.any(Sref[iR][ia, ib]):
                    continue
                B = B0 + R @ cell
                if np.linalg.norm(B - A) >= sites[a].rc + sites[b].rc:
                    continue
                pts, w = _quadrature(A, B, sites[a].rc, sites[b].rc, h, n)
                fa = sites[a].values(pts - A) * w
                fb = sites[b].values(pts - B)
                Sq[iR, ia, ib] = fa @ fb.T
                for c in range(3):
                    X[iR, ia, ib, c] = (fa * pts[:, c]) @ fb.T
    return X, Sq


def rmatrix_file(path, real, cell, cells, no, kind, source, channel_shells, origin):
    """Writes the position matrix file for the R vectors of cells [(R, H, S)].
    channel_shells: the shells of the run (NP, SOC) or of the up and down
    channels (SP, or SOC on two tables with different zetas: two sets of
    orbitals; X and S are then computed per channel, against the S block of
    that channel).
    Returns the largest |S(quadrature) - S(YAeHMOP)|."""
    Rs = [tuple(int(x) for x in R) for R, _, _ in cells]
    Xs, dS = [], 0.0
    for c, shells in enumerate(channel_shells):
        sites = []
        for ia, (sym, _) in enumerate(real):
            sh = shells.get(f"{sym.upper()}({ia + 1})", shells.get(sym.upper()))
            if sh is None:
                raise SystemExit(f"element {sym} not in {origin}")
            sites.append(Atom(sh))
        if sum(len(s.orbs) for s in sites) != no:
            raise SystemExit(f"the orbitals of {origin} differ from the ones of the run")
        blk = slice(c * no, (c + 1) * no)
        Sref = [np.real(S[blk, blk]) for _, _, S in cells]
        X, Sq = position_matrix(real, sites, cell, Rs, Sref)
        Xs.append(X)
        dS = max(dS, max(np.abs(Sq[i] - Sref[i]).max() for i in range(len(Rs))))
    if kind == "SOC" and len(Xs) == 1:
        Xs = Xs * 2                     # same orbitals for both spins
    write_rmatrix(path, Rs, Xs, source)
    return dS


def channel_shells(kind, outs, params):
    """Slater functions for --rmatrix: from --params (one file; two, up and down,
    for SP and for SOC on two tables with different zetas) or from the .out of
    each run; and a description of where from."""
    two = kind == "SP" or outs[0]["soc_two_zetas"]
    if params:
        if len(params) != (2 if two else 1):
            raise SystemExit("--params takes one parameter file (two, up and down, for SP and for "
                             "SOC on two tables with different zetas)")
        return [read_params(p) for p in params], " + ".join(params)
    if kind == "SOC" and two:
        raise SystemExit("SOC on two tables with different zetas: the .out does not list the "
                         "down functions in full; give --params UP DOWN")
    return [o["shells"] for o in outs], "the .out (4 decimals)"


# ------------------------------------------------------------------ output
def write_tb(path, kind, cell, nb, cells):
    with open(path, "w") as f:
        f.write(f"{kind}\n0.000000\n0.000000\n")
        for v in cell:
            f.write(f"{v[0]:.8f}  {v[1]:.8f}  {v[2]:.8f}\n")
        f.write(f"{nb}\n{len(cells)}\n#rcell x   rcell y   rcell z   i   j   ReH   ImH   S\n")
        for R, H, S in cells:
            r = R @ cell
            for i in range(nb):
                for j in range(nb):
                    f.write(f"{r[0]:.8f}   {r[1]:.8f}   {r[2]:.8f}   {i + 1}   {j + 1}   "
                            f"{H[i, j].real:.10f}   {H[i, j].imag:.10f}   {S[i, j]:.10f}\n")


def write_basis(path, labels, spins):
    with open(path, "w") as f:
        f.write("#bindex aspecie ax ay az l m spin\n")
        b = 0
        for spin in spins:
            for sym, p, orb in labels:
                b += 1
                l, m = ORBITAL_LM[orb]
                f.write(f"{b}   {sym}   {p[0]:.8f}   {p[1]:.8f}   {p[2]:.8f}   {l}   {m}   {spin}\n")


def write_rmatrix(path, Rs, Xs, source):
    """seedname_r.dat layout (as siesta2wtb.py --rmatrix): integer R, then m n Re(x) Im(x)
    Re(y) Im(y) Re(z) Im(z), m (home-cell orbital) fastest. Xs: one X per spin block
    (r is spin diagonal; the up-down blocks are zero)."""
    nspin, no = len(Xs), Xs[0].shape[1]
    zero = np.zeros(3)
    with open(path, "w") as f:
        f.write(f"yaehmop2wtb.py: <m,0|r|n,R> (Angstrom) of the Slater basis of {source}\n")
        f.write(f"{nspin * no}\n{len(Rs)}\n")
        for iR, (r1, r2, r3) in enumerate(Rs):
            for s in range(nspin):
                for n_ in range(no):
                    for s2 in range(nspin):
                        for m in range(no):
                            x = Xs[s][iR, m, n_] if s2 == s else zero
                            f.write("%5d%5d%5d%5d%5d" % (r1, r2, r3, s2 * no + m + 1, s * no + n_ + 1)
                                    + "".join("%16.10f%16.10f" % (v, 0.0) for v in x) + "\n")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", help="the .bind input of the YAeHMOP run (outputs next to it)")
    ap.add_argument("-o", "--outdir", default=".", help="output folder [.]")
    ap.add_argument("--name", default="huckel", help="prefix of the TB file [huckel]")
    ap.add_argument("--efermi", type=float, default=None,
                    help="energy zero in eV [Fermi level of the .out]")
    ap.add_argument("--cutoff", type=float, default=1e-10,
                    help="drop R vectors with all |H|, |S| below this [1e-10]")
    ap.add_argument("--rmatrix", action="store_true",
                    help="also write <name>-<TYPE>_r.dat, the position matrix of the basis")
    ap.add_argument("--params", nargs="+", default=None, metavar="PARMS",
                    help="--rmatrix: the YAeHMOP parameter file of the run (SP: the up and the down "
                         "file), for the exponents in full precision [the .out]")
    a = ap.parse_args()

    run, outdir = Path(a.run), Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    info = read_bind(run)
    dim, ov, kpts, cell, real = info["dim"], info["overlaps"], info["kpoints"], info["cell"], info["real"]
    N = mesh_size(kpts, dim, ov)

    stem, ext = (run.stem, run.suffix) if run.suffix else (run.name, "")
    if Path(f"{run}.HAM").exists():
        runs = {"": run}
    else:
        runs = {t: run.with_name(f"{stem}_{t}{ext}") for t in ("up", "down")}
        if not all(Path(f"{r}.HAM").exists() for r in runs.values()):
            raise SystemExit(f"no {run}.HAM (nor _up/_down): run YAeHMOP with 'Dump Hamil' and 'Dump Overlap'")
    outs = {t: read_out(f"{r}.out") for t, r in runs.items()}
    out0 = next(iter(outs.values()))
    kind = "SP" if "up" in runs else "SOC" if out0["soc"] else "NP"

    labels = []
    for ia, (sym, p) in enumerate(real):
        sh = out0["shells"].get(f"{sym.upper()}({ia + 1})", out0["shells"].get(sym.upper()))
        if sh is None:
            raise SystemExit(f"element {sym} not found in the parameters of the .out")
        labels += [(sym, p, orb) for s in "spdf" if s in sh for orb in ORBITALS[s]]
    no = len(labels)
    nb = no if kind == "NP" else 2 * no

    ef_common = None
    if kind == "SP" and a.efermi is None and Path(f"{run}.spin").exists():
        v = re.findall(r"^#Fermi_Energy:\s*(\S+)", Path(f"{run}.spin").read_text(), re.M)
        ef_common = float(v[-1]) if v else None
    print(f"{run}: {kind}, {no} orbitals{' x 2 spins' if kind != 'NP' else ''}, "
          f"mesh {N[0]}x{N[1]}x{N[2]}, overlap radius {ov}")
    blocks = {}
    for t, r in runs.items():
        Hk, Sk = read_dump(f"{r}.HAM"), read_dump(f"{r}.OV")
        if len(Hk) != len(kpts) or Hk.shape[1] != (2 * no if kind == "SOC" else no):
            raise SystemExit(f"{r}.HAM: {len(Hk)} k points of size {Hk.shape[1]}; expected "
                             f"{len(kpts)} of size {2 * no if kind == 'SOC' else no}")
        e0 = a.efermi if a.efermi is not None else ef_common if ef_common is not None else outs[t]["fermi"]
        if e0 is None:
            raise SystemExit(f"{r}.out has no #Fermi_Energy: add 'Average Properties' or give --efermi")
        cells, eH, eS = real_space(Hk, Sk, kpts, dim, ov, e0, a.cutoff)
        if kind != "SOC" and max(np.abs(h.imag).max() for _, h, _ in cells) > 1e-6:
            raise SystemExit(f"{r}: H(R) is not real")
        blocks[t] = cells
        print(f"  {('channel ' + t + ': ') if t else ''}E_F = {e0:.6f} eV (file holds H - E_F S); "
              f"{len(cells)} R vectors; H(k), S(k) rebuilt to {eH:.1e} eV, {eS:.1e}")

    if kind == "SP":
        up, dn = ({tuple(R): (h, s) for R, h, s in blocks[t]} for t in ("up", "down"))
        Z = np.zeros((no, no))
        cells = [(np.array(R), np.block([[up.get(R, (Z, Z))[0], Z], [Z, dn.get(R, (Z, Z))[0]]]),
                  np.block([[up.get(R, (Z, Z))[1], Z], [Z, dn.get(R, (Z, Z))[1]]]))
                 for R in sorted(set(up) | set(dn), key=lambda R: (sum(map(abs, R)), R))]
    else:
        cells = blocks[""]
    tb = outdir / f"{a.name}-{kind}.dat"
    write_tb(tb, kind, cell, nb, cells)
    write_basis(outdir / f"basis_set-{kind}", labels, (0,) if kind == "NP" else (1, -1))
    print(f"written: {tb}, {outdir / f'basis_set-{kind}'}")

    if a.rmatrix:
        rfile = outdir / f"{a.name}-{kind}_r.dat"
        shells, origin = channel_shells(kind, list(outs.values()), a.params)
        dS = rmatrix_file(rfile, real, cell, cells, no, kind, run.name, shells, origin)
        print(f"written: {rfile} (BSE_CENTER_FILE); exponents from {origin}; "
              f"max |S(quadrature) - S(YAeHMOP)| = {dS:.1e}")


if __name__ == "__main__":
    main()
