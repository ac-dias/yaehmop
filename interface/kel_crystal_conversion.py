#!/usr/bin/env python3
"""kel_crystal_conversion.py — crystal structures of DFT codes -> YAeHMOP (bind) input.

Reads the structure of a VASP, FHI-aims, SIESTA, DFTB+ or Quantum ESPRESSO calculation and writes a
periodic YAeHMOP input: "Geometry Crystallographic" with fractional coordinates, the lattice vectors
defined by helper atoms (copies of atom 1 shifted by a, b, c), the "Lattice" block (periodic dimension,
overlap images, vector definitions) and "Crystal Spec" (a, b, c, alpha, beta, gamma). Optional blocks:
THE CONST, Electrons/Charge, Band (path), K points (Gamma-centred mesh), Average Properties.

Single file, standard library only (Python >= 3.7).

Examples
  kel_crystal_conversion.py --format VASP --in POSCAR
  kel_crystal_conversion.py --format AIMS --in geometry.in --out MoS2.bind
  kel_crystal_conversion.py --format SIESTA --in MoS2.fdf --dimension 2 --const 1.75
  kel_crystal_conversion.py --format DFTB --in geo_end.gen
  kel_crystal_conversion.py --format QE --in scf.in --band "G 0 0 0, M 0.5 0 0, K 1/3 1/3 0, G 0 0 0" \\
      --band-points 40 --kmesh 11 11 1 --average-properties

Formats and files
  VASP    POSCAR / CONTCAR (VASP 5 with species line; VASP 4 needs --species), Direct or Cartesian,
          negative scale = volume, Selective dynamics accepted
  AIMS    geometry.in (lattice_vector, atom, atom_frac; Angstrom)
  SIESTA  .fdf (LatticeConstant, LatticeVectors or LatticeParameters, AtomicCoordinatesFormat,
          ChemicalSpeciesLabel, AtomicCoordinatesAndAtomicSpecies), .XV (Bohr) or .STRUCT_OUT
  DFTB    .gen (C/S/F types; C = cluster, refused) or dftb_in.hsd with Geometry = GenFormat {...}
          (inline block or <<< "file.gen")
  QE      pw.x input (&SYSTEM ibrav = 0 with CELL_PARAMETERS, or ibrav 1, 2, 3, 4, 6, 8 with celldm/A;
          ATOMIC_POSITIONS alat, bohr, angstrom or crystal)

Periodicity
  YAeHMOP periodic directions are the first lattice vectors (Lattice 2 = a, b). The periodic dimension
  is detected from the vacuum (an empty slab of at least --vacuum Angstrom across a lattice direction)
  unless --dimension is given; vacuum directions are moved to the end (cyclic permutation). Along the
  non-periodic directions the atoms are made contiguous: a DFT code is periodic in 3D, so an atom of a
  slab may sit at z = 0.999 instead of -0.001 and would be isolated from its layer in a 2D YAeHMOP run
  (disable with --no-contiguous).
"""
from __future__ import annotations
import argparse
import math
import re
import sys
from pathlib import Path

BOHR = 0.529177210903  # Angstrom

ELEMENTS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As "
            "Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb "
            "Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk "
            "Cf Es Fm Md No Lr").split()
SYMBOLS = {e.lower(): e for e in ELEMENTS[1:]}


# ----------------------------------------------------------------------------- small linear algebra
def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def cross(u, v):
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


def norm(u):
    return math.sqrt(dot(u, u))


def det3(m):
    return dot(m[0], cross(m[1], m[2]))


def inv3(m):
    """Inverse of a 3x3 matrix given as rows."""
    d = det3(m)
    if abs(d) < 1e-12:
        raise ValueError("singular lattice (zero cell volume)")
    c = [cross(m[1], m[2]), cross(m[2], m[0]), cross(m[0], m[1])]  # columns of the inverse * d
    return [[c[j][i] / d for j in range(3)] for i in range(3)]


def cart_to_frac(cell, xyz):
    """cell rows = a, b, c; r = f1 a + f2 b + f3 c -> f = r . inv(cell)."""
    inv = inv3(cell)
    return [sum(xyz[k] * inv[k][i] for k in range(3)) for i in range(3)]


def frac_to_cart(cell, f):
    return [sum(f[i] * cell[i][k] for i in range(3)) for k in range(3)]


# ----------------------------------------------------------------------------- helpers
class Structure:
    def __init__(self, cell, symbols, frac, source=""):
        if len(symbols) != len(frac) or not symbols:
            raise ValueError("no atoms read")
        self.cell = [list(map(float, v)) for v in cell]
        self.symbols = list(symbols)
        self.frac = [list(map(float, f)) for f in frac]
        self.source = source

    @classmethod
    def from_cartesian(cls, cell, symbols, xyz, source=""):
        return cls(cell, symbols, [cart_to_frac(cell, r) for r in xyz], source)


def element(label):
    """'Ti', 'Ti_sv', 'TI', 'S1', 'Mo.pbe' -> chemical symbol."""
    m = re.match(r"([A-Za-z]{1,2})", label.strip())
    if not m:
        raise ValueError(f"cannot read a chemical symbol from {label!r}")
    s = m.group(1)
    if s.lower() in SYMBOLS:
        return SYMBOLS[s.lower()]
    if s[0].lower() in SYMBOLS:
        return SYMBOLS[s[0].lower()]
    raise ValueError(f"unknown element {label!r}")


def num(text):
    """Float that accepts Fortran exponents (1.0d0) and fractions (1/3)."""
    t = text.strip().lower().replace("d", "e")
    if "/" in t:
        a, b = t.split("/")
        return float(a) / float(b)
    return float(t)


def strip_comment(line, chars="#!"):
    for c in chars:
        i = line.find(c)
        if i >= 0:
            line = line[:i]
    return line


# ----------------------------------------------------------------------------- VASP
def read_vasp(path, species=None):
    lines = [l.rstrip("\n") for l in Path(path).read_text().splitlines()]
    if len(lines) < 8:
        raise ValueError("POSCAR too short")
    scale = [num(x) for x in lines[1].split()[:3]]
    cell = [[num(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    if len(scale) == 1 and scale[0] < 0:  # negative = target volume
        f = (abs(scale[0]) / abs(det3(cell))) ** (1.0 / 3.0)
        cell = [[f * x for x in v] for v in cell]
    elif len(scale) == 3:
        cell = [[scale[k] * v[k] for k in range(3)] for v in cell]
    else:
        cell = [[scale[0] * x for x in v] for v in cell]
    i = 5
    tokens = lines[i].split()
    if all(re.match(r"^[A-Za-z]", t) for t in tokens):
        names = [element(t.split("/")[0]) for t in tokens]
        i += 1
    else:
        if not species:
            raise ValueError("VASP 4 POSCAR without a species line: give --species, e.g. --species Ti S")
        names = [element(s) for s in species]
    counts = [int(x) for x in lines[i].split()]
    if len(counts) != len(names):
        raise ValueError(f"{len(names)} species but {len(counts)} counts")
    i += 1
    if lines[i].strip()[:1].lower() == "s":  # Selective dynamics
        i += 1
    mode = lines[i].strip()[:1].lower()
    cart = mode in ("c", "k")
    i += 1
    symbols, pos = [], []
    for name, n in zip(names, counts):
        for _ in range(n):
            pos.append([num(x) for x in lines[i].split()[:3]])
            symbols.append(name)
            i += 1
    if cart:
        sc = scale[0] if len(scale) == 1 and scale[0] > 0 else 1.0
        if len(scale) == 1 and scale[0] < 0:
            sc = (abs(scale[0]) / abs(det3([[num(x) for x in lines[k].split()[:3]] for k in (2, 3, 4)]))) ** (1 / 3)
        return Structure.from_cartesian(cell, symbols, [[sc * x for x in r] for r in pos], "VASP")
    return Structure(cell, symbols, pos, "VASP")


# ----------------------------------------------------------------------------- FHI-aims
def read_aims(path):
    cell, symbols, frac, cart = [], [], [], []
    for raw in Path(path).read_text().splitlines():
        t = strip_comment(raw, "#").split()
        if not t:
            continue
        key = t[0].lower()
        if key == "lattice_vector":
            cell.append([num(x) for x in t[1:4]])
        elif key == "atom":
            cart.append(([num(x) for x in t[1:4]], element(t[4])))
        elif key == "atom_frac":
            frac.append(([num(x) for x in t[1:4]], element(t[4])))
    if len(cell) != 3:
        raise ValueError(f"geometry.in has {len(cell)} lattice_vector lines; a crystal needs 3")
    pos, symbols = [], []
    for r, s in cart:
        pos.append(cart_to_frac(cell, r)); symbols.append(s)
    for f, s in frac:
        pos.append(f); symbols.append(s)
    return Structure(cell, symbols, pos, "FHI-aims")


# ----------------------------------------------------------------------------- SIESTA
def _fdf_lines(path, seen=None):
    """fdf lines with %include files expanded and comments removed."""
    path = Path(path)
    seen = seen or set()
    if path.resolve() in seen:
        return []
    seen.add(path.resolve())
    out = []
    for raw in path.read_text().splitlines():
        line = strip_comment(raw, "#!;").strip()
        if not line:
            continue
        m = re.match(r"%include\s+(\S+)", line, re.I)
        if m:
            out += _fdf_lines(path.parent / m.group(1), seen)
        else:
            out.append(line)
    return out


def _fdf_key(k):
    return re.sub(r"[._-]", "", k).lower()


def _fdf_value(lines, key):
    for l in lines:
        t = l.split(None, 1)
        if t and _fdf_key(t[0]) == _fdf_key(key):
            return t[1].strip() if len(t) > 1 else ""
    return None


def _fdf_block(lines, name):
    for i, l in enumerate(lines):
        m = re.match(r"%block\s+(\S+)", l, re.I)
        if m and _fdf_key(m.group(1)) == _fdf_key(name):
            out = []
            for l2 in lines[i + 1:]:
                if re.match(r"%endblock", l2, re.I):
                    return out
                out.append(l2.split())
    return None


def _length(value, default_unit="ang"):
    t = value.split()
    v = num(t[0])
    unit = (t[1] if len(t) > 1 else default_unit).lower()
    if unit.startswith("bohr") or unit == "au":
        return v * BOHR
    if unit.startswith("nm"):
        return v * 10.0
    return v


def cell_from_parameters(a, b, c, alpha, beta, gamma):
    """Cell (rows) from lengths and angles (degrees): a along x, b in the xy plane."""
    ca, cb, cg = (math.cos(math.radians(x)) for x in (alpha, beta, gamma))
    sg = math.sin(math.radians(gamma))
    cx = c * cb
    cy = c * (ca - cb * cg) / sg
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0.0))
    return [[a, 0.0, 0.0], [b * cg, b * sg, 0.0], [cx, cy, cz]]


def read_siesta(path):
    p = Path(path)
    suf = p.suffix.lower()
    if suf == ".xv":
        rows = [l.split() for l in p.read_text().splitlines() if l.strip()]
        cell = [[num(x) * BOHR for x in rows[i][:3]] for i in range(3)]
        n = int(rows[3][0])
        symbols = [ELEMENTS[int(rows[4 + i][1])] for i in range(n)]
        xyz = [[num(x) * BOHR for x in rows[4 + i][2:5]] for i in range(n)]
        return Structure.from_cartesian(cell, symbols, xyz, "SIESTA .XV")
    if suf == ".struct_out" or p.name.upper().endswith("STRUCT_OUT"):
        rows = [l.split() for l in p.read_text().splitlines() if l.strip()]
        cell = [[num(x) for x in rows[i][:3]] for i in range(3)]
        n = int(rows[3][0])
        symbols = [ELEMENTS[int(rows[4 + i][1])] for i in range(n)]
        frac = [[num(x) for x in rows[4 + i][2:5]] for i in range(n)]
        return Structure(cell, symbols, frac, "SIESTA .STRUCT_OUT")
    lines = _fdf_lines(p)
    lc = _fdf_value(lines, "LatticeConstant")
    alat = _length(lc) if lc else 1.0
    lv = _fdf_block(lines, "LatticeVectors")
    lp = _fdf_block(lines, "LatticeParameters")
    if lv:
        cell = [[alat * num(x) for x in r[:3]] for r in lv[:3]]
    elif lp:
        a, b, c, al, be, ga = (num(x) for x in lp[0][:6])
        cell = cell_from_parameters(a * alat, b * alat, c * alat, al, be, ga)
    else:
        raise ValueError("fdf without %block LatticeVectors or LatticeParameters: not a crystal")
    species = {}
    for r in _fdf_block(lines, "ChemicalSpeciesLabel") or []:
        species[int(r[0])] = ELEMENTS[int(r[1])] if int(r[1]) > 0 else element(r[2])
    coords = _fdf_block(lines, "AtomicCoordinatesAndAtomicSpecies")
    if not coords:
        raise ValueError("fdf without %block AtomicCoordinatesAndAtomicSpecies")
    fmt = _fdf_key(_fdf_value(lines, "AtomicCoordinatesFormat") or "NotScaledCartesianBohr")
    symbols = [species[int(r[3])] for r in coords]
    xyz = [[num(x) for x in r[:3]] for r in coords]
    if fmt in ("fractional", "scaledbylatticevectors"):
        return Structure(cell, symbols, xyz, "SIESTA .fdf")
    if fmt in ("ang", "notscaledcartesianang"):
        pass
    elif fmt in ("bohr", "notscaledcartesianbohr"):
        xyz = [[x * BOHR for x in r] for r in xyz]
    elif fmt == "scaledcartesian":
        xyz = [[x * alat for x in r] for r in xyz]
    else:
        raise ValueError(f"AtomicCoordinatesFormat {fmt!r} not supported")
    return Structure.from_cartesian(cell, symbols, xyz, "SIESTA .fdf")


# ----------------------------------------------------------------------------- DFTB+
def _parse_gen(text, source):
    rows = [strip_comment(l, "#").split() for l in text.splitlines()]
    rows = [r for r in rows if r]
    n, kind = int(rows[0][0]), rows[0][1].upper()
    if kind == "C":
        raise ValueError("gen type C is a cluster (no lattice): YAeHMOP periodic input needs S or F")
    names = [element(x) for x in rows[1]]
    symbols, pos = [], []
    for r in rows[2:2 + n]:
        symbols.append(names[int(r[1]) - 1])
        pos.append([num(x) for x in r[2:5]])
    origin = [num(x) for x in rows[2 + n][:3]]
    cell = [[num(x) for x in rows[3 + n + i][:3]] for i in range(3)]
    if kind == "F":
        return Structure(cell, symbols, pos, source)
    return Structure.from_cartesian(cell, symbols, [[r[k] - origin[k] for k in range(3)] for r in pos], source)


def read_dftb(path):
    p = Path(path)
    text = p.read_text()
    if p.suffix.lower() == ".gen":
        return _parse_gen(text, "DFTB+ .gen")
    m = re.search(r"Geometry\s*=\s*GenFormat\s*\{(.*?)\}", text, re.S | re.I)
    if not m:
        raise ValueError("no 'Geometry = GenFormat { ... }' in the HSD file (give the .gen file instead)")
    body = m.group(1)
    inc = re.search(r'<<<\s*"?([^"\s]+)"?', body)
    if inc:
        return _parse_gen((p.parent / inc.group(1)).read_text(), "DFTB+ .gen")
    return _parse_gen(body, "DFTB+ .hsd")


# ----------------------------------------------------------------------------- Quantum ESPRESSO
def _qe_namelist(text):
    """Values of all namelists (&CONTROL ... /), one per line or several on a line."""
    vals = {}
    body = "\n".join(strip_comment(l, "!") for l in text.splitlines())
    # a namelist ends at a '/' that closes a line (not one inside a quoted path such as './pseudo/')
    for m in re.finditer(r"&(\w+)(.*?)(?:^|\s|,)/[ \t]*$", body, re.S | re.M):
        for item in re.findall(r"([\w()]+)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)", m.group(2)):
            vals[item[0].strip().lower()] = item[1].strip().strip("'\"")
    return vals


def _qe_card(text, name):
    m = re.search(rf"^\s*{name}\b\s*[\{{(]?\s*([A-Za-z_]*)\s*[\}})]?\s*$", text, re.M | re.I)
    if not m:
        return None, []
    rows = []
    for l in text[m.end():].splitlines():
        t = strip_comment(l, "!#").split()
        if not t:
            if rows:
                break
            continue
        if re.match(r"^[A-Z_]{4,}$", t[0]):
            break
        rows.append(t)
    return m.group(1).lower(), rows


def _qe_cell(nl, card_unit, card_rows):
    def get(k):
        return num(nl[k]) if k in nl else None
    ibrav = int(get("ibrav") or 0)
    alat = None
    if get("celldm(1)"):
        alat = get("celldm(1)") * BOHR
    elif get("a"):
        alat = get("a")
    if ibrav == 0:
        if len(card_rows) < 3:
            raise ValueError("ibrav = 0 needs CELL_PARAMETERS")
        v = [[num(x) for x in r[:3]] for r in card_rows[:3]]
        unit = card_unit or ("alat" if alat else "bohr")
        f = {"angstrom": 1.0, "bohr": BOHR}.get(unit, alat if unit == "alat" else None)
        if f is None:
            raise ValueError("CELL_PARAMETERS alat without celldm(1) or A")
        return [[f * x for x in r] for r in v], alat
    if alat is None:
        raise ValueError(f"ibrav = {ibrav} needs celldm(1) or A")
    a = alat
    ba = get("celldm(2)") or (get("b") / a if get("b") else None)
    ca = get("celldm(3)") or (get("c") / a if get("c") else None)
    if ibrav == 1:
        cell = [[a, 0, 0], [0, a, 0], [0, 0, a]]
    elif ibrav == 2:
        cell = [[-a / 2, 0, a / 2], [0, a / 2, a / 2], [-a / 2, a / 2, 0]]
    elif ibrav == 3:
        cell = [[a / 2, a / 2, a / 2], [-a / 2, a / 2, a / 2], [-a / 2, -a / 2, a / 2]]
    elif ibrav == 4:
        cell = [[a, 0, 0], [-a / 2, a * math.sqrt(3) / 2, 0], [0, 0, a * ca]]
    elif ibrav == 6:
        cell = [[a, 0, 0], [0, a, 0], [0, 0, a * ca]]
    elif ibrav == 8:
        cell = [[a, 0, 0], [0, a * ba, 0], [0, 0, a * ca]]
    else:
        raise ValueError(f"ibrav = {ibrav} not supported: use ibrav = 0 with CELL_PARAMETERS")
    return cell, alat


def read_qe(path):
    text = Path(path).read_text()
    nl = _qe_namelist(text)
    cu, crow = _qe_card(text, "CELL_PARAMETERS")
    cell, alat = _qe_cell(nl, cu, crow)
    pu, prow = _qe_card(text, "ATOMIC_POSITIONS")
    nat = int(num(nl["nat"])) if "nat" in nl else len(prow)
    prow = prow[:nat]
    if not prow:
        raise ValueError("no ATOMIC_POSITIONS")
    symbols = [element(r[0]) for r in prow]
    pos = [[num(x) for x in r[1:4]] for r in prow]
    unit = pu or "alat"
    if unit == "crystal":
        return Structure(cell, symbols, pos, "Quantum ESPRESSO")
    f = {"angstrom": 1.0, "bohr": BOHR, "alat": alat}.get(unit)
    if f is None:
        raise ValueError(f"ATOMIC_POSITIONS {unit} not supported")
    return Structure.from_cartesian(cell, symbols, [[f * x for x in r] for r in pos], "Quantum ESPRESSO")


READERS = {"vasp": read_vasp, "aims": read_aims, "siesta": read_siesta, "dftb": read_dftb, "qe": read_qe}
ALIASES = {"poscar": "vasp", "contcar": "vasp", "fhi-aims": "aims", "fhiaims": "aims", "dftb+": "dftb",
           "gen": "dftb", "espresso": "qe", "quantum-espresso": "qe", "pw": "qe", "pwscf": "qe"}


# ----------------------------------------------------------------------------- periodicity
def heights(cell):
    """Distance between opposite faces of the cell along each lattice direction."""
    v = abs(det3(cell))
    return [v / norm(cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]


def largest_gap(values):
    """Largest empty interval of fractional coordinates (mod 1): (size, start of the occupied block)."""
    z = sorted(x % 1.0 for x in values)
    gaps = [(z[i + 1] - z[i], z[i + 1]) for i in range(len(z) - 1)] + [(z[0] + 1.0 - z[-1], z[0])]
    return max(gaps)


def vacuum_axes(st, vacuum):
    h = heights(st.cell)
    return [ax for ax in range(3) if largest_gap([f[ax] for f in st.frac])[0] * h[ax] >= vacuum]


def permute_axes(st, order):
    st.cell = [st.cell[i] for i in order]
    st.frac = [[f[i] for i in order] for f in st.frac]


def make_contiguous(st, axes):
    """Atoms of a slab/wire wrapped across the vacuum are brought back next to their block."""
    moved = 0
    for ax in axes:
        _, start = largest_gap([f[ax] for f in st.frac])
        mean_old = sum(f[ax] for f in st.frac) / len(st.frac)
        new = [((f[ax] - start) % 1.0) + start for f in st.frac]
        shift = round(mean_old - sum(new) / len(new))
        for f, z in zip(st.frac, new):
            z += shift
            if abs(z - f[ax]) > 1e-9:
                moved += 1
            f[ax] = z
    return moved


def lengths_angles(cell):
    a, b, c = (norm(v) for v in cell)
    ang = lambda u, v: math.degrees(math.acos(max(-1.0, min(1.0, dot(u, v) / (norm(u) * norm(v))))))
    return (a, b, c), (ang(cell[1], cell[2]), ang(cell[0], cell[2]), ang(cell[0], cell[1]))


# ----------------------------------------------------------------------------- writer
def parse_band(text):
    """'G 0 0 0, M 0.5 0 0, K 1/3 1/3 0' -> [(label, [k1, k2, k3]), ...]."""
    out = []
    for item in text.split(","):
        t = item.split()
        if not t:
            continue
        if len(t) != 4:
            raise ValueError(f"band point {item.strip()!r}: use LABEL k1 k2 k3")
        out.append((re.sub(r"[^A-Za-z0-9]", "", t[0])[:3] or "K", [num(x) for x in t[1:]]))
    if len(out) < 2:
        raise ValueError("a band path needs at least two points")
    return out


def write_bind(st, dim, overlaps, title, const=None, electrons=None, charge=None, band=None,
               band_points=40, kmesh=None, average=False, notes=()):
    n = len(st.symbols)
    (a, b, c), (al, be, ga) = lengths_angles(st.cell)
    L = [f"; {line}" for line in notes]
    L += [title, "", "Geometry Crystallographic", str(n + dim)]
    for i, (s, f) in enumerate(zip(st.symbols, st.frac), 1):
        L.append(f"{i:4d} {s:<2s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    for j in range(dim):  # helper atoms: atom 1 + lattice vector j
        f = list(st.frac[0]); f[j] += 1.0
        L.append(f"{n + j + 1:4d} {st.symbols[0]:<2s} {f[0]: .12f} {f[1]: .12f} {f[2]: .12f}")
    L += ["", "Lattice", str(dim), " ".join(str(overlaps if j < dim else 0) for j in range(3))]
    L += [f"1 {n + j + 1}" for j in range(dim)]
    L += ["", "Crystal Spec", f"{a:.12f} {b:.12f} {c:.12f}", f"{al:.12f} {be:.12f} {ga:.12f}"]
    if const is not None:
        L += ["", f"THE CONST {const:.12f}"]
    if electrons is not None:
        L += ["", "Electrons", f"{electrons:g}"]
    if charge is not None:
        L += ["", "Charge", f"{charge:g}"]
    if band:
        L += ["", "Band", str(band_points), str(len(band))]
        L += [f"{lab:<4s} {k[0]: .12f} {k[1]: .12f} {k[2]: .12f}" for lab, k in band]
    if average:
        L += ["", "Average Properties"]
    if kmesh:
        mesh = [m if j < dim else 1 for j, m in enumerate(kmesh)]
        pts = []
        for i in range(mesh[0]):
            for j in range(mesh[1]):
                for k in range(mesh[2]):
                    pts.append((i / mesh[0], j / mesh[1], k / mesh[2]))
        L += ["", "K points", str(len(pts))]
        L += [f"{p[0]:.15f} {p[1]:.15f} {p[2]:.15f} 1.0" for p in pts]
    return "\n".join(L) + "\n"


# ----------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="kel_crystal_conversion.py",
        description="Convert a crystal structure of VASP, FHI-aims, SIESTA, DFTB+ or Quantum ESPRESSO "
                    "into a periodic YAeHMOP (bind) input.",
        epilog="Examples" + __doc__.split("Examples", 1)[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--format", "-f", required=True, type=str.lower,
                    help="input format: VASP, AIMS, SIESTA, DFTB, QE (case-insensitive)")
    ap.add_argument("--in", "-i", dest="infile", required=True, type=Path, help="input structure file")
    ap.add_argument("--out", "-o", type=Path, default=None,
                    help="output YAeHMOP input (default: <input name>.bind; '-' writes to the screen)")
    ap.add_argument("--dimension", "-d", type=int, choices=(1, 2, 3), default=None,
                    help="periodic dimension (default: detected from the vacuum)")
    ap.add_argument("--vacuum", type=float, default=6.0,
                    help="empty thickness (Angstrom) that marks a non-periodic direction (default 6)")
    ap.add_argument("--overlaps", type=int, default=5,
                    help="overlap images along each periodic direction, the 'Lattice' block (default 5)")
    ap.add_argument("--no-contiguous", action="store_true",
                    help="keep atoms wrapped across the vacuum as they are")
    ap.add_argument("--species", nargs="+", help="species names for a VASP 4 POSCAR without species line")
    ap.add_argument("--title", default=None, help="title line of the YAeHMOP input")
    ap.add_argument("--const", type=float, default=None, help="Wolfsberg-Helmholz K (THE CONST)")
    ap.add_argument("--electrons", type=float, default=None, help="number of valence electrons (Electrons)")
    ap.add_argument("--charge", type=float, default=None, help="total charge (Charge)")
    ap.add_argument("--band", default=None,
                    help='band path in fractional reciprocal coordinates, e.g. "G 0 0 0, M 0.5 0 0, K 1/3 1/3 0, G 0 0 0"')
    ap.add_argument("--band-points", type=int, default=40, help="points per band line (default 40)")
    ap.add_argument("--kmesh", type=int, nargs=3, metavar=("N1", "N2", "N3"),
                    help="Gamma-centred k mesh written as 'K points' (non-periodic directions use 1)")
    ap.add_argument("--average-properties", action="store_true", help="add 'Average Properties'")
    args = ap.parse_args(argv)

    fmt = ALIASES.get(args.format, args.format)
    if fmt not in READERS:
        ap.error(f"unknown format {args.format!r}; use VASP, AIMS, SIESTA, DFTB or QE")
    if not args.infile.is_file():
        ap.error(f"input file not found: {args.infile}")
    if args.overlaps < 1:
        ap.error("--overlaps must be >= 1")
    try:
        st = read_vasp(args.infile, args.species) if fmt == "vasp" else READERS[fmt](args.infile)
        notes = [f"YAeHMOP input written by kel_crystal_conversion.py from {args.infile.name} ({st.source})"]
        vac = vacuum_axes(st, args.vacuum)
        if args.dimension is None:
            dim = 3 - len(vac)
            if dim < 1:
                raise ValueError("vacuum along all three directions: a molecule, not a crystal")
            how = "detected from the vacuum"
        else:
            dim = args.dimension
            how = "given with --dimension"
            if len(vac) != 3 - dim:
                notes.append(f"warning: {len(vac)} vacuum direction(s) found, but dimension {dim} was given")
        nonper = [ax for ax in range(3) if ax in vac][:3 - dim] if len(vac) >= 3 - dim else []
        if dim < 3 and not nonper:
            nonper = list(range(dim, 3))
        order = [ax for ax in range(3) if ax not in nonper] + nonper
        if order != [0, 1, 2]:
            for k in range(3):  # prefer a cyclic permutation (keeps the handedness)
                cyc = [(k + i) % 3 for i in range(3)]
                if set(cyc[dim:]) == set(nonper):
                    order = cyc
                    break
            permute_axes(st, order)
            notes.append(f"lattice vectors reordered {order} so that the periodic ones come first")
        notes.append(f"periodic dimension {dim} ({how}); overlap images {args.overlaps}")
        if dim < 3 and not args.no_contiguous:
            moved = make_contiguous(st, list(range(dim, 3)))
            if moved:
                notes.append(f"{moved} atom coordinate(s) moved across the vacuum to keep the slab/wire contiguous")
        notes.append(f"the last {dim} atom(s) only define the lattice vectors (1 -> {len(st.symbols) + 1}"
                     + (f", 1 -> {len(st.symbols) + 2}" if dim > 1 else "")
                     + (f", 1 -> {len(st.symbols) + 3}" if dim > 2 else "") + ")")
        band = parse_band(args.band) if args.band else None
        formula = {}
        for s in st.symbols:
            formula[s] = formula.get(s, 0) + 1
        title = args.title or "Generic EHT: " + "".join(f"{s}{n if n > 1 else ''}" for s, n in formula.items())
        text = write_bind(st, dim, args.overlaps, title, args.const, args.electrons, args.charge, band,
                          args.band_points, args.kmesh, args.average_properties, notes)
    except (ValueError, KeyError, IndexError) as exc:
        ap.error(f"{args.infile}: {exc}")
    if args.out is not None and str(args.out) == "-":
        sys.stdout.write(text)
        return
    out = args.out or args.infile.with_name(args.infile.name + ".bind")
    out.write_text(text)
    print(f"written {out}: {len(st.symbols)} atoms, dimension {dim}")


if __name__ == "__main__":
    main()
