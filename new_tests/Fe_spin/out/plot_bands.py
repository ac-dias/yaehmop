#!/usr/bin/env python3
"""Band structure figures from YAeHMOP (yaehmop-dev) outputs.

Run it inside an out/ folder:

    python3 plot_bands.py [--emin -6] [--emax 4] [--dir .]

The kind of run is found from the files in the folder:
  * spin polarized: <name>_up.<ext>.band and <name>_down.<ext>.band
    -> bands_spin.png (spin up, spin down, both)
  * SOC: <name>.band with a <name>.SOC file next to it; if nosoc/<name>.band
    exists, the bands without SOC are drawn too -> bands_soc.png

Energies are relative to the Fermi level of each run (#Fermi_Energy in the
.out; for an insulator YAeHMOP puts it at the top of the valence band).
The k-path distances use the reciprocal lattice from the 'Crystal Spec' of the
input file (the copy of the input in the folder).
"""
import argparse, re
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_band(path):
    lines = Path(path).read_text().splitlines()
    i = next(n for n, l in enumerate(lines) if "Special points were done with" in l)
    nsp = int(lines[i].split()[0])
    nst = int(next(l for l in lines if "orbitals in the unit cell" in l).split()[0])
    j = next(n for n, l in enumerate(lines) if "orbitals in the unit cell" in l)
    specials = []
    for l in lines[j + 1:j + 1 + nsp]:
        p = l.split()
        specials.append((p[0], np.array([float(x) for x in p[1:4]])))
    kp, ev, cur = [], [], []
    for l in lines[next(n for n, l in enumerate(lines) if "Begin band data" in l) + 1:]:
        if l.startswith("#END"):
            break
        if l.startswith("; K point"):
            if cur:
                ev.append(cur)
            cur = []
            kp.append([float(x) for x in l.split(":")[1].split()])
        elif l.strip():
            cur.append(float(l))
    ev.append(cur)
    ev = np.array(ev)
    assert ev.shape[1] == nst, "number of bands does not match the header"
    return np.array(kp), ev, specials


def fermi(out_file):
    vals = re.findall(r"#Fermi_Energy:\s*(\S+)", Path(out_file).read_text())
    return float(vals[-1]) if vals else None


def reciprocal(input_file):
    """Reciprocal vectors (rows, 1/A, with 2 pi) from the Crystal Spec section."""
    lines = [l for l in Path(input_file).read_text().splitlines() if not l.startswith(";")]
    i = next(n for n, l in enumerate(lines) if l.strip().lower() == "crystal spec")
    a, b, c = (float(x) for x in lines[i + 1].split()[:3])
    al, be, ga = (np.radians(float(x)) for x in lines[i + 2].split()[:3])
    A = np.zeros((3, 3))
    A[0] = [a, 0, 0]
    A[1] = [b * np.cos(ga), b * np.sin(ga), 0]
    cx = c * np.cos(be)
    cy = c * (np.cos(al) - np.cos(be) * np.cos(ga)) / np.sin(ga)
    A[2] = [cx, cy, np.sqrt(c * c - cx * cx - cy * cy)]
    return 2 * np.pi * np.linalg.inv(A).T


def path_axis(kp, specials, B):
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(kp, axis=0) @ B, axis=1))]
    ticks, labels = [], []
    for name, k in specials:
        idx = np.where(np.linalg.norm(kp - k, axis=1) < 1e-5)[0]
        for n in idx:
            if not any(abs(d[n] - t) < 1e-8 for t in ticks):
                ticks.append(d[n]); labels.append("Γ" if name.upper() in ("G", "GAMMA") else name)
    order = np.argsort(ticks)
    return d, [ticks[o] for o in order], [labels[o] for o in order]


def decorate(ax, d, ticks, labels, emin, emax, title):
    for t in ticks:
        ax.axvline(t, color="0.75", lw=0.6, zorder=0)
    ax.axhline(0, color="0.4", lw=0.7, ls=":", zorder=0)
    ax.set_xticks(ticks); ax.set_xticklabels(labels)
    ax.set_xlim(d[0], d[-1]); ax.set_ylim(emin, emax)
    ax.set_title(title, fontsize=10)


def plot_spin(folder, up, dn, emin, emax):
    name = up.name.replace("_up", "").replace(".band", "")
    inp = folder / name if (folder / name).exists() else next(folder.glob("*.bind"))
    B = reciprocal(inp)
    data = {}
    for tag, f in (("up", up), ("down", dn)):
        kp, ev, sp = read_band(f)
        ef = fermi(f.with_suffix(".out"))
        data[tag] = (kp, ev - ef, sp, ef)
    d, ticks, labels = path_axis(data["up"][0], data["up"][2], B)
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.8), sharey=True)
    col = {"up": "#1f5fbf", "down": "#c8372d"}
    for ax, tag in zip(axs[:2], ("up", "down")):
        ev = data[tag][1]
        ax.plot(d, ev, color=col[tag], lw=1.2)
        decorate(ax, d, ticks, labels, emin, emax,
                 f"spin {tag}  (E_F = {data[tag][3]:.3f} eV)")
    for tag, ls in (("up", "-"), ("down", "--")):
        ev = data[tag][1]
        for n in range(ev.shape[1]):
            axs[2].plot(d, ev[:, n], color=col[tag], lw=1.1, ls=ls,
                        label=f"spin {tag}" if n == 0 else None)
    decorate(axs[2], d, ticks, labels, emin, emax, "both channels")
    axs[2].legend(loc="lower right", fontsize=8, frameon=False)
    axs[0].set_ylabel("E − E_F of each channel (eV)")
    fig.suptitle(f"{name.split('.')[0]}: spin-polarized EHT (yaehmop-dev)", fontsize=11)
    fig.tight_layout()
    out = folder / "bands_spin.png"
    fig.savefig(out, dpi=200)
    print("wrote", out)


def plot_soc(folder, band, emin, emax):
    name = band.name[:-len(".band")]
    B = reciprocal(folder / name)
    kp, ev, sp = read_band(band)
    ef = fermi(folder / (name + ".out"))
    d, ticks, labels = path_axis(kp, sp, B)
    fig, ax = plt.subplots(figsize=(6, 5))
    nosoc = folder / "nosoc" / band.name
    if nosoc.exists():
        kp0, ev0, _ = read_band(nosoc)
        ef0 = fermi(nosoc.with_suffix(".out"))
        for n in range(ev0.shape[1]):
            ax.plot(d, ev0[:, n] - ef0, color="0.55", lw=1.0, ls="--",
                    label=f"without SOC (E_F = {ef0:.3f} eV)" if n == 0 else None)
    for n in range(ev.shape[1]):
        ax.plot(d, ev[:, n] - ef, color="#c8372d", lw=1.2,
                label=f"with SOC (E_F = {ef:.3f} eV)" if n == 0 else None)
    decorate(ax, d, ticks, labels, emin, emax, f"{name.split('.')[0]}: EHT + SOC (yaehmop-dev)")
    ax.set_ylabel("E − E_F (eV)")
    ax.legend(loc="lower right", fontsize=8, frameon=True)
    fig.tight_layout()
    out = folder / "bands_soc.png"
    fig.savefig(out, dpi=200)
    print("wrote", out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default=Path(__file__).resolve().parent, type=Path)
    ap.add_argument("--emin", type=float, default=-6.0)
    ap.add_argument("--emax", type=float, default=4.0)
    a = ap.parse_args()
    folder = a.dir.resolve()
    ups = sorted(folder.glob("*_up.*band"))
    if ups:
        up = ups[0]
        dn = folder / up.name.replace("_up.", "_down.")
        plot_spin(folder, up, dn, a.emin, a.emax)
        return
    for band in sorted(folder.glob("*.band")):
        if (folder / (band.name[:-len(".band")] + ".SOC")).exists():
            plot_soc(folder, band, a.emin, a.emax)
            return
    raise SystemExit("no spin-polarized or SOC band file found in " + str(folder))


if __name__ == "__main__":
    main()
