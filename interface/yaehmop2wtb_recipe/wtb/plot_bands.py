#!/usr/bin/env python3
"""Plot the band structure computed by wtb.x (out/bands.dat), colored by <s_z>.

    python3 plot_bands.py            # -> bands.png

out/bands.dat has one block per band ("#band  n"), with the columns
k-path, energy (eV, relative to the VBM), <sx>, <sy>, <sz>, orbital weight.
out/KLABELS.dat has the k-path coordinate of each high-symmetry point, and the
labels come from kpoints.dat. Needs numpy and matplotlib.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# bands: split the file at the "#band" lines
bands, cur = [], []
for line in open("out/bands.dat"):
    if line.lstrip().startswith("#"):
        if cur:
            bands.append(np.array(cur))
        cur = []
        continue
    if line.split():
        cur.append([float(x) for x in line.split()])
if cur:
    bands.append(np.array(cur))

# high-symmetry points: positions from out/KLABELS.dat, names from kpoints.dat
ticks = np.loadtxt("out/KLABELS.dat")
lines = [l.split() for l in open("kpoints.dat").read().splitlines()[2:] if l.split()]
names = [lines[0][3]] + [l[3] for l in lines[1::2]]
names = ["Γ" if n.upper() in ("G", "GAMMA") else n for n in names]

fig, ax = plt.subplots(figsize=(5, 6))
for b in bands:
    ax.plot(b[:, 0], b[:, 1], color="0.8", lw=0.8, zorder=1)
    sc = ax.scatter(b[:, 0], b[:, 1], c=b[:, 4], cmap="coolwarm", vmin=-1, vmax=1, s=4, zorder=2)
for t in ticks:
    ax.axvline(t, color="k", lw=0.5)
ax.axhline(0, color="k", ls=":", lw=0.6)
ax.set_xticks(ticks)
ax.set_xticklabels(names)
ax.set_xlim(ticks[0], ticks[-1])
ax.set_ylim(-3, 4)
ax.set_ylabel(r"$E - E_\mathrm{VBM}$ (eV)")
ax.set_title("MoS$_2$ with SOC: YAeHMOP → WanTiBEXOS", fontsize=10)
fig.colorbar(sc, ax=ax, label=r"$\langle s_z \rangle$")
fig.tight_layout()
fig.savefig("bands.png", dpi=150)
print("bands.png written")
