"""
Skript zur visualisierung der Wirbelstärke mit Skala
"""


import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from src.simulation import VortexSimulation

warnings.filterwarnings("ignore")

FIG = "bilder"
CKPT = "ckpt"
os.makedirs(FIG, exist_ok=True)
os.makedirs(CKPT, exist_ok=True)

dt_gitter = {160: 0.033, 256: 0.012}

plt.rcParams["font.size"] = 9
plt.rcParams["figure.dpi"] = 130

T = 12.2
t_periode = [180.0 - T, 180.0 - 2 * T / 3, 180.0 - T / 3]


def panel(ax, omega, x, y, label, vmax=1.5):
    mesh = ax.pcolormesh(x, y, omega, cmap="RdBu_r", shading="auto",
                         vmin=-vmax, vmax=vmax, rasterized=True)
    zyl = plt.Circle((0, 0), 1.0, facecolor="0.8", edgecolor="k",
                     linewidth=0.7, zorder=5)
    ax.add_patch(zyl)
    ax.set_aspect("equal")
    ax.set_xlim(-4, 22)
    ax.set_ylim(-5, 5)
    ax.set_ylabel(r"$y/a$")
    ax.text(0.985, 0.90, label, transform=ax.transAxes,
            ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6",
                      alpha=0.9))
    return mesh


def save(fig, name):
    fig.savefig(f"{FIG}/{name}.pdf", bbox_inches="tight")
    fig.savefig(f"{FIG}/{name}.png", bbox_inches="tight")
    plt.close(fig)
    print("  ->", name, flush=True)


def gitter(g):
    sim = VortexSimulation(100, g, g, dt=dt_gitter[g])
    r = np.exp(sim.grid_xi)
    return r * np.cos(sim.grid_theta), r * np.sin(sim.grid_theta)


def feld(Re, g, t):
    return np.load(f"{CKPT}/field_Re{Re}_{g}_t{t:.1f}.npy")


def rechnen(Re, g, ziele, budget):
    dt = dt_gitter[g]
    ck = f"{CKPT}/Re{Re}_{g}.npz"
    sim = VortexSimulation(Re, g, g, dt=dt)

    if os.path.exists(ck):
        z = np.load(ck)
        sim.omega = z["omega"]
        sim.psi = z["psi"]
        sim.u_xi = z["u_xi"]
        sim.u_theta = z["u_theta"]
        done = int(z["done"])
    else:
        done = 0

    n_ges = int(180 // dt)
    stop = min(n_ges, done + budget)

    for k in range(done + 1, stop + 1):
        sim.omega = sim._step()
        t = k * dt
        for tz in ziele:
            if abs(t - tz) < dt / 2:
                np.save(f"{CKPT}/field_Re{Re}_{g}_t{tz:.1f}.npy", sim.omega)

    np.savez(ck, omega=sim.omega, psi=sim.psi, u_xi=sim.u_xi,
             u_theta=sim.u_theta, done=stop)

    print(f"Re={Re}, {g}x{g}: Schritt {stop} von {n_ges}, t = {stop * dt:.1f}",
          flush=True)
    if not np.all(np.isfinite(sim.omega)):
        print("  Achtung, NaN oder Inf im Feld", flush=True)

    return stop >= n_ges


def abb_wirbelstrasse():
    x, y = gitter(256)

    fig, axes = plt.subplots(3, 1, figsize=(7.4, 5.6), sharex=True)
    for ax, t in zip(axes, t_periode):
        mesh = panel(ax, feld(100, 256, t), x, y, rf"$t = {t:.1f}$")
    axes[-1].set_xlabel(r"$x/a$")
    fig.suptitle(r"Wirbelstärke $\omega$, $Re = 100$, Gitter $256\times256$, "
                 rf"eine Ablöseperiode $T \approx {T:.1f}$", y=0.93)
    cb = fig.colorbar(mesh, ax=axes, pad=0.015, fraction=0.03)
    cb.set_label(r"$\omega$")
    save(fig, "wirbelstrasse_periode_Re100")

    # Re = 100 gegen Re = 150, jeweils am Ende
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 4.0), sharex=True)
    for ax, Re in zip(axes, (100, 150)):
        mesh = panel(ax, feld(Re, 256, 180.0), x, y,
                     rf"$Re = {Re}$,  $t = 180$")
    axes[-1].set_xlabel(r"$x/a$")
    fig.suptitle(r"Einfluss der Reynolds-Zahl, Gitter $256\times256$", y=0.96)
    cb = fig.colorbar(mesh, ax=axes, pad=0.015, fraction=0.03)
    cb.set_label(r"$\omega$")
    save(fig, "vergleich_Re100_Re150")


def abb_entstehung():
    g = 160
    dt = dt_gitter[g]
    zeiten = [20.0, 60.0, 120.0, 180.0]

    sim = VortexSimulation(100, g, g, dt=dt)
    x, y = gitter(g)

    fig, axes = plt.subplots(4, 1, figsize=(7.4, 7.4), sharex=True)
    done = 0
    for ax, t in zip(axes, zeiten):
        n = int(round(t / dt))
        for _ in range(n - done):
            sim.omega = sim._step()
        done = n
        mesh = panel(ax, sim.omega, x, y, rf"$t = {t:.0f}$")
    axes[-1].set_xlabel(r"$x/a$")
    fig.suptitle(r"Entstehung der Wirbelstraße, $Re = 100$, "
                 r"Gitter $160\times160$", y=0.94)
    cb = fig.colorbar(mesh, ax=axes, pad=0.015, fraction=0.03)
    cb.set_label(r"$\omega$")
    save(fig, "entstehung_Re100")


if __name__ == "__main__":
    modus = sys.argv[1] if len(sys.argv) > 1 else "alles"

    if modus == "alles":
        for Re in (100, 150):
            if Re == 100:
                ziele = t_periode + [180.0]
            else:
                ziele = [180.0]
            rechnen(Re, 256, ziele, 10**9)
            z = np.load(f"{CKPT}/Re{Re}_256.npz")
            np.save(f"{CKPT}/field_Re{Re}_256_t180.0.npy", z["omega"])
        abb_wirbelstrasse()
        abb_entstehung()
        print("fertig, alles in bilder/", flush=True)

    elif modus == "dev":
        abb_entstehung()

    elif modus == "plot":
        abb_wirbelstrasse()

    else:
        Re = int(sys.argv[2])
        ziele = t_periode + [180.0] if Re == 100 else [180.0]
        rechnen(Re, 256, ziele, int(modus))
