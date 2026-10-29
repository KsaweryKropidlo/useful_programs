"""Generate figures for the oblique projectile motion report."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


OUT = Path(__file__).parent / "files" / "rzut_ukosny"
OUT.mkdir(parents=True, exist_ok=True)

G = 9.81
V0 = 36.0
ALPHA_DEG = 12.0
ALPHA = np.deg2rad(ALPHA_DEG)
MASS = 0.145
RADIUS = 0.0365
AREA = np.pi * RADIUS**2
RHO = 1.225
CD = 0.47
DT = 0.001

BLUE, ORANGE, GREEN, PURPLE, INK = (
    "#2463A6", "#D47722", "#168578", "#8455A4", "#233044"
)
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.2,
    "lines.linewidth": 2.2,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})


def acceleration(vx, vy, drag_coefficient=CD):
    speed = np.hypot(vx, vy)
    drag_scale = RHO * drag_coefficient * AREA / (2 * MASS)
    return -drag_scale * speed * vx, -G - drag_scale * speed * vy


def derivative(state, drag_coefficient=CD):
    x, y, vx, vy = state
    ax, ay = acceleration(vx, vy, drag_coefficient)
    return np.array([vx, vy, ax, ay])


def rk4_step(state, dt, drag_coefficient=CD):
    k1 = derivative(state, drag_coefficient)
    k2 = derivative(state + 0.5 * dt * k1, drag_coefficient)
    k3 = derivative(state + 0.5 * dt * k2, drag_coefficient)
    k4 = derivative(state + dt * k3, drag_coefficient)
    return state + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6


def simulate(drag_coefficient=CD, dt=DT):
    state = np.array([
        0.0, 0.0, V0 * np.cos(ALPHA), V0 * np.sin(ALPHA)
    ])
    times = [0.0]
    states = [state.copy()]
    time = 0.0

    while time < 10.0:
        next_state = rk4_step(state, dt, drag_coefficient)
        next_time = time + dt
        if next_state[1] < 0:
            fraction = state[1] / (state[1] - next_state[1])
            landing_state = state + fraction * (next_state - state)
            times.append(time + fraction * dt)
            states.append(landing_state)
            break
        times.append(next_time)
        states.append(next_state.copy())
        state = next_state
        time = next_time
    else:
        raise RuntimeError("The projectile did not return to ground level.")

    return np.asarray(times), np.asarray(states)


def save(fig, filename):
    fig.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main():
    vacuum_time = 2 * V0 * np.sin(ALPHA) / G
    vacuum_range = V0**2 * np.sin(2 * ALPHA) / G
    vacuum_height = (V0 * np.sin(ALPHA))**2 / (2 * G)
    vacuum_peak_time = V0 * np.sin(ALPHA) / G
    vacuum_vx = V0 * np.cos(ALPHA)
    vacuum_vy = V0 * np.sin(ALPHA)

    t_drag, s_drag = simulate(CD)
    t_vac, s_vac = simulate(0.0)
    t_double, s_double = simulate(2 * CD)

    x, y, vx, vy = s_drag.T
    speed = np.hypot(vx, vy)
    drag_force = 0.5 * RHO * CD * AREA * speed**2
    kinetic = 0.5 * MASS * speed**2
    potential = MASS * G * y
    mechanical = kinetic + potential

    fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
    ax.plot(s_vac[:, 0], s_vac[:, 1], color=INK, ls="--",
            label="Bez oporu")
    ax.plot(x, y, color=BLUE, label=f"Opór kwadratowy, $C_D={CD:g}$")
    ax.plot(s_double[:, 0], s_double[:, 1], color=ORANGE,
            label=f"Opór kwadratowy, $C_D={2 * CD:g}$")
    ax.scatter([0], [0], color=GREEN, zorder=4, label="Punkt wyrzutu")
    ax.set(xlabel="Położenie poziome x [m]",
           ylabel="Wysokość y [m]",
           title="Tor ruchu: wpływ oporu powietrza")
    ax.set_ylim(bottom=0)
    ax.legend()
    save(fig, "tory_ruchu.png")

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                             layout="constrained")
    axes[0].plot(t_vac, s_vac[:, 0], color=INK, ls="--",
                 label="Bez oporu")
    axes[0].plot(t_drag, x, color=BLUE, label="Z oporem")
    axes[0].set_ylabel("x [m]")
    axes[0].set_title("Położenie w funkcji czasu")
    axes[0].legend()
    axes[1].plot(t_vac, s_vac[:, 1], color=INK, ls="--",
                 label="Bez oporu")
    axes[1].plot(t_drag, y, color=BLUE, label="Z oporem")
    axes[1].set(xlabel="Czas t [s]", ylabel="y [m]")
    axes[1].legend()
    save(fig, "polozenie_czas.png")

    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                             layout="constrained")
    axes[0].plot(t_drag, vx, color=BLUE, label="$v_x$")
    axes[0].plot(t_drag, vy, color=ORANGE, label="$v_y$")
    axes[0].axhline(0, color=INK, lw=0.8)
    axes[0].set_ylabel("Prędkość [m/s]")
    axes[0].set_title("Składowe prędkości przy oporze powietrza")
    axes[0].legend()
    axes[1].plot(t_drag, speed, color=PURPLE, label="Szybkość")
    axes[1].set(xlabel="Czas t [s]", ylabel="v [m/s]")
    axes[1].legend()
    save(fig, "predkosc_czas.png")

    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(t_vac, 0.5 * MASS * np.sum(s_vac[:, 2:]**2, axis=1),
            color=ORANGE, label="Kinetyczna bez oporu")
    ax.plot(t_vac, MASS * G * s_vac[:, 1],
            color=GREEN, label="Potencjalna bez oporu")
    ax.plot(t_drag, kinetic, color=ORANGE, ls="--",
            label="Kinetyczna z oporem")
    ax.plot(t_drag, potential, color=GREEN, ls="--",
            label="Potencjalna z oporem")
    ax.plot(t_drag, mechanical, color=BLUE, lw=2.5,
            label="Mechaniczna z oporem")
    ax.set(xlabel="Czas t [s]", ylabel="Energia [J]",
           title="Energia mechaniczna i praca oporu")
    ax.legend(ncol=2, fontsize=9)
    save(fig, "energia.png")

    fig, ax = plt.subplots(figsize=(9, 4.5), layout="constrained")
    ax.plot(t_drag, drag_force, color=PURPLE)
    ax.set(xlabel="Czas t [s]", ylabel="Wartość siły oporu $F_D$ [N]",
           title="Siła oporu zależy od kwadratu szybkości")
    save(fig, "sila_oporu.png")

    print(f"v0 = {V0:.3f} m/s; alpha = {ALPHA_DEG:.3f} deg "
          f"= {ALPHA:.6f} rad")
    print(f"v0x = {vacuum_vx:.6f} m/s; v0y = {vacuum_vy:.6f} m/s")
    print(f"Vacuum: flight = {vacuum_time:.6f} s, "
          f"peak time = {vacuum_peak_time:.6f} s, "
          f"height = {vacuum_height:.6f} m, "
          f"range = {vacuum_range:.6f} m")
    print(f"Drag: flight = {t_drag[-1]:.6f} s, "
          f"height = {np.max(y):.6f} m, range = {x[-1]:.6f} m, "
          f"landing speed = {speed[-1]:.6f} m/s")
    print(f"2*Cd: flight = {t_double[-1]:.6f} s, "
          f"height = {np.max(s_double[:, 1]):.6f} m, "
          f"range = {s_double[-1, 0]:.6f} m")
    print(f"Area = {AREA:.8f} m^2; "
          f"drag coefficient per mass = "
          f"{RHO * CD * AREA / (2 * MASS):.8f} 1/m")
    print(f"Mechanical energy lost = "
          f"{mechanical[0] - mechanical[-1]:.6f} J")
    print(f"Figures saved to {OUT}")


if __name__ == "__main__":
    main()
