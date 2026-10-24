#%% imports and plot settings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.animation import FuncAnimation, PillowWriter
from scipy.signal import find_peaks


plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 12,
    "figure.titlesize": 18
})

#####################%% TWO BODY CIRCULAR ORBIT #####################

#%% plotting two body orbits
N = 2
dt = 0.001

x, y, vx, vy = np.loadtxt("data/data_two_body.txt", unpack=True)

fig, ax = plt.subplots(figsize=(7, 7))

end = 100

ax.plot(x[0:end:2], y[0:end:2],
        color="blue", label="Particle 1 trajectory")
ax.plot(x[1:end:2], y[1:end:2],
        color="red", label="Particle 2 trajectory")

ax.scatter(x[:2], y[:2],
           c=["blue", "red"], s=100,
           edgecolor="black", zorder=5)

ax.quiver(x[:2], y[:2], vx[:2], vy[:2],
          color=["blue", "red"],
          angles="xy", scale_units="xy",
          scale=2, width=0.007, zorder=6)

for i, color in [(0, "blue"), (1, "red")]:
    ax.annotate(
        r"$\vec{v}_0$",
        xy=(x[i] + vx[i]/4, y[i] + vy[i]/4),
        xytext=(-12, 0) if i == 0 else (12, 0),
        textcoords="offset points",
        color=color,
        fontsize=16,
        ha="right" if i == 0 else "left",
        va="center",
        zorder=7
    )

for i, color in [(0, "blue"), (1, "red")]:
    ax.annotate(
        "",
        xy=(x[end-2+i], y[end-2+i]),
        xytext=(x[end-20+i], y[end-20+i]),
        arrowprops=dict(arrowstyle="->", color=color,
                        lw=2, mutation_scale=18)
    )

theta = np.linspace(0, 2*np.pi, 500)
ax.plot(0.5*np.cos(theta), 0.5*np.sin(theta),
        "k--", alpha=0.5, label="Analytical orbit")

ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Two-body circular orbits, R=0.5")
ax.grid(True, alpha=0.3)
ax.legend()

ax.set_xlim(-0.65, 0.70)
ax.set_ylim(-0.65, 0.65)

plt.tight_layout()

plt.show()

#%% find orbital period

distance = 100
max_iteration = len(x)//N
counter=1

while(distance>0.0001):

    distance = np.sqrt((x[counter*2]-x[0])**2 + (y[counter*2]-y[0])**2)
    counter = counter+1
    if(counter>=max_iteration):
        print("Particle did not return within the tolerance.")
        break

else:
    print("Particle returned at step", counter - 1)
    print("Numerical period:", (counter - 1)*dt)   

#%%##########################%% 2D N BODY ###########################

#%% plot initial distribution

N = 10
dt = 0.001

data = np.loadtxt("data/data_2d.txt").reshape(-1, N, 4)

x = data[:, :, 0]
y = data[:, :, 1]
vx = data[:, :, 2]
vy = data[:, :, 3]

t = 0

fig, ax = plt.subplots(figsize=(8, 8))

circle = plt.Circle((0, 0), 1, color="steelblue",
                    fill=False, linestyle="--", linewidth=2)
ax.add_patch(circle)

ax.scatter(x[t], y[t], color="blue", s=50,
           zorder=3, label="Particles")

ax.quiver(x[t], y[t], vx[t], vy[t],
          color="red", angles="xy", scale_units="xy",
          scale=2, width=0.005, label="Velocity")

ax.set_aspect("equal")
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)

ax.set_xlabel("x", fontsize=16)
ax.set_ylabel("y", fontsize=16)
ax.set_title(f"t = {t * dt:.3f}", fontsize=18)

ax.tick_params(axis="both", labelsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=12)

plt.tight_layout()

plt.show()

#%% plot three timestamps of body motion

N = 1000
N_plot = 100
dt = 0.01

x, y, vx, vy = np.loadtxt("data/eps/eps0d001/data_2d.txt", unpack=True)

timesteps = [0, 25, 50]

rng = np.random.default_rng(42)
particle_indices = rng.choice(N, size=N_plot, replace=False)

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, t in zip(axes, timesteps):

    start = t * N

    circle = plt.Circle((0, 0), 1, color="steelblue",
                        fill=False, linestyle="--", linewidth=2)
    ax.add_patch(circle)

    idx = start + particle_indices

    radius = np.sqrt(x[idx]**2 + y[idx]**2)

    inside = radius <= 1
    outside = radius > 1

    ax.scatter(x[idx][inside], y[idx][inside],
               color="blue", s=50, zorder=3, label="Inside R=1")

    ax.scatter(x[idx][outside], y[idx][outside],
               color="orange", s=50, zorder=3, label="Outside R=1")

    ax.quiver(x[idx], y[idx], vx[idx], vy[idx],
              color="red", angles="xy", scale_units="xy",
              scale=2, width=0.005, label="Velocity")

    ax.set_aspect("equal")
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)

    ax.set_xlabel("x", fontsize=16)
    ax.set_ylabel("y", fontsize=16)
    ax.set_title(f"t = {t * dt:.3f}", fontsize=18)

    ax.tick_params(axis="both", labelsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig("plots/2d_3plot_1000N_select10_colour_eps0d0001.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% plotting close encounter

N = 100
dt = 0.01
eps = 0.01
close_distance = 0.5
window_before = 30
window_after = 30

data = np.loadtxt("data/dt/dt0d01/data_2d.txt").reshape(-1, N, 4)
x, y, vx, vy = [data[:, :, i] for i in range(4)]

ax = np.gradient(vx, dt, axis=0)
ay = np.gradient(vy, dt, axis=0)
acc = np.sqrt(ax**2 + ay**2)

encounters = []

for i in range(N):
    for j in range(i + 1, N):
        r = np.sqrt((x[:, i] - x[:, j])**2 + (y[:, i] - y[:, j])**2)
        minima, _ = find_peaks(-r)

        for t in minima:
            if r[t] >= close_distance or t < window_before or t + window_after >= len(x):
                continue

            score = np.max(acc[t-window_before:t+window_after+1, i]) + \
                    np.max(acc[t-window_before:t+window_after+1, j])

            encounters.append((score, t, i, j))

if not encounters:
    raise ValueError("No complete close encounters found.")

score, t_close, p1, p2 = max(encounters)

start = t_close - window_before
end = t_close + window_after + 1
times = np.arange(start, end)*dt

r = np.sqrt((x[:, p1] - x[:, p2])**2 + (y[:, p1] - y[:, p2])**2)

fig, ax1 = plt.subplots(figsize=(9, 5))
ax2 = ax1.twinx()

ax1.plot(times, r[start:end], color="purple", linewidth=2, label="Separation")
ax2.plot(times, acc[start:end, p1], color="blue", linewidth=2, label=f"Particle {p1+1}")
ax2.plot(times, acc[start:end, p2], color="red", linewidth=2, label=f"Particle {p2+1}")

ax1.axvline(t_close*dt, color="black", linestyle="--",
            linewidth=1.5, label="Closest approach")

ax1.set_xlabel("Time", fontsize=16)
ax1.set_ylabel(r"$|\vec{X}_i-\vec{X}_j|$", fontsize=16, color="purple")
ax2.set_ylabel(r"$|\vec{a}_i|$", fontsize=16)
ax1.set_title("Particle separation and acceleration during close encounter", fontsize=18)

ax1.tick_params(axis="both", labelsize=14)
ax2.tick_params(axis="y", labelsize=14)
ax1.grid(True, alpha=0.3)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=12)

plt.tight_layout()
plt.show()

#%% close encounter simulation with different softening parameters

N = 100
t = np.linspace(-0.5, 0.5, 1000)

r_min = 0.02
v = 1.0
r = np.sqrt(r_min**2 + (v*t)**2)

a_0 = 1 / (N * r**2)
a_001 = r / (N * (r**2 + 0.01**2)**1.5)
a_01 = r / (N * (r**2 + 0.1**2)**1.5)

fig, ax = plt.subplots(figsize=(9, 5))

plt.plot(t, a_0, label=r"$\epsilon=0$", linewidth=2)
plt.plot(t, a_001, label=r"$\epsilon=0.01$", linewidth=2)
plt.plot(t, a_01, label=r"$\epsilon=0.1$", linewidth=2)

plt.axvline(0, color="black", linestyle="--",
            label="Closest approach")

plt.xlabel("Time relative to closest approach")
plt.ylabel(r"$|\mathbf{a}_{i\leftarrow j}|$")
plt.title("Gravitational acceleration during a close encounter")
plt.legend()
plt.grid(alpha=0.25)
plt.tight_layout()
plt.savefig("plots/close_encounter_softening.pdf", bbox_inches="tight")
plt.show()

#%% animation

data = np.loadtxt("data/data_2d.txt").reshape(-1, N, 4)

x = data[:, :, 0]
y = data[:, :, 1]
vx = data[:, :, 2]
vy = data[:, :, 3]

frames = np.arange(0, len(data), 10)

fig, ax = plt.subplots(figsize=(8, 8))

circle = plt.Circle((0, 0), 1, color="steelblue",
                    fill=False, linestyle="--", linewidth=2)
ax.add_patch(circle)

particles = ax.scatter([], [], color="blue", s=50,
                       zorder=3, label="Particles")

arrows = ax.quiver(x[0], y[0], vx[0], vy[0],
                   color="red", angles="xy", scale_units="xy",
                   scale=2, width=0.005, label="Velocity")

ax.set_aspect("equal")
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)
ax.set_xlabel("x", fontsize=16)
ax.set_ylabel("y", fontsize=16)
ax.tick_params(axis="both", labelsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=14)

def update(t):

    particles.set_offsets(np.column_stack((x[t], y[t])))

    arrows.set_offsets(np.column_stack((x[t], y[t])))
    arrows.set_UVC(vx[t], vy[t])

    ax.set_title(f"Particle motion: t = {t*dt:.3f}", fontsize=18)

    return particles, arrows

animation = FuncAnimation(
    fig, update, frames=frames,
    interval=50, blit=False
)

animation.save("Nbody_2d.gif", writer=PillowWriter(fps=20))

plt.show()

#%% plot energy

dt = 0.001

Ep, Ek, Lx, Ly, Lz, Vbx, Vby = np.loadtxt("data/quantity_2d.txt", unpack=True)
time = np.arange(len(Ep)) * dt
Etotal = Ep + Ek

plt.figure(figsize=(9, 6))
plt.plot(time, Ek, color="purple", linewidth=2, label=r"$E_k$")
plt.plot(time, Ep, color="darkorange", linewidth=2, label=r"$E_p$")
plt.plot(time, Etotal, color="darkblue", linewidth=2, label=r"$E_T$")

plt.xlabel("Time")
plt.ylabel("Energy (dimensionless)")
plt.title("Energy of 2D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N1000_energy_2d.pdf", bbox_inches="tight")
plt.show()

#%% energy deviation
E = Ek + Ep

relative_error = (E - E[0]) / abs(E[0])

plt.figure(figsize=(9, 5))
plt.plot(time, relative_error)
plt.xlabel("Time (dimensionless)")
plt.ylabel("Relative energy error")
plt.title("Conservation of total mechanical energy")
plt.grid(alpha=0.3)
plt.savefig("plots/N1000_energy_deviation_2d.pdf", bbox_inches="tight")
plt.show()
print("Maximum relative energy error:",
      np.max(np.abs(relative_error)))

#%% energy for epsilon values 
dt = 0.001

eps_values = {0.001: "0d001", 0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.001: "darkviolet", 0.01: "darkred", 0.1: "darkgreen", 1.0: "darkblue"}
linestyles = {0.001: "-", 0.01: "-", 0.1: "-", 1.0: "-"}

plt.figure(figsize=(9, 6))

for eps, name in eps_values.items():
    data = np.loadtxt(f"data/eps/eps{name}/quantity_2d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    plt.plot(time, Etotal, color=colors[eps], linestyle=linestyles[eps],
             linewidth=2, label=rf"$\epsilon={eps}$")

plt.xlabel("Time")
plt.ylabel("Total energy (dimensionless)")
plt.title("Total energy of 2D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()

plt.show()

#%% energy deviation for eps 
dt = 0.001

eps_values = {0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.01: "darkred", 0.1: "darkgreen", 1.0: "darkblue"}
linestyles = {0.01: "-", 0.1: "-", 1.0: "-"}

plt.figure(figsize=(8, 5))

for eps, name in eps_values.items():
    data = np.loadtxt(f"data/eps/eps{name}/quantity_2d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    relative_error = np.abs((Etotal - Etotal[0]) / Etotal[0])

    plt.plot(time, relative_error, color=colors[eps],
             linestyle=linestyles[eps], linewidth=2,
             label=rf"$\epsilon={eps}$")

    print(f"epsilon = {eps}: Maximum relative energy error = {np.max(relative_error):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Absolute relative energy error $|\Delta E/E_0|$")
plt.title("Conservation of total mechanical energy - 2D N-body")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()

plt.show()

#%% energy comparison for different time steps 

dt_values = {0.1: "dt0d1", 0.01: "dt0d01", 0.001: "dt0d001"}
colors = {0.1: "darkblue", 0.01: "darkred", 0.001: "darkgreen"}

plt.figure(figsize=(9, 6))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_2d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    plt.plot(time, Etotal, color=colors[dt], linewidth=2,
             label=rf"$\Delta t={dt}$")

plt.xlabel("Time")
plt.ylabel("Total energy (dimensionless)")
plt.title("Total energy of 2D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_dt_comparison_2d.pdf", bbox_inches="tight")
plt.show()

#%% energy deviation for different time steps 

plt.figure(figsize=(8, 5))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_2d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    relative_error = np.abs((Etotal - Etotal[0]) / Etotal[0])

    plt.plot(time, relative_error, color=colors[dt], linewidth=2,
             label=rf"$\Delta t={dt}$")

    print(f"dt = {dt}: Maximum relative energy error = {np.max(relative_error):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Absolute relative energy error $|\Delta E/E_0|$")
plt.title("Conservation of total mechanical energy - 2D N-body")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_deviation_dt_comparison_2d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum
plt.figure(figsize=(9, 5))
plt.scatter(time, Lz, label=r'$L_z$', color='darkgreen', s=1)
plt.xlabel("Time (dimensionless)")
plt.ylabel("Angular Momentum (dimensionless)")
plt.title("Angular Momentum of 2D N-body system")
plt.grid(alpha=0.3)
plt.legend()
plt.savefig("plots/N100_angular_momentum_2d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum relative error
L_error = (Lz - Lz[0]) / abs(Lz[0])

plt.figure(figsize=(9, 5))
plt.plot(time, L_error, color="green", linewidth=1)

plt.xlabel("Time (dimensionless)")
plt.ylabel(r"Relative angular momentum error $\delta L_z$")
plt.title("Conservation of angular momentum — 2D N-body")
plt.grid(alpha=0.3)
plt.savefig("plots/N100_angular_momentum_rleative_error_2d.pdf", bbox_inches="tight")
plt.show()

print("Initial angular momentum:", Lz[0])
print("Maximum relative error:", np.max(np.abs(L_error)))

#%% barycentric v
dt = 0.001

Ep, Ek, Lx, Ly, Lz, Vbar_x, Vbar_y = np.loadtxt(
    "data/quantity_2d.txt", unpack=True
)

time = np.arange(len(Vbar_x)) * dt

Vbar = np.sqrt(Vbar_x**2 + Vbar_y**2)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(time, Vbar_x, label=r"$v_{\mathrm{bar},x}$")
axes[0].plot(time, Vbar_y, label=r"$v_{\mathrm{bar},y}$")

axes[0].set_xlabel("Time (dimensionless)")
axes[0].set_ylabel("Barycentric velocity")
axes[0].set_title("Barycentric velocity components - 2D")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].plot(time, Vbar, color="green")

axes[1].set_xlabel("Time (dimensionless)")
axes[1].set_ylabel(r"$|\mathbf{v}_{\mathrm{bar}}|$")
axes[1].set_title("Magnitude of barycentric velocity - 2D")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("plots/N100_barycentric_velocity_2d.pdf", bbox_inches="tight")
plt.show()

print("Maximum barycentric speed:", np.max(Vbar))

#%% plot backwards data

N = 10
dt = 0.001

forward = np.loadtxt("data/data_2d.txt").reshape(-1, N, 4)
backward = np.loadtxt("data/data_2d_backwards.txt").reshape(-1, N, 4)

plt.figure(figsize=(8, 8))

for i in range(N):
    plt.plot(forward[:,i,0], forward[:,i,1],
             color="darkblue", linewidth=1, alpha=0.7,
             label="Forward" if i == 0 else None)

    plt.plot(backward[:,i,0], backward[:,i,1],
             color="darkorange", linestyle="--", linewidth=1,
             alpha=0.7, label="Backward" if i == 0 else None)

plt.scatter(forward[0,:,0], forward[0,:,1],
            color="green", s=30, label="Initial positions", zorder=3)

plt.scatter(backward[-1,:,0], backward[-1,:,1],
            color="red", marker="x", s=35,
            label="Recovered positions", zorder=4)

plt.xlabel("x (dimensionless)")
plt.ylabel("y (dimensionless)")
plt.title("RK4 forward and backward trajectories - 2D")
plt.axis("equal")
plt.grid(alpha=0.3)
plt.legend()
plt.show()

#%% forward and backward animation 

N = 10
dt = 0.001

forward = np.loadtxt("data/data_2d.txt").reshape(-1, N, 4)
backward = np.loadtxt("data/data_2d_backwards.txt").reshape(-1, N, 4)

data = np.concatenate((forward, backward[1:]), axis=0)

x = data[:, :, 0]
y = data[:, :, 1]
vx = data[:, :, 2]
vy = data[:, :, 3]

n_forward = len(forward)

frames = np.arange(0, len(data), 10)

fig, ax = plt.subplots(figsize=(8, 8))

circle = plt.Circle((0, 0), 1, color="steelblue",
                    fill=False, linestyle="--", linewidth=2)
ax.add_patch(circle)

particles = ax.scatter([], [], color="blue", s=50,
                       zorder=3, label="Particles")

initial_positions = ax.scatter(
    forward[0, :, 0],
    forward[0, :, 1],
    color="limegreen",
    marker="x",
    s=100,
    linewidths=2,
    zorder=5,
    label="Initial positions"
)

arrows = ax.quiver(x[0], y[0], vx[0], vy[0],
                   color="red", angles="xy", scale_units="xy",
                   scale=2, width=0.005, label="Velocity")

ax.set_aspect("equal")
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)
ax.set_xlabel("x", fontsize=16)
ax.set_ylabel("y", fontsize=16)
ax.tick_params(axis="both", labelsize=14)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=14)

def update(t):

    particles.set_offsets(np.column_stack((x[t], y[t])))

    arrows.set_offsets(np.column_stack((x[t], y[t])))
    arrows.set_UVC(vx[t], vy[t])

    if t < n_forward:
        direction = "Forward"
        current_time = t * dt
    else:
        direction = "Backward"
        current_time = (2*n_forward - 2 - t) * dt

    ax.set_title(
        f"{direction} integration: t = {current_time:.3f}",
        fontsize=18
    )

    return particles, arrows

animation = FuncAnimation(
    fig, update, frames=frames,
    interval=50, blit=False
)

animation.save("Nbody_2d_forward_backward_rk4.gif",
               writer=PillowWriter(fps=20))

plt.show()

#%% ########################### 3D N BODY ###########################

#%% plot initial condition

N = 10

x, y, z, vx, vy, vz = np.loadtxt("data/data_3d.txt", unpack=True)

fig = plt.figure(figsize=(10, 9), dpi=130)
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2*np.pi, 70)
v = np.linspace(0, np.pi, 40)

X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))

ax.plot_surface(X, Y, Z, color='lightskyblue', alpha=0.09, edgecolor='none', shade=False, label="Sphere R=1")
ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5, color='steelblue', alpha=0.13, linewidth=0.45)

ax.scatter(x[:N], y[:N], z[:N], color='blue', s=38, depthshade=False, zorder=5, label="Particles")

ax.quiver(x[:N], y[:N], z[:N], vx[:N], vy[:N], vz[:N], color='red', length=0.28,
          normalize=True, arrow_length_ratio=0.25, linewidth=1.6, label="Velocity")

ax.plot(0,0,[-1,1],'k',alpha=0.8, linewidth=0.5, label="Axis of rotations")

ax.set_box_aspect([1, 1, 1])
ax.view_init(elev=30, azim=-45)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
ax.set_title(f"t=0.0000")
ax.legend(fontsize=12)
plt.savefig("plots/N10_initial_3d.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% plot 3d data time stamps

N = 1000

fig = plt.figure(figsize=(10, 9), dpi=130)
ax = fig.add_subplot(111, projection='3d')

u = np.linspace(0, 2*np.pi, 70)
v = np.linspace(0, np.pi, 40)

X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))

ax.plot_surface(X, Y, Z, color='lightskyblue', alpha=0.09, edgecolor='none', shade=False, label="Sphere R=1")
ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5, color='steelblue', alpha=0.13, linewidth=0.45)

N = 10
start = 1000
end = start+N
ax.scatter(x[0:10], y[0:10], z[0:10], color='blue', s=38, depthshade=False, zorder=5,label=f"i={1}")
ax.scatter(x[start:end], y[start:end], z[start:end], color='m', s=38, depthshade=False, zorder=5,label=f"i={start}")

start = 5000
end = start+N
ax.scatter(x[start:end], y[start:end], z[start:end], color='r', s=38, depthshade=False, zorder=5, label=f"i={start}")

ax.set_box_aspect([1, 1, 1])
ax.view_init(elev=30, azim=-45)
plt.legend()

plt.show()

#%% plot some particles only over 3 time steps

N = 1000
N_plot = 100
dt = 0.1

x, y, z, vx, vy, vz = np.loadtxt("data/data_3d.txt", unpack=True)

timesteps = [0, 2, 5]

r3 = np.sqrt(
    x.reshape(-1, N)**2 +
    y.reshape(-1, N)**2 +
    z.reshape(-1, N)**2
)

rng = np.random.default_rng(42)
particle_indices = np.argsort(r3[50])[-N_plot:]

fig = plt.figure(figsize=(18, 6), dpi=130)

u = np.linspace(0, 2*np.pi, 70)
v = np.linspace(0, np.pi, 40)

X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))

for i, t in enumerate(timesteps):

    ax = fig.add_subplot(1, 3, i+1, projection='3d')

    ax.plot_surface(X, Y, Z, color='lightskyblue', alpha=0.09,
                    edgecolor='none', shade=False)

    ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5,
                      color='steelblue', alpha=0.13, linewidth=0.45)
    idx = t * N + particle_indices

    radius = np.sqrt(x[idx]**2 + y[idx]**2 + z[idx]**2)

    inside = radius <= 1
    outside = radius > 1

    ax.scatter(x[idx][inside], y[idx][inside], z[idx][inside],
            color='blue', s=38, depthshade=False,
            label="Inside R=1")

    ax.scatter(x[idx][outside], y[idx][outside], z[idx][outside],
            color='orange', s=45, depthshade=False,
           label="Outside R=1")

    ax.quiver(x[idx], y[idx], z[idx],
              vx[idx], vy[idx], vz[idx],
              color='red', length=0.28,
              normalize=True, arrow_length_ratio=0.25,
              linewidth=1.6, label="Velocity")

    ax.plot(0, 0, [-1, 1], 'k', alpha=0.8,
            linewidth=0.5, label="Axis of rotation")

    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=30, azim=-45)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_zlim(-1.3, 1.3)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.set_title(f"t={t * dt:.4f}")
    ax.legend(fontsize=10)

plt.tight_layout()
plt.savefig("plots/3d_3plot_1000N_select100_colour.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% trajectory of 1 particle for 3 time steps

N = 100
particle = 90
t_max = 2.0

simulations = {0.1: "dt0d1", 0.01: "dt0d01", 0.001: "dt0d001"}
colors = {0.1: "blue", 0.01: "red", 0.001: "green"}
linestyles = {0.1: "-", 0.01: "--", 0.001: ":"}

data_2d = {}
data_3d = {}

for dt, folder in simulations.items():
    data_2d[dt] = np.loadtxt(f"data/dt/{folder}/data_2d.txt").reshape(-1, N, 4)
    data_3d[dt] = np.loadtxt(f"data/dt/{folder}/data_3d.txt").reshape(-1, N, 6)

fig = plt.figure(figsize=(16, 7), dpi=130)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122, projection="3d")

circle = plt.Circle((0, 0), 1, color="steelblue", fill=False,
                    linestyle="--", linewidth=2, label="Initial radius R=1")
ax1.add_patch(circle)

for dt, data in data_2d.items():
    time = np.arange(len(data)) * dt
    mask = time <= t_max + 1e-10

    x = data[mask, particle, 0]
    y = data[mask, particle, 1]
    print(f"dt={dt}, t_max={t_max}: {mask.sum()} frames")

    ax1.plot(x, y, color=colors[dt], linestyle=linestyles[dt],
             linewidth=2.5, label=fr"$\Delta t={dt}$")
    ax1.scatter(x[0], y[0], color=colors[dt], s=70,
                edgecolor="black", zorder=5)

ax1.set_aspect("equal")
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title(f"2D particle {particle+1} trajectory ($t\\leq {t_max}$)")
ax1.grid(True, alpha=0.3)
ax1.legend()

u = np.linspace(0, 2*np.pi, 70)
v = np.linspace(0, np.pi, 40)

X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))

ax2.plot_surface(X, Y, Z, color="lightskyblue", alpha=0.09,
                 edgecolor="none", shade=False)
ax2.plot_wireframe(X, Y, Z, rstride=5, cstride=5,
                   color="steelblue", alpha=0.13, linewidth=0.45)

coords = []

for dt, data in data_3d.items():
    time = np.arange(len(data)) * dt
    mask = time <= t_max + 1e-10

    x = data[mask, particle, 0]
    y = data[mask, particle, 1]
    z = data[mask, particle, 2]

    ax2.plot(x, y, z, color=colors[dt], linestyle=linestyles[dt],
             linewidth=2.5, label=fr"$\Delta t={dt}$")
    ax2.scatter(x[0], y[0], z[0], color=colors[dt],
                s=60, depthshade=False, edgecolor="black")

    coords.extend([x, y, z])

ax2.plot([0, 0], [0, 0], [-1, 1], "k", alpha=0.8,
         linewidth=0.5, label="Axis of rotation")

limit = max(1.3, max(np.max(np.abs(c)) for c in coords) * 1.1)

ax2.set_xlim(-limit, limit)
ax2.set_ylim(-limit, limit)
ax2.set_zlim(-limit, limit)
ax2.set_box_aspect((1, 1, 1))

ax2.view_init(elev=30, azim=-45)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.set_title(f"3D particle {particle+1} trajectory ($t\\leq {t_max}$)")
ax2.legend()

plt.tight_layout()
plt.savefig("plots/particle_trajectory_2d_3d_dt_comparison.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% compare one particle trajectory for three time steps

N = 1000
particle = 8
t_max = 3.0

simulations = {
    0.1: "dt0d1",
    0.01: "dt0d01",
    0.001: "dt0d001"
}

colors = {
    0.1: "blue",
    0.01: "red",
    0.001: "green"
}

trajectories_2d = {}
trajectories_3d = {}

for dt, folder in simulations.items():

    data2 = np.loadtxt(f"plots/{folder}/data_2d.txt")
    data2 = data2.reshape(-1, N, 4)

    data3 = np.loadtxt(f"plots/{folder}/data_3d.txt")
    data3 = data3.reshape(-1, N, 6)

    n2 = min(len(data2), int(t_max / dt) + 1)
    n3 = min(len(data3), int(t_max / dt) + 1)

    trajectories_2d[dt] = data2[:n2, particle, :2]
    trajectories_3d[dt] = data3[:n3, particle, :3]

#%% expulsion for one time step

N = 1000
dt = 0.1

x2, y2, vx2, vy2 = np.loadtxt(
    "data/data_2d.txt", unpack=True
)

x3, y3, z3, vx3, vy3, vz3 = np.loadtxt(
    "data/data_3d.txt", unpack=True
)

r2 = np.sqrt(x2.reshape(-1, N)**2 +
             y2.reshape(-1, N)**2)

r3 = np.sqrt(x3.reshape(-1, N)**2 +
             y3.reshape(-1, N)**2 +
             z3.reshape(-1, N)**2)

fraction_2d = np.mean(r2 > 1, axis=1)
fraction_3d = np.mean(r3 > 1, axis=1)

t2 = np.arange(len(fraction_2d)) * dt
t3 = np.arange(len(fraction_3d)) * dt

plt.figure(figsize=(8, 5))

plt.plot(t2, fraction_2d, label="2D", linewidth=2)
plt.plot(t3, fraction_3d, label="3D", linewidth=2)

plt.xlabel("Time")
plt.ylabel("Fraction of particles outside R=1")
plt.title("Radial expansion: 2D vs 3D")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("plots/N1000_3d_ejection_rate.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% expulsion for three time steps

N = 100

simulations = {

    0.1:   "dt0d1",
    0.01:  "dt0d01",
    0.001: "dt0d001"
}

colors = {

    0.1:   ("darkblue", "cornflowerblue"),
    0.01:  ("darkred", "lightcoral"),
    0.001: ("darkgreen", "mediumseagreen")
}

plt.figure(figsize=(9, 6))

for dt, folder in simulations.items():

    x2, y2, vx2, vy2 = np.loadtxt(
        f"data/dt/{folder}/data_2d.txt",
        unpack=True
    )

    x3, y3, z3, vx3, vy3, vz3 = np.loadtxt(
        f"data/dt/{folder}/data_3d.txt",
        unpack=True
    )

    r2 = np.sqrt(
        x2.reshape(-1, N)**2 +
        y2.reshape(-1, N)**2
    )

    r3 = np.sqrt(
        x3.reshape(-1, N)**2 +
        y3.reshape(-1, N)**2 +
        z3.reshape(-1, N)**2
    )

    fraction_2d = np.mean(r2 > 1, axis=1)
    fraction_3d = np.mean(r3 > 1, axis=1)

    t2 = np.arange(len(fraction_2d)) * dt
    t3 = np.arange(len(fraction_3d)) * dt

    color_2d, color_3d = colors[dt]

    plt.plot(t2, fraction_2d, color=color_2d, linewidth=2, label=f"2D, dt={dt}")
    plt.plot(t3, fraction_3d, color=color_3d, linewidth=2, label=f"3D, dt={dt}")

plt.xlabel("Time", fontsize=16)
plt.ylabel("Fraction of particles outside R=1", fontsize=16)
plt.title("Radial expansion: 2D vs 3D", fontsize=18)

plt.tick_params(axis="both", labelsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11, ncol=2)

plt.tight_layout()

plt.savefig("plots/N100_ejection_rate_dt_comparison.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% plot energy 3d

dt = 0.01

Ep, Ek, Lx, Ly, Lz, Vbx, Vby, Vbz = np.loadtxt("data/dt/dt0d001/quantity_3d.txt", unpack=True)
time = np.arange(len(Ep)) * dt
Etotal = Ep + Ek

plt.figure(figsize=(9, 6))
plt.plot(time, Ek, color="purple", linewidth=2, label=r"$E_k$")
plt.plot(time, Ep, color="darkorange", linewidth=2, label=r"$E_p$")
plt.plot(time, Etotal, color="darkblue", linewidth=2, label=r"$E_T$")

plt.xlabel("Time")
plt.ylabel("Energy (dimensionless)")
plt.title("Energy of 3D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N1000_energy_3d.pdf", bbox_inches="tight")
plt.show()

#%% energy deviation

relative_error = (Etotal - Etotal[0]) / abs(Etotal[0])

plt.figure(figsize=(9, 5))
plt.plot(time, relative_error)
plt.xlabel("Time")
plt.ylabel("Relative energy error")
plt.title("Conservation of total mechanical energy - 3D N-body")
plt.grid(alpha=0.3)
plt.savefig("plots/N1000_energy_deviation_3d.pdf", bbox_inches="tight")
plt.show()

print("Maximum relative energy error:",
      np.max(np.abs(relative_error)))

#%% energy for epsilon values
dt = 0.001

eps_values = {0.001: "0d001", 0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.001: "darkviolet", 0.01: "darkred", 0.1: "darkgreen", 1.0: "darkblue"}
linestyles = {0.001: "-", 0.01: "-", 0.1: "-", 1.0: "-"}

plt.figure(figsize=(9, 6))

for eps, name in eps_values.items():
    data = np.loadtxt(f"data/eps/eps{name}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    plt.plot(time, Etotal, color=colors[eps], linestyle=linestyles[eps],
             linewidth=2, label=rf"$\epsilon={eps}$")

plt.xlabel("Time")
plt.ylabel("Total energy (dimensionless)")
plt.title("Total energy of 3D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_eps_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% energy deviation for eps
dt = 0.001

eps_values = {0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.01: "darkred", 0.1: "darkgreen", 1.0: "darkblue"}
linestyles = {0.01: "-", 0.1: "-", 1.0: "-"}

plt.figure(figsize=(8, 5))

for eps, name in eps_values.items():
    data = np.loadtxt(f"data/eps/eps{name}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    relative_error = (Etotal - Etotal[0]) / abs(Etotal[0])

    plt.plot(time, relative_error, color=colors[eps],
             linestyle=linestyles[eps], linewidth=2,
             label=rf"$\epsilon={eps}$")

    print(f"epsilon = {eps}: Maximum relative energy error = {np.max(np.abs(relative_error)):.6e}")

plt.yscale('log')
plt.xlabel("Time")
plt.ylabel("Relative energy error")
plt.title("Conservation of total mechanical energy - 3D N-body")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_deviation_eps_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% energy comparison for different time steps

dt_values = {0.1: "dt0d1", 0.01: "dt0d01", 0.001: "dt0d001"}
colors = {0.1: "darkblue", 0.01: "darkred", 0.001: "darkgreen"}

plt.figure(figsize=(9, 6))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    plt.plot(time, Etotal, color=colors[dt], linewidth=2,
             label=rf"$\Delta t={dt}$")

plt.xlabel("Time")
plt.ylabel("Total energy (dimensionless)")
plt.title("Total energy of 3D N-body system")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_dt_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% energy deviation for different time steps

plt.figure(figsize=(8, 5))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    relative_error = np.abs((Etotal - Etotal[0]) / Etotal[0])

    plt.plot(time, relative_error, color=colors[dt], linewidth=2,
             label=rf"$\Delta t={dt}$")

    print(f"dt = {dt}: Maximum relative energy error = {np.max(relative_error):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Absolute relative energy error $|\Delta E/E_0|$")
plt.title("Conservation of total mechanical energy - 3D N-body")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_energy_deviation_dt_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum

plt.figure(figsize=(9, 5))
plt.plot(time, Lx, label=r'$L_x$')
plt.plot(time, Ly, label=r'$L_y$')
plt.plot(time, Lz, label=r'$L_z$')

plt.xlabel("Time (dimensionless)")
plt.ylabel("Angular momentum (dimensionless)")
plt.title("Angular momentum of 3D N-body system")
plt.grid(alpha=0.3)
plt.legend()
plt.show()

#%% angular momentum relative error

L0 = np.sqrt(Lx[0]**2 + Ly[0]**2 + Lz[0]**2)

L_error = np.sqrt(
    (Lx - Lx[0])**2 +
    (Ly - Ly[0])**2 +
    (Lz - Lz[0])**2
) / L0

plt.figure(figsize=(9, 5))
plt.plot(time, L_error, color="green", linewidth=1)

plt.xlabel("Time (dimensionless)")
plt.ylabel("Relative angular momentum error")
plt.title("Conservation of angular momentum - 3D N-body")
plt.grid(alpha=0.3)
plt.show()

print("Initial angular momentum magnitude:", L0)
print("Maximum relative angular momentum error:",
      np.max(L_error))

#%% angular momentum comparison for different time steps

dt_values = {0.1: "dt0d1", 0.01: "dt0d01", 0.001: "dt0d001"}
colors = {0.1: "darkblue", 0.01: "darkred", 0.001: "darkgreen"}
linestyles = {"Lx": "-", "Ly": "--", "Lz": ":"}

plt.figure(figsize=(9, 6))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_3d.txt")
    Lx, Ly, Lz = data[:, 2], data[:, 3], data[:, 4]
    time = np.arange(len(Lx)) * dt

    plt.plot(time, Lx, color=colors[dt], linestyle="-", linewidth=2,
             label=rf"$L_x$, $\Delta t={dt}$")
    plt.plot(time, Ly, color=colors[dt], linestyle="--", linewidth=2,
             label=rf"$L_y$, $\Delta t={dt}$")
    plt.plot(time, Lz, color=colors[dt], linestyle=":", linewidth=2,
             label=rf"$L_z$, $\Delta t={dt}$")

plt.xlabel("Time")
plt.ylabel("Angular momentum (dimensionless)")
plt.title("Angular momentum of 3D N-body system")
plt.grid(alpha=0.3)
plt.legend(ncol=3)

plt.tight_layout()
plt.savefig("plots/N100_angular_momentum_dt_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum relative error for different time steps

plt.figure(figsize=(8, 5))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_3d.txt")
    Lx, Ly, Lz = data[:, 2], data[:, 3], data[:, 4]
    time = np.arange(len(Lx)) * dt

    L0 = np.sqrt(Lx[0]**2 + Ly[0]**2 + Lz[0]**2)

    L_error = np.sqrt(
        (Lx - Lx[0])**2 +
        (Ly - Ly[0])**2 +
        (Lz - Lz[0])**2
    ) / L0

    plt.plot(time, L_error, color=colors[dt], linewidth=2,
             label=rf"$\Delta t={dt}$")

    print(f"dt = {dt}: Maximum relative angular momentum error = {np.max(L_error):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Relative angular momentum error $|\Delta \mathbf{L}|/|\mathbf{L}_0|$")
plt.title("Conservation of angular momentum - 3D N-body")
plt.grid(alpha=0.3)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_angular_momentum_deviation_dt_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% barycentric velocity 

Vbar = np.sqrt(Vbx**2 + Vby**2 + Vbz**2)

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(time, Vbx, label=r"$v_{\mathrm{bar},x}$")
axes[0].plot(time, Vby, label=r"$v_{\mathrm{bar},y}$")
axes[0].plot(time, Vbz, label=r"$v_{\mathrm{bar},z}$")

axes[0].set_xlabel("Time (dimensionless)")
axes[0].set_ylabel("Barycentric velocity")
axes[0].set_title("Barycentric velocity components - 3D")
axes[0].legend()
axes[0].grid(alpha=0.3)

axes[1].plot(time, Vbar, color="green")

axes[1].set_xlabel("Time (dimensionless)")
axes[1].set_ylabel(r"$|\mathbf{v}_{\mathrm{bar}}|$")
axes[1].set_title("Magnitude of barycentric velocity - 3D")
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()

print("Maximum barycentric speed:", np.max(Vbar))

#%% barycentric velocity comparison for different time steps 

dt_values = {0.1: "dt0d1", 0.01: "dt0d01", 0.001: "dt0d001"}
colors = {0.1: "darkblue", 0.01: "darkred", 0.001: "darkgreen"}
linestyles = {"x": "-", "y": "--", "z": ":"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for dt, folder in dt_values.items():
    data = np.loadtxt(f"data/dt/{folder}/quantity_3d.txt")
    Vbx, Vby, Vbz = data[:, 5], data[:, 6], data[:, 7]
    time = np.arange(len(Vbx)) * dt

    Vbar = np.sqrt(Vbx**2 + Vby**2 + Vbz**2)

    axes[0].plot(time, Vbx, color=colors[dt], linestyle="-", linewidth=2,
                 label=rf"$v_{{\mathrm{{bar}},x}}$, $\Delta t={dt}$")
    axes[0].plot(time, Vby, color=colors[dt], linestyle="--", linewidth=2,
                 label=rf"$v_{{\mathrm{{bar}},y}}$, $\Delta t={dt}$")
    axes[0].plot(time, Vbz, color=colors[dt], linestyle=":", linewidth=2,
                 label=rf"$v_{{\mathrm{{bar}},z}}$, $\Delta t={dt}$")

    axes[1].plot(time, Vbar, color=colors[dt], linewidth=2,
                 label=rf"$\Delta t={dt}$")

    print(f"dt = {dt}: Maximum barycentric speed = {np.max(Vbar):.6e}")

axes[0].set_xlabel("Time")
axes[0].set_ylabel("Barycentric velocity (dimensionless)")
axes[0].set_title("Barycentric velocity components - 3D")
axes[0].grid(alpha=0.3)
axes[0].legend(ncol=3, fontsize=9)

axes[1].set_xlabel("Time")
axes[1].set_ylabel(r"$|\mathbf{v}_{\mathrm{bar}}|$")
axes[1].set_title("Magnitude of barycentric velocity - 3D")
axes[1].grid(alpha=0.3)
axes[1].legend()

plt.tight_layout()
plt.savefig("plots/N100_barycentric_velocity_dt_comparison_3d.pdf",
            bbox_inches="tight")
plt.show()

#%% epsilon comparison

N = 100
particle = 67
dt = 0.01
t_max = 3
n_steps = int(t_max / dt) + 1
eps_values = {0.001: "0d001", 0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.001: "blue", 0.01: "red", 0.1: "green", 1.0: "purple"}
linestyles = {0.001: "-", 0.01: "--", 0.1: ":", 1.0: "-."}

fig = plt.figure(figsize=(13, 6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122, projection="3d")

for eps, name in eps_values.items():
    data2 = np.loadtxt(f"data/data_2d_eps{name}.txt").reshape(-1, N, 4)[:n_steps]
    data3 = np.loadtxt(f"data/data_3d_eps{name}.txt").reshape(-1, N, 6)[:n_steps]

    x2, y2 = data2[:, particle-1, 0], data2[:, particle-1, 1]
    x3, y3, z3 = data3[:, particle-1, 0], data3[:, particle-1, 1], data3[:, particle-1, 2]

    ax1.plot(x2, y2, color=colors[eps], linestyle=linestyles[eps], linewidth=2, label=fr"$\epsilon={eps}$")
    ax2.plot(x3, y3, z3, color=colors[eps], linestyle=linestyles[eps], linewidth=2, label=fr"$\epsilon={eps}$")

ax1.scatter(x2[0], y2[0], color="black", s=50, zorder=5)
ax2.scatter(x3[0], y3[0], z3[0], color="black", s=50, depthshade=False)

theta = np.linspace(0, 2*np.pi, 200)
ax1.plot(np.cos(theta), np.sin(theta), color="steelblue", linestyle="--", alpha=0.7, label="Initial radius R=1")

u = np.linspace(0, 2*np.pi, 70)
v = np.linspace(0, np.pi, 40)

X = np.outer(np.cos(u), np.sin(v))
Y = np.outer(np.sin(u), np.sin(v))
Z = np.outer(np.ones_like(u), np.cos(v))

ax2.plot_surface(X, Y, Z, color="lightskyblue", alpha=0.09, edgecolor="none", shade=False)
ax2.plot_wireframe(X, Y, Z, rstride=5, cstride=5, color="steelblue", alpha=0.13, linewidth=0.45)

ax1.set_title(f"2D particle {particle} trajectory", fontsize=16)
ax2.set_title(f"3D particle {particle} trajectory", fontsize=16)

ax1.set_xlabel("x", fontsize=14)
ax1.set_ylabel("y", fontsize=14)
ax2.set_xlabel("x", fontsize=14)
ax2.set_ylabel("y", fontsize=14)
ax2.set_zlabel("z", fontsize=14)

ax2.set_xlim(-1.3, 1.3)
ax2.set_ylim(-1.3, 1.3)
ax2.set_zlim(-1.3, 1.3)
ax2.set_box_aspect((1, 1, 1))

ax1.set_aspect("equal")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=10)
ax2.legend(fontsize=10)

plt.tight_layout()
plt.savefig("plots/epsilon_2d_3d_dt_particle_67.pdf", dpi=300, bbox_inches="tight")
plt.show()

#%% expulsion for three epsilon values

N = 100
dt = 0.01

eps_values = {0.001: "0d001", 0.01: "0d01", 0.1: "0d1", 1.0: "1"}
colors = {0.001: ("darkblue", "cornflowerblue"), 0.01: ("darkred", "lightcoral"), 0.1: ("darkgreen", "mediumseagreen"), 1.0: ("darkviolet", "plum")}

plt.figure(figsize=(9, 6))

for eps, name in eps_values.items():

    x2, y2, vx2, vy2 = np.loadtxt(f"data/data_2d_eps{name}.txt", unpack=True)
    x3, y3, z3, vx3, vy3, vz3 = np.loadtxt(f"data/data_3d_eps{name}.txt", unpack=True)

    r2 = np.sqrt(x2.reshape(-1, N)**2 + y2.reshape(-1, N)**2)
    r3 = np.sqrt(x3.reshape(-1, N)**2 + y3.reshape(-1, N)**2 + z3.reshape(-1, N)**2)

    fraction_2d = np.mean(r2 > 1, axis=1)
    fraction_3d = np.mean(r3 > 1, axis=1)

    t2 = np.arange(len(fraction_2d)) * dt
    t3 = np.arange(len(fraction_3d)) * dt

    color_2d, color_3d = colors[eps]

    plt.plot(t3, fraction_3d, color=color_3d, linewidth=2, label=fr"3D, $\epsilon={eps}$")

plt.xlabel("Time", fontsize=16)
plt.ylabel("Fraction of particles outside R=1", fontsize=16)
plt.title("Influence of gravitational softening on radial expansion - 3D", fontsize=18)

plt.tick_params(axis="both", labelsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11, ncol=2)

plt.tight_layout()
plt.savefig("N100_ejection_rate_eps_comparison_3d.pdf", dpi=300, bbox_inches="tight")
plt.show()


#%% plot trajectories of 3 integration methods

N = 100
N_plot = 10
R = 1.0

methods = ["leapfrog", "rk4", "rk8"]
colors = plt.cm.tab10(np.arange(N_plot))

data = {}

for method in methods:
    filename = f"data/integration_method/{method}/data_3d.txt"
    data[method] = np.loadtxt(filename).reshape(-1, N, 6)

u = np.linspace(0, 2*np.pi, 60)
v = np.linspace(0, np.pi, 40)

sphere_x = R * np.outer(np.cos(u), np.sin(v))
sphere_y = R * np.outer(np.sin(u), np.sin(v))
sphere_z = R * np.outer(np.ones_like(u), np.cos(v))

all_pos = np.concatenate([
    data[m][:, :N_plot, :3].reshape(-1, 3)
    for m in methods
])

mins = np.minimum(all_pos.min(axis=0), -R)
maxs = np.maximum(all_pos.max(axis=0), R)

center = (mins + maxs) / 2
half_range = np.max(maxs - mins) / 2 * 1.05

fig = plt.figure(figsize=(18, 6))

for j, method in enumerate(methods):

    ax = fig.add_subplot(1, 3, j + 1, projection="3d")

    ax.plot_surface(
        sphere_x, sphere_y, sphere_z,
        color="steelblue",
        alpha=0.08,
        edgecolor="none",
        shade=False
    )

    pos = data[method][:, :N_plot, :3]

    for i in range(N_plot):

        x = pos[:, i, 0]
        y = pos[:, i, 1]
        z = pos[:, i, 2]

        ax.plot(
            x, y, z,
            color=colors[i],
            linewidth=1.4,
            alpha=0.9
        )

        ax.scatter(
            x[-1], y[-1], z[-1],
            color=colors[i],
            s=45,
            edgecolor="black",
            linewidth=0.5,
            depthshade=False
        )

    ax.set_title(method.upper(), fontsize=15, pad=15)

    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("y", fontsize=12)
    ax.set_zlabel("z", fontsize=12)

    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)

    ax.set_box_aspect((1, 1, 1), zoom=1.25)

plt.subplots_adjust(left=0.01, right=0.99, wspace=0.02)

plt.show()

#%% compare trajectories numerically
N = 100
dt = 0.001
methods = ["leapfrog", "rk4", "rk8"]

data = {}
for method in methods:
    filename = f"data/integration_method/{method}/data_3d.txt"
    data[method] = np.loadtxt(filename).reshape(-1, N, 6)

nframes = min(len(data[m]) for m in methods)
time = np.arange(nframes) * dt

pos_lf = data["leapfrog"][:nframes, :, :3]
pos_rk4 = data["rk4"][:nframes, :, :3]
pos_rk8 = data["rk8"][:nframes, :, :3]

diff_rk4 = np.linalg.norm(pos_lf - pos_rk4, axis=2)
diff_rk8 = np.linalg.norm(pos_lf - pos_rk8, axis=2)
diff_rk4_rk8 = np.linalg.norm(pos_rk4 - pos_rk8, axis=2)

mean_rk4 = np.mean(diff_rk4, axis=1)
mean_rk8 = np.mean(diff_rk8, axis=1)
mean_rk4_rk8 = np.mean(diff_rk4_rk8, axis=1)

plt.figure(figsize=(8, 5))

plt.plot(time, mean_rk4,
         label="Leapfrog vs RK4",
         linewidth=4, color="darkblue")

plt.plot(time, mean_rk8,
         label="Leapfrog vs RK8",
         linewidth=4, color="darkorange", linestyle="--")

plt.plot(time, mean_rk4_rk8,
         label="RK4 vs RK8",
         linewidth=4, color="darkgreen", linestyle="-.")

plt.xlabel("Time")
plt.ylabel("Mean particle position difference")
plt.title("Difference in Particle Positions")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.show()

print(f"Final mean difference (LF vs RK4): {mean_rk4[-1]:.6e}")
print(f"Final mean difference (LF vs RK8): {mean_rk8[-1]:.6e}")
print(f"Final mean difference (RK4 vs RK8): {mean_rk4_rk8[-1]:.6e}")

print(f"Maximum mean difference (LF vs RK4): {mean_rk4.max():.6e}")
print(f"Maximum mean difference (LF vs RK8): {mean_rk8.max():.6e}")
print(f"Maximum mean difference (RK4 vs RK8): {mean_rk4_rk8.max():.6e}")

#%% compare integration energy
dt = 0.001

methods = {"leapfrog": "Leapfrog", "rk4": "RK4", "rk8": "RK8"}
colors = {"leapfrog": "darkgreen", "rk4": "darkblue", "rk8": "darkorange"}
linestyles = {"leapfrog": "-", "rk4": "-", "rk8": "--"}

plt.figure(figsize=(9, 6))

for method, label in methods.items():
    data = np.loadtxt(f"data/integration_method/{method}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt

    plt.plot(time, Etotal, color=colors[method], linewidth=4, label=label, linestyle=linestyles[method])

plt.xlabel("Time")
plt.ylabel("Total energy (dimensionless)")
plt.title("Total Energy of 3D N-body System")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.ticklabel_format(axis="y", style="plain", useOffset=False)
plt.savefig("plots/N100_energy_method_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% relative energy

dt = 0.001

methods = {"leapfrog": "Leapfrog", "rk4": "RK4", "rk8": "RK8"}
colors = {"leapfrog": "darkgreen", "rk4": "darkblue", "rk8": "darkorange"}
linestyles = {"leapfrog": "-", "rk4": "-", "rk8": "--"}

plt.figure(figsize=(8, 5))

for method, label in methods.items():
    data = np.loadtxt(f"data/integration_method/{method}/quantity_3d.txt")
    Ep, Ek = data[:, 0], data[:, 1]
    Etotal = Ep + Ek
    time = np.arange(len(Etotal)) * dt
    relative_error = (Etotal - Etotal[0]) / abs(Etotal[0])

    plt.plot(time, np.abs(relative_error), color=colors[method],
             linewidth=2, label=label)
    print(f"{label}: Maximum relative energy error = {np.max(np.abs(relative_error)):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel("Absolute relative energy error")
plt.title("Conservation of Total Mechanical Energy - 3D N-body")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("plots/N100_energy_deviation_method_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum

dt = 0.001
methods = {"leapfrog": "Leapfrog", "rk4": "RK4", "rk8": "RK8"}
colors = {"leapfrog": "darkgreen", "rk4": "darkblue", "rk8": "darkorange"}
linestyles = {"Lx": "-", "Ly": "--", "Lz": ":"}

plt.figure(figsize=(9, 6))

for method, label in methods.items():
    data = np.loadtxt(f"data/integration_method/{method}/quantity_3d.txt")
    Lx, Ly, Lz = data[:, 2], data[:, 3], data[:, 4]
    time = np.arange(len(Lx)) * dt
    for L, values in zip(linestyles, [Lx, Ly, Lz]):
        plt.plot(time, values, color=colors[method],
                 linestyle=linestyles[L], linewidth=2,
                 label=rf"${L[0]}_{L[1]}$, {label}")

plt.xlabel("Time")
plt.ylabel("Angular momentum (dimensionless)")
plt.title("Angular Momentum: Integration Method Comparison")
plt.grid(alpha=0.3)
plt.legend(ncol=3)
plt.tight_layout()
plt.savefig("plots/N100_angular_momentum_method_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% angular momentum relative error

dt = 0.001
methods = {"leapfrog": "Leapfrog", "rk4": "RK4", "rk8": "RK8"}
colors = {"leapfrog": "darkgreen", "rk4": "darkblue", "rk8": "darkorange"}

plt.figure(figsize=(8, 5))

for method, label in methods.items():
    data = np.loadtxt(f"data/integration_method/{method}/quantity_3d.txt")
    Lx, Ly, Lz = data[:, 2], data[:, 3], data[:, 4]
    time = np.arange(len(Lx)) * dt
    L0 = np.sqrt(Lx[0]**2 + Ly[0]**2 + Lz[0]**2)
    L_error = np.sqrt((Lx-Lx[0])**2 + (Ly-Ly[0])**2 + (Lz-Lz[0])**2) / L0

    plt.plot(time, L_error, color=colors[method], linewidth=2, label=label)
    print(f"{label}: Maximum relative angular momentum error = {np.max(L_error):.6e}")

plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"Relative angular momentum error $|\Delta \vec{L}|/|\vec{L}_0|$")
plt.title("Conservation of Angular Momentum - 3D N-body")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("plots/N100_angular_momentum_deviation_method_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% barycentric velocity

dt = 0.001
methods = {"leapfrog": "Leapfrog", "rk4": "RK4", "rk8": "RK8"}
colors = {"leapfrog": "darkgreen", "rk4": "darkblue", "rk8": "darkorange"}

linestyles = {"x": "-", "y": "--", "z": ":"}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for method, label in methods.items():
    data = np.loadtxt(f"data/integration_method/{method}/quantity_3d.txt")
    Vbx, Vby, Vbz = data[:, 5], data[:, 6], data[:, 7]
    time = np.arange(len(Vbx)) * dt
    Vbar = np.sqrt(Vbx**2 + Vby**2 + Vbz**2)

    for component, values in zip(linestyles, [Vbx, Vby, Vbz]):
        axes[0].plot(time, values, color=colors[method], linestyle=linestyles[component],
                     linewidth=2, label=rf"$v_{{\mathrm{{bar}},{component}}}$, {label}")
    axes[1].plot(time, Vbar, color=colors[method], linewidth=2, label=label)
    print(f"{label}: Maximum barycentric speed = {np.max(Vbar):.6e}")

axes[0].set_xlabel("Time")
axes[0].set_ylabel("Barycentric velocity (dimensionless)")
axes[0].set_title("Barycentric Velocity Components - 3D")
axes[0].grid(alpha=0.3)
axes[0].legend(ncol=3, fontsize=9)

axes[1].set_xlabel("Time")
axes[1].set_ylabel(r"$|\vec{v}_{\mathrm{bar}}|$")
axes[1].set_title("Magnitude of Barycentric Velocity - 3D")
axes[1].grid(alpha=0.3)
axes[1].legend()

plt.tight_layout()
plt.savefig("plots/N100_barycentric_velocity_method_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% backward error

N = 100
dim = 3

forward = np.loadtxt(
    "data/data_3d.txt"
).reshape(-1, N, 2*dim)

backward = np.loadtxt(
    "data/data_3d_backwards.txt"
).reshape(-1, N, 2*dim)

position_error = np.linalg.norm(
    backward[-1, :, :dim] - forward[0, :, :dim],
    axis=1
)

velocity_error = np.linalg.norm(
    backward[-1, :, dim:] - forward[0, :, dim:],
    axis=1
)

print("Number of steps:", len(forward) - 1)
print("Mean position error:", np.mean(position_error))
print("Maximum position error:", np.max(position_error))
print("Mean velocity error:", np.mean(velocity_error))
print("Maximum velocity error:", np.max(velocity_error))

#%% position errors

steps = np.array([10, 100, 1000, 10000])

leapfrog = np.array([
    1.1969591984239968e-18,
    1.3801025948792657e-17,
    5.108502929774833e-15,
    1.1794785464142628
])

rk4 = np.array([
    7.498257389786673e-13,
    5.450187635402142e-9,
    2.1196632282910492e-6,
    1.0112505256614117
])

rk8 = np.array([
    2.9940891620166473e-16,
    2.5098833064426608e-14,
    1.110252964176191e-11,
    1.0349574212308135
])

leapfrog_soft = 1.6422613416739514e-12

plt.figure(figsize=(9, 5.5))

plt.plot(steps, leapfrog, 'o-', color='steelblue',
         linewidth=2, markersize=7, label='Leapfrog')

plt.plot(steps, rk4, 'o-', color='darkorange',
         linewidth=2, markersize=7, label='RK4')

plt.plot(steps, rk8, 'o-', color='seagreen',
         linewidth=2, markersize=7, label='RK8')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('Number of Integration Steps')
plt.ylabel('Mean Position Recovery Error')
plt.title(r'Time Reversibility ($N=100$, $\Delta t=0.001$)')

plt.grid(True, which='major', alpha=0.3)
plt.grid(True, which='minor', alpha=0.12)
plt.legend()

plt.tight_layout()
plt.savefig("plots/N100_reversibility_comparison_3d.pdf", bbox_inches="tight")
plt.show()

#%% #################### COMPUTATIONAL COMPLEXITY ####################
#%% complexity for leapfrog

N_2d = np.array([10,100,500,1000,2000, 5000])
N_3d = N_2d

T_2d = np.array([0.04057, 0.7051, 14.52, 57.07, 226.8, 1460.9])
T_3d = np.array([0.05330, 0.8213, 16.13, 65.14, 250.7, 1528.96])

coeff_2d = np.polyfit(N_2d, T_2d, 2)
coeff_3d = np.polyfit(N_3d, T_3d, 2)

poly_2d = np.poly1d(coeff_2d)
poly_3d = np.poly1d(coeff_3d)

N_fit = np.linspace(0, np.max(N_2d), 500)

plt.figure(figsize=(8, 5))

plt.scatter(N_2d, T_2d, color="darkblue", label="2D data", zorder=3)
plt.scatter(N_3d, T_3d, color="darkred", label="3D data", zorder=3)

plt.plot(N_fit, poly_2d(N_fit), color="darkblue",
         linestyle="--", label="2D quadratic fit")
plt.plot(N_fit, poly_3d(N_fit), color="darkred",
         linestyle="--", label="3D quadratic fit")

plt.xlabel("Number of particles, $N$")
plt.ylabel("Execution time (s)")
plt.title("Computational Complexity of Leapfrog Integration")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig("plots/computational_complexity.pdf")
plt.show()

print("2D coefficients [a, b, c]:", coeff_2d)
print("3D coefficients [a, b, c]:", coeff_3d)

#%% checking order by getting slope of log log plot

N = N_2d

slope_2d, intercept_2d = np.polyfit(
    np.log10(N[2:]), np.log10(T_2d[2:]), 1
)

slope_3d, intercept_3d = np.polyfit(
    np.log10(N[2:]), np.log10(T_3d[2:]), 1
)

plt.figure(figsize=(8, 5))

plt.scatter(N[2:], T_2d[2:], color="darkblue", label="2D data")
plt.scatter(N[2:], T_3d[2:], color="darkred", label="3D data")

plt.plot(N[2:], 10**intercept_2d * N[2:]**slope_2d,
         "--", color="darkblue",
         label=fr"2D fit: slope = {slope_2d:.3f}")

plt.plot(N[2:], 10**intercept_3d * N[2:]**slope_3d,
         "--", color="darkred",
         label=fr"3D fit: slope = {slope_3d:.3f}")

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Number of particles, $N$")
plt.ylabel("Execution time (s)")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.title("Computational Complexity (Log–Log Scale)")
plt.tight_layout()
plt.savefig("plots/computational_complexity_log.pdf")
plt.show()

#%% complexity for rk

N = np.array([10, 100, 500, 1000, 2000])

rk4_2d = np.array([0.061479, 2.262227, 53.033375,
                   227.460165, 844.786652])

rk4_3d = np.array([0.067863, 2.512471, 58.356186,
                   231.375005, 923.357640])

rk8_2d = np.array([0.112363, 6.404672, 156.996522,
                   630.150703, 2493.946722])

rk8_3d = np.array([0.124942, 7.227244, 172.005377,
                   688.309639, 2734.207511])

datasets = {
    "RK4 - 2D": (rk4_2d, "darkblue"),
    "RK4 - 3D": (rk4_3d, "darkred"),
    "RK8 - 2D": (rk8_2d, "darkgreen"),
    "RK8 - 3D": (rk8_3d, "darkorange")
}

#%% quadratic fit

coefficients = {}
polynomials = {}

for label, (T, color) in datasets.items():
    coefficients[label] = np.polyfit(N, T, 2)
    polynomials[label] = np.poly1d(coefficients[label])

N_fit = np.linspace(0, np.max(N), 500)

plt.figure(figsize=(8, 5))

for label, (T, color) in datasets.items():
    plt.scatter(N, T, color=color,
                label=f"{label} data", zorder=3)

for label, (T, color) in datasets.items():
    plt.plot(N_fit, polynomials[label](N_fit),
             color=color, linestyle="--",
             label=f"{label} quadratic fit")

plt.xlabel("Number of particles, $N$")
plt.ylabel("Execution time (s)")
plt.title("Computational Complexity of RK4 and RK8 Integration")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

plt.savefig("plots/rk_quadratic_fit.pdf")
plt.show()

for label, coeff in coefficients.items():
    print(f"{label} coefficients [a, b, c]:", coeff)

#%% checking order using log-log plot

N_log = N[2:]

slopes = {}
intercepts = {}

for label, (T, color) in datasets.items():
    slopes[label], intercepts[label] = np.polyfit(
        np.log10(N_log), np.log10(T[2:]), 1
    )

plt.figure(figsize=(8, 5))

for label, (T, color) in datasets.items():
    plt.scatter(N_log, T[2:], color=color,
                label=f"{label} data", zorder=3)

for label, (T, color) in datasets.items():
    plt.plot(
        N_log,
        10**intercepts[label] * N_log**slopes[label],
        "--", color=color,
        label=f"{label} fit: slope = {slopes[label]:.3f}"
    )

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Number of particles, $N$")
plt.ylabel("Execution time (s)")
plt.title("Computational Complexity (Log–Log Scale)")
plt.legend()
plt.grid(alpha=0.3, which="both")
plt.tight_layout()

plt.savefig("plots/rk_loglog_fit.pdf")
plt.show()

for label, slope in slopes.items():
    print(f"{label} scaling exponent: {slope:.6f}")
