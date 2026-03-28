import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import json
import numpy as np

with open('sim_data.json', 'r') as f:
    data = json.load(f)

# timing
dt = 0.0001
log_interval = 10
effective_dt = log_interval * dt  # 0.001s per logged frame

# interval frames logged
step: int = 8

# real sim‑time between animation frames
frame_dt = effective_dt * step  # 0.008s
# FuncAnimation interval is in milliseconds
interval_ms = frame_dt * 1000   # 8 ms

# lists for position
positions: list[list[float]] = data['r']
x: list[float] = []
y: list[float] = []
z: list[float] = []

for p in positions:
    x.append(p[0])
    y.append(p[1])
    z.append(p[2])

# list of orientation
psi = data['psi']

def get_direction(theta, phi):
    dx = np.sin(phi)
    dy = -np.sin(theta) * np.cos(phi)
    dz = np.cos(theta) * np.cos(phi)
    return dx, dy, dz

# plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')

pad = 0.1
ax.set_xlim(min(x) - pad, max(x) + pad)
ax.set_ylim(min(y) - pad, max(y) + pad)
ax.set_zlim(min(z) - pad, max(z) + pad)

point, = ax.plot([], [], [], 'ro', markersize=6)        # type: ignore
line_fwd, = ax.plot([], [], [], color='black', linewidth=2)
line_bwd, = ax.plot([], [], [], color='green', linewidth=2)
line_length = max(z) * 0.15
trail, = ax.plot([], [], [], 'b--', linewidth=1, alpha=0.5)

frame_indices = list(range(0, len(x), step))

def update(frame) -> tuple:
    px, py, pz = x[frame], y[frame], z[frame]
    point.set_data_3d([px], [py], [pz])                  # type: ignore
    trail.set_data_3d(x[:frame+1], y[:frame+1], z[:frame+1])  # type: ignore

    theta = psi[frame][0]
    phi = psi[frame][1]
    dx, dy, dz = get_direction(theta, phi)
    mag = np.sqrt(dx**2 + dy**2 + dz**2)
    dx, dy, dz = dx / mag, dy / mag, dz / mag

    line_fwd.set_data_3d([px, px + dx*line_length], [py, py + dy*line_length], [pz, pz + dz*line_length])  # type: ignore
    line_bwd.set_data_3d([px, px - dx*line_length], [py, py - dy*line_length], [pz, pz - dz*line_length])  # type: ignore

    # current sim time for this logged frame
    t = frame * effective_dt
    ax.set_title(f't = {t:.2f} s')

    return (point, trail, line_fwd, line_bwd)

ani = FuncAnimation(fig, update, frames=frame_indices,interval=interval_ms, blit=False)
plt.show()