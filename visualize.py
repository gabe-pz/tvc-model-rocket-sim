import numpy as np
import matplotlib.pyplot as plt #type: ignore
from matplotlib.animation import FuncAnimation #type: ignore
from mpl_toolkits.mplot3d import Axes3D #type: ignore
import json


def load_data(filepath="sim_data.json"):
    with open(filepath, "r") as f:
        data = json.load(f)
    return data


def rocket_body_points(pos, psi, length=0.5):
    """
    Given position [x,y,z] and orientation [theta, phi],
    returns the nose and tail points of the rocket as two 3D points.
    theta = rotation about x (pitch)
    phi = rotation about y (yaw-ish)
    """
    theta, phi = psi

    # rocket axis direction in world frame
    # body z-axis rotated by theta about x, phi about y
    dx = -np.sin(phi) * np.cos(theta)
    dy = np.sin(theta)
    dz = np.cos(theta) * np.cos(phi)

    direction = np.array([dx, dy, dz])
    half = (length / 2) * direction

    nose = np.array(pos) + half
    tail = np.array(pos) - half

    return nose, tail


def animate(data, skip=1):
    """
    Animates rocket orientation and position over time.
    data: dict with keys 'r', 'psi', each a list of [x,y,z] or [theta,phi]
    skip: only draw every N-th frame (speeds up playback)
    """
    positions = data["r"][::skip]
    orientations = data["psi"][::skip]
    n_frames = len(positions)

    # figure out axis limits from trajectory
    pos_arr = np.array(positions)
    pad = 1.0
    x_min, x_max = pos_arr[:, 0].min() - pad, pos_arr[:, 0].max() + pad
    y_min, y_max = pos_arr[:, 1].min() - pad, pos_arr[:, 1].max() + pad
    z_min, z_max = pos_arr[:, 2].min() - pad, pos_arr[:, 2].max() + pad

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # rocket body line and nose marker
    (rocket_line,) = ax.plot([], [], [], "b-", linewidth=4)
    (nose_dot,) = ax.plot([], [], [], "ro", markersize=6)
    (trail_line,) = ax.plot([], [], [], "g--", linewidth=0.8, alpha=0.5)

    # ground plane indicator
    ax.plot(
        [x_min, x_max], [0, 0], [0, 0], "k-", linewidth=0.5, alpha=0.3
    )

    trail_x, trail_y, trail_z = [], [], []

    def init():
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z (up)")
        ax.set_title("TVC Rocket Simulation")
        return rocket_line, nose_dot, trail_line

    def update(frame):
        pos = positions[frame]
        psi = orientations[frame]

        nose, tail = rocket_body_points(pos, psi, length=0.5)

        rocket_line.set_data([tail[0], nose[0]], [tail[1], nose[1]])
        rocket_line.set_3d_properties([tail[2], nose[2]]) #type: ignore

        nose_dot.set_data([nose[0]], [nose[1]])
        nose_dot.set_3d_properties([nose[2]])#type: ignore

        trail_x.append(pos[0])
        trail_y.append(pos[1])
        trail_z.append(pos[2])
        trail_line.set_data(trail_x, trail_y)
        trail_line.set_3d_properties(trail_z)#type: ignore

        ax.set_title(f"TVC Rocket  |  t = {frame * skip * 0.001:.3f}s")

        return rocket_line, nose_dot, trail_line

    anim = FuncAnimation(
        fig,
        update,
        frames=n_frames,
        init_func=init,
        interval=1, 
        blit=False,
    )

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    data = load_data("sim_data.json")
    animate(data, skip=10)