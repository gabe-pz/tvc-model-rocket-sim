import numpy as np
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.animation import FuncAnimation  # type: ignore
from mpl_toolkits.mplot3d import Axes3D  # type: ignore
import json


def load_data(filepath="sim_data.json"):
    with open(filepath, "r") as f:
        data = json.load(f)
    return data


def animate(data, skip=1):
    positions = data["r"][::skip]
    orientations = data["psi"][::skip]
    n_frames = len(positions)

    pos_arr = np.array(positions)
    pad = 2.0
    x_min, x_max = pos_arr[:, 0].min() - pad, pos_arr[:, 0].max() + pad
    y_min, y_max = pos_arr[:, 1].min() - pad, pos_arr[:, 1].max() + pad
    z_min, z_max = pos_arr[:, 2].min() - pad, pos_arr[:, 2].max() + pad

    trajectory_extent = max(x_max - x_min, y_max - y_min, z_max - z_min)
    arrow_length = trajectory_extent * 0.05

    directions = []
    for i in range(n_frames):
        if i < n_frames - 1:
            d = np.array(positions[i + 1]) - np.array(positions[i])
        else:
            d = np.array(positions[i]) - np.array(positions[i - 1])
        norm = np.linalg.norm(d)
        if norm > 1e-12:
            d = d / norm
        directions.append(d)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    #Trail line
    (trail_line,) = ax.plot([], [], [], "g--", linewidth=0.8, alpha=0.5)

    #Ground plane indicator
    ax.plot([x_min, x_max], [0, 0], [0, 0], "k-", linewidth=0.5, alpha=0.3)

    trail_x, trail_y, trail_z = [], [], []
    rocket_arrow = [None]

    def init():
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        ax.set_zlim(z_min, z_max)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z (up)")
        ax.set_title("TVC Rocket Simulation")
        return (trail_line,)

    def update(frame):
        pos = positions[frame]
        d = directions[frame] * arrow_length

        #Tail is half the body length behind the center
        tail = np.array(pos) - d / 2

        #Remove old arrow
        if rocket_arrow[0] is not None:
            rocket_arrow[0].remove()

        #Purple arrow
        rocket_arrow[0] = ax.quiver( #type: ignore 
            tail[0], tail[1], tail[2],
            d[0], d[1], d[2],
            color="purple",
            linewidth=4,
            arrow_length_ratio=0.3,
        )

        #Update trail
        trail_x.append(pos[0])
        trail_y.append(pos[1])
        trail_z.append(pos[2])
        trail_line.set_data(trail_x, trail_y)
        trail_line.set_3d_properties(trail_z)  # type: ignore

        ax.set_title(f"TVC Rocket  |  t = {frame * skip * 0.001:.3f}s")
        return (trail_line,)

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