import numpy as np
import matplotlib.pyplot as plt  # type: ignore
from matplotlib.animation import FuncAnimation  # type: ignore
from mpl_toolkits.mplot3d import Axes3D  # type: ignore
from mpl_toolkits.mplot3d import proj3d  # type: ignore
import json


def load_data(filepath="sim_data.json"):
    # Read the simulation data from a JSON file
    with open(filepath, "r") as f:
        data = json.load(f)
    return data


def animate(data, skip=1):
    # Grab every skip-th frame of position and orientation
    positions = data["r"][::skip]
    orientations = data["psi"][::skip]
    n_frames = len(positions)

    # Convert to numpy array so we can slice columns easily
    pos_arr = np.array(positions)

    # Work out axis limits with some breathing room
    pad = 2.0
    x_min = pos_arr[:, 0].min() - pad
    x_max = pos_arr[:, 0].max() + pad
    y_min = pos_arr[:, 1].min() - pad
    y_max = pos_arr[:, 1].max() + pad
    z_min = pos_arr[:, 2].min() - pad
    z_max = pos_arr[:, 2].max() + pad

    # Arrow size scales with how far the rocket travels
    trajectory_extent = max(x_max - x_min, y_max - y_min, z_max - z_min)
    arrow_length = trajectory_extent * 0.05

    # Pre-compute a unit direction vector for each frame
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

    # Set up the 3D figure
    fig = plt.figure(figsize=(10, 8))  # type: ignore
    ax = fig.add_subplot(111, projection="3d")  # type: ignore

    #Draw a filled green ground plane at z = 0
    ground_x = np.array([x_min, x_max, x_max, x_min])
    ground_y = np.array([y_min, y_min, y_max, y_max])
    ground_z = np.array([0.0, 0.0, 0.0, 0.0])
    verts = [list(zip(ground_x, ground_y, ground_z))]
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # type: ignore
    ground = Poly3DCollection(verts, alpha=0.25, facecolor="green", edgecolor="green")  # type: ignore
    ax.add_collection3d(ground)  # type: ignore

    #Trail line that shows where the rocket has been
    trail_line, = ax.plot([], [], [], "b", linewidth=0.8, alpha=0.5)  # type: ignore

    #Lists to accumulate trail points
    trail_x = []  # type: list
    trail_y = []  # type: list
    trail_z = []  # type: list

    # We'll draw the rocket as an annotation arrow in SCREEN (pixel) space.
    # This means no matter which way the rocket points in 3D, the arrow is
    # always the same pixel length on screen — no perspective foreshortening.
    ARROW_PX = 40  # arrow length in pixels — tweak to taste

    # Mutable container for the annotation so we can remove/redraw it
    rocket_annot = [None]

    def _project_3d_to_2d(point_3d):
        """Convert a 3D data point to 2D display (pixel) coordinates."""
        x2d, y2d, _ = proj3d.proj_transform(
            point_3d[0], point_3d[1], point_3d[2], ax.get_proj()
        )
        return ax.transData.transform((x2d, y2d))

    def init():
        # Set axis ranges and labels once
        ax.set_xlim(x_min, x_max)  # type: ignore
        ax.set_ylim(y_min, y_max)  # type: ignore
        ax.set_zlim(z_min, z_max)  # type: ignore
        ax.set_xlabel("X")  # type: ignore
        ax.set_ylabel("Y")  # type: ignore
        ax.set_zlabel("Z (up)")  # type: ignore
        ax.set_title("TVC Rocket Simulation")  # type: ignore
        return (trail_line,)

    def update(frame):
        pos = np.array(positions[frame])
        d = directions[frame]

        # Project position and a point slightly ahead into screen space
        screen_pos = _project_3d_to_2d(pos)
        screen_ahead = _project_3d_to_2d(pos + d * arrow_length)

        # Get the 2D direction on screen, then normalise to fixed pixel length
        screen_dir = screen_ahead - screen_pos
        screen_norm = np.linalg.norm(screen_dir)
        if screen_norm > 1e-6:
            screen_dir = screen_dir / screen_norm
        else:
            screen_dir = np.array([0.0, 1.0])

        # Tail and tip in pixel coordinates, centred on the projected position
        tail_px = screen_pos - screen_dir * (ARROW_PX / 2)
        tip_px  = screen_pos + screen_dir * (ARROW_PX / 2)

        # Remove previous annotation
        if rocket_annot[0] is not None:
            rocket_annot[0].remove()

        # Draw a fixed-pixel-size arrow using annotate in display coords
        rocket_annot[0] = ax.annotate( #type: ignore
            "",
            xy=tip_px,#type: ignore
            xytext=tail_px,#type: ignore
            xycoords="figure pixels",
            textcoords="figure pixels",
            arrowprops=dict(
                arrowstyle="-|>",
                color="purple",
                lw=2.5,
                mutation_scale=15,
            ),
        )

        # Append current position to the trail
        trail_x.append(pos[0])
        trail_y.append(pos[1])
        trail_z.append(pos[2])
        trail_line.set_data(trail_x, trail_y)  # type: ignore
        trail_line.set_3d_properties(trail_z)  # type: ignore

        # Show elapsed simulation time in the title
        ax.set_title(f"TVC Rocket  |  t = {frame * skip * 0.001:.3f}s")  # type: ignore
        return (trail_line,)

    # Build and display the animation
    anim = FuncAnimation(  # type: ignore
        fig,
        update,
        frames=n_frames,
        init_func=init,
        interval=1,
        blit=False,
    )

    plt.tight_layout()  # type: ignore
    plt.show()  # type: ignore


if __name__ == "__main__":
    data = load_data("sim_data.json")
    animate(data, skip=10)