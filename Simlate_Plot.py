import Parameters_Structures as ps
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib.animation import FuncAnimation

import matplotlib.pyplot as plt

# Simulate the diffusion for 30 seconds
dt = 0.05
num_steps = int(10 / dt)
my_cell = ps.cell(1.0, 1.0, 1.0, 200, 200)

# Set up the figure and 3D axis
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the ellipsoid
def plot_ellipsoid():
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)

    x = ps.ellipsoid_axis_a * np.outer(np.cos(u), np.sin(v))
    y = ps.ellipsoid_axis_b * np.outer(np.sin(u), np.sin(v))
    z = ps.ellipsoid_axis_c * np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='g', alpha=0.1)

colors = {'A': 'r', 'B': 'b'}

# Function to animate the diffusion
def animate(step):
    ax.cla()
    print(f'Step {step}')
    plot_ellipsoid()

    my_cell.update_status(dt)
    for j in my_cell.mem_molecules:
        x = ps.ellipsoid_axis_a * np.sin(j.theta) * np.cos(j.phi)
        y = ps.ellipsoid_axis_b * np.sin(j.theta) * np.sin(j.phi)
        z = ps.ellipsoid_axis_c * np.cos(j.theta)
        ax.scatter(x, y, z, color=colors[j.type])

ax.set_box_aspect([1,1,1])
ax.view_init(elev=10, azim=100)
# Create the animation
ani = FuncAnimation(fig, animate, frames=num_steps, interval=dt, repeat=False)

# plt.show()
# Save the animation as a video
ani.save('/Users/wuxiaoyu/Documents/GitHub/PolaSim/simulation2.mp4', writer='ffmpeg', fps=20)