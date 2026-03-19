import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
import json 
import numpy as np 


with open('sim_data.json', 'r') as f:
    data = json.load(f) 

#interval frames logged
step: int = 8

#lists for postion
positions: list[list[float]] = data['r'] 
x: list[float] = []
y: list[float] = []
z: list[float] = []
for p in positions:
    x.append(p[0]) 
    y.append(p[1]) 
    z.append(p[2]) 

#list of orintation
psi = data['psi'] 

#get direction of vector
def get_direction(theta, phi):
    # Rotate around X by theta, then around Y by phi
    dx = np.sin(phi)
    dy = -np.sin(theta) * np.cos(phi)
    dz = np.cos(theta) * np.cos(phi)
    return dx, dy, dz

#plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d') 

ax.set_title('Flight Path') 
ax.set_xlabel('x[m]') 
ax.set_ylabel('y[m]')
ax.set_zlabel('z[m]')

#set limits for axes
pad = 0.1
ax.set_xlim(min(x) - pad, max(x) + pad)
ax.set_ylim(min(y) - pad, max(y) + pad)
ax.set_zlim(min(z) - pad, max(z) + pad)

#create cg point, that gets traced out
point, = ax.plot([], [], [], 'ro', markersize=6)#type: ignore

#create arrow that points from cg
arrow = ax.quiver(0, 0, 0, 0, 0, 0, color='blue', length=2)  # type: ignore

trail, = ax.plot([], [], [], 'g--', linewidth=1, alpha=0.5) 
arrow_size: float = max(z) * 0.1 
def update(frame):
    point.set_data_3d([x[frame]], [y[frame]], [z[frame]]) #type: ignore
    
    #feed trail the postion data from start up to current frame 
    trail.set_data_3d(x[:frame+1], y[:frame+1], z[:frame+1]) #type: ignore

    global arrow
    arrow.remove()
    
    theta = psi[frame][0] 
    phi = psi[frame][1] 
    
    dx, dy, dz = get_direction(theta, phi)
    
    arrow = ax.quiver(x[frame], y[frame], z[frame], dx, dy, dz,
                       color='blue', length=arrow_size, normalize=True)#type: ignore 
    
    return point, arrow


ani = FuncAnimation(fig, update, frames=range(0, len(x), step), interval=1)
plt.show()
