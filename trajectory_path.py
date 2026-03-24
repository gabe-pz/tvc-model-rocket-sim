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

#create line that passes through cg
line_fwd, = ax.plot([], [], [], color='black', linewidth=2) 
line_bwd, = ax.plot([], [], [], color='green', linewidth=2) 
line_length = max(z) * 0.15

#trajectory line
trail, = ax.plot([], [], [], 'b--', linewidth=1, alpha=0.5) 

def update(frame) -> tuple:
    #get current point of cg
    px = x[frame]
    py = y[frame] 
    pz = z[frame] 

    point.set_data_3d([px], [py], [pz]) #type: ignore
    
    #feed trail the postion data from start up to current frame 
    trail.set_data_3d(x[:frame+1], y[:frame+1], z[:frame+1]) #type: ignore
    
    theta = psi[frame][0] 
    phi = psi[frame][1] 
    
    #create unit vector for current point
    (dx, dy, dz) = get_direction(theta, phi) 
    mag = np.sqrt(dx**2 + dy**2 + dz**2)
    dx = dx / mag
    dy = dy / mag
    dz = dz / mag

    line_fwd.set_data_3d([px, px + dx*line_length], [py, py + dy*line_length], [pz, pz + dz*line_length]) #type: ignore
    line_bwd.set_data_3d([px, px - dx*line_length], [py, py - dy*line_length], [pz, pz - dz*line_length]) #type: ignore
    
    return (point, trail, line_fwd, line_bwd)

ani = FuncAnimation(fig, update, frames=range(0, len(x), step), interval=1)
plt.show()
