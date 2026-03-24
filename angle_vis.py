import matplotlib.pyplot as plt
import json, math


with open('sim_data.json', 'r') as f:
    data = json.load(f) 


psi = data['psi']
theta: list[float] = []
phi: list[float] = [] 
for p in psi:
    theta.append(math.degrees(p[0])) 
    phi.append(math.degrees(p[1])) 

#log time
dt = 0.0001 
log_interval = 10
effective_dt = log_interval*dt

time: list[float] = []
for i in range(len(theta)):
    time.append(i * effective_dt)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

ax1.plot(time, theta, 'b-', label='Theta (actual)')
ax1.axhline(y=0, color='r', linestyle='--', label='Setpoint')
ax1.set_ylabel('Theta (deg)')
ax1.legend()
ax1.grid(True)

ax2.plot(time, phi, 'b-', label='Phi (actual)')
ax2.axhline(y=0, color='r', linestyle='--', label='Setpoint')
ax2.set_ylabel('Phi (deg)')
ax2.set_xlabel('Time (s)')
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.show()