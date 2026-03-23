import numpy as np  
import math, json, os, random 
import rocket_math as rm

#clear file before new run
if os.path.isfile("sim_data.json"):
    os.remove("sim_data.json")

#rocket constants
mass: float = 1.01074445935 
distance_to_thrust_vector: float = 0.6477
center_of_pressure: float = 0.08775954
center_of_gravity: float = 0.4059174

M_arm_thrust_b: list[float] = [0.0, 0.0, center_of_gravity-distance_to_thrust_vector] 
M_arm_aero_force_b: list[float] = [0.0, 0.0, center_of_gravity-center_of_pressure] 
I_xx: float = 0.0249899588
I_yy: float = 0.0249868814

#sim constants
sim_time: float = 20
burn_time: float = 3.5
dt: float  = 0.0001 
log_interval: int = 10

#physical constants
g: float = 9.81
drag_coef: float = 0.291
reference_area: float = 0.00456
rho: float = 1.187 

#wind 
u: float = 7.0 
I_u: float = 0.08 
sigma_u: float = I_u*u 
alpha_coef: float = 5.0/3.0 

a_constants: list[float] = [0.0, 0.0, 0.0]
a_constants[0] = 1.0

for k in range(1, 3): 
    a_constants[k] = (k-1-alpha_coef/2.0)*a_constants[k-1]/k

pink_noise_std: float = 2.252

x_prev_1: float = 0.0 
x_prev_2: float = 0.0

t_total: float = sim_time        # Simulate 20 seconds
dt_turb: float = 0.05        # One sample every 0.05s (20 Hz)
n_samples = int(t_total / dt_turb) + 1

wind = []

for n in range(n_samples):
    w_n = np.random.randn()

    x_n = w_n - a_constants[1] * x_prev_1 - a_constants[2] * x_prev_2
    x_normalized = x_n / pink_noise_std

    U_n = u + sigma_u * x_normalized

    wind.append(U_n)
    x_prev_2 = x_prev_1
    x_prev_1 = x_n

turb_times = np.arange(0, t_total + dt_turb, dt_turb)[:n_samples]
sim_times = np.arange(0, t_total, dt) 
wind_at_sim = np.interp(sim_times, turb_times, wind) 

#thrust vector function
def f_thrust_b(alpha: float, beta: float, t: float) ->list:
    f_thrust_mag: float = 14.4 
    if(t <= burn_time):
        return [f_thrust_mag*math.sin(alpha), f_thrust_mag*math.sin(beta), f_thrust_mag*math.cos(alpha)*math.cos(beta)]
    else:
        return [0.0, 0.0, 0.0]


def main() -> None:
    #init log stuff
    log_r: list[list[float]] = []
    log_psi: list[list[float]] = []
    step: float = 0

    #init postion and its derivatives
    r: list[float] = [0.0, 0.0, 0.0]
    v: list[float] = [0.0, 0.0, 0.0]
    a: list[float] = [0.0, 0.0, 0.0] 

    #init rotation stuff
    q: list[float] = [1.0, 0.0, 0.0, 0.0] 
    q_dot: list[float] = [0.0, 0.0, 0.0, 0.0] 
    psi: list[float] = [0.0, 0.0] # psi = (theta, phi),  

    #init angular acceleration and velocity  
    omega_b: list[float] = [0.0, 0.0, 0.0] 
    alpha_b: list[float]= [0.0, 0.0, 0.0] 
    omega_b_q: list[float] = [0.0, 0.0, 0.0, 0.0] 

    #init force and torque
    F_w: list[float] = [0.0, 0.0, 0.0] 
    torque_b = [] 

    #inital start values for angle of tvc
    alpha: float = 0.0
    beta: float = 0.0

    wind_speed = 0.0
    wind_direction = [1.0, 0.25, 0.3] 
    wind_vector: list[float] = [wind_speed*wind_direction[0], wind_speed*wind_direction[1], wind_speed*wind_direction[2]] 
    v_rel_fluid: list[float] = [0.0, 0.0, 0.0] 
    F_d_w: list[float] = [0.0, 0.0, 0.0] 
    F_d_b: list[float] = [0.0, 0.0, 0.0] 
    v_rel_mag: float = 0.0

    for i, t in enumerate(sim_times):
        #1. Force and Position

        wind_speed = wind_at_sim[i]
        wind_vector[0] = wind_speed*wind_direction[0] 
        wind_vector[1] = wind_speed*wind_direction[1] 
        wind_vector[2] = wind_speed*wind_direction[2] 
        

        v_rel_fluid[0] = v[0] - wind_vector[0] 
        v_rel_fluid[1] = v[1] - wind_vector[1] 
        v_rel_fluid[2] = v[2] - wind_vector[2] 
        v_rel_mag = math.sqrt(v_rel_fluid[0]**2+v_rel_fluid[1]**2+v_rel_fluid[2]**2) 

        F_d_w[0] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid[0]*v_rel_mag)
        F_d_w[1] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid[1]*v_rel_mag)
        F_d_w[2] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid[2]*v_rel_mag) 

        F_d_b = rm.rotate_v_b(q, F_d_w)
        
        #Force of thrust, from body to world
        F_thrust_b = f_thrust_b(math.radians(alpha), math.radians(beta), t) 
        F_thrust_w = rm.rotate_v_w(q, F_thrust_b)   
        
        #Sum of forces in world
        F_w = [F_thrust_w[0]+F_d_w[0], F_thrust_w[1]+F_d_w[1], F_thrust_w[2]-mass*g+F_d_w[2]] 

        
        #Compute accleration
        a[0] = F_w[0] / mass
        a[1] = F_w[1] / mass 
        a[2] = F_w[2] / mass

        #Numerical integration for v and r
        v[0] += dt*a[0]  
        v[1] += dt*a[1] 
        v[2] += dt*a[2]  
        r[0] += dt*v[0]   
        r[1] += dt*v[1] 
        r[2] += dt*v[2]

        #2. Torque and Rotation 

        #Torque on the rocket
        torque_thrust_b = np.cross(M_arm_thrust_b, F_thrust_b) 
        torque_d_b = np.cross(M_arm_aero_force_b, F_d_b) 

        #sum of torque on rocket
        torque_b = [torque_thrust_b[0]+torque_d_b[0], torque_thrust_b[1]+torque_d_b[1], torque_thrust_b[2]+torque_d_b[2]]  

        #Compute angular acceleration
        alpha_b[0] = torque_b[0] / I_xx
        alpha_b[1] = torque_b[1] / I_yy

        #Numerical integration for omega
        omega_b[0] += dt*alpha_b[0]
        omega_b[1] += dt*alpha_b[1] 

        #Compute rate of change of quaternion
        omega_b_q = rm.vec_to_pure_q(omega_b)
        q_dot = rm.multiply_q_p(q, omega_b_q)
        q_dot[0] = 1/2*q_dot[0] 
        q_dot[1] = 1/2*q_dot[1] 
        q_dot[2] = 1/2*q_dot[2] 
        q_dot[3] = 1/2*q_dot[3] 

        #Numerical integraion for q
        q[0] += dt*q_dot[0]
        q[1] += dt*q_dot[1]
        q[2] += dt*q_dot[2]
        q[3] += dt*q_dot[3]

        #Renormalize q to account for numerical drift
        q = rm.normalize_q(q) 
        
        #Convert quaternion to euler angles 
        psi = rm.q_to_euler(q) 

        #Check if reach ground
        if r[2] <= 0 and t > 0.1: 
            break

        #Log data 
        if step % log_interval == 0:
            log_r.append([r[0], r[1], r[2]]) 
            log_psi.append([psi[0], psi[1]]) 
        step += 1 

    #Final print log
    print(f"Final position: x={r[0]}, y={r[1]}, z={r[2]}")  
    print(f"Final velocity: vx={v[0]}, vy={v[1]}, vz={v[2]}")
    print(f"Final angles:   theta={math.degrees(psi[0])} deg, phi={math.degrees(psi[1])} deg")    

    #Log flight data to json
    data = {"r": log_r, "psi": log_psi} 
    with open("sim_data.json", "w") as f:
        json.dump(data, f) 
    
if __name__ == "__main__":
    main()