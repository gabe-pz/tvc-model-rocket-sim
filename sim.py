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
norm_coef: float = 2.0

reference_area: float = 0.00456
rho: float = 1.187 

#wind generation
u: float = 5.0 
I_u: float = 0.05
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
dt_turb: float = 0.05            # One sample every 0.05s (20 Hz)
n_samples = int(t_total / dt_turb) + 1

wind = []

for _ in range(n_samples):
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
#end wind generation 

#thrust vector function
def f_thrust_b(alpha: float, beta: float, t: float) ->list:
    f_thrust_mag: float = 15
    if(t <= burn_time):
        return [f_thrust_mag*math.sin(alpha), f_thrust_mag*math.sin(beta), f_thrust_mag*math.cos(alpha)*math.cos(beta)]
    else:
        return [0.0, 0.0, 0.0]

def clamp(value: float, min_val: float, max_val: float) -> float:
    return min(max(value, min_val), max_val) 

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

    #init wind 
    wind_speed: float = 0.0
    wind_direction: list[float] = [1.0/math.sqrt(2), 1.0/math.sqrt(2), 0.0]   
    wind_velocity: list[float] = [wind_speed*wind_direction[0], wind_speed*wind_direction[1], wind_speed*wind_direction[2]] 

    #init drag 
    v_rel_fluid_w: list[float] = [0.0, 0.0, 0.0] 
    v_rel_fluid_b: list[float] = [0.0, 0.0, 0.0] 
    v_rel_fluid_mag: float = 0.0

    alpha_x_aoa: float = 0.0 
    alpha_y_aoa: float = 0.0

    F_drag_w: list[float] = [0.0, 0.0, 0.0] 
    F_drag_b: list[float] = [0.0, 0.0, 0.0] 
    
    F_norm_w: list[float] = [0.0, 0.0, 0.0]
    F_norm_b: list[float] = [0.0, 0.0, 0.0]

    error_x: float = 0.0
    error_y: float = 0.0 

    kp_x: float = 0.75
    kp_y: float = 0.75 
    for i, t in enumerate(sim_times):
        #1. Force and Position
        
        error_x = 0 - math.degrees(psi[0]) 
        error_y = 0 - math.degrees(psi[1]) 

        alpha = clamp((kp_x*error_x), -5, 5)
        beta = clamp((kp_y*error_y), -5, 5) 
        #determine wind
        wind_speed = wind_at_sim[i] 
        wind_velocity[0] = wind_speed*wind_direction[0] 
        wind_velocity[1] = wind_speed*wind_direction[1] 
        wind_velocity[2] = wind_speed*wind_direction[2] 
        
        #determine veloctiy of rocket relative to fluid in world frame and then body frame
        v_rel_fluid_w[0] = v[0] - wind_velocity[0] 
        v_rel_fluid_w[1] = v[1] - wind_velocity[1] 
        v_rel_fluid_w[2] = v[2] - wind_velocity[2] 
        v_rel_fluid_b = rm.rotate_v_b(q, v_rel_fluid_w) 

        v_rel_fluid_mag = math.sqrt(v_rel_fluid_w[0]**2+v_rel_fluid_w[1]**2+v_rel_fluid_w[2]**2) 

        #determine the angles of attacks
        alpha_x_aoa = math.atan2(v_rel_fluid_b[0], v_rel_fluid_b[2])  
        alpha_y_aoa  = math.atan2(v_rel_fluid_b[1], v_rel_fluid_b[2]) 

        #determine drag force in world frame 
        F_drag_w[0] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid_mag*v_rel_fluid_w[0]) 
        F_drag_w[1] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid_mag*v_rel_fluid_w[1]) 
        F_drag_w[2] = -0.5*rho*drag_coef*reference_area*(v_rel_fluid_mag*v_rel_fluid_w[2]) 
        #convert drag force to body frame
        F_drag_b = rm.rotate_v_b(q, F_drag_w) 
        

        #determine normal force in body frame
        F_norm_b[0] = 0.125*rho*norm_coef*reference_area*(v_rel_fluid_mag**2)*alpha_x_aoa 
        F_norm_b[1] = 0.125*rho*norm_coef*reference_area*(v_rel_fluid_mag**2)*alpha_y_aoa 
        #convert normal force to world frame
        F_norm_w = rm.rotate_v_w(q, F_norm_b) 

        #Force of thrust, from body to world
        F_thrust_b = f_thrust_b(math.radians(alpha), math.radians(beta), t) 
        F_thrust_w = rm.rotate_v_w(q, F_thrust_b)   
        
        #Sum of forces in world
        F_w = [F_thrust_w[0]+F_drag_w[0]+F_norm_w[0], F_thrust_w[1]+F_drag_w[1]+F_norm_w[1], F_thrust_w[2]-mass*g+F_drag_w[2]+F_norm_w[2]]  

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
        torque_drag_b = np.cross(M_arm_aero_force_b, F_drag_b) 
        torque_norm_b = np.cross(M_arm_aero_force_b, F_norm_b) 


        #sum of torque on rocket
        torque_b = [torque_thrust_b[0]+torque_drag_b[0]+torque_norm_b[0], torque_thrust_b[1]+torque_drag_b[1]+torque_norm_b[1], torque_thrust_b[2]+torque_drag_b[2]+torque_norm_b[2]]   

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
        if r[2] <= 0 and t > 0.5: 
            print("*****FLIGHT STATS*****")
            print(f"Total flight time: {t:.2f}s")
            print(f"Max altitude: {(max(log_r, key=lambda p: p[2])[2]):.2f}m")
            
            break

        #Log data 
        if step % log_interval == 0:
            log_r.append([r[0], r[1], r[2]]) 
            log_psi.append([psi[0], psi[1]]) 
        step += 1 

    #Final print log
    print()
    print("*****Final r and ψ*****") 
    print(f"Final r: x={r[0]:.2f}, y={r[1]:.2f}, z={r[2]:.2f}")  
    print(f"Final ψ:   theta={math.degrees(psi[0]):.2f} deg, phi={math.degrees(psi[1]):.2f} deg")    

    #Log flight data to json
    data = {"r": log_r, "psi": log_psi} 
    with open("sim_data.json", "w") as f:
        json.dump(data, f) 
    
if __name__ == "__main__":
    main()