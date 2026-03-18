import numpy as np  
import random, math, json


#constants
mass: float = 0.81
g: float = 9.81

M_arm_thrust_b: tuple[float, float, float] = (0.0, 0.0, -0.2) 
radius: float = 0.0762     
length: float = 0.635  

sim_time: float = 15.5
burn_time: float = 3.2
dt: float  = 0.0001 

#** Calculate accurate value ** 
moment_inertia: float = (1/3)*mass*length**2 + (1/4)*mass*radius** 2

def f_thrust_b(alpha: float, beta: float, t: float) -> tuple[float, float, float]:
    f_thrust_mag: float = random.uniform(13.5, 14.5) 
    if(t <= burn_time):
        return (f_thrust_mag*math.sin(alpha), f_thrust_mag*math.sin(beta), f_thrust_mag*math.cos(alpha)*math.cos(beta)) 
    else:
        return (0.0, 0.0, 0.0)

#takes in rocket orintation, thrust vector in rocket frame, anand returns the thrust vector in the world frame
def f_thrust_w(psi: list[float], f_thrust_b: tuple[float, float, float]) -> tuple[float, float, float]: 
    F_thrust_b_x = f_thrust_b[0] 
    F_thrust_b_y = f_thrust_b[1]
    F_thrust_b_z = f_thrust_b[2]
    
    theta = psi[0] 
    phi = psi[1]
    
    return (F_thrust_b_x*math.cos(phi)+F_thrust_b_y*math.sin(theta)*math.sin(phi)+F_thrust_b_z*math.cos(theta)*math.sin(phi), 
            F_thrust_b_y*math.cos(theta)-F_thrust_b_z*math.sin(theta), 
            -F_thrust_b_x*math.sin(phi)+F_thrust_b_y*math.sin(theta)*math.cos(phi)+F_thrust_b_z*math.cos(theta)*math.cos(phi)) 


#takes in the rocket orintation, and transforms a vector r in rocket frame, to a vector r in the world frame
def vec_b_to_w(psi: list[float], vec_b: tuple) -> tuple:
    x = vec_b[0] 
    y = vec_b[1] 
    z = vec_b[2]
    
    theta = psi[0]
    phi = psi[1]
    
    return (math.cos(phi)*x  + math.sin(phi)*math.sin(theta)*y+ math.sin(phi)*math.cos(theta)*z, 
            math.cos(theta)*y - math.sin(theta)*z, 
            -math.sin(phi)*x  + math.cos(phi) *math.sin(theta)*y  + math.cos(phi) *math.cos(theta)*z) 

def vec_to_pure_q(v: list) -> list:
    return [0, v[0], v[1], v[2]]

def normalize_q(q: list) -> list:
    mag_q: float = math.sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2) 
    unit_q: list = [q[0]/mag_q, q[1]/mag_q, q[2]/mag_q, q[3]/mag_q] 

    return unit_q

def multiply_q_p(q: list, p: list) -> list:
    v: list = [0.0, 0.0, 0.0, 0.0]  

    v[0] = (q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3])  
    v[1] = (q[0]*p[1] + q[1]*p[0] + q[2]*p[3] - q[3]*p[2]) 
    v[2] = (q[0]*p[2] - q[1]*p[3] + q[2]*p[0] + q[3]*p[1]) 
    v[3] = (q[0]*p[3] + q[1]*p[2] - q[2]*p[1] - q[3]*p[0])  
    
    return v

def main() -> None:
    #init logs 
    log_r = []
    log_psi = []
    log_interval = 10
    step = 0

    #init postion and its derivatives
    r: list[float] = [0.0, 0.0, 0.0]
    v: list[float] = [0.0, 0.0, 0.0]
    a: list[float] = [0.0, 0.0, 0.0] 

    #init rotation 
    q: list[float] = [1.0, 0.0, 0.0, 0.0] 
    q_dot: list[float] = [0.0, 0.0, 0.0, 0.0] 
    
    theta: float = 0.0
    phi: float = 0.0 
    psi: list[float] = [theta, phi]  

    #angular velocity and acceleration
    omega_b: list[float] = [0.0, 0.0, 0.0] 
    alpha_b: list[float]= [0.0, 0.0, 0.0] 
    omega_b_q: list[float] = [0.0, 0.0, 0.0, 0.0] 

    #init force and torque
    F_w: list[float] = [0.0, 0.0, 0.0] 
    torque_b = [] 

    alpha: float = 0.1
    beta: float = 0
    
    for t in np.arange(0, sim_time, dt):
        #thrust vector of rocket in body frame to world frame
        F_thrust_b = f_thrust_b(math.radians(alpha), math.radians(beta), t) 
        F_thrust_w = f_thrust_w(psi, F_thrust_b)  

        #net force on rocket, in world frame
        F_w = [F_thrust_w[0], F_thrust_w[1], F_thrust_w[2]-mass*g]

        a[0] = F_w[0] / mass
        a[1] = F_w[1] / mass 
        a[2] = F_w[2] / mass

        v[0] += dt*a[0]  
        v[1] += dt*a[1] 
        v[2] += dt*a[2]  
        
        r[0] += dt*v[0]   
        r[1] += dt*v[1] 
        r[2] += dt*v[2]

        #compute torque on rocket, due to thrust in rocket frame
        torque_b = np.cross(M_arm_thrust_b, F_thrust_b)

        alpha_b[0] = torque_b[0] / moment_inertia
        alpha_b[1] = torque_b[1] / moment_inertia

        omega_b[0] += dt*alpha_b[0]
        omega_b[1] += dt*alpha_b[1] 

        #compute rate of change of quaternion
        omega_b_q = vec_to_pure_q(omega_b)
        q_dot = multiply_q_p(q, omega_b_q)
        q_dot[0] = 1/2*q_dot[0] 
        q_dot[1] = 1/2*q_dot[1] 
        q_dot[2] = 1/2*q_dot[2] 
        q_dot[3] = 1/2*q_dot[3] 

        #update orintation 
        q[0] += dt*q_dot[0]
        q[1] += dt*q_dot[1]
        q[2] += dt*q_dot[2]
        q[3] += dt*q_dot[3]

        q = normalize_q(q) 
        
        w = q[0] 
        x = q[1] 
        y = q[2] 
        z = q[3]
        theta = math.asin(max(-1.0, min(1.0, 2*(w*x - z*y))))
        phi = math.atan2(2*(w*y + x*z), 1 - 2*(x**2 + y**2))
        psi = [theta, phi]

        #logging stuff
        if step % log_interval == 0:
            log_r.append([r[0], r[1], r[2]]) 
            log_psi.append([psi[0], psi[1]]) 
        step += 1 

    #Log data to json
    print(f"Final position: x={r[0]:.4f}, y={r[1]:.4f}, z={r[2]:.4f}")
    print(f"Final velocity: vx={v[0]:.4f}, vy={v[1]:.4f}, vz={v[2]:.4f}")
    print(f"Final angles:   theta={math.degrees(theta):.4f} deg, phi={math.degrees(phi):.4f} deg")
    data = {"r": log_r, "psi": log_psi} 
    with open("sim_data.json", "w") as f:
        json.dump(data, f) 
    
if __name__ == "__main__":
    main()