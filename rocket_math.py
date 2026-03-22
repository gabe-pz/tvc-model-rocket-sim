import math

#takes in the rocket orintation, and transforms a vector r in rocket frame, to a vector r in the world frame
def rotate_v_w(q: list[float], v: list[float]) -> list[float]:
    v_q: list[float] = vec_to_pure_q(v) 
    v_1: list[float] = multiply_q_p(q, v_q) 
    q_conjugate: list[float] = conjugate_q(q) 
    v_world_q: list[float] = multiply_q_p(v_1, q_conjugate) 
    v_world: list[float] = [v_world_q[1], v_world_q[2], v_world_q[3]]

    return v_world

#takes in the rocket orintation, and transforms a vector r in world frame, to a vector r in the rocket frame
def rotate_v_b(q: list[float], v: list[float]) -> list[float]:
    v_q: list[float] = vec_to_pure_q(v) 
    q_conjugate: list[float] = conjugate_q(q) 
    v_1: list[float] = multiply_q_p(q_conjugate, v_q)
    v_body_q: list[float] = multiply_q_p(v_1, q) 
    v_body: list[float] = [v_body_q[1], v_body_q[2], v_body_q[3]]

    return v_body
def vec_to_pure_q(v: list) -> list:
    return [0, v[0], v[1], v[2]]

def normalize_q(q: list) -> list:
    mag_q: float = math.sqrt(q[0]**2+q[1]**2+q[2]**2+q[3]**2) 
    unit_q: list = [q[0]/mag_q, q[1]/mag_q, q[2]/mag_q, q[3]/mag_q] 

    return unit_q

def conjugate_q(q: list) -> list:
    return [q[0], -q[1], -q[2], -q[3]] 

def multiply_q_p(q: list, p: list) -> list:
    v: list = [0.0, 0.0, 0.0, 0.0]  
    v[0] = (q[0]*p[0] - q[1]*p[1] - q[2]*p[2] - q[3]*p[3])  
    v[1] = (q[0]*p[1] + q[1]*p[0] + q[2]*p[3] - q[3]*p[2]) 
    v[2] = (q[0]*p[2] - q[1]*p[3] + q[2]*p[0] + q[3]*p[1]) 
    v[3] = (q[0]*p[3] + q[1]*p[2] - q[2]*p[1] + q[3]*p[0])  

    return v

def q_to_euler(q: list) -> list: 
    theta = math.atan2(2*(q[0]*q[1]+q[2]*q[3]), 1-2*(q[1]**2+q[2]**2))
    phi = math.asin(2*(q[0]*q[2] - q[1]*q[3]))
    return [theta, phi]