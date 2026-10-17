import matplotlib.pyplot as plt 
import numpy as np
"""We aim to calculate the force between the Earth and the Sun"""
 # Define the global variables 
G = 6.6743e-11 # m**3/kg*s**2
M = 1.99e30 # mass of Sun in kg 
m = 5.97e24# mass of Earth in kg 



# Define the force function

def F_vector(d, theta):
    F_magnitude = G*M*m / d**2
    Fx = -F_magnitude*np.cos(theta)
    Fy = -F_magnitude*np.sin(theta)
    return Fx, Fy

# Defining distance and angles

mean_distance = 149.6e9
angles = np.linspace(0*np.pi, 2*np.pi,10)

# We send input in the function

Fy_list = []
angle_list = []

for angle in angles:
    Fx, Fy = F_vector(mean_distance, angle)
    Fy_list.append(Fy)
    angle_list.append(angle)
    print(f"theta = {angle:.2f} rad: Fx = {Fx:.2e} N Fy = {Fy:.2e} N")


plt.plot(angle_list,Fy_list,)
plt.xlabel("angle / radians")
plt.ylabel("force / N")
plt.grid()
plt.show()






































