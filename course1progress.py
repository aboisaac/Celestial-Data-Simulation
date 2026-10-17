import numpy as np 
import matplotlib.pyplot as plt

# A spaceship travels to Moon in a linear motion
# With what speed could the spaceship reach the Moon in 3 days
# Velocity equation v = d/t becomes 

d_mean = 3.84e8 # meter
t = 3 * 24 * 60*60 # second

v = d_mean / t # m/s

print(f"{v} m/s")
print()
print()

# What if the speed is t = 2 days, 1 day, 12 hours or 9 hours?
# We define the t's

t1 = 2* 24* 60*60 # hour
t2 = 1 * 24 * 60*60 
t3 = 12 * 60*60 
t4 = 9 * 60*60

# We make a list

times = [t1,t2,t3,t4]
velocities = [] # Store the output-values

# We make a foor-loop

for time in times:
    v = d_mean / time
    velocities.append(v)
    print(f"{v} m/s")
    

# We plot the values

plt.plot(times,velocities, label = "velocities for specific time")
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")
plt.title("Velocity change with time")
plt.grid()
plt.legend()
plt.show()

# Conlusion: Shorter time-travel requires higher velocities



