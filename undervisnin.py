import  numpy as np
import matplotlib.pyplot as plt 

d_mean = 3.84e5 # 3.84*10**5 , km
t = 3*24 # hour 

v = d_mean/t

print(f"The velocity is {v} km/hour")

# We want to test for 2 days, 1 day, 12 hours, and 9 hours

t1 = 2*24 # hours
t2 = 24 
t3 = 12
t4 = 9

times = [t1, t2,t3,t4] # we make a list for input values
velocities = [] #  we make a list for outpur values

for time in times:
    v  = d_mean / time
    velocities.append(v)
    print(f"{v} km/hour")



plt.plot(times,velocities, label = "velocities for specific time")
plt.xlabel("t [km]")
plt.ylabel("v [km/hour]")
plt.grid()
plt.legend()
plt.show()