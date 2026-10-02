import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, odeint
np.set_printoptions(threshold=np.inf)
#I put this line to see the change.
#we define our function for integrating
def fun(t,y,potential):
   
    x1 = y[0]
    x2 = y[1]
    x3 = y[2]
    lamb = y[3]
    
    
    if potential == "powerlaw":
        gamma = 2
    else:
        gamma = 1
    dx1dN = -3*x1 + np.sqrt(6)/2 * lamb*x2**2 + 1/2 * x1 *(3+3*x1**2-3*x2**2+x3**2)
    dx2dN = -np.sqrt(6)/2 * lamb*x1*x2 + 1/2 * x2 *(3+3*x1**2-3*x2**2+x3**2)
    dx3dN = -2*x3 + 1/2 * x3 *(3+3*x1**2-3*x2**2+x3**2)
    dlambdN = -np.sqrt(6)*lamb**2*(gamma-1)*x1
    solution = [dx1dN,dx2dN,dx3dN,dlambdN]
    return solution




#we define our start values
z_start = 2e7
z_end = 0
N_start = np.log(1/(1+z_start))
N_end = np.log(1/(1+z_end))


#we define our values for the two potentials
y0 = [5e-5, 1e-8, 0.9999, 1e9]
y00 = [0,5e-13,0.9999,3/2]

N = np.linspace(N_start,N_end,1000)

#by changing to another string and by changing the y0 to y00  in sol, we make plots and calculations for the two potentials.
potential = "powerlaw"
sol = solve_ivp(lambda N, y : fun(N,y,potential), [N_start, N_end], y0,t_eval =N , rtol = 1e-8, atol = 1e-8)

x1, x2, x3, lamb = sol.y 
N = sol.t 

# we define the densities
Ω_fi = x1**2 + x2**2
Ω_r = x3**2 
Ω_m = 1-x1**2-x2**2-x3**2


z = np.exp(-N)-1

#we define the omega that is important for further calculations
ω = (0.5*(x1*np.sqrt(6))**2 - x2**2 *3/1)/ ((x1*np.sqrt(6))**2+ x2**2 *3/1)






#we plot the densities with respect to the redshift
plt.plot(z,Ω_fi,label ="Ω_fi",color = "blue")
plt.plot(z,Ω_r, label ="Ω_r",color = "green")
plt.plot(z,Ω_m,label ="Ω_m", color = "red")
plt.plot(z,ω , label = "ω_ϕ")
plt.xscale("log")
plt.legend()
plt.title("The density parameters  and The EoS parameter ")
plt.xlabel("redshift z")
plt.ylabel("density")
plt.grid()
plt.show()



#Here we aim to plot the Hubble parameter
H_over_H_0 =np.sqrt( Ω_m[-1]*(1+z)**3 + Ω_r[-1]*(1+z)**4 + Ω_fi[-1]* np.exp(3*(1+ω)))
H_cdm = np.sqrt(0.3*(1+z)**3 + (1-0.3))

plt.plot(z,H_cdm,color = "cyan",label = "ΛCDM")
plt.plot(z,H_over_H_0, color = "purple", label = "Model exponential")
plt.ylabel("H/H_0")
plt.xlabel("z")
plt.xscale("log")
plt.yscale("log")
plt.title("Hubble parameter ")
plt.legend()
plt.grid()
plt.show()

#we calculate the Hubble parameter
H_0  = 1000 *100*0.72 / (1*3.08e22)
c = 299792458
t_0 = 1/(H_0*(1+z)*(np.sqrt( Ω_m[-1]*(1+z)**3 + Ω_r[-1]*(1+z)**4 + Ω_fi[-1]* np.exp(3*(1+ω)))**(1/2))) / 31556926/1e9
print(np.max(t_0))

t_a = 1/ (H_0*(1+z)*np.sqrt(0.3*(1+z)**3 + (1-0.3))**(1/2))  /31556926/1e9
print(np.max(t_a))

dimension_less =  1/((1+z)*(np.sqrt( Ω_m[-1]*(1+z)**3 + Ω_r[-1]*(1+z)**4 + Ω_fi[-1]* np.exp(3*(1+ω)))**(1/2)))
dimension_less1 =  1/ ((1+z)*np.sqrt(0.3*(1+z)**3 + (1-0.3))**(1/2))
print(np.max(dimension_less), np.max(dimension_less1))

#we calculate the dimensionless luminosity distance and plot it
#we define a new z variable
z = np.linspace(0,2,len(H_over_H_0))
distance = (1+z)/H_over_H_0
plt.plot(z,distance,color = "black")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("redshift z")
plt.ylabel("luminosity distance")
plt.title("Luminosity distance exponential potential")
plt.grid()
plt.show()

print(np.max(distance))


