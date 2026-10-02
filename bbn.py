# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 17:02:24 2026

@author: salma
"""

# Imports:
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from astropy import constants as const

from scipy.interpolate import CubicSpline


# Own imports:
#from astropy import units 

from reaction_rates import ReactionRates  # Class with equations for each reaction
from new_cosmoproject2 import BackgroundBBN# Class with background equations




           

class BBN:
    def __init__(self, NR_interacting_species=2, Omega_b0=0.05, N_eff=3) -> None:
        """
        Initialize the BBN calculations including particles up to a given number,
        NR_interacting_species. This is a number between 2 and 8, determining
        which particles to be considered in the model. The particle species are
        n, p, D, T, He3, He4, Li7, Be7. The class also initializes the reaction
        rates and background cosmology, with the option to change the baryon
        density and effective number of neutrinos.

        Parameters
        ----------
        NR_interacting_species : int, optional
            The number of particle species to be considered, by default 2

        Keyword Arguments:
            Omega_b0 : float, optional
                The relative baryon density, by default 0.05
            N_eff : float, optional
                The effective number of neutrinos, by default 3
        """
        self.NR_species = NR_interacting_species
        ast_msg = "NR_interacting_species must be between 2 and 8"
        assert self.NR_species <= 8 and self.NR_species >= 2, ast_msg

        # Define the particle species order and their mass numbers:
        self.species_labels = ["n", "p", "D", "T", "He3", "He4", "Li7", "Be7"]
        self.mass_number = [1, 1, 2, 3, 3, 4, 7, 7]

        # Initialize the reaction rates and background cosmology:
        self.RR = ReactionRates()
        self.new_cosmoproject2 = BackgroundBBN(Omega_b0=Omega_b0, N_eff=N_eff)
        
        self.Y_p_list=[]
        self.Y_D_list=[]
        
        self.Y_T =   []
        self.Y_He3 = []
        self.Y_He4 = []
        self.Y_Li7 = []
        self.Y_Be7 = []

        self.Y_p_array = np.array(self.Y_p_list)
        self.Y_D_array = np.array(self.Y_D_list)
        

    def get_ODE(self, lnT, Y):
        """
        Compute the right hand sides of the ODEs for Boltzmanns equation
        for particle species in Y, up to the number of interacting species.

        Parameters
        ----------
        lnT : float
            The natural log of temperature.
        Y : np.ndarray
            Array of relative number densities at current lnT for each species.
            With all 8 included species Y takes the form:
            n, p, D, T, He3, He4, Li7, Be7 = Y

        Returns
        -------
        np.ndarray
            The differential for each species following the shape of Y.
            E.g. dY[0] corresponds to dY_n, dY[1] to dY_p, dY[2] to dY_D, etc.
            Initialized to zero for all species and add to the change for each species.
        """
        # Set current temperature, Hubble parameter and baryon density:
        T = np.exp(lnT)
        T9 = T / 1e9
        Hubble = self.new_cosmoproject2.get_Hubble(T)
        rho_b = self.new_cosmoproject2.get_rho_b(T)

        # Initialize the differential for each species following the shape of Y,
        # i.e. dY[0] corresponds to dY_n. Initialized to zero for all species.

        # Example of how to include each interactions:
        # ~ (n <-> p) (a 1-3)
        # Relative number densities of neutrons and protons:
        dY = np.zeros_like(Y)
        Y_n, Y_p = Y[0], Y[1]
        Y_D  = Y[2]
        
        
        Y_T, Y_He3, Y_He4, Y_Li7, Y_Be7 = Y[3], Y[4], Y[5], Y[6], Y[7]
        rate_n, rate_p = self.RR.get_weak_force(T9)  # Get the weak force rates
        self.Y_p_list.append(Y_p)
        self.Y_D_list.append(Y_T)
        # for the other reactions

        
       
        # The change to the particles on the left hand side of the reaction equation:
        change_LHS = Y_p * rate_p - Y_n * rate_n
        # Update the changes to each side:
            
        dY[0] += change_LHS  # neutron
        dY[1] -= change_LHS  # proton (opposite sign)

        # Include deuterium
        if self.NR_species > 2:
            Y_D = Y[2]  # Relative number density of deuterium
            # (n+p <-> D+gamma) (b.1)
            Y_np = Y_n * Y_p
            rate_np, rate_D = self.RR.get_np_to_D(T9, rho_b)
            change_LHS =  Y_np * rate_np - Y_D * rate_D


            # Update changes to each side:
            dY[0] -= change_LHS  # neutron
            dY[1] -= change_LHS  # proton
            dY[2] += change_LHS  # deuterium (opposite sign)
          
        # Include tritium
        if self.NR_species > 3:
            Y_T = Y[3]  # Relative number density of tritium
            # (n+D <-> T+gamma) (b.3)
            Y_nD = Y_n * Y_D
            rate_nD, rate_T = self.RR.get_nD_to_T(T9, rho_b)
            change_LHS = Y_nD*rate_nD - Y_T* rate_T

            dY[0] -= change_LHS  # neutron
            dY[2] -= change_LHS  # deuterium
            dY[3] += change_LHS  # tritium



            # TWO EQUAL PARTICLES -- MIND THE FACTORS OF 2
            # (D+D <-> p + T) (b.8)
            #equation 8
            Y_DD = Y_D * Y_D
            Y_pT = Y_p * Y_T
            rate_DD, rate_pT = self.RR.get_DD_to_pT(T9, rho_b)
            # Different relative changes to each side:
            change_LHS = 2 * Y_pT * rate_pT - Y_DD * rate_DD
            change_RHS = 0.5 * Y_DD * rate_DD - Y_pT * rate_pT
            dY[2] += change_LHS  #deutrium 
            dY[1] += change_RHS  #proton
            dY[3] += change_RHS  #tritium




            




            
            


        if self.NR_species > 4:  # Include He3
            YpD = Y_p*Y_D
            YH3 = Y_He3
            rate_pD, rate_He3 = self.RR.pD_to_He3(T9,rho_b)
            change_LHS = YpD*rate_pD - YH3*rate_He3

            dY[1] -= change_LHS  # proton
            dY[2] -= change_LHS  # deutrium
            dY[4] += change_LHS  # He3

            YnHe3 = Y_n*Y_He3
            YpT = Y_p*Y_T
            rate_nHe3, rate_pT= self.RR.nHe3_to_pT(T9,rho_b)
            change_LHS = YnHe3*rate_nHe3 - YpT*rate_pT

            dY[1] += change_LHS  #proton
            dY[0] -= change_LHS  # neutron
            dY[3] += change_LHS  # tritium
            dY[4] -= change_LHS  # Helium 3

            #equation 7
            YDD = Y_D*Y_D
            YnHe3 = Y_n*Y_He3
            rate_DD, rate_nHe3= self.RR.get_DD_to_nHe3(T9,rho_b)
            change_LHS = 2*rate_nHe3*YnHe3-rate_DD*YDD
            change_RHS = 0.5*rate_DD*YDD - rate_nHe3*YnHe3

            dY[2] += change_LHS  #deutrium 
            dY[4] += change_RHS  #He3
            dY[0] += change_RHS  #neutron



            


        if self.NR_species > 5:  # Include He4
            YpT = Y_p*Y_T
            YHe4 = Y_He4
            rate_pT, rate_He4= self.RR.get_pT_to_He4(T9,rho_b)
            change_LHS = YpT*rate_pT - YHe4*rate_He4
            dY[1] -= change_LHS  # proton
            dY[3] -= change_LHS  # tritum
            dY[5] += change_LHS  # helium 4

            YnHe3 = Y_n*Y_He3
            YHe4 = Y_He4
            rate_nHe3, rate_He4= self.RR.get_nHe3_to_He4(T9,rho_b)
            change_LHS = rate_nHe3*YnHe3 - rate_He4*YHe4
            dY[0] -= change_LHS  # neutron
            dY[4] -= change_LHS  # helium 3
            dY[5] += change_LHS  #helium 4


            #equation 9
            YDD = Y_D*Y_D
            YHe4 = Y_He4 
            rate_DD, rate_He4= self.RR.get_DD_to_He4(T9,rho_b)
            change_LHS =2*rate_He4*YHe4 - rate_DD*YDD 
            change_RHS = 0.5*rate_DD*YDD - rate_He4*YHe4

            dY[2] += change_LHS  # Deutrium
            dY[5] += change_RHS  # He4


            #equation 10 
            YDHe3 = Y_D*Y_He3 
            YHe4p = Y_He4 *Y_p
            rate_DHe3, rate_He4p = self.RR.get_DHe3_to_He4p(T9,rho_b)
            change_LHS = rate_DHe3*YDHe3 -rate_He4p*YHe4p

            dY[1] += change_LHS # proton
            dY[2] -= change_LHS  # deutrium
            dY[4] -= change_LHS  # helium 3
            dY[5] += change_LHS  #helium 4





            #equation 11
            YDT = Y_D*Y_T
            YHe4n = Y_He4*Y_n
            rate_DT, rate_He4n = self.RR.get_DT_to_He4n(T9,rho_b)
            change_LHS = rate_DT*YDT-rate_He4n*YHe4n
            dY[0] += change_LHS # neutron
            dY[2] -= change_LHS  # deutrium
            dY[3] -= change_LHS  # titrium
            dY[5] += change_LHS  #helium 4


            
            #equation 15
            YHe3T = Y_He3*Y_T
            YHe4D = Y_He4*Y_D
            rate_He3T, rate_He4D = self.RR.get_He3T_to_He4D(T9,rho_b)
            change_LHS = rate_He3T*YHe3T-rate_He4D*YHe4D
            dY[2] += change_LHS # deutrium
            dY[3] -= change_LHS  # titrium
            dY[4] -= change_LHS  # helium 3
            dY[5] += change_LHS  #helium 4




        if self.NR_species > 6:  # Include Li7
            #equation 17
            YTHe4  = Y_T*Y_He4
            YLi7 = Y_Li7
            rate_THe4, rate_Li7 = self.RR.get_THe4_to_Li7(T9,rho_b)
            change_LHS = rate_THe4*YTHe4-rate_Li7*YLi7

            dY[3] -= change_LHS  #tittrium
            dY[5] -= change_LHS  #helium4
            dY[6] += change_LHS  # lithium

            #equation 20
            YpLi7 = Y_p*Y_Li7
            YHe4He4= Y_He4*Y_He4
            rate_pLi7, rate_He4He4 = self.RR.get_pLi7_to_He4He4(T9,rho_b)
            change_LHS = 0.5*rate_He4He4*YHe4He4 - rate_pLi7*YpLi7
            change_RHS = -rate_He4He4*YHe4He4 + 2*rate_pLi7*YpLi7

            dY[1] += change_LHS  # proton
            dY[6] += change_LHS  # lithium 
            dY[5] += change_RHS  # He4 


            
        if self.NR_species > 7:  # Include Be7
           #equation 16
           YHe3He4 = Y_He3*Y_He4
           YBe7 = Y_Be7
           rate_He3He4, rate_Be7 = self.RR.get_He3He4_to_Be7(T9,rho_b)
           change_LHS = rate_He3He4*YHe3He4-rate_Be7*YBe7
            
           dY[4] -= change_LHS  #helium 3
           dY[5] -= change_LHS  #helium 4
           dY[7] += change_LHS  #Be7
            #equation 18
           YnBe7 = Y_n*Y_Be7
           YpLi7 = Y_p*Y_Li7
           rate_nBe7, rate_pLi7 = self.RR.get_nBe7_to_pLi7(T9,rho_b)
           change_LHS = rate_nBe7*YnBe7-rate_pLi7*YpLi7
           dY[0] -= change_LHS  # neutron
           dY[1] += change_LHS  # proton
           dY[6] += change_LHS  # lithium
           dY[7] -= change_LHS  #beryllium7

        #equation 21
                  
           YnBe7 = Y_n*Y_Be7
           YHe4He4 = Y_He4*Y_He4

           rate_nBe7, rate_He4He4 = self.RR.get_nBe7_to_He4He4(T9,rho_b)
           change_LHS =0.5*rate_He4He4*YHe4He4-rate_nBe7*YnBe7
           change_RHS = 2*rate_nBe7*YnBe7-rate_He4He4*YHe4He4

           dY[0] += change_LHS  # neutron
           dY[7] += change_LHS  # Be7
           dY[5] += change_RHS  # He4
           
       
           
        
        
        
        # Each reaction equation are missing a factor(-1/Hubble)
        dY = -dY / Hubble
        return dY

  
    def get_IC(self, T_init):
        """
        Defines the initial condition array used in solve_BBN.
        The only nonzero values are for neutrons and protons, but the shape
        of Y_init is determined by the number of particles included.

        Parameters
        ----------
        T_init : float
            The initial temperature.

        Returns
        -------
        np.ndarray
            The initial condition array for the relative number densities of each species.
            
        """
        
        # Initialize the array with zeros for all included species:
        Y_init = np.zeros(self.NR_species)
        # Calculate the initial number densities for neutrons and protons:
        #print(const.m_n.decompose("kg"))
        #Yn_init = (1 + np.exp(const.m_n.decompose("kg")-const.m_p.decompose("kg"))*const.c.decompose("kg")**2/(const.k_B.decompose("kg")*100e9))**-1
        #Yp_init = 1 - Yn_init
        
        #Yn_init = (1 + np.exp((1.6749275e-27-1.67262193e-27)*299792458**2 / (1.380649e-23*100e9)))**-1
        Yn_init = (1 + np.exp((const.m_n.cgs.value-const.m_p.cgs.value)*const.c.cgs.value**2 / (const.k_B.cgs.value*100e9)))**-1
        Yp_init = 1 - Yn_init
        
        # Assign the initial number densities to the correct indices in Y_init:
        Y_init[0] = Yn_init
        Y_init[1] = Yp_init
        return Y_init

    def solve_BBN(self, T_init=100e9, T_end=0.01e9, NR_sol_points=1001):
        """
        Solves the BBN-system for a given range of temperature values

        Keyword Arguments:
            T_init {float} -- the initial temperature (default: {100e9})
            T_end {float} -- the final temperature (default: {0.01e9})
        """
        
        T_init = 100e9
        T_end = 0.01e9
        
        # Get initial conditions for the solver:
        Y_init = self.get_IC(T_init)

        # Define the logarithmic temperature range for the solver
        # as our equations are defined over ln(T):
        lnT_range = [np.log(T_init), np.log(T_end)]
        lnT_eval = np.linspace(lnT_range[0], lnT_range[-1], NR_sol_points)
    
        # Solve the ODE-system using scipy.solve_ivp:
        sol = solve_ivp(
            self.get_ODE,
            t_span=lnT_range,
            y0=Y_init,
            t_eval=lnT_eval,
            method="Radau",
            rtol=1e-12,
            atol=1e-12,
        )

        # Store the results:
        self.T = np.exp(sol.t)
        self.Y_i = sol.y


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    # Example use:
    # Initiate the system including 3 species (neutrons, protons and deuterium):
    bbn = BBN(NR_interacting_species=8)
    bbn.solve_BBN(T_end=0.1e9)  # solve the system until end temperature 0.1*10^9 K

    # Plot the mass fraction for each species:
    fig, ax = plt.subplots()
    for i, y in enumerate(bbn.Y_i):
        ax.loglog(bbn.T, bbn.mass_number[i] * y, label=bbn.species_labels[i])
        plt.xlim(10**11,10**7)
        plt.ylim(10**-11,10**-1)

        #plt.xlim(10**11,10**8)
        #plt.ylim(10**-3,10**0)
    #ax.plot( bbn.T,(1 + np.exp((const.m_n.cgs.value-const.m_p.cgs.value)*const.c.cgs.value**2 / (const.k_B.cgs.value*bbn.T)))**-1, ':', color = "blue", label = "Sigma A_iY_i")
    #ax.plot( bbn.T,1-((1 + np.exp((const.m_n.cgs.value-const.m_p.cgs.value)*const.c.cgs.value**2 / (const.k_B.cgs.value*bbn.T)))**-1),':', color = "orange", label = "Sigma A_iY_i")
    plt.grid()
    plt.xlabel("T[K]")
    plt.ylabel("Mass fraction [A_iY_i]")
    ax.legend()
    plt.show()





"""

Y_p = []
Y_D = []
Y_T = []
Y_He3 = []
Y_He4 = []
Y_Li7 = []
Y_Be7 = []


Y_p_a = np.array(Y_p)#iv.Y_p_array
Y_D_a = np.array(Y_D)#iv.Y_D_array
#Y_T_a = np.array(iv.Y_T)
#Y_He3_a = np.array(iv.Y_He3)
#Y_He4_a = np.array(iv.Y_He4)
#Y_Li7_a = np.array(iv.Y_Li7)
#Y_Be7_a = np.array(iv.Y_Be7)


for i in range(len(omega_b0)):

    bbn = BBN(iv.NR_species, Omega_b0=omega_b0[i])
    bbn.solve_BBN(T_end=0.1e9)
    Y_p.append(bbn.Y_i[1,-1])
    Y_D.append(bbn.Y_i[2,-1])
    #Y_T.append(bbn.Y_i[3,-1])
    #Y_He3.append(bbn.Y_i[4,-1])
    #Y_He4.append(bbn.Y_i[5,-1])
    #Y_Li7.append(bbn.Y_i[6,-1])
    #Y_Be7.append(bbn.Y_i[7,-1])


omega_log  = np.log(omega_b0)
cs_p = Cubicspline(omega_log,np.log(Y_p_a) )
cs_D = Cubicspline(omega_log, np.log(Y_D_a) )
#cs_T = Cubicspline(omega_log, np.log(Y_T_a) )
#cs_He3 = Cubicspline(omega_log,np.log(Y_He3_a) )
#cs_He4 = Cubicspline(omega_log,np.log(Y_He4_a) )
#cs_Li7 = Cubicspline(omega_log,np.log(Y_Li7_a)) 
#cs_Be7 = Cubicspline(omega_log,np.log(Y_Be7_a) )



plt.plot(omega_log, cs_D/cs_p)
#plt.plot(omega_log, cs_He3/cs_p)
#plt.plot(omega_log, cs_Li7/cs_p)
plt.show()

"""








