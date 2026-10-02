# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 16:17:42 2026

@author: salma
"""
import numpy as np
from astropy import constants as constant
from astropy import units as unit

class BackgroundBBN:
    def __init__(self, Omega_b0 = 0.05, N_eff = 3.0):
        """
        Initialize the background cosmology defined by input parameters.

        Parameters
        ----------
        Omega_b0 : float, optional
            The relative baryon density, by default 0.05
        N_eff : float, optional
            The effective number of neutrinos, by default 3.0
        """

        # Define constants and parameters:
        self.Omega_b0 = Omega_b0
        self.N_eff = N_eff
        self.H0 = 100*0.7 /3.09e19
        #self.H0 =((100*0.7*unit.km)/(unit.s*unit.Mpc)).cgs.value
        self.rho_crit0 = 9.2e-27/1000 #kgm**-3 #3*self.H0/(8*np.pi*constant.G.cgs.value)#
        self.Omega_r0 = 10e-4
        self.kb = constant.k_B.cgs.value  #1.380649e-23# m2 kg s-2 K-1
        self.mp = constant.m_p.cgs.value  #1.6726219e-27 #kg
        #self.mH = constant.m_H.cgs.value  #1.6735575e-27
        self.mn = constant.m_n.cgs.value  #1.6749e-27 
        self.pi = np.pi
        self.h_planck = constant.h.cgs.value #Js
        self.h_bar = constant.hbar.cgs.value #1.0545718e-34 #Js
        self.c = constant.c.cgs.value #m/s
        self.pi = np.pi
        self.e = np.exp(1)
        self.Q = self.c*(self.mn-self.mp)
        
        
        
    def get_a(self, T):
        """
        Compute the scale factor as a function of temperature.

        Parameters
        ----------
        T : float or np.ndarray
            The temperature in Kelvin.

        Returns
        -------
        float or np.ndarray
            The scale factor.
        """
        T0 = 2.725
        a = (T0/T)
        return a

    def get_Hubble(self, T):
        """
        Compute the Hubble parameter as a function of temperature.

        Parameters
        ----------
        T : float or np.ndarray
            The temperature in Kelvin.

        Returns
        -------
        float or np.ndarray
            The Hubble parameter in 1/s.
        """
        H = self.H0*np.sqrt(self.Omega_r0)*self.get_a(T)**(-2)
        return H

    def get_rho_b(self, T):
        """
        Compute the energy density of baryons as a function of temperature.

        Parameters
        ----------
        T : float or np.ndarray
            The temperature in Kelvin.

        Returns
        -------
        float or np.ndarray
            The energy density of baryons in g/cm^3.
            
        """
         
        
        #nb = 2.404/self.pi**2 * (self.kb*T / (self.h_bar*self.c))**3 *2.7e-8*0.05*0.7**2 
         
        #rho_b = nb* self.mp
        
        rho_b = self.Omega_b0*self.rho_crit0*(self.get_a(T)**-3)
         
        return rho_b
        
ic = BackgroundBBN()

print(ic.get_rho_b(2.725))


