import numpy as np
import astropy.constants as c

# Physical constants
Grav = c.G.cgs.value
k_b = c.k_B.cgs.value
m_p = c.m_p.cgs.value

# Mass constants
MS = c.M_sun.cgs.value
MJ = c.M_jup.cgs.value
ME = c.M_earth.cgs.value

# Length constants
au = c.au.cgs.value 
RE = c.R_earth.cgs.value

# Time constants
day = 24 * 3600
year = 365.25 * day

class Star:
    '''
    A class that represents a star

    Attributes
    ----------
    M : stellar mass (float) [g]
    P : period of rotation (float) [s]
    
    Methods
    ----------
    OmegaK: Keplerian orbital frequency
    
    '''
    def __init__(self, Mass = MS, Period = 7 * day):
        self.M = Mass
        self.P = Period

    def OmegaK(self, r): 
        return np.sqrt(Grav * self. M / r ** 3)    

class Planet:
    rhoConst =  {'ice' : 0.95, 'gas' : 1.326, 'rock' : 5.6, 'iron' : 10} 
    
    '''
    A class to represent a planet. 
    
    Attributes
    ----------
    M       : Total mass (float) [g]
    Z       : solid mass fraction (float) [none]
    fsolid  : mass fractions of solid components (dict) [none]
    fgas    : mass fractions of gaseous components (dict) [none]
    R       : Total radius (float) [cm]
    rho     : average (float) [g / cm^3]
    '''

    def __init__(self, M : float, Z : float, fsolid : dict):
        self.fsolid = fsolid
        self.Z = Z
        self.M = M
    
    
    @property
    def rho(self):
        return self._rho
    @rho.getter
    def rho(self):
        return 1 / sum(f / self.rhoConst[el] for el, f in self.fsolid.items())
    @property
    def R(self):
        return self._R
    @R.getter
    def R(self):
        return (3 * self.M / (4 * np.pi * self.rho)) ** (1 / 3)


