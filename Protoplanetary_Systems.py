from Systems import Star, Planet
import numpy as np
from astropy import constants as c
import OL18 #to calculate accretion efficiency
from PP import pebble_predictor
import matplotlib.pyplot as plt
from scipy import interpolate

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

def eta(rgrid, pg, rhog, OmegaK):
    rInt = np.zeros(np.size(rgrid) + 1)
    rInt[0] = 1.5 * rgrid[0] - 0.5 * rgrid[1]
    rInt[1:-1] = 0.5 * (rgrid[1:] + rgrid[:-1])
    rInt[-1] = 1.5 * rgrid[-1] - 0.5 * rgrid[-2]
    dr = rInt[1:] - rInt[:-1]
    pgInt = np.interp(rInt, rgrid, pg)
    return (pgInt[1:] - pgInt[:-1]) / dr[:] / (2. * rhog * OmegaK ** 2. * rgrid)





class Disc:

    '''
    A class used to represent the protoplanetary disc of a star
    
    Attributes
    ----------
    S       : star central star of the system 
    mdisk   : total gas mass of disc (float) [g]
    rgrid   : radial grid (ndarray) [cm] 
    rout    : critical radius of disc (float) [cm]
    Rsl     : snowline radius (float) [cm]
    alpha   : turbulence parameter (float) [none]
    beta    : temperature power law index (float) [none]
    zeta    : gas surface density power law index (float) [none]
    tgrid   : temporal grid (ndarray) [s]  
    
    Methods
    ----------
    T(r)            : gas temperature at r (float) [K]
    OmegaK(r)       : Keplerian frequency at r (float) [rad/s]
    cs(r)           : sound speed at r (float) [cm/s]
    hgas(r)         : dimensionless gas scale height at r (float) [none]
    Miso(r)         : pebble isolation mass at r (float) [g]
    SigmaGas(r)     : gas surface density at r (float) [g / cm^2]
    SigmaDust(r)    : dust surface density at r (float) [g / cm^2]
    rhog(r)         : midplane gas density at r (float) [g / cm^3]
    pg(r)           : gas pressure at r (float) [g / cm s^2]
    nu(r)           : viscocity at r (float) [g / cm s]
    Mdotg(t)        : gas flux onto the star after t (float) [g / s]
    vfrag(r)        : fragmentation velocity at r (float) [cm / s]
    solid_f(r)      : fractions of solids in the disc at r (float) [none]
    '''

    default_parameters = {'Z0' : 0.02, 'alpha' : 1e-4, 'beta' : 1 / 2, 'zeta' : 1, 'rout' : 30 * au, 'fInitial' : {'ice' : 1 / 2, 'iron' : 1 / 3}, 'vfragInitial': {'dust' : 1e2, 'ice' : 1e3}, 'rhop' : 1.25}

    def __init__(self, Star : Star = Star(), **pars) -> None:
        self.S = Star

        self.mdisk = 0.2 * Star.M

        if not 'rgrid' in pars:
            Rco = (Grav * Star.M * Star.P ** 2 / (3 * np.pi ** 2)) ** (1 / 3) # Corotation radius of the disc. Adapted from equation 3.127 of Philip J. Armitage 2017
            Rout = 1e3 * au
            Nr = 300
            self.rgrid = np.logspace(np.log10(Rco), np.log10(Rout), Nr)
        else: self.rgrid = pars['rgrid']
        if not 'tgrid' in pars:
            endtime = 1e7 * year
            Nt = 1000
            self.tgrid = np.logspace(0 * np.log10(year), np.log10(endtime), Nt)
        else: self.rgrid = pars['tgrid']
        for par in self.default_parameters:
            if par in pars: setattr(self, par, pars[par])
            else: setattr(self, par, self.default_parameters[par])
        
        self.Rsl = (150 / self.T(au)) ** (-1 / self.beta) * au # Snowline radius of the disc. Calculated using the snowline temperature and temperature profile from Table 1 of Gijs D. Mulders et al. 2021

        self.st_array, self.flux_array = pebble_predictor(
            rgrid = self.rgrid, 
            tgrid = self.tgrid,
            Mstar = self.S.M,
            SigmaGas = self.SigmaGas(self.rgrid),
            T = self.T(self.rgrid),
            SigmaDust = np.vectorize(self.SigmaDust)(self.rgrid),
            alpha = self.alpha,
            vfrag = np.vectorize(self.vfrag)(self.rgrid),
            rhop = self.rhop                        
        )
        self.eta_array = eta(
            rgrid = self.rgrid, 
            pg = self.pg(self.rgrid), 
            rhog = self.rhog(self.rgrid), 
            OmegaK = self.S.OmegaK(self.rgrid)
        )
        
        self.flux = interpolate.RegularGridInterpolator((self.tgrid, self.rgrid), self.flux_array)
        self.st = interpolate.RegularGridInterpolator((self.tgrid, self.rgrid), self.st_array)
        self.eta = lambda r: np.interp(r, self.rgrid, self.eta_array)


    def T(self, r) -> float: 
        '''Temperature of the disc in the midplane at radius r. Adapted from Table 1 of Gijs D. Mulders et al. 2021'''
        return 280 * (r / au) ** (-self.beta) * (self.S.M / MS) ** 0.5
    def cs(self, r) -> float: 
        '''Sound speed in the disc at radius r. Adapted from equation 7 of Jesper Nielsen et al. 2023'''
        return np.sqrt(k_b * self.T(r) / (2.3 * m_p))
    def hgas(self, r) -> float: 
        '''Aspect ratio of gas disc at radius r. Adapted from pebble_predictor'''
        return self.cs(r) / (self.S.OmegaK(r) * r)
    def Miso(self, r) -> float: 
        '''Pebble isolation mass in the disc at radius r. Adapted from equation 6 of Joanna Drążkowska et al. 2021 and Gijs D. Mulders et al. 2021'''
        return 40 * ME * self.S.M / MS * (self.hgas(r) / 0.05) ** 3
    def SigmaGas(self, r) -> float: 
        '''Gas surface density in the disc at radius r. Adapted from equation 1 of Joanna Drążkowska et al. 2021'''
        return self.mdisk / (2. * np.pi * self.rout ** 2.) * (r / self.rout) ** (-self.zeta) * np.exp(-1. * (r / self.rout))
    def SigmaDust(self, r) -> float: 
        '''Dust surface density in the disc at radius r. Adapted from equation 2 of Joanna Drążkowska et al. 2021. Modified to include the influence of a snowline'''
        return (self.fInitial['ice'] + self.fDust(r)['ice']) * self.Z0 * self.SigmaGas(r)
    def rhog(self, r) -> float: 
        return self.SigmaGas(r) * self.S.OmegaK(r) / (np.sqrt(2. * np.pi)*self.cs(r))
    def pg(self, r): return self.rhog(r) * self.cs(r) ** 2
    def nu(self, r): return self.alpha * self.cs(r) * self.hgas(r) * r
    def vfrag(self, r) -> float:
        '''fragmentation velocity of pebbles in the disc at radius r. Step function to accomodate influence of snowline'''
        return self.vfragInitial['dust'] if r <= self.Rsl else self.vfragInitial['ice']
    def fDust(self, r) -> dict:
        '''mass fractions of pebble population in the disc at radius r. Step functions to accomodate influence of snowline'''
        ice_frac = self.fInitial['ice'] if r > self.Rsl else 0
        return {
            'ice' : ice_frac,
            'rock' : (1 - self.fInitial['iron']) * (1 - ice_frac),
            'iron' : self.fInitial['iron'] * (1 - ice_frac)
        }
    def Mdotg(self, r) -> float:
        '''Mass flux of gas at radius r. Adapted from equations 12 and 16 of Jesper Nielsen et al. 2023'''
        ur = - 3 / 2 * self.alpha * self.cs(r) * self.hgas(r)
        return - 2 * np.pi * r * ur * self.SigmaGas(r)


class Protoplanet(Planet):
    '''A class used to represent a developing protoplanet. Extends Planet.
    
    Attributes
    ----------
    M                   : total mass (float) [g]
    Z                   : solid mass fraction (float) [none]
    fsolid              : mass fractions of solid components (dict) [none]
    R                   : core radius (float) [cm]
    rho                 : core density (float) [g / cm^3]
    disc                : protoplanetary disc (Disc)
    a                   : semimajor axis (float) [cm]
    Mdot                : mass flux (float) [g / s]
    vmig                : migration velocity (float) [cm / s]
    pebble accretion    : ability to accrete pebbles (bool)
    tau_growth          : growth timescale (float) [s]
    tau_mig             : migration timescale (float) [s]
    
    '''
    def __init__(self, disc : Disc = Disc(), a : float = 1 * au, M: float = 0.001 * ME, t: float = year):
        self.t_points = []
        self.M_points = []
        self.a_points = []
        self.Z_points = []
        self.f_points = {el : [] for el in self.rhoConst}
        self.disc = disc
        super().__init__(disc.S, a, M, 1, disc.fDust(a))
        
        self.t = t
        

        
    @property
    def M(self):
        return self._M
    @M.setter
    def M(self, M):
        self.pebble_accretion = M < self.disc.Miso(self.a)
        if hasattr(self,'M'):
            if self.pebble_accretion:
                for el in self.fsolid:
                    self.fsolid[el] = (self.fsolid[el] * self.M + self.disc.fDust(self.a)[el] * (M - self.M)) / M

            else:
                self.Z = self.Z * self.M / M
        self.Z_points.append(self.Z)

        for el, f in self.fsolid.items():
            self.f_points[el].append(f)

        self._M = M
        self.M_points.append(M) 

    @property
    def a(self):
        return self._a
    @a.setter
    def a(self, a):
        rgrid = self.disc.rgrid
        self.migration = rgrid[1] < a < rgrid[-1]
        a = np.clip(a, rgrid[1], rgrid[-1])
        self._a = a
        self.a_points.append(a)

    @property
    def t(self):
        return self._t
    @t.setter
    def t(self, t):
        tgrid = self.disc.tgrid
        self.simulation = t < tgrid[-1]
        t = np.clip(t, tgrid[0], tgrid[-1])
        self._t = t
        self.t_points.append(t)
        
    
    @property
    def Mdot(self):
        return self._Mdot
    
    @Mdot.getter
    def Mdot(self):
        '''
        mass flux of material onto the planet. Pebble accretion using `pebble-predictor` and `epsilon.tar.gz`. Gas accretion adapted from equations 32, 33, and 34 of Jesper Nielsen et al. 2023
        
        '''
        disc : Disc = self.disc
        if self.pebble_accretion:
            eta = disc.eta(self.a)
            st = disc.st((self.t, self.a))
            flux = disc.flux((self.t, self.a))

            rp = self.R / self.a
            qp = self.M / disc.S.M

            self.eff = OL18.epsilon_general(tau = st, qp = qp, eta = -1 * eta, hgas = disc.hgas(self.a), alphaz = disc.alpha, Rp = rp)

            Mdot = self.eff * flux
        
        else:
            Mdotg = disc.Mdotg(self.a)
            MdotKH = 1e-5 * (ME / year) * (self.M / (10 * ME)) ** 4 * 0.05
            MdotHill = 1.5 * 1e-3 * (ME / year) * (disc.hgas(self.a) / 0.05) ** 4 * (self.M / (10 * ME)) ** (4 / 3) * 0.01 / disc.alpha * (Mdotg / (1e-8 * MS / year))
            MdotHill /= (1 + (self.M / (2.3 * disc.Miso(self.a))) ** 2)

            Mdot = min(MdotKH, MdotHill, Mdotg)
        
        self.tau_growth = self.M / Mdot

        return Mdot
    
    @property
    def vmig(self):
        return self._vmig

    @vmig.getter
    def vmig(self):
        '''
        migration velocity of planet. Adapted from equations 9, 10, and 11 of Jesper Nielsen et al. 2023
        
        '''
        if self.migration:
            disc = self.disc
            k_mig = 2 * (1.36 + 0.62 * disc.beta * 0.43 * disc.zeta)
            v_mig = -k_mig * self.M * disc.SigmaGas(self.a) * self.a ** 2 / ((disc.S.M * disc.hgas(self.a)) ** 2) * self.a * disc.S.OmegaK(self.a)
            v_mig /= (1 + (self.M /(2.3 * disc.Miso(self.a))) ** 2)
            
            self.tau_mig = np.abs(self.a / v_mig)
        
        else: 
            v_mig = 0
        
        return v_mig
    
    def plot(self):
        f, ax = plt.subplots(2,2, figsize = (14,6), sharex = True)

        if self.a_points[0] > self.disc.Rsl and self.a <= self.disc.Rsl:
            RSLi = self.a_points.index(min(self.a_points, key = lambda a: abs(a - self.disc.Rsl)))
            ax[0][0].axvline(x = self.t_points[RSLi] / year, linestyle = ':', color = 'tab:cyan')
            ax[1][0].axvline(x = self.t_points[RSLi] / year, linestyle = ':', color = 'tab:cyan')
            ax[0][1].axhline(y = self.disc.Rsl / au, linestyle = ':', color = 'tab:cyan')
        
        if self.Z < 1:
            Misoi = len(self.Z_points) - list(reversed(self.Z_points)).index(1)
            ax[0][0].plot(self.t_points[Misoi] / year, self.M_points[Misoi] / ME, 'o', color = 'black')
            ax[0][0].plot(np.array(self.t_points) / year, self.disc.Miso(np.array(self.a_points)) / ME, color = 'grey', label = '$M_\\text{total}$', linestyle = '--', alpha = 0.6)
        
        ax[0][0].loglog(np.array(self.t_points) / year, np.array(self.M_points) / ME, color = 'black', label = '$M_\\text{total}$')


        colors = ['tab:red', 'tab:grey', 'tab:blue']
        colors.reverse()


        ax[1][0].plot(np.array(self.t_points) / year, 1 - np.array(self.Z_points), '--', label = 'gas', color = 'black')

        for el, c in zip(self.fsolid, colors):
            ax[1][0].plot(np.array(self.t_points) / year, np.array(self.Z_points) * np.array(self.f_points[el]), color = c, label = el)

        ax[1][0].set_xscale('log')

        ax[0][0].set_ylabel('Mass [$M_\\oplus$]')
        ax[1][0].set_xlabel('time [yr]')
        ax[1][1].set_xlabel('time [yr]')
        ax[1][0].set_ylabel('Mass fraction')
        ax[1][0].legend()
        ax[0][1].set_ylabel('orbital radius [au]')
        ax[0][1].set_xlim(right = 1e7)

        ax[0][1].loglog(np.array(self.t_points) / year, np.array(self.a_points) / au, color = 'white')

        xmin, xmax = ax[0][1].get_xlim()
        ymin, ymax = ax[0][1].get_ylim()

        xmini, xmaxi = (self.disc.tgrid/year).searchsorted(xmin), (self.disc.tgrid/year).searchsorted(xmax)
        ymini, ymaxi = (self.disc.rgrid/au).searchsorted(ymin), (self.disc.rgrid/au).searchsorted(ymax)

        x = ax[0][1].imshow(self.disc.flux_array[xmini:xmaxi, ymini:ymaxi] / ME * year, aspect = 'auto', extent = (xmin, xmax, ymin, ymax), norm = 'log')

        f.colorbar(x, label = 'Pebble mass flux [$M_\\oplus$ yr$^{-1}$]')




