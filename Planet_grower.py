import Protoplanetary_Systems as dq
import numpy as np

def planet_grower(p : dq.Protoplanet):
    disc : dq.Disc = p.disc

    tgrid = disc.tgrid
    rgrid = disc.rgrid

    while p.simulation:
        Mdot = p.Mdot

        if p.migration:
            vmig = p.vmig
            tstep = 0.01 * min(p.tau_growth, p.tau_mig, 10 * p.t)
        else:
            tstep = 0.01 * min(p.tau_growth, 10 * p.t)
        
        p.a += tstep * vmig
        p.M += tstep * Mdot
        p.t += tstep

    return p

