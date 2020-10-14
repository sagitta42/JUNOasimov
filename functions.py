import numpy as np

## to enable backslash e.g. in \text
# import matplotlib as mpl
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

# -------------------------------------------------------------

# masses in MeV
m_n = 939.57 # neutron
m_p = 938 # proton
m_e = 0.511 # positron

# -------------------------------------------------------------
# 
# # color coding for mass ordering
# COLORS = {'NO': 'b', 'IO': 'r'}

# -------------------------------------------------------------


def antinu_flux(E):
    '''
    Reactor antineutrino flux, E in any units
    '''

    print ('Calculating reactor flux...')

    ## f: relative fission contribution (actually varies over time)
    fU235 = 0.58
    fPu239 = 0.3
    fU238 = 0.07
    fPu241 = 0.05

    return fU235 * np.exp(0.870 - 0.16*E - 0.091*E**2)\
            + fPu239 * np.exp(0.896 - 0.239*E - 0.0981*E**2)\
            + fU238 * np.exp(0.976 - 0.162*E - 0.079*E**2)\
            + fPu241 * np.exp(0.793 - 0.08*E - 0.1085*E**2)


def sigmaIBD(E):
    '''
    Cross section of neutrino capture via IBD
    E in MeV
    '''

    print ('Calculating cross section...')

    # total energy of positron
    Ee = E - (m_n - m_p)

    # Ee^2 = (pe c)^2 + (me c^2)^2
    pe = np.sqrt(Ee**2 - m_e**2) # momentum of positron

    return 0.0952 * Ee * pe * 1e-42 #cm^2


def gaus(x, mu, sigma, a=1):
    ''' Gauss function '''
    return a / sigma / np.sqrt(2*np.pi) * np.exp( -0.5 * (x-mu)**2 / sigma**2 )


def smear_gaus(fpoints, x, reso):
    '''
    Smear a function at each point in x with gaus
    fpoints: y points of the function
    x: corresponding x points
    reso: energy resolution at 1 MeV (from 0.0 to 1.0)

    '''
    # empty array
    res = np.zeros(len(x))

    # sum gaussian distributions along the original point for each point
    for i in range(len(x)):
        # current x point
        xp = x[i]
        #  mean = x point, amplitude -> y point at xp
        # sigma(E) / E = reso / sqrt(E)
        res += gaus(x, mu=xp, sigma=reso*np.sqrt(xp), a=fpoints[i])

    return res
