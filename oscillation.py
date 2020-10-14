import numpy as np
import pandas as pd
import sys

## different sets of parameters
PARAMS = {}
## NuFIT 4.1 (July 2019) with SK atm data
PARAMS['NuFIT'] = {}
# simplified symmetric sigma
PARAMS['NuFIT']['th12'] = np.array([33.82, 0.78]) # deg
PARAMS['NuFIT']['th23'] = np.array([48.6, 1.4]) # deg
PARAMS['NuFIT']['th13'] = np.array([8.6, 0.13]) # deg
PARAMS['NuFIT']['dm2sol'] = np.array([7.39, 0.21])*1e-5 # ev^2
PARAMS['NuFIT']['dm2atm'] = np.array([2.528, 0.031])*1e-3 # eV^2

# PARAMS['sin2theta12'] = [0.310, 0.013]
# PARAMS['sin2theta12'] = [0.310, 0.013]

## LEGEND proposal plot
PARAMS['LEGEND'] = {}
PARAMS['LEGEND']['s2th12'] = np.array([0.297, 0])
PARAMS['LEGEND']['s2th23'] = np.array([0, 0]) # not relevant for the LEGEND plot
PARAMS['LEGEND']['s2th13'] = np.array([0.0216, 0])
PARAMS['LEGEND']['dm2sol'] = np.array([73.7*1e-6, 0])# eV^2
PARAMS['LEGEND']['dm2atm'] = np.array([2540*1e-6, 0.]) # eV^2



class Parameters():
    def __init__(self,source='NuFIT'):
        '''
        Set of oscillation parameters based on measurements
        source: source of values
                'NuFIT': July 2019 NuFIT results (http://www.nu-fit.org/?q=node/211)
                'LEGEND': values reported on the plot in the LEGEND proposal
        '''

        ### -------------------------------------------------

        ## get values
        values = PARAMS[source]

        ## init table: columns -> parameters, row 0 - mean, row 2 - sigma
        self.params = pd.DataFrame()

        ## NuFIT are angles in degrees
        if source == 'NuFIT':
            for theta in ['th12','th23','th13']:
                self.params[theta] = values[theta]  / 180. * np.pi # in rad

        ## LEGEND plot are in sin^2 (theta_ij)
        elif source == 'LEGEND':
            for theta in ['th12','th13', 'th23']:
                self.params[theta] = np.arcsin(np.sqrt(values['s2'+theta])) # in rad

        for dm in ['dm2sol', 'dm2atm']:
            self.params[dm] = values[dm] # eV^2

        ### -------------------------------------------------

        ## transpose: columns -> Mean, Sigma, rows -> parameters
        self.params = self.params.transpose()
        # add column names
        self.params.rename(columns = {0: 'mean', 1: 'sigma'}, inplace=True)

        print ('Oscillation parameters based on', source)
        print ('---------------')
        print (self.params)
        print ('---------------')


    def generate_input(self, N):
        '''
        N > 0: generate N sets of oscillation parameters from a Gaussian distribution
            with (mu,sigma) given by the measurement
        N = 0: return default parameters (mu,sigma)
        '''

        if N == 0: return self.default_params()

        ## init table: column -> i (0 to N), row -> parameters
        inp = pd.DataFrame()

        for i in range(N):
            # random generation (pd.Series with index as parameters)
            rnd = self.generate_params()
            # add column
            inp[i] = rnd

        ## transpose
        inp = inp.transpose()

        print ('Generated {0} inputs'.format(N))
        print ('---------------')
        print (inp)
        print ('---------------')

        return inp


    def default_params(self):
        ''' Default mean measured values rather than random '''
        res = self.params.transpose()
        res = res.drop('sigma', axis=0)
        res = res.reset_index()
        res = res.drop('index', axis=1)
        res = res.loc[0] # N = 0 returns a table with one row
        return res


    def generate_params(self):
        ''' Generate a random value from a Gaussian distribution for each parameter '''

        ## small hepler function
        def gauss(par):
            return np.random.normal(par['mean'], par['sigma'])

        return self.params.apply(gauss, axis=1)





class Model():
    def __init__(self, th12, th23, th13, dm2sol, dm2atm, alpha=0, beta=0):
        '''
        th12, th23, th13: mixing angles in radians
        dm2_sol: solar mass splitting in eV^2
        dm2_atm: atmospheric mass splitting in eV^2
        alpha, beta: Majorana phases in radians
        '''

        self.th12 = th12
        self.th23 = th23
        self.th13 = th13
        self.dm2sol = dm2sol
        self.dm2atm = dm2atm
        self.alpha = alpha
        self.beta = beta

        # useful for variable Pee
        # self.params = {'th12': self.th12, 'th23': self.th23, 'th13': self.th13,
        #                 'dm2sol': self.dm2sol, 'dm2atm': self.dm2atm}


    def get_mBB(self, m_lightest, mo):
        '''
        array, string -> array

        Calculates and returns effective Majorana mass (mBB) as a function of
        the mass of the lightest neutrino (m_lightest in eV)
        mo ['NO'|'IO']: mass ordering (assuming dm2atm is absolute value)
        '''

        ## masses depending on the mass ordering
        if mo == 'NO':
            m1 = m_lightest
            m2 = np.sqrt(m1**2 + self.dm2sol)
            m3 = np.sqrt(m2**2 + self.dm2atm)
        elif mo == 'IO':
            m3 = m_lightest
            m1 = np.sqrt(m3**2 + self.dm2atm)
            m2 = np.sqrt(m1**2 + self.dm2sol)

        ## calculation of mbb

        ## mbb = | c12^2 c13^2 m1 + s12^2 c13^2 m2 e^ialpha + s13^2 m3 e^ibeta
        ## let a = c12^2 c13^2 m1
        ## b = s12^2 c13^2 m2
        ## c = s13^2 m3
        ## -> mbb = | a + b e^ialpha + c e^ibeta |
        ## mbb^2 = a^2 + b^2 +c^2 + 2bc cos(alpha -beta) + 2ab(cos(alpha) + cos(beta))

        c13sq = np.cos(self.th13)**2 # cos^2(theta_13)

        a = np.cos(self.th12)**2 * c13sq * m1
        b = np.sin(self.th12)**2 * c13sq * m2
        c = np.sin(self.th13)**2 * m3

        ## square of Majorana mass
        mbb2 = a**2 + b**2 + c**2\
                + 2*b*c*np.cos(self.alpha - self.beta)\
                + 2*a*( b*np.cos(self.alpha) + c*np.cos(self.beta) )

        ## if we switch phases around
        # mbb2 = a**2 + b**2 + c**2\
        #         + 2*a*b*np.cos(self.alpha - self.beta)\
        #         + 2*c*( a*np.cos(self.alpha) + b*np.cos(self.beta) )

        return np.sqrt(mbb2)


    def get_pee(self, L, E, mo=None):
        '''
        float, array, string -> array

        Calculates and returns electron antinu survival probability as a function of energy

        E: energy in MeV
        L: distance in km
        mo: mass ordering ['NO'|'IO'|None]
            If mass ordering is not given (i.e. None is given), value of dm2atm is taken
            directly. Otherwise, dm2atm is taken as absolute value, and the factor +-
            is decided based on given mo
        '''

        # if MO is given, we take the abs dm2atm and a factor
        # if not, dm2atm can be negative
        factor = 1 if mo == None else {'NO': 1, 'IO':-1}[mo]
        dm2atm = factor * self.dm2atm


        # formula dm^2 L/ 4E -> 1.27 dm^2[eV^2] L[km]/E[GeV]
        Delta21 = 1.27 * self.dm2sol * L / (E/1000.) # L [km] E [MeV]
        # absolute value by default, non-absolute if dm2atm is given as negative
        Delta31 = 1.27 * dm2atm * L / (E/1000.) #* factor


        s22th13 = np.sin(2*self.th13)**2
        s2th12 = np.sin(self.th12)**2
        s2Delta21 = np.sin(Delta21)**2

        pee = 1 - np.cos(self.th13)**4 * np.sin(2*self.th12)**2 * s2Delta21 \
            - s22th13 * np.sin(Delta31)**2\
            - s2th12 * s22th13 * s2Delta21  * np.cos(2*Delta31)\
            + s2th12 / 2 * s22th13 * np.sin(2*Delta21) * np.sin(2*Delta31)

        return pee
