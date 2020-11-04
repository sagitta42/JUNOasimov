import sys

# oscillation module (class Model)
sys.path.append('/home/sagitta/mymodules')
from oscillation import *

sys.path.append('/home/sagitta/mymodules/myplot')
from myplot import *

from copy import deepcopy

## to enable backslash e.g. in \text
# import matplotlib as mpl
# mpl.rcParams['text.usetex'] = True
# mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command

# ----------------------------------------------------------------------------------------------

# JUNO baseline
Ljuno = 53 # km

# number of target protons
target_p = 1.5 * 1e33

# energy points in MeV
# Epoints = np.arange(0.1,10,5000)
Epoints = np.linspace(0.1,10,6000)
# Epoints = np.linspace(0.1,10,5000)

# masses in MeV
m_n = 939.57 # neutron
m_p = 938 # proton
m_e = 0.511 # positron
# threshold of IBD in MeV
thr_ibd = m_n - m_p + m_e

# ----------------------------------------------------------------------------------------------

# oscillation parameters
par = Parameters(source='NuFIT')
osc = par.generate_input(N=0) # N = 0 -> take default osc parameters
# osc = osc_pars.loc[0] # N = 0 returns a table with one row

# create a model based on those parameters
model_default = Model(osc['th12'], osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])

# ----------------------------------------------------------------------------------------------

LABELS = {
'unosc': 'no oscillation',
'osc0': r'$\theta_{13} = 0$',
'oscNO': 'NO',
'oscIO': 'IO',
'osc': 'custom'
}

COLORS = {
'unosc': 'k', 'osc0': 'g', 'oscNO': 'b', 'oscIO': 'r',
'osc': 'k',
'data': 'k'
}

# ----------------------------------------------------------------------------------------------


def plot_osc_spectra(E,L,reso=0):
    '''
    Plot JUNO spectrum
    1) unoscillated
    2) oscillated with theta13 = 0
    3) oscillated with NO
    3) oscillated with IO

    smear: energy resolution in percentage
            if non-zero, Gaussian smear will be applied
            otherwise, no smearing is applied
    '''

    print ('###### JUNO spectra')

    spec = Spectra(E,L)

    # get reactor spectrum
    spec.get_spectrum()
    # apply oscillation with NO and IO, and with theta13 = 0 for comparison
    spec.oscillate()
    # apply det resp
    # spec.det_response(1100)
    # smear if resolution is given
    if reso:
        spec.smear(reso)

    # get spectra PDFs
    # spec.get_pdfs()
    # return

    # convert to histo
    # spec.histogramize(Nbins=250,sample=False,Nev=100000)

    # remove unosc
    spec.spectra.pop('unosc')
    spec.spectra.pop('osc0')
    mp = Plot((10,8))

    # plot all resulting spectra
    # spec.plot_histos(mp)
    spec.plot_spectra(mp)

    fname = 'juno_histos_theta13'
    # fname = 'juno_spectrum_theta13'
    if reso:
        mp.fig.suptitle('Resolution {0} %'.format(reso))
        fname += '_smear{0}'.format(reso)

    mp.figure(fname + '.png')


# ----------------------------------------------------------------------------------------------


class Spectra():
    ''' Class that contains methods of constructing spectra, PDFs and histograms
        for unoscillated and oscillated scenarios.
        Methods also include smearing and sampling based on PDF.
    '''

    def __init__(self, E, L):
        ''' E: energy in MeV, L: distance in km '''
        # take energy above IBD threshold
        self.Epoints = E[E >= thr_ibd]
        # energy points: start with MeV, later can be in p.e.
        self.E = self.Epoints
        self.Npe = False # we start with MeV
        # bin edges after constructing histograms will be saved here
        self.Ebinned = []

        self.L = L

        # spectra (continuous function)
        self.spectra = {}
        # PDFs (continuous function norm to 1)
        self.pdfs = {}
        # binned histograms
        self.histos = {}


    def get_spectrum(self):
        # JUNO spectrum (unoscillated)
        self.spectra['unosc'] = target_p * sigmaIBD(self.Epoints) * antinu_flux(self.Epoints)


    def oscillate(self,motypes=['NO','IO'],model=model_default):
        ''' Create oscillated spectrum based on the unoscillated one
                obtained with get_spectrum().
            motypes [list]: mass ordering. Possible values in the list: 'NO', 'IO', None.
                By default both NO and IO are constructed.
                If None is given, instead of the default value for dm2atm,
                the given one will be taken.
            mod [Model]: Model object which calculates Pee for the energy points
                By default model with standard mean measured parameters is taken
        '''

        # oscillation when theta13 = 0
        mod0 = deepcopy(model)
        mod0.th13 = 0

        pee0 = mod0.get_pee(L=self.L, E=self.Epoints, mo='NO') # MO doesn't matter for th13 = 0
        self.spectra['osc0'] = pee0*self.spectra['unosc']

        # oscillated flux for given mass orderings
        for mo in motypes:
            pee = model.get_pee(L=self.L, E=self.Epoints, mo=mo)
            self.spectra['osc' + mo if mo != None else 'osc'] = pee*self.spectra['unosc']


    def det_response(self, LY, sp=None):
        '''
        Convert MeV to Npe
        LY: value of p.e./MeV @ 1 MeV
        '''
        spectra = self.spectra if sp == None else [sp]

        for sp in spectra:
            self.spectra[sp] = self.spectra[sp] * LY # p.e.

        self.E = self.Epoints * LY # p.e.
        self.Npe = True


    def smear(self, reso=-1, sp=None):
        '''
        Convolve each spectrum with Gaus, The smear_gaus() function is contained
            in functions.py
        reso: energy resolution in %
        [!] not implemented: normalizing to the same number of events
        Is used when constructing PDFs, but there it doesn't matter
            cause they get normalized to 1 anyway, but for actual spectra
            it is needed
        '''

        # print 'Smearing...'
        # OINK
        # histos = self.histos if sp == None else [sp]
        spectra = self.spectra if sp == None else [sp]
        # OINK
        for sp in spectra:
        # for sp in histos:
            # print ('...', sp)
            # OINK
            # self.histos[sp] = smear_gaus(self.histos[sp], self.Ebinned, reso/100.)
            res = -1 if self.Npe else reso/100. # -1 given to Gaus means use Npe
            self.spectra[sp] = smear_gaus(self.spectra[sp], self.E, res)


    def plot_spectra(self,mp):
        ''' Plot continuous spectra on a given Myplot object mp'''
        for sp in self.spectra:
            mp.ax.plot(self.E, self.spectra[sp], color=COLORS[sp], label=LABELS[sp])

        mp.ax.set_xlabel('Npe' if self.Npe else 'E [MeV]')
        mp.legend()
        mp.pretty()


    def plot_histos(self,mp,order=None):
        '''
        Plot histograms on a given Myplot object mp
        order [list]: order in which to plot (if None, default is taken)
        '''

        order = self.histos if order == None else order
        for sp in order:
            print (sp)
            mp.ax.step(self.Ebinned, self.histos[sp], color=COLORS[sp], label=LABELS[sp])

        xlabel = 'Npe' if self.Npe else 'E [MeV]'
        mp.ax.set_xlabel(xlabel)
        mp.legend()
        mp.pretty(large=3)


    def get_pdfs(self, sp=None):
        ''' Construct and save PDFs of given spectrum.
            If no spectrum is given, PDFs are constructed for all present spectra'''

        spectra = self.spectra if sp == None else [sp]
        for sp in spectra:
            self.pdfs[sp] = self.spectra[sp] / self.spectra[sp].sum()
            ## save as csv
            # fname = sp + '.csv'
            # pd.DataFrame({'x': self.E, 'y': pdf }).to_csv(fname,\
            #     header=False, index=False, sep = ' ')
            # print '-->', fname


    def histogramize(self,Nbins,Nev,sample,sp=None):
        '''
        Create histogram for the given spectrum sp with Nbins bins
            normalized to Nev based on the PDF. If no spectrum name is given,
            histograms are constructed for all present spectra

        If sample is True, the histo is created by sampling Nev events based on
            the probabilities given by the PDF. The resulting bin content is int,
            resembling data.
            Otherwise, the PDF is binned and scaled by Nev. Resulting bin content
                is float (used as reference for fitting).
        '''

        pdfs = self.pdfs if sp == None else [sp]

        for sp in pdfs:
            if sample:
                # sample Nev values based on the PDF
                values = np.random.choice(self.E, size=Nev, p=self.pdfs[sp])
                # construct a histo
                self.histos[sp], bin_edges = np.histogram(values, bins=Nbins, range=[min(self.E), max(self.E)])
            else:
                # bin the PDF and scale by Nev
                self.histos[sp], bin_edges = np.histogram(self.E, bins=Nbins, range=[min(self.E), max(self.E)], weights=self.pdfs[sp]*Nev)

        # new binned energy points
        self.Ebinned = bin_edges[:-1]


# -------------------------------------------------------------

# energy resolution function
def energy_reso(E, a, c):
    return np.sqrt(a**2 * E + c**2)

def gaus(x, mu, sigma, a=1):
    ''' Gauss function '''
    return a / sigma / np.sqrt(2*np.pi) * np.exp( -0.5 * (x-mu)**2 / sigma**2 )


def smear_gaus(fpoints, x, reso):
    '''
    Smear a function at each point in x with gaus
    fpoints: y points of the function
    x: corresponding x points (MeV or Npe)
    reso: energy resolution at 1 MeV (in %) for the case of MeV.
        For the case of Npe, reso=-1 should be given

    '''
    # empty array
    res = np.zeros(len(x))

    # sum gaussian distributions along the original point for each point
    for i in range(len(x)):
        # current x point
        xp = x[i]
        #  mean = x point, amplitude -> y point at xp
        # sigma(E) / E = reso / sqrt(E)
        sigma = np.sqrt(xp) if reso == -1 else reso*np.sqrt(xp)
        res += gaus(x, mu=xp*1., sigma=sigma, a=fpoints[i])

    return res


# -------------------------------------------------------------

def antinu_flux(E):
    '''
    Reactor antineutrino flux, E in any units
    '''

    # print ('Calculating reactor flux...')

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

    # print ('Calculating cross section...')

    # total energy of positron
    Ee = E - (m_n - m_p)

    # Ee^2 = (pe c)^2 + (me c^2)^2
    pe = np.sqrt(Ee**2 - m_e**2) # momentum of positron

    return 0.0952 * Ee * pe * 1e-42 #cm^2

# -------------------------------------------------------------


def plot_flux_sigma(E):
    '''
    Subplot 1: reactor flux and IBD cross section (second y axis)
    Subplot 2: resulting unoscillated flux
    '''

    print ('###### Reactor flux, IBD cross section and unoscillated flux')


    mp = Plot((10,8),2, sharex=True)

    ## reactor flux
    flux = antinu_flux(E)
    mp.axes[0].plot(E, flux)

    ## IBD cross section
    # energy above IBD threshold
    Eibd = E[E >= thr_ibd]
    sigma = sigmaIBD(Eibd)

    # plot on second Y axis
    mp.add_axis(col='r')
    # ax 1 is the 2nd subplot, ax 2 is the 2nd Y axis of subplot 1
    mp.axes[2].plot(Eibd, sigma, color='r')

    ## measured flux
    # flux above IBD (to match with sigma)
    # flux = flux[len(E) - len(Eibd):]
    flux = flux[E >= thr_ibd]
    fluxJUNO = flux * sigma * target_p
    mp.axes[1].plot(Eibd, fluxJUNO) # axis 1 is the second subplot

    ## prettification
    mp.axes[0].set_ylabel(r'Reactor flux [cm${}^{-2}$]')
    mp.axes[2].set_ylabel(r'IBD cross section [cm${}^2$]', color='r')
    mp.axes[1].set_xlabel('E [MeV]')
    mp.axes[1].set_ylabel('Spectrum at JUNO')

    mp.pretty()
    mp.figure('flux_sigma_spectrum_vs_E.png')


def plot_P_vs_E(E,L,mod=model_default):
    ''' Plot survival probability as a function of energy for given baseline '''

    print ('###### Electron antinu survival probability')


    mp = Plot((10,8))

    # get_pee() takes energy in MeV
    mp.ax.plot(E, mod.get_pee(L=L, E=E, mo='NO'), label='NO', color='b')
    mp.ax.plot(E, mod.get_pee(L=L, E=E, mo='IO'), label='IO', color='r')


    mp.ax.set_xlabel('E [MeV]')
    mp.ax.set_ylabel(r'$P(\bar\nu_e \to \bar\nu_e)$')
    mp.fig.suptitle('Survival probability at L = {0} km'.format(L))

    mp.legend(title='Mass ordering', ncol=2)
    mp.pretty()
    mp.figure('P_vs_E.png')



if __name__ == '__main__':

    # plot_P_vs_E(Epoints, Ljuno)
    # plot_flux_sigma(Epoints)

    ## JUNO spectrum without smearing
    # plot_osc_spectra(Epoints, Ljuno)

    ## save JUNO spectra PDFs for 3% resolution
    plot_osc_spectra(Epoints, Ljuno, 3)


    ## plot JUNO spectrum for different resolutions
    # for reso in [0.1, 1, 2, 3, 4, 5]:
    #     print 'Resolution {0} %'.format(reso)
    #     plot_osc_spectra(Epoints, Ljuno, reso) # smear
    #     print '--------'
