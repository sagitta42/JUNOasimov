## Example
# python juno_fit.py Nev=100000 fit_method=chi2 data_method=asimov save

### Needs to be run with python3 because iminuit in python2 has a bug
from juno_osc import *
import math
import time
import os

from iminuit import Minuit

import itertools

###########################
## ~~~ Initial setup ~~~ ##
###########################

defaults = {
'fit_method': 'chi2',
'Nbins': 200,
'data_method': 'asimov'
}

all_params = ['th12', 'th13', 'a', 'b']
# all_params = ['mu', 'sig', 'a', 'b']

#free_params = ['th12', 'th13', 'a', 'b']
free_params = ['th12', 'a', 'b']
#free_params = ['th13', 'a']
#free_params = ['th12', 'a']
# free_params = ['mu', 'sig', 'eres', 'b']
# free_params = ['mu', 'sig', 'a', 'b']
# free_params = ['mu', 'eres', 'b']
# free_params = ['mu', 'eres']
fixed_params = np.setdiff1d(all_params, free_params)

# ------------------------------------------------------------------

def user_dict():
    '''
    parsing user input
    Example:
    python juno_fit.py Nev=10000 fit_method=chi2 data_method=asimov Nbins=200 save
    Returns dict
    {'Nev': '1000', 'fit_method': 'chi2', 'data_method': asimov, 'Nbins':200}
    If not given, variables are set to default values defined in dictionary default at the bottom
    '''
    # config = open(sys.argv[1])
    # lines = config.readlines()
    lines = sys.argv[1:]

    # read user input
    mydict = {}
    for line in lines:
        if not '=' in line: continue
        var,val = line.split('=')
        mydict[var] = val.strip()

    # assign defaults
    for key in defaults:
        if not key in mydict:
            mydict[key] = defaults[key]

    return mydict

# ------------------------------------------------------------------

# Argument check upon starting
if len(sys.argv) == 1:
    print ('Available variables')
    print ('\t Nev: number of events in the simulated/Asimov sample')
    print ('\t fit_method: chi2 or lkl (default chi2)')
    print ('\t data_method: sample, asimov or gaus (default sample)')
    print ('\t Nbins: number of bins (default 250)')
    sys.exit(1)

user_input = user_dict()    

# ------------------------------------------------------------------

# global variables

Nev = int(user_input['Nev'])
fit_method = user_input['fit_method']
data_method = user_input['data_method']


# JUNO baseline
Ljuno = 53 # km

# MO = 'IO'
MO = 'NO'

# ~2 - 10 MeV, bin = 30 keV @ 1 MeV -> 250 bins
Nbins = int(user_input['Nbins'])

LY = 1100 # p.e./MeV @ 1 MeV
# Ereso = 1/np.sqrt(LY)*100 # %, assumed resolution
#apar = 1/np.sqrt(LY) # MeV^(1/2)
apar = 2.9/100 # MeV^{1/2}, 2.9%
bpar = 0.8/100 # no units, 0.8%

# ------------------------------------------------------------------
# standard oscillation model is already initialized in juno_osc.py as model_default

Epoints = np.linspace(0.1,10,6000)
# Epoints = np.linspace(0.1,10,10000)
# Epoints = np.linspace(0.1,10,5000)

# JUNO spectra for E and L that will define "data"
spec0 = Spectra(Epoints,Ljuno)
# after this the energy will be changed to be above IBD threshold

## debug mode for fitting a gaussian
mu_gaus = Epoints[int(len(Epoints)/2)]
sig_gaus = np.sqrt(mu_gaus) / 2

if data_method == "gaus":
    MO = 'Gaus'
    print("Debug mode: gaus")
    spec0.make_gaus(mu_gaus, sig_gaus)
else:
    print ('### Get unoscillated reactor spectrum')
    spec0.get_spectrum()
    print ('### Get standard oscillated {0} spectrum'.format(MO))
    spec0.oscillate([MO])

# print ('### Apply detector resopnse')
# spec0.det_response(LY, 'osc'+MO)

print ('### Smear it')
spec0.smear(a=apar, b=bpar, sp='osc' + MO)
# get PDFs of both
print ('### Get corresponding PDFs (norm to 1)')
spec0.get_pdfs()

# data points for likelihood calculation (will be constructed in each fit iteration)
data = np.array([])

######################
## ~~~ Main fit ~~~ ##
######################

# original values of variables
orig = {'a': apar, 'b': bpar, 'th12': osc['th12'], 'th13': osc['th13'],
        'mu': mu_gaus, 'sig': sig_gaus}

AXIS_LABEL = {'a': r'a [MeV$^{1/2}$]', 'b': 'b',
            'th12': r'$\theta_{12}$ [rad]', 'th13': r'$\theta_{13}$ [rad]',
            'mu': r'$\mu$ [MeV]', 'sig': r'$\sigma$ [MeV]'}


def juno_fit(fit_method, sample='sample', plot=False):
    '''
    fit_method (string): method of fitting
        'chi2': minimize chi square ( function chisquare() )
         Other option: 'lkl' -> removed for now

    sample (string): way of constructing toy data
        'sample': sample Nev events based on the PDF
        'asimov': Asimov dataset (scale PDF)
        'gaus': debug mode, fit a Gaussian (mu and sigma defined on top)
        
    plot (bool): True -> plot fit result, profile likelihood scan and contour.
                False -> return fit results
        
    Free variables defined on top
    '''
    
    global data

    print ('Mass ordering:', MO)

    print ('Original parameters:')
    for par in free_params:
        print('\t', par, ':', orig[par])

    # construct histograms with given number of bins assuming a number of events
    print ('### Construct {0} data histogram for {1} events based on PDF'.format(MO,Nev))
    # True -> sample -> get events randomly according to PDF rather than bin the PDF
    sample_bool = {'sample': True, 'asimov': False, 'gaus': False}[sample]
    spec0.histogramize(Nbins, Nev, sample_bool, 'osc' + MO)
    # now this is our generated data
    data = spec0.histos['osc' + MO] #np array

    # -----------------------------------------------------------

    print('Fit method:', {'lkl': 'likelihood', 'chi2': 'Chi^2', 'chi2mx': 'Chi^2 (mx)'}[fit_method])
    func = {'lkl': None, 'chi2': chisquare, 'chi2mx': None}[fit_method]
    errdef = 0.5 if fit_method == 'lkl' else 1.0

    ## JUNO mode
    mn = Minuit(func, th12=osc['th12'], th13=osc['th13'], a=apar, b=bpar, errordef=errdef)
        # error_a=0.001, error_b = 0.01, error_mu = 0.01, error_sig = 0.01,\
        # limit_mu=(0,10), limit_sig=(0,10), limit_a=(0,1), limit_b=(0,1))

    ## debug mode (Gauss)
    # mn = Minuit(func, a=apar, b=bpar, mu=mu_gaus, sig=sig_gaus, errordef=errdef,\
    #     error_a=0.001, error_b = 0.01, error_mu = 0.01, error_sig = 0.01,\
    #     limit_mu=(0,10), limit_sig=(0,10), limit_a=(0,1), limit_b=(0,1))

    # fix parameters that are not free
    for par in fixed_params:
        mn.fixed[par] = True

    print(mn.get_param_states())

    print ('Minimizing...')
    start_time = time.time()
    mn.migrad()
    print ('Running time:', time.time() - start_time, 's')
    print ('Result:', mn.values)
    if 'chi2' in fit_method:
        print('Chi2/Ndof:', mn.fval / len(data[data > 0]))

    # common name for all future plots
    fname = 'result_{0}_{1}_{2}_Nev{3}_smeared_Nbins{4}_{5}'.format(data_method, fit_method, MO, Nev, Nbins, '-'.join(free_params))

    # -----------------------------------------------------------

    # return result
    if not plot:
        return mn.values, mn.fval

    # -----------------------------------------------------------

    # apply these values

    spec = Spectra(Epoints,Ljuno)

    if data_method == "gaus":
        print("### Make Gaus spectrum")
        spec.make_gaus(mn.values['mu'], mn.values['sig'])
    else:
        print ('### Get unoscillated reactor spectrum')
        spec.get_spectrum()
        model_res = Model(mn.values['th12'], osc['th23'], mn.values['th13'], osc['dm2sol'], osc['dm2atm'])
        spec.oscillate(motypes=[MO], model=model_res)

    ## apply detector response
    # LYres = 1./(mn.values['eres']/100)**2
    # spec.det_response(LYres)#, 'osc'+MO)
    # smear
    spec.smear(a=mn.values['a'], b=mn.values['b'], sp='osc'+MO)
    # get resulting PDF
    spec.get_pdfs('osc'+MO)
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc'+MO)

    # add data here for faster plotting
    spec.histos['data'] = data

    # plot the fit
    mp = Plot((10,8))

    LABELS['data'] = fit_label(MO, orig)
    LABELS['osc'+MO] = fit_label('Result', mn.values)
    spec.plot_histos(mp, ['data', 'osc' + MO])

    if 'chi2' in fit_method:
        txt = r'$\chi^2/N_{dof}$' + ' = {0}'.format( round(mn.fval / len(data[data > 0]), 2) )
        mp.ax.text(0.05, 0.9, txt, transform=mp.ax.transAxes, fontsize=15)

    mp.fig.suptitle('Result based on a dataset of {0} events ({1})'.format(Nev, fit_method))
    mp.figure(fname + '.png')

    del mp

    # -----------------------------------------------------------

    # print ('Profile likelihood scan...')
    # mp = Plot((10,8))
    # start_time = time.time()
    # mn.draw_profile("th12", bins=20)
    # print ('Running time:', time.time() - start_time, 's')
    # mp.ax.set_xlabel(r'$\theta_{12}$')
    # mp.ax.set_ylabel('-lnL')
    # mp.pretty()
    # mp.figure(fname + '_scan.png')
    # del mp

    # -----------------------------------------------------------

    #print('Minos...')
    #start_time = time.time()
    #mn.minos(sigma=2)
    #print ('Running time:', time.time() - start_time, 's')
    
    bins_cont = 15
    print('-> Contours with {0} bins'.format(bins_cont))
    # run contours for every pair of variables
    for vars in itertools.combinations(free_params, 2):
    # for vars in [('mu', 'a')]:
        print(vars)
        try:
            mp = Plot((10,8))
            
            start_time = time.time()
            # minimize wrt other parameters
            #mn.draw_mncontour(vars[0], vars[1], numpoints=bins_cont, nsigma=2)
            # keep other parameters fixed (eq to mncontour for 2 free vars)
            mn.draw_contour(vars[0], vars[1], bins=bins_cont, bound=2)
            print ('Running time:', time.time() - start_time, 's')
            # best result of Minuit
            mp.ax.plot(mn.values[vars[0]], mn.values[vars[1]], 'ko')
            # original input
            mp.ax.plot(orig[vars[0]], orig[vars[1]], 'go')

            mp.ax.set_xlabel(AXIS_LABEL[vars[0]])
            mp.ax.set_ylabel(AXIS_LABEL[vars[1]])
            mp.pretty()
            cname = 'contour' + fname.split('result')[1] + '_' + 'VS'.join(vars)
            mp.figure(cname + '.png')
        except:
            print("Minuit fails to do the contour")



###################################
## ~~~ Functions and helpers ~~~ ##
###################################

def chisquare(th12, th13, a, b):
    '''
    Chi square calculation based on JUNO spectrum
    Variables: theta12, theta13 (osc params) and a, b (detector response params)
    '''
## Debug mode -> Gaus    
# def chisquare(mu, sig, a, b):
    spec = Spectra(Epoints,Ljuno)

    ## Debug mode
    # spec.make_gaus(mu, sig)

    ## JUNO mode
    spec.get_spectrum()
    # model with current parameters
    model_current = Model(th12, osc['th23'], th13, osc['dm2sol'], osc['dm2atm'])
    # apply oscillation to original spectrum with given model
    spec.oscillate(motypes=[MO], model=model_current)

    # apply detector response
    # spec.det_response(1./(eres/100.)**2, 'osc')

    # smear
    spec.smear(a=a, b=b, sp='osc'+MO)
    # obtain corresponding PDF
    spec.get_pdfs('osc'+MO)
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc'+MO)
    # now this histo is used as the expected numbers
    exp = spec.histos['osc'+MO]

    # chisquare of each energy point (bin)
    chi2 = 0
    for i in range(len(data)):
        if data[i] == 0:
            continue
        chi2 += ( (data[i]  - exp[i]) / np.sqrt(data[i]) )**2

    return chi2


def fit_label(txt, dct):
    '''
    Generate fit label
    
    txt (string): text for title of the legend (e.g. "Result" or "Original")
    dct (dict): dictionary of the format {variable: value}
    
    Example:
    fit_label("Result", mn.values)
    -> in the legend
    Result
    th12=...
    a=...
    b=...
    
    The values shown are for free fit variables set on top of the code
    '''
    
    label = txt + '\n'

    for par in free_params:
        label += par + ' = {0}\n'.format(round(dct[par], 2))

    return label



if __name__ == '__main__':

    juno_fit(fit_method, data_method, plot=True)


