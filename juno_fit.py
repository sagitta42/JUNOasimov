### Needs to be run with python3 because iminuit in python2 has a bug

from juno_osc import *
import math
import time

from iminuit import Minuit

# number of events in the dataset
if len(sys.argv) == 1:
    print('Syntax: python juno_fit.py Nev [fit_method]')
    print('Nev: number of events in the Asimov dataset (int)')
    print('fit_method: lkl or chi2 (default lkl)')
    sys.exit(1)

Nev = int(sys.argv[1])
fit_method = sys.argv[2] if len(sys.argv) > 2 else 'lkl'
# Nev = 1000

#Epoints = np.linspace(0.1,10,10000)
Epoints = np.linspace(0.1,10,5000)

# JUNO baseline
Ljuno = 53 # km

# MO = 'IO'
MO = 'NO'

# number of bins
Nbins = 250
# Nbins = 500

Ereso = 3 # %, assumed resolution

# standard oscillation model is already initialized in juno_osc.py as model_default

# JUNO spectra for E and L
spec = Spectra(Epoints,Ljuno)
# after this the energy will be changed to be above IBD threshold

print ('### Get unoscillated reactor spectrum')
spec.get_spectrum()
print ('### Get standard oscillated {0} spectrum'.format(MO))
spec.oscillate([MO])
print ('### Smear it')
spec.smear(Ereso, 'osc' + MO)
# get PDFs of both
print ('### Get corresponding PDFs (norm to 1)')
spec.get_pdfs()

# construct histograms with given number of bins assuming a number of events
print ('### Construct {0} data histogram for {1} events based on PDF'.format(MO,Nev))
# True -> sample Nev
spec.histogramize(Nbins, Nev, True, 'osc' + MO)

# data points for likelihood calculation
data = spec.histos['osc' + MO]


def juno_fit(fit_method):
    '''
    Fit of Asimov dataset
    Due to smearing each fit takes around 2 minutes
    Current setting: theta_12 and dm^2_21 are free variables
    '''

    print ('Mass ordering:', MO)

    print ('Original parameters:')
    print('\t Theta12:', osc['th12'], 'rad')
    print('\t Energy resolution @ 1 MeV:', Ereso, '%')
    #, osc['dm2sol'])

    # -----------------------------------------------------------
    print('Fit method:', {'lkl': 'likelihood', 'chi2': 'Chi^2'}[fit_method])
    func = {'lkl': negloglkl, 'chi2': chisquare}[fit_method]
    mn = Minuit(func, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(chisquare, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(negloglkl, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(negloglkl, th12=osc['th12'], dm21=osc['dm2sol'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    print ('Minimizing...')
    start_time = time.time()
    mn.migrad()
    print ('Running time:', time.time() - start_time, 's')
    print ('Result:', mn.values)
    if fit_method == 'chi2':
        print('Chi2/Ndof:', mn.fval / len(data))

    # common for all future plots
    fname = 'result_{0}_{1}_Nev{2}_smeared_Nbins{3}'.format(fit_method, MO, Nev, Nbins)


    # -----------------------------------------------------------

    # apply these values
    model_res = Model(mn.values['th12'], osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])
    # model_res = Model(mn.values['th12'], osc['th23'], osc['th13'], mn.values['dm21'], osc['dm2atm'])
    # construct oscillated spectrum
    spec.oscillate(motypes=[None], model=model_res) # here we need the original spectrum back
    # smear
    spec.smear(mn.values['eres'], 'osc')
    # spec.smear(3, 'osc')
    # get resulting PDF
    spec.get_pdfs('osc')
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc')

    # plot the fit
    mp = Plot((10,8))

    LABELS['osc' + MO] = fit_label(MO, osc)
    LABELS['osc'] = fit_label('Result', mn.values)
    spec.plot_histos(mp, ['osc' + MO, 'osc'])
    # spec.plot_spectra(mp)

    if fit_method == 'chi2':
        txt = r'$\chi^2/N_{dof}$' + ' = {0}'.format( round(mn.fval / len(data), 2) )
        mp.ax.text(0.05, 0.9, txt, transform=mp.ax.transAxes, fontsize=15)

    mp.fig.suptitle('Result based on a dataset of {0} events ({1})'.format(Nev, fit_method))
    mp.figure(fname + '.png')
    #
    # del mp

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

    # print('Contour...')
    # mp = Plot((10,8))
    # start_time = time.time()
    # mn.draw_mncontour("th12", "dm21", numpoints=10, nsigma=2)
    # print ('Running time:', time.time() - start_time, 's')
    # mp.ax.plot(mn.values['th12'], mn.values['dm21'], 'rx')
    # mp.ax.set_xlabel(r'$\theta_{12}$')
    # mp.ax.set_ylabel(r'$\Delta m^2_{21}$')
    # mp.pretty(stretch='float')
    # mp.figure(fname + '_contour.png')

#    # save results
#    df.to_csv(fname + '.csv', index=False)
#    print '-->', fname + '.csv'



    # likelihood with original values
#    lkl_asimov = negloglkl([ osc['th12'], osc['dm2sol'] ])

    # is actually -ln lambda A
#    lambdaA = lkl_fit / lkl_asimov

#    qmu = 2*lambdaA






def negloglkl(th12,eres):
# def negloglkl(th12,dm21,eres):
    '''
    Variables:
        oscillation parameters th12 and dm^2_21
        Energy resolution @ 1 MeV
    '''

    # model with current parameters
    dm21 = osc['dm2sol']
    model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    # apply oscillation to original spectrum with given model
    spec.oscillate(motypes=[None], model=model_current)
    # smear
    spec.smear(eres, 'osc')
    # obtain corresponding PDF
    spec.get_pdfs('osc')
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc')
    # now this histo is used as the expected numbers
    exp = spec.histos['osc']

    # log likelihood of each energy point (bin)
    loglkl = 0
    for i in range(len(data)):
        if data[i] == 0 or exp[i] == 0:
            continue
        loglkl += data[i] * math.log(exp[i]) - exp[i] - math.log(np.math.factorial(data[i]))
    # loglkl = data * array_log(pdf) - pdf - array_log(array_fac(data))
    # replace NaN with zero since it shouldn't count in the sum
    # loglkl = np.nan_to_num(loglkl)
    return -loglkl
    # return -sum(loglkl)


def chisquare(th12,eres):
    # model with current parameters
    dm21 = osc['dm2sol']
    model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    # apply oscillation to original spectrum with given model
    spec.oscillate(motypes=[None], model=model_current)
    # smear
    spec.smear(eres, 'osc')
    # obtain corresponding PDF
    spec.get_pdfs('osc')
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc')
    # now this histo is used as the expected numbers
    exp = spec.histos['osc']

    # chisquare of each energy point (bin)
    chi2 = 0
    for i in range(len(data)):
        if data[i] == 0:
            continue
        chi2 += ( (data[i]  - exp[i]) / np.sqrt(data[i]) )**2
    # loglkl = data * array_log(pdf) - pdf - array_log(array_fac(data))
    # replace NaN with zero since it shouldn't count in the sum
    # loglkl = np.nan_to_num(loglkl)
    return chi2



def chisquare_mx(th12,eres):
    # model with current parameters
    dm21 = osc['dm2sol']
    model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    # apply oscillation to original spectrum with given model
    spec.oscillate(motypes=[None], model=model_current)
    # smear
    spec.smear(eres, 'osc')
    # obtain corresponding PDF
    spec.get_pdfs('osc')
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc')
    # now this histo is used as the expected numbers
    exp = spec.histos['osc']

    

    # chisquare of each energy point (bin)
    chi2 = 0
    for i in range(len(data)):
        if data[i] == 0:
            continue
        chi2 += ( (data[i]  - exp[i]) / np.sqrt(data[i]) )**2
    # loglkl = data * array_log(pdf) - pdf - array_log(array_fac(data))
    # replace NaN with zero since it shouldn't count in the sum
    # loglkl = np.nan_to_num(loglkl)
    return chi2


def fit_label(txt, dct):
    label = txt + '\n'
    # for th in ['']:
    # label += '$\\theta_{{{0}}}$'.format(th) + ' = {0}\n'.format(round(dct['th{0}'.format(th)], 6))
    label += r'$\theta_{12}$' + ' = {0}\n'.format(round(dct['th12'], 6))
    # dmlabel = 'dm2sol' if 'dm2sol' in dct else 'dm21'
    # label += r'$\Delta m^2_{21}$' + ' = {0}\n'.format(round(dct[dmlabel], 6))
    eres = round(dct['eres'],2) if 'eres' in dct else Ereso
    label += 'Ereso = {0}%'.format(eres)
    # label = r'{}'.format(label)
    return label


juno_fit(fit_method)



# ----------------------------------------------------------------------------------------------



# def array_fac(arr):
#     ''' Array factorial '''
#     print arr
#     return np.array([np.math.factorial(x) for x in arr])
#     # return np.array([np.NaN if x == np.NaN else np.math.factorial(x) for x in arr])
#
#
# def array_log(arr):
#     '''
#     Log of array. Numpy cannot get log of large numbers
#     so have to use math.log
#     '''
#
#     return np.array([math.log(x) for x in arr])
