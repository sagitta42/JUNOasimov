### Needs to be run with python3 because iminuit in python2 has a bug
from juno_osc import *
import math
import time
import os

from iminuit import Minuit

###########################
## ~~~ Initial setup ~~~ ##
###########################

defaults = {
'fit_method': 'chi2',
'Nbins': 250,
'data_method': 'sample'
}

def user_dict():
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


if len(sys.argv) == 1:
    print ('Available variables')
    print ('\t Nev: number of events in the simulated/Asimov sample')
    print ('\t fit_method: chi2 or lkl (default chi2)')
    print ('\t data_method: sample or Asimov (default sample)')
    print ('\t Nbins: number of bins (default 250)')
    # print('Syntax: python juno_fit.py Nev [fit_method]')
    # print('Nev: number of events in the Asimov dataset (int)')
    # print('fit_method: lkl or chi2 (default lkl)')
    sys.exit(1)

user_input = user_dict()

Nev = int(user_input['Nev'])
# Nev = int(sys.argv[1])
fit_method = user_input['fit_method']
data_method = user_input['data_method']
# fit_method = sys.argv[2] if len(sys.argv) > 2 else 'lkl'

# Epoints = np.linspace(0.1,10,10000)
Epoints = np.linspace(0.1,10,6000)
# Epoints = np.linspace(0.1,10,5000)

# JUNO baseline
Ljuno = 53 # km

# MO = 'IO'
MO = 'NO'

# ~2 - 10 MeV, bin = 30 keV @ 1 MeV -> 250 bins
Nbins = int(user_input['Nbins'])
# Nbins = 250
# number of bins
# Nbins = 2500
# Nbins = 1500
# Nbins = 1000
# Nbins = 500
# Nbins = 400
# Nbins = 300 # -> 33 keV bin

# Ereso = 3 # %, assumed resolution
LY = 1100 # p.e./MeV @ 1 MeV
Ereso = 1/np.sqrt(LY)*100 # %, assumed resolution

# ------------------------------------------------------------------


# standard oscillation model is already initialized in juno_osc.py as model_default

# JUNO spectra for E and L
spec0 = Spectra(Epoints,Ljuno)
# after this the energy will be changed to be above IBD threshold

print ('### Get unoscillated reactor spectrum')
spec0.get_spectrum()
print ('### Get standard oscillated {0} spectrum'.format(MO))
spec0.oscillate([MO])
print ('### Apply detector resopnse')
spec0.det_response(LY, 'osc'+MO)
# OINK
print ('### Smear it')
spec0.smear(sp='osc' + MO)
# spec0.smear(Ereso, 'osc' + MO)
# get PDFs of both
print ('### Get corresponding PDFs (norm to 1)')
spec0.get_pdfs()

# construct histograms with given number of bins assuming a number of events
# print ('### Construct {0} data histogram for {1} events based on PDF'.format(MO,Nev))
# # True -> sample Nev
# spec.histogramize(Nbins, Nev, True, 'osc' + MO)
# spec0.histogramize(Nbins, Nev, True, 'osc' + MO)
# mp = Plot((10,8))
# spec0.plot_spectra(mp)
# # spec0.plot_histos(mp)
# mp.figure()
# sys.exit()

# data points for likelihood calculation
data = np.array([])
# data = spec.histos['osc' + MO] #np array

######################
## ~~~ Main fit ~~~ ##
######################

def juno_fit(fit_method, sample='sample', plot=False):
    '''
    fit_method (string): method of fitting
        'chi2': minimize chi square ( function chisquare() )
         'lkl': minimize negative log likelihood ( function negloglkl() )

    sample (string): way of constructing toy data
        'sample': sample Nev events based on the PDF
        'asimov': Asimov dataset (scale PDF)

    Note: pay attention to free variables, currently the implementation is not flexible
        (chi2 and lkl functions have to be redefined)
    '''
    global data

    print ('Mass ordering:', MO)

    print ('Original parameters:')
    print('\t Theta13:', osc['th13'], 'rad')
    # print('\t Theta12:', osc['th12'], 'rad')
    print('\t Energy resolution @ 1 MeV:', Ereso, '%')
    #, osc['dm2sol'])

    # construct histograms with given number of bins assuming a number of events
    print ('### Construct {0} data histogram for {1} events based on PDF'.format(MO,Nev))
    # True -> sample -> get events randomly according to PDF rather than bin the PDF
    sample_bool = {'sample': True, 'asimov': False}[sample]
    spec0.histogramize(Nbins, Nev, sample_bool, 'osc' + MO)
    # spec0.histogramize(Nbins, Nev, True, 'osc' + MO)
    data = spec0.histos['osc' + MO] #np array


    # -----------------------------------------------------------
    print('Fit method:', {'lkl': 'likelihood', 'chi2': 'Chi^2', 'chi2mx': 'Chi^2 (mx)'}[fit_method])
    func = {'lkl': negloglkl, 'chi2': chisquare, 'chi2mx': chisquare_mx}[fit_method]
    # fit only Ereso
    errdef = 0.5 if fit_method == 'lkl' else 1.0
    mn = Minuit(func, th13=osc['th13'], eres=Ereso, errordef=errdef)
    # mn = Minuit(func, th12=osc['th12'], eres=Ereso, errordef=errdef)
    # mn = Minuit(func, eres=Ereso, errordef=errdef)
    # mn = Minuit(func, eres=Ereso, errordef=errdef, limit_eres=(0,10)) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(func, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(chisquare, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(negloglkl, th12=osc['th12'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    # mn = Minuit(negloglkl, th12=osc['th12'], dm21=osc['dm2sol'], eres=Ereso, errordef=0.5) # 0.5 for negative log likelihood function, 1 for least squares
    print ('Minimizing...')
    start_time = time.time()
    mn.migrad()
    print ('Running time:', time.time() - start_time, 's')
    print ('Result:', mn.values)
    if 'chi2' in fit_method:
        print('Chi2/Ndof:', mn.fval / len(data[data > 0]))

    # common for all future plots
    fname = 'result_{0}_{1}_{2}_Nev{3}_smeared_Nbins{4}'.format(data_method, fit_method, MO, Nev, Nbins)

    # -----------------------------------------------------------

    # return result
    if not plot:
        return mn.values, mn.fval


    # -----------------------------------------------------------

    # apply these values
    # fit only Ereso
    # model_res = Model(osc['th12'], osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])
    # model_res = Model(mn.values['th12'], osc['th23'], osc['th13'], mn.values['dm21'], osc['dm2atm'])
    # construct oscillated spectrum
    spec = Spectra(Epoints,Ljuno)
    print ('### Get unoscillated reactor spectrum')
    spec.get_spectrum()
    model_res = Model(osc['th12'], osc['th23'], mn.values['th13'], osc['dm2sol'], osc['dm2atm'])
    # model_res = Model(mn.values['th12'], osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])
    spec.oscillate(motypes=['NO'], model=model_res)
    # spec.oscillate(motypes=['NO'])#, model=model_default) # here we need the original spectrum back
    # spec.oscillate(motypes=[None], model=model_default) # here we need the original spectrum back
    # spec.oscillate(motypes=[None], model=model_res) # here we need the original spectrum back
    # apply detector response
    LYres = 1./(mn.values['eres']/100)**2
    spec.det_response(LYres)#, 'osc'+MO)
    # smear
    # spec.smear(Ereso, 'osc'+MO)
    spec.smear(sp='osc'+MO)
    # spec.smear(mn.values['eres'], 'osc'+MO)
    # spec.smear(3, 'osc')
    # get resulting PDF
    spec.get_pdfs('osc'+MO)
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc'+MO)

    # add data here for faster plotting
    spec.histos['data'] = data

    # plot the fit
    mp = Plot((10,8))

    LABELS['data'] = fit_label(MO, {'eres':Ereso, 'th13': osc['th13']})
    # LABELS['data'] = fit_label(MO, {'eres':Ereso, 'th12': osc['th12']})
    # LABELS['data'] = fit_label(MO, {'eres':Ereso})
    LABELS['osc'+MO] = fit_label('Result', mn.values)
    spec.plot_histos(mp, ['data', 'osc' + MO])
    # spec.plot_spectra(mp)

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

    print('Contour...')
    mp = Plot((10,8))
    start_time = time.time()

    ## keep other params fixed
    # draw_contour(self, x, y, bins=50, bound=2, **deprecated_kwargs
    ## minimize wrt other params
    # draw_mncontour(self, x, y, nsigma=2, numpoints=100)
    ## equivalent in case of only 2 params

    mn.draw_contour("eres", "th13", bins=10, bound=2)
    # mn.draw_mncontour("eres", "th13", numpoints=10, nsigma=2)
    # mn.draw_mncontour("th12", "dm21", numpoints=10, nsigma=2)
    print ('Running time:', time.time() - start_time, 's')
    mp.ax.plot(mn.values['eres'], mn.values['th13'], 'rx')
    # mp.ax.plot(mn.values['th12'], mn.values['dm21'], 'rx')
    mp.ax.set_xlabel('Energy resolution (%)')
    # mp.ax.set_xlabel(r'$\theta_{12}$')
    mp.ax.set_ylabel(r'$\theta_{13}$ (rad)')
    # mp.ax.set_ylabel(r'$\Delta m^2_{21}$')
    # mp.pretty(stretch='float')
    mp.pretty()
    cname = 'contour' + fname.split('result')[1]
    mp.figure(cname + '.png')

#    # save results
#    df.to_csv(fname + '.csv', index=False)
#    print '-->', fname + '.csv'



    # likelihood with original values
#    lkl_asimov = negloglkl([ osc['th12'], osc['dm2sol'] ])

    # is actually -ln lambda A
#    lambdaA = lkl_fit / lkl_asimov

#    qmu = 2*lambdaA

########################################
## ~~~ Perform fit multiple times ~~~ ##
########################################

def multiple_fits(Nfits):
    '''
    Perform multiple fits and plot the distribution of results
    '''

    folder = 'results/'
    figname = "results_{0}_{1}_Nev{2}_Nbins{3}_Nfits{4}".format(fit_method, MO, Nev, Nbins, Nfits)
    globname = "results_{0}_{1}_Nev{2}_Nbins{3}.csv".format(fit_method, MO, Nev, Nbins)

    if not os.path.exists(folder + globname):
        df = pd.DataFrame(columns=['Ereso', 'th13', 'chi2ndof'])
        # df = pd.DataFrame(columns=['Ereso', 'th12', 'chi2ndof'])
        # df = pd.DataFrame(columns=['res', 'chi2ndof'])
        df.to_csv(folder + globname, index=False)

    ## read total collection of results
    if Nfits == 0:
        print('[Reading global results]')
        df = pd.read_csv(folder + globname)

    ## read the results if saved
    elif os.path.exists(folder + figname + '.csv'):
        print('[Reading saved results]')
        df = pd.read_csv(folder + figname + '.csv')

    ## perform the fits and collect the results
    else:
        print('Performing fits...')
        # results = []
        results = {'Ereso':[], 'th13':[], 'chi2ndof':[]}
        # results = {'Ereso':[], 'th12':[], 'chi2ndof':[]}
        # chi2_res = []

        for i in range(Nfits):
            print('~~~~~~~~~~~~~~~~~~~~~~~~', i+1)
            res, chi2res = juno_fit(fit_method, plot=False)
            results['Ereso'].append(res['eres'])
            results['th13'].append(res['th13'])
            # results['th12'].append(res['th12'])
            # results.append(res['eres'])
            results['chi2ndof'].append(chi2res/len(data[data > 0]))
            # chi2_res.append(chi2res/len(data[data > 0]))


        ## save results
        df = pd.DataFrame(results)
        # df = pd.DataFrame({'res': results, 'chi2ndof': chi2_res})
        dfglob = pd.read_csv(folder + globname)
        dfglob = pd.concat([dfglob, df], ignore_index=True)
        print(dfglob)
        df.to_csv(folder + figname + '.csv', index=False)
        print('-->', folder + figname + '.csv')
        dfglob.to_csv(folder + globname, index=False)
        print('-->', folder + globname)


    ## plot the results
    # nplots = 2 if 'chi2' in fit_method else 1
    # figsize = (10,8) if 'chi2' in fit_method else (10,4)
    nplots = 'paramspace'
    figsize = (12,12)
    mp = Plot(figsize, nplots, sharex=True)

    # ----------------------------------------

    # mp.axes[0].hist(df['res'], bins=30, histtype='step', linewidth=1.5)
    # # mp.axes[0].set_xlabel('Energy resolution (%)')
    # # mp.axes[1].hist(chi2_res, bins=40, histtype='step', linewidth=1.5)
    # if 'chi2' in fit_method:
    #     mp.axes[1].plot(df['res'], df['chi2ndof'], '.', ms=10)
    #     mp.axes[1].set_ylabel(r'$\chi^2/N_{dof}$')
    #     mp.axes[1].set_xlabel('Energy resolution (%)')
    # # mp.axes[1].set_xlabel(r'$\chi^2/N_{dof}$')

    # ----------------------------------------

    n_bins = 20
    mp.axes[0].scatter(df['Ereso'], df['th13'])#, markersize=10)
    mp.axes[0].scatter([Ereso], [osc['th13']], color='r')
    mp.axes[1].hist(df['Ereso'], bins=n_bins, histtype='step')
    mp.axes[1].axvline(Ereso, color='r', linestyle='--')
    mp.axes[2].hist(df['th13'], bins=n_bins, histtype='step', orientation='horizontal')
    mp.axes[2].axhline(osc['th13'], color='r', linestyle='--')#, orientation='horizontal')
    mp.axes[0].set_xlabel('Energy resolution (%)')
    mp.axes[0].set_ylabel(r'$\theta_{13}$ (rad)')

    R = np.corrcoef(df['Ereso'], df['th13'])
    R = R[0][1]

    mp.axes[0].text(0.05, 0.95, 'R = {0}'.format(round(R,2)), transform=mp.axes[0].transAxes, fontsize=13)


    # ----------------------------------------


    nfits = len(df) if Nfits == 0 else Nfits
    mp.axes[1].set_title('Results from {0} fits on {1} events'.format(nfits, Nev))
    # mp.fig.suptitle('Results from {0} fits on {1} events'.format(nfits, Nev))
    mp.pretty(large=0)

    mp.figure(figname + '.png')





###################################
## ~~~ Functions and helpers ~~~ ##
###################################


# def negloglkl(eres):
def negloglkl(th12,eres):
# def negloglkl(th12,dm21,eres):
    '''
    Variables:
        oscillation parameters th12 and dm^2_21
        Energy resolution @ 1 MeV
    '''

    # model with current parameters
    # dm21 = osc['dm2sol']
    # model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    spec = Spectra(Epoints,Ljuno)
    spec.get_spectrum()

    # apply oscillation to original spectrum with given model
    model_current = Model(th12, osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])
    spec.oscillate(motypes=[None], model=model_current)
    # spec.oscillate(motypes=[None])
    # apply detector response
    spec.det_response(1./(eres/100.)**2, 'osc')
    # spec.det_response(LY, 'osc')
    # smear
    spec.smear(sp='osc')
    # spec.smear(eres, 'osc')
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


# def chisquare(eres):
def chisquare(th13,eres):
# def chisquare(th12,eres):
    # model with current parameters
    # dm21 = osc['dm2sol']
    # model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    spec = Spectra(Epoints,Ljuno)
    spec.get_spectrum()

    # apply oscillation to original spectrum with given model
    model_current = Model(osc['th12'], osc['th23'], th13, osc['dm2sol'], osc['dm2atm'])
    # model_current = Model(th12, osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])
    spec.oscillate(motypes=[None], model=model_current)
    # spec.oscillate(motypes=[None])
    # apply detector response
    spec.det_response(1./(eres/100.)**2, 'osc')
    # spec.det_response(LY, 'osc')
    # smear
    spec.smear(sp='osc')
    # spec.smear(eres, 'osc')
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
    # print ('chi2:', chi2)
    return chi2



def chisquare_mx(th12,eres):
    # model with current parameters
    dm21 = osc['dm2sol']
    model_current = Model(th12, osc['th23'], osc['th13'], dm21, osc['dm2atm'])

    # apply oscillation to original spectrum with given model
    spec.oscillate(motypes=[None], model=model_current)
    # apply detector response
    # spec.det_response(LY, 'osc')
    # smear
    spec.smear(eres, 'osc')
    # obtain corresponding PDF
    spec.get_pdfs('osc')
    # construct a histogram based on it
    # False -> norm by Nev
    spec.histogramize(Nbins, Nev, False, 'osc')
    # now this histo is used as the expected numbers
    exp = spec.histos['osc']

    # ignore zero bins
    data0 = data[data > 0]
    exp0 = exp[data > 0]
    # vector data - exp
    diff = data0 - exp0
    # sigma^2 matrix
    V = np.diag(data0)

    chi2 = np.transpose(diff).dot(np.linalg.inv(V)).dot(diff)
    return chi2


def cov_mx(M):
    ''' Calculate covariance matrix based on M samples '''

    eres = 3 # %
    Nbins = len(data)
    vres = np.zeros((Nbins, Nbins))

    for i in range(M):
        model_current = Model(osc['th12'], osc['th23'], osc['th13'], osc['dm2sol'], osc['dm2atm'])

        # apply oscillation to original spectrum with given model
        spec.oscillate(motypes=[None], model=model_current)
        # apply detector response
        # spec.det_response(LY, 'osc')
        # smear
        spec.smear(eres, 'osc')
        # obtain corresponding PDF
        spec.get_pdfs('osc')
        # construct a histogram based on it (sample)
        # False -> norm by Nev
        spec.histogramize(Nbins, Nev, False, 'osc')
        # now this histo is used as the expected numbers
        exp = spec.histos['osc']

        # use "data" as nominal
        for j in range(Nbins):
            for k in range(Nbins):
                if j == k: continue
                vres[j][k] += (data[j] - exp[j]) * (data[k] - exp[k])

    vres = vres / M
    return vres




def fit_label(txt, dct):
    label = txt + '\n'
    # for th in ['']:
    th = 13
    label += '$\\theta_{{{0}}}$'.format(th) + ' = {0}\n'.format(round(dct['th{0}'.format(th)], 6))
    # label += r'$\theta_{12}$' + ' = {0}\n'.format(round(dct['th12'], 6))
    # dmlabel = 'dm2sol' if 'dm2sol' in dct else 'dm21'
    # label += r'$\Delta m^2_{21}$' + ' = {0}\n'.format(round(dct[dmlabel], 6))
    eres = round(dct['eres'],2)
    label += 'Ereso = {0}%'.format(eres)
    # label = r'{}'.format(label)
    return label


def compare_pdfs():
    pdfs = {}

    # resos = [0.1]
    resos = [3]
    # resos = [0.1, 3, 30]

    spec = Spectra(Epoints,Ljuno)

    # get reactor spectrum
    spec.get_spectrum()
    for reso in resos:
        # construct oscillated spectrum
        spec.oscillate()#motypes=['NO'])#, model=model_default) # here we need the original spectrum back
        # spec.oscillate(motypes=[None], model=model_res) # here we need the original spectrum back
        # apply detector response
        # spec.det_response(LY, 'oscNO')
        # OINK
        spec.smear(reso, 'oscNO')
        # get resulting PDF
        spec.get_pdfs('oscNO')
        # construct a histogram based on it
        # False -> norm by Nev
        # spec.histogramize(Nbins, Nev, False, 'oscNO')
        # smear histo OINK
        # spec.smear(reso, 'osc')
        # pdfs[reso] = spec.histos['oscNO']

        # pdfs[reso] = spec.spectra['oscNO']
        pdfs[reso] = spec.pdfs['oscNO']


    # plot the fit
    mp = Plot((10,8))

    for reso in resos:
        mp.ax.plot(spec.E, pdfs[reso], label='{0}%'.format(reso), linewidth=7)
        # mp.ax.step(spec.Ebinned, pdfs[reso], label='{0}%'.format(reso))

    # mp.legend(out=False, ncol=1)
    mp.pretty()
    mp.figure('compare_pdfs.png')







if __name__ == '__main__':

    # perform one fit and plot (uncomment plotting)
    juno_fit(fit_method, data_method, plot=True)

    # perform multiple fits and plot the distribution of results
    # multiple_fits(int(sys.argv[3]))

    # compare 3% and 30% reso PDFs
    # compare_pdfs()



    # print(cov_mx(100))



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
