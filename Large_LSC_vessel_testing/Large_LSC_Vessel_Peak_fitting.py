import scipy.stats
from scipy import special
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
import numpy as np
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL
import pickle
from numba_stats import truncnorm, truncexpon, norm, expon
from uncertainties import ufloat

### limit on LSC energy and limits on NaI coincidence energy
x_limit = (100, 3000)
NaI_limit = [[[0], [2, [1700, 2600]]], [[0], [4, [1300, 1900]]]]

class truncErfc_gen(rv_continuous):
    ''' Class for a truncated complementary error function
    a and b are bounds
    m is shift
    s is scale
    Shift and scale do not affect bounds
    Normalization of pdf did not work for scipy 1.5.2, but does for scipy 1.10.1
    '''

    def _argcheck(self, a, b, m, s):
        return (a < b) and (s > 0)

    def _get_support(self, a, b, m, s):
        return a, b

    def _basIntegErfc(self, a, b):  # Integral between a and b, neither shifted nor scaled: int_ab erfc(x) dx
        if True:
            return (np.exp(-a * a) - np.exp(-b * b)) / np.sqrt(np.pi) + b * special.erfc(b) - a * special.erfc(
                a)  # This analytical method is faster than numerical integration
        else:
            return integrate.quad(special.erfc, a, b)[0]  # Slow normalization

    def _integErfc(self, a, b, m,
                   s):  # Integral of shifted and scaled function but with fixed bounds: int_ab erfc((x-m)/s) dx
        return s * self._basIntegErfc((a - m) / s, (b - m) / s)

    def _pdf(self, x, a, b, m, s):
        return special.erfc((x - m) / s) / self._integErfc(a, b, m, s)  # Faster normalization

    def _cdf(self, x, a, b, m, s):
        return self._integErfc(a, x, m, s) / self._integErfc(a, b, m, s)

truncErfc = truncErfc_gen(name='truncErfc', momtype=1)

### take the coincidence with conditions on NaI energies
def coincidence2(df, Ch_coincidence):
    dfCo = len(Ch_coincidence) * [0]
    k = 1
    if len(Ch_coincidence) < 2:
        return df
    else:
        dfCo[0] = df[df['E' + str(Ch_coincidence[0][0]) + '_AU'] > 0]
        for el in Ch_coincidence[1:]:
            if len(el) == 1:
                dfCo[k] = dfCo[k - 1][dfCo[k - 1]['E' + str(el[0]) + '_AU'] > 0]
                k += 1
            else:
                condition_energy = (dfCo[k - 1]['E' + str(el[0]) + '_AU'] > el[1][0]) & (
                            dfCo[k - 1]['E' + str(el[0]) + '_AU'] < el[1][1])
                dfCo[k] = dfCo[k - 1].loc[condition_energy, :].copy()
                k += 1
        dfCoin = dfCo[-1]
        return dfCoin

### the form of the return is forced by iminuit module
def erf_gauss (x, x_limit0, x_limit1, n_erf, n_gauss, loc, scale, sigma, mu):

    return (n_erf+n_gauss, n_erf * truncErfc.pdf(x, x_limit0, x_limit1, loc, scale) + n_gauss * truncnorm.pdf(x, x_limit0, x_limit1, loc=mu, scale=sigma))

def erf_gauss_fit(dfl, csv, name, plot = True):
    #############################################################################################################################################
    ###  This function is the main function that is used to fit the data. Not fully sure how everything works right now but I will add to this 
    ### 
    ###  dfl: Data frame containing the data that was used. 
    ###  csv: Blank array that will be used to store information about the fit. This is done through using a dictionary.
    ###  name: Name of the file that is used to fit the peak. 
    ###  Plot: Turns on and off the plotting function. 
    ###
    #############################################################################################################################################
    ret = []
    chlst = [0, 2, 4]
    ch = chlst[0]
    otherChLst = [c for c in chlst if not (c == ch)]
    limit_coinc = NaI_limit

    for i, col in zip(range(0, 2), ['r', 'b']) :
        x_limit = (100, 3500)
        dfc = coincidence2(dfl, limit_coinc[i])

        ### apply the energy window
        dfc = dfc[dfc['E' + str(ch) + '_AU'] > x_limit[0]]
        dfc = dfc[dfc['E' + str(ch) + '_AU'] < x_limit[1]]

        ### plot the spectrum histogram
        nBins = 50
        bins, step = np.linspace(x_limit[0], x_limit[1], nBins + 1, retstep=True)
        n1, x1, useless = plt.hist(dfc['E' + str(ch) + '_AU'], bins=bins, density=False, log=False, histtype='step', lw=3, color="k",
                 label='C' + str(ch) + '-' + 'C' + str(otherChLst[i]))
        plt.close()

        ### iminuit at work (fine tuning can be needed for initial value for convergence)
        c = ExtendedUnbinnedNLL(dfc['E' + str(ch) + '_AU'], erf_gauss)
        m = Minuit(c, x_limit0 = x_limit[0], x_limit1 = x_limit[1], n_erf=1500, n_gauss=2000, loc = 1500, scale=500, sigma=200, mu=1500)
        m.limits["n_gauss", "n_erf", "loc", "sigma", "mu"] = (0, None)
        m.limits["scale"] = (1, None)
        m.simplex()
        m.migrad()
        m.hesse()

        fig, ax = plt.subplots(1, 1)

        x_limit = (max(0, m.values['mu'] - 4 * m.values['sigma']), min(4000, m.values['mu'] + 4 * m.values['sigma']))
        dfc = dfc[dfc['E' + str(ch) + '_AU'] > x_limit[0]]
        dfc = dfc[dfc['E' + str(ch) + '_AU'] < x_limit[1]]
        bins, step = np.linspace(x_limit[0], x_limit[1], nBins + 1, retstep=True)
        n1, x1, useless = ax.hist(dfc['E' + str(ch) + '_AU'], bins=bins, density=False, log=False, histtype='step',
                                  lw=3, color="k", label='C' + str(ch) + '-' + 'C' + str(otherChLst[i]))
        c = ExtendedUnbinnedNLL(dfc['E' + str(ch) + '_AU'], erf_gauss)
        m = Minuit(c, x_limit0=x_limit[0], x_limit1=x_limit[1], n_erf=900, n_gauss=2500, loc=m.values['loc'], scale=m.values['scale'], sigma=m.values['sigma'],
                   mu=m.values['mu'])
        m.limits["n_gauss", "n_erf", "loc", "sigma", "mu"] = (0, None)
        m.limits["scale"] = (1, None)
        m.fixed['x_limit0'] = True
        m.fixed['x_limit1'] = True
        m.simplex()
        m.migrad()
        m.hesse()

        xx = np.linspace(x_limit[0], x_limit[1], 1000)
        yy = step * erf_gauss(xx, x_limit[0], x_limit[1], m.values['n_erf'], m.values['n_gauss'], m.values['loc'],
                              m.values['scale'],
                              m.values['sigma'], m.values['mu'])[1]
        ax.plot(xx, yy, label="fit", color=col, lw=3)
        ax.set_xlabel("energy (AU)")
        ax.set_ylabel("count per bins")
        ax.legend(["data", "fit"])
        ax.set_title(name + "_CH" + str(otherChLst[i]) + ".png")

        ### chisquare calculation
        bins_count, u = np.histogram(dfc['E' + str(ch) + '_AU'], bins=bins)
        xxx = 0.5 * (bins[1:] + bins[:-1])
        yyy = step * erf_gauss(xxx, x_limit[0], x_limit[1], m.values['n_erf'], m.values['n_gauss'], m.values['loc'], m.values['scale'], m.values['sigma'],
                               m.values['mu'])[1]

        delta2 = np.sum([(bbins - yyyy) ** 2 / yyyy for bbins, yyyy in zip(bins_count, yyy) if yyyy != 0])
        df = nBins - m.npar
        p_value = chi2.sf(delta2, df)

        ### printing the parameters and additional information
        textstr = '\n'.join((r'n_gauss' + ' = ' + str(ufloat(m.values['n_gauss'], m.errors['n_gauss'])),
                             r'n_erf' + ' = ' + str(ufloat(m.values['n_erf'], m.errors['n_erf'])),
                             r'loc' + ' = ' + str(ufloat(m.values['loc'], m.errors['loc'])),
                             r'scale' + ' = ' + str(ufloat(m.values['scale'], m.errors['scale'])),
                             r'$\sigma$' + ' = ' + str(ufloat(m.values['sigma'], m.errors['sigma'])),
                             r'FWHM' + ' = ' + str(
                                 ufloat(2 * np.sqrt(2 * np.log(2)) * m.values['sigma'],
                                        2 * np.sqrt(2 * np.log(2)) * m.errors['sigma']))
                             ,
                             r'$\mu$' + ' = ' + str(ufloat(m.values['mu'], m.errors['mu'])),
                             r'$\chi^2$' + ' = ' + str("%.1f" % delta2),
                             r'freedom' + ' = ' + str(nBins - 5),
                             r'p_value' + ' = ' + str("%.1e" % p_value)))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=9, verticalalignment='top', bbox=props)

        for ext in ['pdf', 'png'] :
            plt.savefig("Figures\\LSC_campaign\\Thresholded\\" + ext + "\\" + name + "_CH" + str(otherChLst[i]) + '.' +  ext)
        if plot :
            plt.show()
        plt.close()
        if i == 0 :
            n_gauss = m.values['n_gauss'],
            err_n_gauss = m.errors['n_gauss']
            peak = m.values['mu']
            err_peak = m.errors['mu']
            sig_save = m.values['sigma']
            err_sig_save = m.errors['sigma']
            FWHM = 2 * np.sqrt(2 * np.log(2)) * m.values['sigma']
            err_FWHM = 2 * np.sqrt(2 * np.log(2)) * m.errors['sigma']
            conv = m.valid
            pval = p_value
    if name.split('_')[1][-2:] != 'mL':
        print(name.split('_')[1][-2:] != "mL")
        salt = name.split('_')[1]
    else:
        salt = 'No salt'

    dens = 0.
    for i in range(len(name.split('_'))) :

        if name.split('_')[i][-2:] == 'mL':
            if dens == 0 :
                dens = 1
            dens *= float(name.split('_')[i][:-2])
        if name.split('_')[i][-1] == 'M' :
            if dens == 0 :
                dens = 1
            dens *= float(name.split('_')[i][:-1])
    dens *= 1.95
    csv.append({'Name' : name,
                'LSC' : name.split('_')[0],
                'Salt' : salt,
                'density' : dens,
                'Temp_mean' : 0,
                'Temp_strd_dev' : 0,
                'gauss_count_1' : n_gauss,
                'gauss_count_2' : m.values['n_gauss'],
                'err_gauss_count_1' : err_n_gauss,
                'err_gauss_count_2' : m.errors['n_gauss'],
                'sigma_1' : sig_save,
                'sigma_2' : m.values['sigma'],
                'err_sigma_1' : err_sig_save,
                'err_sigma_2' : m.errors['sigma'],
                'Peak_value_1' : peak,
                'Peak_value_2' : m.values['mu'],
                'err_Peak_value_1' : err_peak,
                'err_Peak_value_2' : m.errors['mu'],
                'FWHM_1' : FWHM,
                'FWHM_2' : 2 * np.sqrt(2 * np.log(2)) * m.values['sigma'],
                'err_FWHM_1' : err_FWHM,
                'err_FWHM_2' : 2 * np.sqrt(2 * np.log(2)) * m.errors['sigma'],
                'p_value_1' : pval*100,
                'p_value_2' : p_value*100,
                'convergence_1' : conv,
                'convergence_2' : m.valid
                })
    return()