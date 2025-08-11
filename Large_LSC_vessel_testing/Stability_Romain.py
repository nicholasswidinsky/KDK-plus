from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from scipy.stats import *
from scipy import special
import numpy as np
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL
from numba_stats import truncnorm
from uncertainties import ufloat
import os

def stab_analysis(dfc, slice_number, slice_total, name, log):
    if not os.path.exists('\\Boulot\\Pycharm\\Projects\\KDK+\\Figures\\Stability\\' + name):
        os.makedirs('\\Boulot\\Pycharm\\Projects\\KDK+\\Figures\\Stability\\' + name)
    directory = ('\\Boulot\\Pycharm\\Projects\\KDK+\\Figures\\Stability\\' + name + '\\' +
                 str(slice_number + 1) + '_out_of_' + str(slice_total))
    if not os.path.exists(directory):
        os.makedirs(directory)
    print(directory)
    dfc = dfc.copy()
    nbins_energy = 100
    nbins_time = 350
    bins_energy = np.linspace(50, 4000, nbins_energy + 1)
    bins_time = np.linspace(-300, 300, nbins_time + 1)
    dfc['deltaT'] = [(k - i) / 1000 for k, i in zip(dfc['T0_ps'], dfc['T4_ps'])]

    ### automated coincidence timing cut
    curve, useless, useless = plt.hist(dfc['deltaT'], bins=bins_time, log=log, density=False,
                                       histtype='step', lw=3, color='k')
    test_peaks, useless = find_peaks(curve, prominence=1000)
    k = 0
    while curve[k] < 0.01 * max(curve):
        k += 1
    time_cut_min = bins_time[k]
    while curve[k] > 0.01 * max(curve):
        k += 1
    time_cut_max = bins_time[k]
    plt.close()
    dfc_time = dfc[(dfc['deltaT'] > time_cut_min) & (dfc['deltaT'] < time_cut_max)]
    bins_time_reduced = np.linspace(time_cut_min, time_cut_max, nbins_time)

    ### automated LYSO energy cut
    energy_cut = dfc_time[(dfc_time['E0_AU'] > 1300) & (dfc_time['E0_AU'] < 2300)]
    curve, useless, useless = plt.hist(energy_cut['E4_AU'], bins=bins_energy, density=False, log=log, histtype='step', lw=3,
                                       color='k')
    plt.close()
    peaks, useless = find_peaks(curve, prominence=max(curve) / 2)
    k = 0
    while curve[k] < 0.02 * max(curve):
        k += 1
    energy_cut_min = bins_energy[k]
    energy_cut_max = bins_energy[peaks][0] + 0.8 * (bins_energy[peaks][0] - bins_energy[k])
    plt.close()
    dfc_energy = dfc_time[(dfc_time['E4_AU'] > energy_cut_min) & (dfc_time['E4_AU'] < energy_cut_max)]

    ### gaussian peak fit
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
    x_lim0 = 1000
    x_lim1 = 2500
    n_events = len(dfc_energy[(dfc_energy['E0_AU'] < x_lim1) & (dfc_energy['E0_AU'] > x_lim0)])

    def erf_gauss(x, x_limit0, x_limit1, n_erf, n_gauss, loc, scale, sigma, mu):
        return (n_erf + n_gauss,
                n_erf * truncErfc.pdf(x, x_limit0, x_limit1, loc, scale) +
                n_gauss * truncnorm.pdf(x, x_limit0, x_limit1, loc=mu, scale=sigma))

    ### incertitude evenements totaux ?
    c = ExtendedUnbinnedNLL(dfc_energy['E0_AU'], erf_gauss)
    m = Minuit(c, x_limit0=x_lim0, x_limit1=x_lim1, n_erf=0.6 * n_events, n_gauss=0.4*n_events, loc=1800, scale=500, sigma=200,
               mu=1800)
    m.limits["loc", "sigma", "mu"] = (0, None)
    m.limits['n_erf', 'n_gauss'] = (0, n_events)
    m.limits["scale"] = (1, None)
    m.fixed['x_limit0'] = True
    m.fixed['x_limit1'] = True
    m.simplex()
    m.migrad()
    m.hesse()
    #print(m.params)
    nBins_fit = 100
    bins, step = np.linspace(m.values['x_limit0'], m.values['x_limit1'], nBins_fit + 1, retstep=True)
    xx = np.linspace(m.values['x_limit0'], m.values['x_limit1'], 1000)
    yy = step * erf_gauss(xx, m.values['x_limit0'], m.values['x_limit1'], m.values['n_erf'], m.values['n_gauss'], m.values['loc'],
                          m.values['scale'], m.values['sigma'], m.values['mu'])[1]
    plt.close()

    bins_count, u = np.histogram(dfc_energy['E0_AU'], bins=bins)
    xxx = 0.5 * (bins[1:] + bins[:-1])
    yyy = step * erf_gauss(xxx, m.values['x_limit0'], m.values['x_limit1'], m.values['n_erf'], m.values['n_gauss'], m.values['loc'],
                           m.values['scale'], m.values['sigma'],
                           m.values['mu'])[1]

    delta2 = np.sum([(bbins - yyyy) ** 2 / yyyy for bbins, yyyy in zip(bins_count, yyy) if yyyy != 0])
    df = nBins_fit - m.npar
    p_value = chi2.sf(delta2, df)

    plt.hist(dfc['deltaT'], bins=bins_time, log=True, density=False, histtype='step', lw=3, color='k')
    plt.vlines(time_cut_min, 0, 10000, colors='r', linestyles='dashed')
    plt.vlines(time_cut_max, 0, 10000, colors='r', linestyles='dashed')
    plt.xlabel('deltaT (ns)')
    plt.ylabel('count per bins')
    plt.title('Coincidence deltaT between LSC/LYSO')
    plt.savefig(directory + '\\coincidence_deltaT_full.png')
    plt.savefig(directory + '\\coincidence_deltaT_full.pdf')
    plt.close()

    plt.hist(dfc['deltaT'], bins=bins_time_reduced, log=False, density=False, histtype='step', lw=3, color='k')
    plt.xlabel('deltaT (ns)')
    plt.ylabel('count per bins')
    plt.title('Coincidence deltaT between LSC/LYSO')
    plt.savefig(directory + '\\coincidence_deltaT.png')
    plt.savefig(directory + '\\coincidence_deltaT.pdf')
    plt.close()

    plt.hist(dfc_energy['E0_AU'], bins=bins, density=False, log=log, lw=3, histtype='step', color='k')
    plt.plot(xx, yy, 'r')
    plt.xlabel('charge integral of LSC (AU)')
    plt.ylabel('count per bins')
    plt.title('Fit for the Gaussian part of the spectrum')
    plt.legend(['data spectrum', 'fit'])
    textstr = '\n'.join((r'n_gauss' + ' = ' + str(ufloat(m.values['n_gauss'], m.errors['n_gauss'])),
                         r'$\sigma$' + ' = ' + str(ufloat(m.values['sigma'], m.errors['sigma'])),
                         r'FWHM' + ' = ' + str(
                             ufloat(2 * np.sqrt(2 * np.log(2)) * m.values['sigma'],
                                    2 * np.sqrt(2 * np.log(2)) * m.errors['sigma'])),
                         r'$\mu$' + ' = ' + str(ufloat(m.values['mu'], m.errors['mu'])),
                         r'$\chi^2$' + ' = ' + str("%.1f" % delta2),
                         r'freedom' + ' = ' + str(nBins_fit - m.npar),
                         r'p_value' + ' = ' + str("%.1e" % p_value)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.gca().text(0.10, 0.05, textstr, transform=plt.gca().transAxes, fontsize=9, verticalalignment='bottom', bbox=props)
    plt.savefig(directory + '\\gaussian_fit.png')
    plt.savefig(directory + '\\gaussian_fit.pdf')
    plt.close()

    for data, legend, title in zip([dfc, dfc_time, dfc_energy], ['2D_original_spectrum', '2D_coincidence_deltaT_cut', '2D_LYSO_energy_cut'], ['raw', 'deltaT cut', 'LYSO energy']) :
        plt.hist2d(data['E0_AU'], data['E4_AU'], bins = bins_energy, norm = 'log')
        plt.xlabel('charge integral of LSC (AU)')
        plt.ylabel('charge integral of LYSO (AU)')
        plt.title(title + ' 2D histogram between LSC/LYSO')
        plt.savefig(directory + '\\' + legend + '.png')
        plt.savefig(directory + '\\' + legend + '.pdf')
        plt.close()

    for k, legend in zip([0, 4], ['LSC', 'LYSO']):
        plt.hist(dfc['E' + str(k) + '_AU'], bins=bins_energy, density=False, log=log, lw=3, histtype='step', color='k')
        plt.xlabel('charge integral (AU)')
        plt.ylabel('count per bins')
        plt.title('raw data spectrum of ' + legend)
        plt.legend(['raw'])
        plt.savefig(directory + '\\' + legend + '_original_spectrum.png')
        plt.savefig(directory + '\\' + legend + '_original_spectrum.pdf')

        plt.hist(dfc_time['E' + str(k) + '_AU'], bins=bins_energy, density=False, log=log, lw=3, histtype='step', color='b')
        plt.title('deltaT cut data spectrum of ' + legend)
        plt.legend(['raw', 'deltaT cut'])
        plt.savefig(directory + '\\' + legend + '_coincidence_deltaT_cut.png')
        plt.savefig(directory + '\\' + legend + '_coincidence_deltaT_cut.pdf')

        plt.hist(dfc_energy['E' + str(k) + '_AU'], bins=bins_energy, density=False, log=log, lw=3, histtype='step', color='r')
        plt.title('LYSO energy cut data spectrum of ' + legend)
        plt.legend(['raw', 'deltaT_cut', 'energy LYSO cut'])
        plt.savefig(directory + '\\' + legend + '_LYSO_energy_cut.png')
        plt.savefig(directory + '\\' + legend + '_LYSO_energy_cut.pdf')
        plt.close()

        plt.hist2d(dfc['E' + str(k) + '_AU'], dfc['deltaT'], bins=(bins_energy, bins_time_reduced), norm='log')
        plt.xlabel('charge integral (AU)')
        plt.ylabel('deltaT (ns)')
        plt.title('deltaT VS charge integral of ' + legend)
        plt.savefig(directory + '\\' + legend + '_coincidence_deltaT.png')
        plt.savefig(directory + '\\' + legend + '_coincidence_deltaT.pdf')
        plt.close()
    return ufloat(m.values['n_gauss'], m.errors['n_gauss']), ufloat(m.values['sigma'], m.errors['sigma']), ufloat(m.values['mu'], m.errors['mu']), p_value
