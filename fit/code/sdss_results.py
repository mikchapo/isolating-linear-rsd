"""Plot SDSS results."""

from scipy.integrate import odeint
from scipy.integrate import quad

import matplotlib.pyplot as plt
import numpy as np
import os


def growth_factor(y, a, omegam0):
    """Calculate growth factor."""
    delta, delta_prime = y
    H2_H02 = omegam0 / (a**3.) + (1. - omegam0)
    dyda = [delta_prime, ((omegam0 / 2 / (a**3.) / H2_H02 - 1.) * delta_prime +
                          omegam0 * delta / 2. / (a**4.) / H2_H02) * 3 / a]
    return dyda


def fsigma8_numerical(zs, sigma8=0.811, omegam0=0.315):
    """Calculate fsig8 numerically."""
    scale_factors = np.flip(1. / (1. + zs))
    y0 = [scale_factors[0], 1.]
    deltas = odeint(growth_factor, y0, scale_factors, args=(omegam0,))
    fsig8s = np.flip(sigma8 * scale_factors * deltas[:, 1] / deltas[-1, 0])
    return fsig8s


def omegam(z, omegam0=0.315):
    """Calculate Omega_m."""
    return omegam0 * (z + 1.)**3. / (omegam0 * (z + 1.)**3. + 1. - omegam0)


def fsigma8_approximate(z, sigma8=0.811, omegam0=0.315):
    """Calculate fsig8 approximately."""
    ell = 0.83 + 0.19 / (omegam0**0.8)
    beta = 0.98 + 0.01 / (omegam0**1.2)
    gamma = 0.64 + 0.09 / (omegam0**0.4)
    print("lambda:", ell, "beta:", beta, "gamma:", gamma)
    return ell * sigma8 * omegam(z, omegam0=omegam0)**(gamma) / (1. + z)**beta


def D_integrand(a, omegam, omegal):
    """Compute growth factor integrand."""
    return ((omegam / (a**3.) / (omegam / (a**3.) + omegal))**0.55 - 1.) / a


def calc_D(z, omegam, omegal):
    """Calculate growth factor."""
    a = 1. / (1. + z)
    integral = quad(D_integrand, 0., a, args=(omegam, omegal))
    return a*np.exp(integral[0])


# cosmology functions
def cosgrowRhoDE(z=1, w0=-1, wprime=0, rhoDE=1):
    """Alex function."""
    return (rhoDE * ((1 / (1 + z))**(-(3 + 3 * w0 + 6 * wprime))) *
            np.exp(-6 * wprime * (1 - 1 / (1 + z))))


def Einv(z, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime):
    """Alex function."""
    return 1./np.sqrt(OmegaR * (1 + z)**4 + OmegaM * (1 + z)**3 +
                      OmegaK * (1 + z)**2 +
                      OmegaL*cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1))


def Einva3(a, OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime):
    """Alex function."""
    return 1/((a**3) * (np.sqrt(OmegaR * a**(-4) + OmegaM * a**(-3) +
                                OmegaK * a**(-2) + OmegaL *
                                cosgrowRhoDE(z=1/a-1, w0=w0, wprime=wprime,
                                             rhoDE=1)))**3)


def growthRate(z, H0, OmegaM, OmegaL, OmegaR, w0, wprime, Sigma8):
    """Alex function."""
    OmegaK = 1 - OmegaM - OmegaL - OmegaR
    OmegaSum = (OmegaR * (1 + z)**4 + OmegaM * (1 + z)**3 + OmegaK * (1 + z)**2
                + OmegaL * cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1))
    Hz = H0 * np.sqrt(OmegaSum)
    OmegaRatz = (OmegaR * (1 + z)**4) / OmegaSum
    OmegaMatz = (OmegaM * (1 + z)**3) / OmegaSum
    OmegaLatz = (OmegaL * cosgrowRhoDE(z=z, w0=w0, wprime=wprime, rhoDE=1) /
                 OmegaSum)
    OmegaKatz = (OmegaK * (1 + z)**2) / OmegaSum
    Factor = ((5 * OmegaM / 2) * (Hz / H0) * (1 + z) *
              quad(Einva3, 0, 1/(1+z),
                   args=(OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime))[0])
    # return Factor
    Factor0 = ((5 * OmegaM / 2) *
               quad(Einva3, 0, 1,
                    args=(OmegaM, OmegaL, OmegaR, OmegaK, w0, wprime))[0])
    Sigma8atz = Sigma8 * (Factor/Factor0)/(1 + z)
    Rate = Sigma8atz * (-1 - OmegaMatz / 2 + OmegaLatz + (5 * OmegaMatz) /
                        (2 * Factor))
    return Rate


def plot_planck():
    """Alex code."""
    # load planck chain
    chainNum = 4
    for n in range(1, chainNum+1):
        # this chain sent in email as well
        c = np.loadtxt('../data/base_plikHM_TTTEEE_lowl_lowE_%d.txt' % (n))
        if n == 1:
            chain = c
        else:
            chain = np.append(chain, c, axis=0)

    # change your redshift range here
    z = np.linspace(0, 1.6, 20)

    # calculate fsigma8 at each change point, if you use a differnt chain you
    # have to chains the chain[c][i] to the correct indexes (see fsigma8 func
    # above for what should be what)
    x = []
    y = []
    for i in z:
        yl = []
        x.append(i)
        for c in range(len(chain)):
            # every 500th chain point, otherwise this takes forever but this
            # can vary and won't change your results outside of noise as it's
            # take every 500th point not subsampled at random
            if c % 500 == 0:
                yl.append(growthRate(i, chain[c][29], chain[c][31],
                                     chain[c][30], 0, -1, 0, chain[c][35]))
        y.append(yl)

    # gets the one and 2 sigma regions at each z picked
    yl = []
    yh = []
    yl2 = []
    yh2 = []
    ym = []
    for i in y:
        mean = np.average(i)
        std = np.std(i)
        yh.append(mean+std)
        yl.append(mean-std)
        yh2.append(mean+2*std)
        yl2.append(mean-2*std)
        ym.append(mean)

    # plotting
    plt.plot(x, ym, c='k', label=r'Planck 2018 + $\Lambda$CDM')
    plt.fill_between(x, yl, yh, alpha=0.5, color='grey')
    plt.fill_between(x, yl2, yh2, alpha=0.3, color='grey')


def paper_plot(thesis=False, paper_2022=False, defence=False):
    """Plot for paper."""
    zs = np.linspace(0., 100., 100001)
    fsig8s = fsigma8_numerical(zs)
    fsig8s_approx = fsigma8_approximate(zs)

    # print("z:", zs[737], "fsig8:", fsig8s[737])
    # print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s[700])

    z_uchuu = 0.7
    sigma8_0_uchuu = 0.8159
    omegam0_uchuu = 0.3089

    fsig8s_uchuu = fsigma8_numerical(zs, sigma8=sigma8_0_uchuu,
                                     omegam0=omegam0_uchuu)

    print("z Uchuu:", zs[570], "fsig8 Uchuu:", fsig8s_uchuu[570])
    print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s_uchuu[700])
    sigma8_070_uchuu = (sigma8_0_uchuu *
                        calc_D(z_uchuu, omegam0_uchuu, 1. - omegam0_uchuu) /
                        calc_D(0., omegam0_uchuu, 1. - omegam0_uchuu))
    print("sigma8(z=0.70) Uchuu:", sigma8_070_uchuu, "f(z=0.70) Uchuu:",
          fsig8s_uchuu[700] / sigma8_070_uchuu)

    # params[2] = (params[2] * calc_D(z_eboss, params[0], 1. - params[0]) /
    #              calc_D(z_aemulus, params[0], 1. - params[0]))
    # 0.75160436
    # fsig8s_uchuu = fsigma8_numerical(zs, sigma8=0.75160436, omegam0=0.3367394)
    # fsig8s_uchuu_aemulus = fsigma8_numerical(zs, sigma8=0.80191536,
    #                                          omegam0=0.3367394)

    # gamma_f_uchuu = 1.1133033

    # print("gamma_f:", gamma_f_uchuu)
    # print("fsig8 uchuu:", fsig8s_uchuu[570])
    # print("fsig8 uchuu aemulus:", fsig8s_uchuu_aemulus[700])
    # print("gamma_f*fsig8 uchuu:", gamma_f_uchuu*fsig8s_uchuu[570])

    mgs = (0.15, 0.53, 0.19)
    lange_lz = (0.25, 0.471, 0.024)
    lowz = (0.38, 0.497, 0.039)
    lange_hz = (0.40, 0.431, 0.025)
    cmass = (0.51, 0.458, 0.035)
    reid = (0.57, 0.450, 0.011)
    lrg = (0.698, 0.473, 0.044)
    # Redshift is offset to prevent overlap, best fit is an average of two redshift
    # models
    lrg_only = (0.737, 0.433, 0.066)
    elg = (0.85, 0.315, 0.095)
    qso = (1.48, 0.462, 0.045)
    chapman22 = (0.727, 0.365, 0.025)  # Full fit
    chapman22_ls = (0.727, 0.408, 0.038)
    lrg_ss = (0.757, 0.364, 0.042)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    zhai1 = (0.26, 0.404, 0.03)
    zhai2 = (0.41, 0.444, 0.025)
    zhai3 = (0.56, 0.385, 0.019)
    yuan22 = (0.52, 0.444, 0.016)

    print("Tension: %f-sigma" % ((fsig8s[737] - lrg_ss[1]) / lrg_ss[2]))
    print("Tension: %f-sigma Approx." % ((fsig8s_approx[737] - lrg_ss[1]) /
                                         lrg_ss[2]))

    print(fsig8s[700], fsig8s_approx[700])
    print(fsig8s[737], fsig8s_approx[737])

    marker_scaling = 1.5

    if thesis:
        fig_width = 6.375
        aspect_ratio = 3. / 4.

    elif defence and paper_2022:
        fig_width = 8.
        aspect_ratio = 4.5 / 9.5

    elif defence:
        fig_width = 8.
        aspect_ratio = 5. / 9.5

    else:
        fig_width = 10.
        aspect_ratio = 4.5 / 9.5

    plt.figure(figsize=(fig_width, fig_width * aspect_ratio), dpi=300)
    # fig_width = 6.
    # plt.figure(figsize=(fig_width, fig_width / 6. * 4.5), dpi=300)
    # plt.plot(zs, fsig8s, "k-")
    plot_planck()
    # plt.plot(zs, fsig8s_approx, "r-")
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="bs", label="SDSS MGS",
                 markersize=4*marker_scaling)
    plt.errorbar(lowz[0], lowz[1], yerr=lowz[2], fmt="bo", label="BOSS Galaxy",
                 markersize=4*marker_scaling)
    plt.errorbar(cmass[0], cmass[1], yerr=cmass[2], fmt="bo",
                 markersize=4*marker_scaling)
    plt.errorbar(lrg[0], lrg[1], yerr=lrg[2], fmt="b^", label="CMASS+eBOSS LRG",
                 markersize=5*marker_scaling)
    plt.errorbar(lrg_only[0], lrg_only[1], yerr=lrg_only[2], fmt="bv",
                 label="Large-scale eBOSS LRG", markersize=5*marker_scaling)
    plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="bd", label="eBOSS ELG",
                 markersize=5*marker_scaling)
    plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="bx", label="eBOSS QSO",
                 mew=2*marker_scaling, markersize=5*marker_scaling)
    plt.errorbar(lange_lz[0], lange_lz[1], yerr=lange_lz[2], fmt="go",
                 label="Lange et al. 2022", markersize=4*marker_scaling,
                 mfc='white')
    plt.errorbar(lange_hz[0], lange_hz[1], yerr=lange_hz[2], fmt="go",
                 markersize=4*marker_scaling,
                 mfc='white')
    if not paper_2022:
        plt.errorbar(zhai1[0], zhai1[1], yerr=zhai1[2], fmt="mo",
                     label="Zhai et al. 2023", markersize=4*marker_scaling,
                     mfc='white')
        plt.errorbar(zhai2[0], zhai2[1], yerr=zhai2[2], fmt="mo",
                     markersize=4*marker_scaling,
                     mfc='white')
        plt.errorbar(zhai3[0], zhai3[1], yerr=zhai3[2], fmt="mo",
                     markersize=4*marker_scaling,
                     mfc='white')
        # plt.errorbar(wigglez1[0], wigglez1[1], yerr=wigglez1[2], fmt="m*",
        #              label="WiggleZ", markersize=7*marker_scaling)
        # plt.errorbar(wigglez2[0], wigglez2[1], yerr=wigglez2[2], fmt="m*",
        #              markersize=7*marker_scaling)
        # plt.errorbar(wigglez3[0], wigglez3[1], yerr=wigglez3[2], fmt="m*",
        #              markersize=7*marker_scaling)
        # plt.errorbar(wigglez4[0], wigglez4[1], yerr=wigglez4[2], fmt="m*",
        #              markersize=7*marker_scaling)
        plt.errorbar(yuan22[0], yuan22[1], yerr=yuan22[2], fmt="co",
                     label="Yuan et al. 2022", markersize=4*marker_scaling,
                     mfc='white')
    plt.errorbar(reid[0], reid[1], yerr=reid[2], fmt="yo",
                 label="Reid et al. 2014", markersize=4*marker_scaling,
                 mfc='white')

    if thesis:
        if paper_2022:
            plt.errorbar(lrg_ss[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="This Chapter", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)

        else:
            plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="Chapter 3",
                         markersize=5*marker_scaling, mfc='white')
            plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
                         label="This Chapter", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)

    elif defence:
        if paper_2022:
            plt.errorbar(lrg_ss[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="Chapman et al. 2022", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)

        else:
            plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="Chapman et al. 2022",
                         markersize=5*marker_scaling, mfc='white')
            plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
                         label="Chapman et al. 2023", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)

    else:
        if paper_2022:
            plt.errorbar(lrg_ss[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="This Work", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)

        else:
            plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2],
                         fmt="rv", label="Chapman et al. 2022",
                         markersize=5*marker_scaling, mfc='white')
            plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
                         label="This Work", markersize=6*marker_scaling,
                         elinewidth=2, capsize=4, capthick=2)
    # plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
    #              label="Chapman et al. in prep", markersize=5*marker_scaling)
    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 1.6)

    if thesis:
        if not paper_2022:
            plt.ylim(0.2, 0.80)

    elif defence and not paper_2022:
        plt.ylim(0.25, 0.75)

    else:
        plt.ylim(0.25, 0.7)

    plt.legend(ncol=2)

    if defence:
        if paper_2022:
            plt.savefig("/Users/mike/Documents/administrative/thesis/defence/"
                        "images/sdss_results_2022.png")

        else:
            plt.savefig("/Users/mike/Documents/administrative/thesis/defence/"
                        "images/sdss_results_2023.png")


    elif paper_2022:
        plt.savefig("/Users/mike/Documents/research/RSD/output/thesis_plots/"
                    "f-discussion/sdss_results_thesis.png")

    elif thesis:
        plt.savefig("../output/plots/sdss_results_thesis.png")

    else:
        plt.savefig("../output/plots/sdss_results.jpg")
    # plt.savefig("../../../../administrative/postdocs/sdss_results.jpg")
    # plt.savefig("../../conferences/cosmology_from_home/sdss_results.jpg")


def rsd_plot(defence=False):
    """Plot for RSD section of thesis."""
    dFGS = (0.067, 0.423, 0.055)
    mgs = (0.15, 0.53, 0.19)
    gama1 = (0.18, 0.36, 0.09)
    gama2 = (0.38, 0.44, 0.06)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    boss1 = (0.385, 0.497, (0.039**2. + 0.024**2.)**0.5)
    boss2 = (0.51, 0.458, (0.035**2. + 0.015**2.)**0.5)
    boss3 = (0.61, 0.436, (0.034**2. + 0.009**2.)**0.5)
    vipers1 = (0.6, 0.55, 0.12)
    vipers2 = (0.86, 0.40, 0.11)
    # eboss_lrg = (0.698, 0.473, 0.044)
    # eboss_elg = (0.85, 0.315, 0.095)
    # eboss_qso = (1.48, 0.462, 0.045)

    marker_scaling = 1.5
    if defence:
        fig_width = 4.5
    else:
        fig_width = 6.375
    aspect_ratio = 3. / 4.

    plt.figure(figsize=(fig_width, fig_width * aspect_ratio), dpi=300)
    plot_planck()
    plt.errorbar(dFGS[0], dFGS[1], yerr=dFGS[2], fmt="b^", label="6DFGS",
                 markersize=5*marker_scaling)
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="gs", label="SDSS MGS",
                 markersize=4*marker_scaling)
    plt.errorbar(gama1[0], gama1[1], yerr=gama1[2], fmt="rv", label="GAMA",
                 markersize=5*marker_scaling)
    plt.errorbar(gama2[0], gama2[1], yerr=gama2[2], fmt="rv",
                 markersize=5*marker_scaling)
    plt.errorbar(wigglez1[0], wigglez1[1], yerr=wigglez1[2], fmt="c*",
                 label="WiggleZ", markersize=7*marker_scaling)
    plt.errorbar(wigglez2[0], wigglez2[1], yerr=wigglez2[2], fmt="c*",
                 markersize=7*marker_scaling)
    plt.errorbar(wigglez3[0], wigglez3[1], yerr=wigglez3[2], fmt="c*",
                 markersize=7*marker_scaling)
    plt.errorbar(wigglez4[0], wigglez4[1], yerr=wigglez4[2], fmt="c*",
                 markersize=7*marker_scaling)
    plt.errorbar(boss1[0], boss1[1], yerr=boss1[2], fmt="mo",
                 label="BOSS", markersize=4*marker_scaling)
    plt.errorbar(boss2[0], boss2[1], yerr=boss2[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(boss3[0], boss3[1], yerr=boss3[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(vipers1[0], vipers1[1], yerr=vipers1[2], label="VIPERS",
                 fmt="yd", markersize=5*marker_scaling)
    plt.errorbar(vipers2[0], vipers2[1], yerr=vipers2[2],
                 fmt="yd", markersize=5*marker_scaling)
    # plt.errorbar(lrg[0], lrg[1], yerr=lrg[2], fmt="b^",
    #              label="CMASS+eBOSS LRG", markersize=5*marker_scaling)
    # plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="bv", label="eBOSS ELG",
    #              markersize=5*marker_scaling)
    # plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="b>", label="eBOSS QSO",
    #              mew=2*marker_scaling, markersize=5*marker_scaling)

    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 1.)
    plt.ylim(0.2, 0.85)
    plt.legend(ncol=2)

    if defence:
        plt.savefig("/Users/mike/Documents/administrative/thesis/defence/"
                    "images/rsd_results.png")

    else:
        plt.savefig("/Users/mike/Documents/administrative/thesis/plots/"
                    "rsd_results.png")


def pd_app_plot():
    zs = np.linspace(0., 100., 100001)
    fsig8s = fsigma8_numerical(zs)
    fsig8s_approx = fsigma8_approximate(zs)

    # print("z:", zs[737], "fsig8:", fsig8s[737])
    # print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s[700])

    z_uchuu = 0.7
    sigma8_0_uchuu = 0.8159
    omegam0_uchuu = 0.3089

    fsig8s_uchuu = fsigma8_numerical(zs, sigma8=sigma8_0_uchuu,
                                     omegam0=omegam0_uchuu)

    print("z Uchuu:", zs[570], "fsig8 Uchuu:", fsig8s_uchuu[570])
    print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s_uchuu[700])
    sigma8_070_uchuu = (sigma8_0_uchuu *
                        calc_D(z_uchuu, omegam0_uchuu, 1. - omegam0_uchuu) /
                        calc_D(0., omegam0_uchuu, 1. - omegam0_uchuu))
    print("sigma8(z=0.70) Uchuu:", sigma8_070_uchuu, "f(z=0.70) Uchuu:",
          fsig8s_uchuu[700] / sigma8_070_uchuu)

    # params[2] = (params[2] * calc_D(z_eboss, params[0], 1. - params[0]) /
    #              calc_D(z_aemulus, params[0], 1. - params[0]))
    # 0.75160436
    # fsig8s_uchuu = fsigma8_numerical(zs, sigma8=0.75160436, omegam0=0.3367394)
    # fsig8s_uchuu_aemulus = fsigma8_numerical(zs, sigma8=0.80191536,
    #                                          omegam0=0.3367394)

    # gamma_f_uchuu = 1.1133033

    # print("gamma_f:", gamma_f_uchuu)
    # print("fsig8 uchuu:", fsig8s_uchuu[570])
    # print("fsig8 uchuu aemulus:", fsig8s_uchuu_aemulus[700])
    # print("gamma_f*fsig8 uchuu:", gamma_f_uchuu*fsig8s_uchuu[570])

    mgs = (0.15, 0.53, 0.19)
    lange_lz = (0.25, 0.471, 0.024)
    lowz = (0.38, 0.497, 0.039)
    lange_hz = (0.40, 0.431, 0.025)
    cmass = (0.51, 0.458, 0.035)
    reid = (0.57, 0.450, 0.011)
    lrg = (0.698, 0.473, 0.044)
    # Redshift is offset to prevent overlap, best fit is an average of two redshift
    # models
    lrg_only = (0.703, 0.433, 0.066)
    elg = (0.85, 0.315, 0.095)
    qso = (1.48, 0.462, 0.045)
    chapman22 = (0.737, 0.365, 0.025)  # Full fit
    chapman22_ls = (0.720, 0.408, 0.038)
    lrg_ss = (0.737, 0.364, 0.042)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    zhai1 = (0.26, 0.404, 0.03)
    zhai2 = (0.41, 0.444, 0.025)
    zhai3 = (0.56, 0.385, 0.019)
    yuan22 = (0.52, 0.444, 0.016)

    print("Tension: %f-sigma" % ((fsig8s[737] - lrg_ss[1]) / lrg_ss[2]))
    print("Tension: %f-sigma Approx." % ((fsig8s_approx[737] - lrg_ss[1]) /
                                         lrg_ss[2]))

    print(fsig8s[700], fsig8s_approx[700])
    print(fsig8s[737], fsig8s_approx[737])

    marker_scaling = 1.5

    text_frac = 0.8
    fig_width = 7. * text_frac
    plt.figure(figsize=(fig_width, fig_width/4.*3.), dpi=300)
    # fig_width = 6.
    # plt.figure(figsize=(fig_width, fig_width / 6. * 4.5), dpi=300)
    # plt.plot(zs, fsig8s, "k-")
    plot_planck()
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="bs",
                 label="Howlett et al. 2015",
                 markersize=4*marker_scaling)
    plt.errorbar(lowz[0], lowz[1], yerr=lowz[2], fmt="bo",
                 label="Alam et al. 2017",
                 markersize=4*marker_scaling)
    plt.errorbar(cmass[0], cmass[1], yerr=cmass[2], fmt="bo",
                 markersize=4*marker_scaling)
    plt.errorbar(lrg_only[0], lrg_only[1], yerr=lrg_only[2], fmt="bv",
                 label="Bautista et al. 2021",
                 markersize=5*marker_scaling)
    plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="bd",
                 label="de Mattia et al. 2021",
                 markersize=5*marker_scaling)
    # plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="bx", label="eBOSS QSO",
    #              mew=2*marker_scaling, markersize=5*marker_scaling)
    plt.errorbar(lange_lz[0], lange_lz[1], yerr=lange_lz[2], fmt="go",
                 label="Lange et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(lange_hz[0], lange_hz[1], yerr=lange_hz[2], fmt="go",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai1[0], zhai1[1], yerr=zhai1[2], fmt="mo",
                 label="Zhai et al. 2023", markersize=4*marker_scaling)
    plt.errorbar(zhai2[0], zhai2[1], yerr=zhai2[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai3[0], zhai3[1], yerr=zhai3[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(yuan22[0], yuan22[1], yerr=yuan22[2], fmt="co",
                 label="Yuan et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(reid[0], reid[1], yerr=reid[2], fmt="yo",
                 label="Reid et al. 2014", markersize=4*marker_scaling)

    plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2], fmt="rv",
                 label="Chapman et al. 2022", markersize=5*marker_scaling,
                 fillstyle='none')
    plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
                 label="Chapman et al. in prep", markersize=5*marker_scaling)
    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 0.9)
    plt.ylim(0.25, 0.7)
    plt.legend(ncol=2)
    plt.savefig("sdss_results.jpg")


def presentation_plot_1():
    zs = np.linspace(0., 100., 100001)
    fsig8s = fsigma8_numerical(zs)
    fsig8s_approx = fsigma8_approximate(zs)

    # print("z:", zs[737], "fsig8:", fsig8s[737])
    # print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s[700])

    z_uchuu = 0.7
    sigma8_0_uchuu = 0.8159
    omegam0_uchuu = 0.3089

    fsig8s_uchuu = fsigma8_numerical(zs, sigma8=sigma8_0_uchuu,
                                     omegam0=omegam0_uchuu)

    print("z Uchuu:", zs[570], "fsig8 Uchuu:", fsig8s_uchuu[570])
    print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s_uchuu[700])
    sigma8_070_uchuu = (sigma8_0_uchuu *
                        calc_D(z_uchuu, omegam0_uchuu, 1. - omegam0_uchuu) /
                        calc_D(0., omegam0_uchuu, 1. - omegam0_uchuu))
    print("sigma8(z=0.70) Uchuu:", sigma8_070_uchuu, "f(z=0.70) Uchuu:",
          fsig8s_uchuu[700] / sigma8_070_uchuu)

    # params[2] = (params[2] * calc_D(z_eboss, params[0], 1. - params[0]) /
    #              calc_D(z_aemulus, params[0], 1. - params[0]))
    # 0.75160436
    # fsig8s_uchuu = fsigma8_numerical(zs, sigma8=0.75160436, omegam0=0.3367394)
    # fsig8s_uchuu_aemulus = fsigma8_numerical(zs, sigma8=0.80191536,
    #                                          omegam0=0.3367394)

    # gamma_f_uchuu = 1.1133033

    # print("gamma_f:", gamma_f_uchuu)
    # print("fsig8 uchuu:", fsig8s_uchuu[570])
    # print("fsig8 uchuu aemulus:", fsig8s_uchuu_aemulus[700])
    # print("gamma_f*fsig8 uchuu:", gamma_f_uchuu*fsig8s_uchuu[570])

    mgs = (0.15, 0.53, 0.19)
    lange_lz = (0.25, 0.471, 0.024)
    lowz = (0.38, 0.497, 0.039)
    lange_hz = (0.40, 0.431, 0.025)
    cmass = (0.51, 0.458, 0.035)
    reid = (0.57, 0.450, 0.011)
    lrg = (0.698, 0.473, 0.044)
    # Redshift is offset to prevent overlap, best fit is an average of two redshift
    # models
    lrg_only = (0.703, 0.433, 0.066)
    elg = (0.85, 0.315, 0.095)
    qso = (1.48, 0.462, 0.045)
    chapman22 = (0.737, 0.365, 0.025)  # Full fit
    chapman22_ls = (0.720, 0.408, 0.038)
    lrg_ss = (0.737, 0.364, 0.042)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    zhai1 = (0.26, 0.404, 0.03)
    zhai2 = (0.41, 0.444, 0.025)
    zhai3 = (0.56, 0.385, 0.019)
    yuan22 = (0.52, 0.444, 0.016)

    print("Tension: %f-sigma" % ((fsig8s[737] - lrg_ss[1]) / lrg_ss[2]))
    print("Tension: %f-sigma Approx." % ((fsig8s_approx[737] - lrg_ss[1]) /
                                         lrg_ss[2]))

    print(fsig8s[700], fsig8s_approx[700])
    print(fsig8s[737], fsig8s_approx[737])

    marker_scaling = 1.5

    text_frac = 0.8
    fig_width = 7. * text_frac
    plt.figure(figsize=(fig_width, fig_width/4.*3.), dpi=300)
    # fig_width = 6.
    # plt.figure(figsize=(fig_width, fig_width / 6. * 4.5), dpi=300)
    # plt.plot(zs, fsig8s, "k-")
    plot_planck()
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="bs", label="SDSS MGS",
                 markersize=4*marker_scaling)
    # plt.errorbar(wigglez1[0], wigglez1[1], yerr=wigglez1[2], fmt="m*",
    #              label="WiggleZ", markersize=7*marker_scaling)
    # plt.errorbar(wigglez2[0], wigglez2[1], yerr=wigglez2[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez3[0], wigglez3[1], yerr=wigglez3[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez4[0], wigglez4[1], yerr=wigglez4[2], fmt="m*",
    #              markersize=7*marker_scaling)
    plt.errorbar(lowz[0], lowz[1], yerr=lowz[2], fmt="go", label="BOSS Galaxy",
                 markersize=4*marker_scaling)
    plt.errorbar(cmass[0], cmass[1], yerr=cmass[2], fmt="go",
                 markersize=4*marker_scaling)
    plt.errorbar(lrg[0], lrg[1], yerr=lrg[2], fmt="r^", label="eBOSS LRG",
                 markersize=5*marker_scaling)
    plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="md", label="eBOSS ELG",
                 markersize=5*marker_scaling)
    plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="yx", label="eBOSS QSO",
                 mew=2*marker_scaling, markersize=5*marker_scaling)
    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 1.6)
    plt.ylim(0.25, 0.7)
    plt.legend(ncol=2)
    plt.savefig("{}/Documents/research/presentations/Berkeley_2022-10-25/"
                "plots/large_scale_fsig8_results"
                ".jpg".format(os.path.expanduser('~')))


def presentation_plot_2():
    zs = np.linspace(0., 100., 100001)
    fsig8s = fsigma8_numerical(zs)
    fsig8s_approx = fsigma8_approximate(zs)

    # print("z:", zs[737], "fsig8:", fsig8s[737])
    # print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s[700])

    z_uchuu = 0.7
    sigma8_0_uchuu = 0.8159
    omegam0_uchuu = 0.3089

    fsig8s_uchuu = fsigma8_numerical(zs, sigma8=sigma8_0_uchuu,
                                     omegam0=omegam0_uchuu)

    print("z Uchuu:", zs[570], "fsig8 Uchuu:", fsig8s_uchuu[570])
    print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s_uchuu[700])
    sigma8_070_uchuu = (sigma8_0_uchuu *
                        calc_D(z_uchuu, omegam0_uchuu, 1. - omegam0_uchuu) /
                        calc_D(0., omegam0_uchuu, 1. - omegam0_uchuu))
    print("sigma8(z=0.70) Uchuu:", sigma8_070_uchuu, "f(z=0.70) Uchuu:",
          fsig8s_uchuu[700] / sigma8_070_uchuu)

    # params[2] = (params[2] * calc_D(z_eboss, params[0], 1. - params[0]) /
    #              calc_D(z_aemulus, params[0], 1. - params[0]))
    # 0.75160436
    # fsig8s_uchuu = fsigma8_numerical(zs, sigma8=0.75160436, omegam0=0.3367394)
    # fsig8s_uchuu_aemulus = fsigma8_numerical(zs, sigma8=0.80191536,
    #                                          omegam0=0.3367394)

    # gamma_f_uchuu = 1.1133033

    # print("gamma_f:", gamma_f_uchuu)
    # print("fsig8 uchuu:", fsig8s_uchuu[570])
    # print("fsig8 uchuu aemulus:", fsig8s_uchuu_aemulus[700])
    # print("gamma_f*fsig8 uchuu:", gamma_f_uchuu*fsig8s_uchuu[570])

    mgs = (0.15, 0.53, 0.19)
    lange_lz = (0.25, 0.471, 0.024)
    lowz = (0.38, 0.497, 0.039)
    lange_hz = (0.40, 0.431, 0.025)
    cmass = (0.51, 0.458, 0.035)
    reid = (0.57, 0.450, 0.011)
    lrg = (0.698, 0.473, 0.044)
    # Redshift is offset to prevent overlap, best fit is an average of two redshift
    # models
    lrg_only = (0.720, 0.433, 0.066)
    elg = (0.85, 0.315, 0.095)
    qso = (1.48, 0.462, 0.045)
    chapman22 = (0.737, 0.365, 0.025)  # Full fit
    chapman22_ls = (0.737, 0.408, 0.038)
    lrg_ss = (0.737, 0.364, 0.042)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    zhai1 = (0.26, 0.404, 0.03)
    zhai2 = (0.41, 0.444, 0.025)
    zhai3 = (0.56, 0.385, 0.019)
    yuan22 = (0.52, 0.444, 0.016)

    print("Tension: %f-sigma" % ((fsig8s[737] - lrg_ss[1]) / lrg_ss[2]))
    print("Tension: %f-sigma Approx." % ((fsig8s_approx[737] - lrg_ss[1]) /
                                         lrg_ss[2]))

    print(fsig8s[700], fsig8s_approx[700])
    print(fsig8s[737], fsig8s_approx[737])

    marker_scaling = 1.5

    fig_width = 10.
    plt.figure(figsize=(fig_width, fig_width / 9.5 * 4.5), dpi=300)
    # fig_width = 6.
    # plt.figure(figsize=(fig_width, fig_width / 6. * 4.5), dpi=300)
    # plt.plot(zs, fsig8s, "k-")
    plot_planck()
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="bs", label="SDSS MGS",
                 markersize=4*marker_scaling)
    plt.errorbar(lowz[0], lowz[1], yerr=lowz[2], fmt="bo", label="BOSS Galaxy",
                 markersize=4*marker_scaling)
    plt.errorbar(cmass[0], cmass[1], yerr=cmass[2], fmt="bo",
                 markersize=4*marker_scaling)
    plt.errorbar(lrg[0], lrg[1], yerr=lrg[2], fmt="b^", label="CMASS+eBOSS LRG",
                 markersize=5*marker_scaling)
    plt.errorbar(lrg_only[0], lrg_only[1], yerr=lrg_only[2], fmt="bv",
                 label="Large-scale eBOSS LRG", markersize=5*marker_scaling)
    plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="bd", label="eBOSS ELG",
                 markersize=5*marker_scaling)
    plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="bx", label="eBOSS QSO",
                 mew=2*marker_scaling, markersize=5*marker_scaling)
    plt.errorbar(lange_lz[0], lange_lz[1], yerr=lange_lz[2], fmt="go",
                 label="Lange et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(lange_hz[0], lange_hz[1], yerr=lange_hz[2], fmt="go",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai1[0], zhai1[1], yerr=zhai1[2], fmt="mo",
                 label="Zhai et al. 2023", markersize=4*marker_scaling)
    plt.errorbar(zhai2[0], zhai2[1], yerr=zhai2[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai3[0], zhai3[1], yerr=zhai3[2], fmt="mo",
                 markersize=4*marker_scaling)
    # plt.errorbar(wigglez1[0], wigglez1[1], yerr=wigglez1[2], fmt="m*",
    #              label="WiggleZ", markersize=7*marker_scaling)
    # plt.errorbar(wigglez2[0], wigglez2[1], yerr=wigglez2[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez3[0], wigglez3[1], yerr=wigglez3[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez4[0], wigglez4[1], yerr=wigglez4[2], fmt="m*",
    #              markersize=7*marker_scaling)
    plt.errorbar(yuan22[0], yuan22[1], yerr=yuan22[2], fmt="co",
                 label="Yuan et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(reid[0], reid[1], yerr=reid[2], fmt="yo",
                 label="Reid et al. 2014", markersize=4*marker_scaling)

    plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2], fmt="rv",
                 label="Chapman et al. 2022", markersize=5*marker_scaling)
    # plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
    #              label="This Work", markersize=5*marker_scaling)
    # plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
    #              label="Chapman et al. in prep", markersize=5*marker_scaling)
    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 1.6)
    plt.ylim(0.25, 0.7)
    plt.legend(ncol=2)
    plt.savefig("{}/Documents/research/presentations/Berkeley_2022-10-25/"
                "plots/sdss_fsig8_results"
                ".jpg".format(os.path.expanduser('~')))


def presentation_plot_3():
    zs = np.linspace(0., 100., 100001)
    fsig8s = fsigma8_numerical(zs)
    fsig8s_approx = fsigma8_approximate(zs)

    # print("z:", zs[737], "fsig8:", fsig8s[737])
    # print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s[700])

    z_uchuu = 0.7
    sigma8_0_uchuu = 0.8159
    omegam0_uchuu = 0.3089

    fsig8s_uchuu = fsigma8_numerical(zs, sigma8=sigma8_0_uchuu,
                                     omegam0=omegam0_uchuu)

    print("z Uchuu:", zs[570], "fsig8 Uchuu:", fsig8s_uchuu[570])
    print("z Uchuu:", zs[700], "fsig8 Uchuu:", fsig8s_uchuu[700])
    sigma8_070_uchuu = (sigma8_0_uchuu *
                        calc_D(z_uchuu, omegam0_uchuu, 1. - omegam0_uchuu) /
                        calc_D(0., omegam0_uchuu, 1. - omegam0_uchuu))
    print("sigma8(z=0.70) Uchuu:", sigma8_070_uchuu, "f(z=0.70) Uchuu:",
          fsig8s_uchuu[700] / sigma8_070_uchuu)

    # params[2] = (params[2] * calc_D(z_eboss, params[0], 1. - params[0]) /
    #              calc_D(z_aemulus, params[0], 1. - params[0]))
    # 0.75160436
    # fsig8s_uchuu = fsigma8_numerical(zs, sigma8=0.75160436, omegam0=0.3367394)
    # fsig8s_uchuu_aemulus = fsigma8_numerical(zs, sigma8=0.80191536,
    #                                          omegam0=0.3367394)

    # gamma_f_uchuu = 1.1133033

    # print("gamma_f:", gamma_f_uchuu)
    # print("fsig8 uchuu:", fsig8s_uchuu[570])
    # print("fsig8 uchuu aemulus:", fsig8s_uchuu_aemulus[700])
    # print("gamma_f*fsig8 uchuu:", gamma_f_uchuu*fsig8s_uchuu[570])

    mgs = (0.15, 0.53, 0.19)
    lange_lz = (0.25, 0.471, 0.024)
    lowz = (0.38, 0.497, 0.039)
    lange_hz = (0.40, 0.431, 0.025)
    cmass = (0.51, 0.458, 0.035)
    reid = (0.57, 0.450, 0.011)
    lrg = (0.698, 0.473, 0.044)
    # Redshift is offset to prevent overlap, best fit is an average of two redshift
    # models
    lrg_only = (0.720, 0.433, 0.066)
    elg = (0.85, 0.315, 0.095)
    qso = (1.48, 0.462, 0.045)
    chapman22 = (0.737, 0.365, 0.025)  # Full fit
    chapman22_ls = (0.737, 0.408, 0.038)
    lrg_ss = (0.754, 0.364, 0.042)
    wigglez1 = (0.22, 0.42, 0.07)
    wigglez2 = (0.41, 0.45, 0.04)
    wigglez3 = (0.6, 0.43, 0.04)
    wigglez4 = (0.78, 0.38, 0.04)
    zhai1 = (0.26, 0.404, 0.03)
    zhai2 = (0.41, 0.444, 0.025)
    zhai3 = (0.56, 0.385, 0.019)
    yuan22 = (0.52, 0.444, 0.016)

    print("Tension: %f-sigma" % ((fsig8s[737] - lrg_ss[1]) / lrg_ss[2]))
    print("Tension: %f-sigma Approx." % ((fsig8s_approx[737] - lrg_ss[1]) /
                                         lrg_ss[2]))

    print(fsig8s[700], fsig8s_approx[700])
    print(fsig8s[737], fsig8s_approx[737])

    marker_scaling = 1.5

    fig_width = 10.
    plt.figure(figsize=(fig_width, fig_width / 9.5 * 4.5), dpi=300)
    # fig_width = 6.
    # plt.figure(figsize=(fig_width, fig_width / 6. * 4.5), dpi=300)
    # plt.plot(zs, fsig8s, "k-")
    plot_planck()
    plt.errorbar(mgs[0], mgs[1], yerr=mgs[2], fmt="bs", label="SDSS MGS",
                 markersize=4*marker_scaling)
    plt.errorbar(lowz[0], lowz[1], yerr=lowz[2], fmt="bo", label="BOSS Galaxy",
                 markersize=4*marker_scaling)
    plt.errorbar(cmass[0], cmass[1], yerr=cmass[2], fmt="bo",
                 markersize=4*marker_scaling)
    plt.errorbar(lrg[0], lrg[1], yerr=lrg[2], fmt="b^", label="CMASS+eBOSS LRG",
                 markersize=5*marker_scaling)
    plt.errorbar(lrg_only[0], lrg_only[1], yerr=lrg_only[2], fmt="bv",
                 label="Large-scale eBOSS LRG", markersize=5*marker_scaling)
    plt.errorbar(elg[0], elg[1], yerr=elg[2], fmt="bd", label="eBOSS ELG",
                 markersize=5*marker_scaling)
    plt.errorbar(qso[0], qso[1], yerr=qso[2], fmt="bx", label="eBOSS QSO",
                 mew=2*marker_scaling, markersize=5*marker_scaling)
    plt.errorbar(lange_lz[0], lange_lz[1], yerr=lange_lz[2], fmt="go",
                 label="Lange et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(lange_hz[0], lange_hz[1], yerr=lange_hz[2], fmt="go",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai1[0], zhai1[1], yerr=zhai1[2], fmt="mo",
                 label="Zhai et al. 2023", markersize=4*marker_scaling)
    plt.errorbar(zhai2[0], zhai2[1], yerr=zhai2[2], fmt="mo",
                 markersize=4*marker_scaling)
    plt.errorbar(zhai3[0], zhai3[1], yerr=zhai3[2], fmt="mo",
                 markersize=4*marker_scaling)
    # plt.errorbar(wigglez1[0], wigglez1[1], yerr=wigglez1[2], fmt="m*",
    #              label="WiggleZ", markersize=7*marker_scaling)
    # plt.errorbar(wigglez2[0], wigglez2[1], yerr=wigglez2[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez3[0], wigglez3[1], yerr=wigglez3[2], fmt="m*",
    #              markersize=7*marker_scaling)
    # plt.errorbar(wigglez4[0], wigglez4[1], yerr=wigglez4[2], fmt="m*",
    #              markersize=7*marker_scaling)
    plt.errorbar(yuan22[0], yuan22[1], yerr=yuan22[2], fmt="co",
                 label="Yuan et al. 2022", markersize=4*marker_scaling)
    plt.errorbar(reid[0], reid[1], yerr=reid[2], fmt="yo",
                 label="Reid et al. 2014", markersize=4*marker_scaling)

    plt.errorbar(chapman22_ls[0], chapman22_ls[1], yerr=chapman22_ls[2], fmt="rv",
                 label="Chapman et al. 2022", markersize=5*marker_scaling,
                 fillstyle='none')
    plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
                 label="Chapman et al. In Prep.", markersize=5*marker_scaling)
    # plt.errorbar(lrg_ss[0], lrg_ss[1], yerr=lrg_ss[2], fmt="rv",
    #              label="Chapman et al. in prep", markersize=5*marker_scaling)
    plt.xlabel("z")
    plt.ylabel(r'$f\sigma_8$')
    plt.xlim(0., 1.6)
    plt.ylim(0.25, 0.7)
    plt.legend(ncol=2)
    plt.savefig("{}/Documents/research/presentations/Berkeley_2022-10-25/"
                "plots/new_sdss_fsig8_results"
                ".jpg".format(os.path.expanduser('~')))


rsd_plot(defence=True)
paper_plot(thesis=False, paper_2022=False, defence=True)
paper_plot(thesis=False, paper_2022=True, defence=True)
