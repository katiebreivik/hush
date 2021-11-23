import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from legwork.visualisation import plot_sensitivity_curve
from legwork import psd


def plot_LISAcurves(model):
    def func(x, a, b, c, d, e):
        return a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4

    def cosmic_confusion(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        lisa_psd_no_conf = psd.power_spectral_density(
            f, confusion_noise=None, t_obs=4 * u.yr
        )
        conf = 10 ** func(
            x=np.log10(f.value), 
            a=popt[0], 
            b=popt[1], 
            c=popt[2], 
            d=popt[3], 
            e=popt[4]
        ) * t_obs.to(u.s)

        psd_plus_conf = conf + lisa_psd_no_conf
        return psd_plus_conf.to(u.Hz ** (-1))

    resolved = pd.read_hdf("../data/results.hdf", key="resolved_DWDs_{}_{}".format(model, "fiducial"))
    popt = pd.read_hdf("../data/results.hdf", key="conf_fit_DWDs_{}_{}".format(model, "fiducial"))
    popt = popt.values.flatten()
    
    resolved_HeHe = resolved.loc[(resolved.kstar_1 == 10) & (resolved.kstar_2 == 10)]
    resolved_COHe = resolved.loc[((resolved.kstar_1 == 11) & (resolved.kstar_2 == 10)) |
                                 ((resolved.kstar_2 == 11) & (resolved.kstar_1 == 10))]
    resolved_COCO = resolved.loc[(resolved.kstar_1 == 11) & (resolved.kstar_2 == 11)]
    resolved_ONeX = resolved.loc[((resolved.kstar_1 == 12) & (resolved.kstar_2.isin([10,11,12]))) |
                                 ((resolved.kstar_2 == 12) & (resolved.kstar_1.isin([10,11,12])))]
    
    t_obs = 4 * u.yr
    
    
    psd_conf = psd.power_spectral_density(
        f=np.linspace(1e-4, 1e-1, 1000000) * u.Hz, 
        instrument="custom", 
        custom_psd=cosmic_confusion, 
        t_obs=t_obs, 
        L=None, 
        approximate_R=True, 
        confusion_noise=None
    )
    
    
    Heasd = ((1/4 * t_obs)**(1/2) * resolved_HeHe.h_0.values).to(u.Hz**(-1/2))
    COasd = ((1/4 * t_obs)**(1/2) * resolved_COCO.h_0.values).to(u.Hz**(-1/2))
    COHeasd = ((1/4 * t_obs)**(1/2) * resolved_COHe.h_0.values).to(u.Hz**(-1/2))
    ONeasd = ((1/4 * t_obs)**(1/2) * resolved_ONeX.h_0.values).to(u.Hz**(-1/2))

    fig, ax = plt.subplots(1, 4, figsize=(25, 5))
    ax[0].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True)
    ax[0].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[0].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[0].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[0].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:tomato red', label='He + He', rasterized=True)
    lgnd = ax[0].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False, fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    ax[0].text(0.1, 3e-17, model+', SNR > 7: {}'.format(len(Heasd)), fontsize=24, horizontalalignment='right')
    
    ax[2].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True)
    ax[2].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[2].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[2].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[2].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:pink', label='CO + CO', rasterized=True)
    lgnd = ax[2].legend(loc='lower left',  ncol=4, borderaxespad=0, frameon=False, fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    ax[2].text(0.1, 3e-17, model+', SNR > 7: {}'.format(len(COasd)),  fontsize=24, horizontalalignment='right')
    
    ax[1].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True)
    ax[1].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[1].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[1].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[1].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:blurple', label='CO + He', rasterized=True)
    lgnd = ax[1].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False,
                        fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    ax[1].text(0.1, 3e-17, model+', SNR > 7: {}'.format(len(COHeasd)), fontsize=24,
           horizontalalignment='right')
    
    

    ax[3].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True)
    ax[3].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[3].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[3].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[3].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light blue', label='ONe + X')
    lgnd = ax[3].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False,
                        fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    ax[3].text(0.1, 3e-17, model+', SNR > 7: {}'.format(len(ONeasd)), fontsize=24,
               horizontalalignment='right')
    
    for i in range(4):
        ax[i].set_yscale('log')
        ax[i].set_xscale('log')
        ax[i].tick_params(labelsize=20)
        ax[i].set_xlabel(r'f$_{\rm{GW}}$ [Hz]', size=24)
    ax[0].set_ylabel(r'ASD [Hz$^{-1/2}$]', size=24)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)
    plt.savefig("LISA_SNR_{}.pdf".format(model), dpi=100)

    return


for model in ["F50", "FZ"]:
    plot_LISAcurves(model=model)
