import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from legwork.visualisation import plot_sensitivity_curve
from legwork import psd
import dutils
import paths


def plot_LISAcurves(var, model, label, ii):
    
    def func(x, a, b, c, d, e):
        return a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4
    
    if model == 'fiducial_Z':
        cosmic_confusion = dutils.cosmic_confusion_var_fiducial
    elif model == 'alpha25_Z':
        cosmic_confusion = dutils.cosmic_confusion_var_alpha25
    elif model == 'alpha5_Z':
        cosmic_confusion = dutils.cosmic_confusion_var_alpha5
    elif model == 'q3_Z':
        cosmic_confusion = dutils.cosmic_confusion_var_q3
    elif model == 'fiducial':
        cosmic_confusion = dutils.cosmic_confusion_50_fiducial
    elif model == 'alpha25':
        cosmic_confusion = dutils.cosmic_confusion_50_alpha25
    elif model == 'alpha5':
        cosmic_confusion = dutils.cosmic_confusion_50_alpha5
    elif model == 'q3':
        cosmic_confusion = dutils.cosmic_confusion_50_q3
            
    resolved = pd.read_hdf(paths.data / "results.hdf", key="resolved_DWDs_{}_{}".format(var, model))
    
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
    
    
    Heasd = ((t_obs)**(1/2) * resolved_HeHe.h_0.values).to(u.Hz**(-1/2))
    COasd = ((t_obs)**(1/2) * resolved_COCO.h_0.values).to(u.Hz**(-1/2))
    COHeasd = ((t_obs)**(1/2) * resolved_COHe.h_0.values).to(u.Hz**(-1/2))
    ONeasd = ((t_obs)**(1/2) * resolved_ONeX.h_0.values).to(u.Hz**(-1/2))

    
    ax[ii, 0].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True, lw=2)
    ax[ii, 0].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 0].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 0].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 0].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:tomato red', label='He + He', rasterized=True)
    lgnd = ax[ii, 0].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False, fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    if ("Z" not in model) & ("fiducial" not in model):
        ax[ii, 0].text(0.1, 2e-16, label+', '+var+', SNR > 7: {}'.format(len(Heasd)), fontsize=24, horizontalalignment='right')
    else:
        ax[ii, 0].text(0.1, 2e-16, var+', SNR > 7: {}'.format(len(Heasd)), fontsize=24, horizontalalignment='right')
    
    
    ax[ii, 2].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True, lw=2)
    ax[ii, 2].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 2].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 2].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 2].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:pink', label='CO + CO', rasterized=True)
    lgnd = ax[ii, 2].legend(loc='lower left',  ncol=4, borderaxespad=0, frameon=False, fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    if ("Z" not in model) & ("fiducial" not in model):
        ax[ii, 2].text(0.1, 2e-16, label+', '+var+', SNR > 7: {}'.format(len(COasd)),  fontsize=24, horizontalalignment='right')
    else:
        ax[ii, 2].text(0.1, 2e-16, var+', SNR > 7: {}'.format(len(COasd)),  fontsize=24, horizontalalignment='right')
    
    ax[ii, 1].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True, lw=2)
    ax[ii, 1].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 1].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 1].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 1].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:blurple', label='CO + He', rasterized=True)
    lgnd = ax[ii, 1].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False,
                        fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    if ("Z" not in model) & ("fiducial" not in model):
        ax[ii, 1].text(0.1, 2e-16, label+', '+var+', SNR > 7: {}'.format(len(COHeasd)), fontsize=24,
           horizontalalignment='right')
    else:
        ax[ii, 1].text(0.1, 2e-16, var+', SNR > 7: {}'.format(len(COHeasd)), fontsize=24,
           horizontalalignment='right')
   
    

    ax[ii, 3].plot(np.linspace(1e-4, 1e-1, 1000000), psd_conf**0.5, c='black', rasterized=True, lw=2)
    ax[ii, 3].scatter(resolved_HeHe.f_gw, Heasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 3].scatter(resolved_COCO.f_gw, COasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 3].scatter(resolved_COHe.f_gw, COHeasd, zorder=10, color='xkcd:light grey', rasterized=True)
    ax[ii, 3].scatter(resolved_ONeX.f_gw, ONeasd, zorder=10, color='xkcd:light blue', label='ONe + X')
    lgnd = ax[ii, 3].legend(loc='lower left', ncol=4, borderaxespad=0, frameon=False,
                        fontsize=22, markerscale=2.5, handletextpad=0.15)
    lgnd.legendHandles[0]._sizes = [100]
    if ("Z" not in model) & ("fiducial" not in model):
        ax[ii, 3].text(0.1, 2e-16, label+', '+var+', SNR > 7: {}'.format(len(ONeasd)), fontsize=24,
               horizontalalignment='right')
    else:
        ax[ii, 3].text(0.1, 2e-16, var+', SNR > 7: {}'.format(len(ONeasd)), fontsize=24,
               horizontalalignment='right')
    
    for i in range(4):
        ax[ii, i].set_yscale('log')
        ax[ii, i].set_xscale('log')
        ax[ii, i].tick_params(labelsize=20)
        if ii == 1:
            ax[ii, i].set_xlabel(r'f$_{\rm{GW}}$ [Hz]', size=24)
        ax[ii, i].set_ylim(top=5e-16)
    ax[ii, 0].set_ylabel(r'ASD [Hz$^{-1/2}$]', size=24)

    return resolved_HeHe, resolved_COHe, resolved_COCO, resolved_ONeX



models = ["fiducial", "q3", "alpha25", "alpha5"]
model_labels = ["fiducial", "$q3$", "$\\alpha 25$", "$\\alpha 5$"]
for i, model in enumerate(models):
    
    fig, ax = plt.subplots(2, 4, figsize=(25, 10))

    plot_LISAcurves(var="F50", model=model, label=model_labels[i], ii=0)
    plot_LISAcurves(var="FZ", model=model+"_Z", label=model_labels[i], ii=1)
    
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)
    plt.savefig(paths.figures / "LISA_SNR_{}.pdf".format(model), dpi=100)
