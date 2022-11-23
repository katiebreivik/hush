import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import paths

FIREmin = 0.00015
FIREmax = 13.346
Z_sun = 0.02
num = 30

plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rcParams["font.size"] = 14


met_bins = np.logspace(np.log10(FIREmin), np.log10(FIREmax), num) * Z_sun
    
models = ["fiducial",  "alpha25", "alpha5", "q3"]
model_labels = ["fiducial",  "$\\alpha 25$", "$\\alpha 5$", "$q3$"]

fig, ax = plt.subplots(4, 4, figsize=(16, 14))

for m, model in enumerate(models):
    model_Z = model + "_Z"
    numsFZ = pd.read_hdf(paths.data / "results.hdf", key="numLISA_30bins_{}_{}".format("FZ", model_Z))
    numsF50 = pd.read_hdf(paths.data / "results.hdf", key="numLISA_30bins_{}_{}".format("F50", model))

    Henums = numsFZ.He.values
    COHenums = numsFZ.COHe.values
    COnums = numsFZ.CO.values
    ONenums = numsFZ.ONe.values

    Henums05 = numsF50.He.values
    COHenums05 = numsF50.COHe.values
    COnums05 = numsF50.CO.values
    ONenums05 = numsF50.ONe.values

    ax[m, 0].plot(
        np.log10(met_bins[1:]/Z_sun), 
        Henums/1e5, 
        drawstyle='steps-mid', 
        color='xkcd:tomato red', 
        lw=3, 
        label='FZ',
        rasterized=True
    )

    ax[m, 0].plot(
        np.log10(met_bins[1:]/Z_sun), 
        Henums05/1e5,       
        drawstyle='steps-mid', 
        color='xkcd:tomato red', 
        ls='--', 
        lw=3, 
        label='F50', 
        rasterized=True
    )

    
    ax[m, 1].plot(
        np.log10(met_bins[1:]/Z_sun), 
        COHenums/1e5, 
        drawstyle='steps-mid',     
        color='xkcd:blurple', 
        lw=3, 
        label='FZ', 
        rasterized=True
    )
    ax[m, 1].plot(
        np.log10(met_bins[1:]/Z_sun), 
        COHenums05/1e5, 
        drawstyle='steps-mid', 
        color='xkcd:blurple', 
        ls='--', 
        lw=3, 
        label='F50', 
        rasterized=True
    )


    ax[m, 2].plot(
        np.log10(met_bins[1:]/Z_sun), 
        COnums/1e5, 
        drawstyle='steps-mid', 
        color='xkcd:pink', 
        lw=3, 
        label='FZ', 
        rasterized=True
    )

    ax[m, 2].plot(
        np.log10(met_bins[1:]/Z_sun), 
        COnums05/1e5, 
        drawstyle='steps-mid', 
        color='xkcd:pink', 
        ls='--', 
        lw=3, 
        label='F50', 
        rasterized=True
    )


    ax[m, 3].plot(
        np.log10(met_bins[1:]/Z_sun), 
        ONenums/1e5, 
        drawstyle='steps-mid', 
        color='xkcd:light blue', 
        lw=3, 
        label='FZ', 
        rasterized=True
    )

    ax[m, 3].plot(
        np.log10(met_bins[1:]/Z_sun), 
        ONenums05/1e5, 
        drawstyle='steps-mid',
        color='xkcd:light blue', 
        ls='--', 
        lw=3, 
        label='F50', 
        rasterized=True
    )

    
    for i in range(4):
        if (m > 0):
            ax[m, i].text(0.05, 0.875, model_labels[m], fontsize=18, transform=ax[m, i].transAxes)
        if m == 0:
            ax[m, 0].text(0.05, 0.875, 'He + He', fontsize=18, transform=ax[m, 0].transAxes)
            ax[m, 1].text(0.05, 0.875, 'CO + He', fontsize=18, transform=ax[m, 1].transAxes)
            ax[m, 2].text(0.05, 0.875, 'CO + CO', fontsize=18, transform=ax[m, 2].transAxes)
            ax[m, 3].text(0.05, 0.875, 'ONe + X', fontsize=18, transform=ax[m, 3].transAxes)
            ax[m, i].text(0.05, 0.75, model_labels[m], fontsize=18, transform=ax[m, i].transAxes)
        if m == 3:
            ax[m, i].set_xlabel('Log$_{10}$(Z/Z$_\odot$)', fontsize=18)
        ax[m, i].set_xticks([-3, -2, -1, 0, 1.])
        if m == 0:
            ax[m, i].legend(
                loc='lower left', 
                bbox_to_anchor=(-0.02, 1.01), 
                ncol=2, 
                borderaxespad=0, 
                frameon=False, 
                fontsize=15
            )
        ax[m, i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[m, i].yaxis.set_minor_locator(AutoMinorLocator(5))
        ax[m, i].tick_params(labelsize=15)

    ax[m, 0].set_ylabel(r'N$_{f_{\rm{GW}} \geq 10^{-4} \rm{Hz}}$ (Z) [10$^5$]', fontsize=18)
    
    if "fiducial" in model:
        ax[m, 0].set_yticks(np.arange(0.0, 2.5, 0.5))
        ax[m, 0].set_yticklabels(['0.00', '0.50', '1.00', '1.50', '2.00'])
        ax[m, 0].set_ylim(top=2.05)
        
        ax[m, 1].set_yticks(np.arange(0, 20, 4))
        ax[m, 1].set_ylim(top=17.25)
        ax[m, 1].set_yticklabels(np.arange(0, 20, 4).astype(float).astype(str))
        
        ax[m, 2].set_yticks(np.arange(0.0, 3.0, 0.75))
        #ax[m, 2].set_yticklabels(['0.00', '0.50', '1.00', '1.50', '2.00', '2.50'])
        
        ax[m, 3].set_yticks(np.arange(0., 1.0, 0.2))
        ax[m, 3].set_yticklabels(['0.00', '0.20', '0.40', '0.60', '0.80'])
        ax[m, 3].set_ylim(top=0.82)
    
    if "alpha25" in model:
        ax[m, 0].set_yticks(np.arange(0.0, 0.015, 0.003))
        ax[m, 0].set_ylim(top=0.0123)
        
        ax[m, 1].set_yticks(np.arange(0.0, 1.0, 0.2))
        ax[m, 1].set_yticklabels(['0.00', '0.20', '0.40', '0.60', '0.80'])
        ax[m, 1].set_ylim(top=0.83)
        
        ax[m, 2].set_yticks(np.arange(0.0, 0.4, 0.1))
        ax[m, 2].set_yticklabels(['0.00', '0.10', '0.20', '0.30'])
        
        ax[m, 3].yaxis.set_minor_locator(AutoMinorLocator(4))
        
    if "alpha5" in model:
        ax[m, 0].set_yticks([0, 10, 20, 30])
        ax[m, 0].set_yticklabels(['0.0', '10.0', '20.0', '30.0'])
        
        ax[m, 1].set_yticks([0, 20, 40, 60])
        ax[m, 1].set_yticklabels(['0.0', '20.0', '40.0', '60.0'])
        
        ax[m, 2].set_yticks(np.arange(0, 30, 7.5))
        #ax[m, 2].set_yticklabels(['0.0', '5.0', '10.0', 
        #                          '15.0', '20.0', '25.0'])
        ax[m, 2].yaxis.set_minor_locator(AutoMinorLocator(6))
        ax[m, 2].set_ylim(top=26)
        
        ax[m, 3].set_yticks(np.arange(0.0, 1.3, 0.3))
        ax[m, 3].set_yticklabels(['0.00', '0.30', '0.60', '0.90', '1.20'])
        ax[m, 3].yaxis.set_minor_locator(AutoMinorLocator(4))
        
    if "q3" in model:
        ax[m, 0].set_yticks(np.arange(0, 25, 5))
        ax[m, 0].set_yticklabels(['0.0', '5.0', '10.0', '15.0', '20.0'])
        
        ax[m, 1].set_yticks(np.arange(0, 60, 15))
        ax[m, 1].set_yticklabels(['0.0', '15.0', '30.0', '45.0'])
        ax[m, 1].set_ylim(top=50)
        
        ax[m, 2].set_yticks(np.arange(0, 30, 7.5))
        #ax[m, 2].set_yticklabels(['0.0', '5.0', '10.0', 
        #'15.0', '20.0', '25.0'])
        ax[m, 2].set_ylim(top=25)
        
        ax[m, 3].set_yticks(np.arange(0, 0.8, 0.25))
        ax[m, 3].yaxis.set_minor_locator(AutoMinorLocator(6))
        ax[m, 3].set_ylim(top=0.81)

plt.tight_layout()
plt.subplots_adjust(wspace=0.25)
plt.savefig(paths.figures / "lisa_nums.pdf".format(model), dpi=100)
