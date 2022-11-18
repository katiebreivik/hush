from matplotlib.ticker import AutoMinorLocator
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import paths
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

models = ["fiducial", "alpha25", "alpha5", "q3"]
model_labels = ["fiducial", "$\\alpha 25$", "$\\alpha 5$", "$q3$"]
metmodel = "F50"

FIREmin = 0.00015
FIREmax = 13.346
Z_sun = 0.02
num = 30
met_bins = np.logspace(np.log10(FIREmin), np.log10(FIREmax), num)
met_mids = (met_bins[1:] + met_bins[:-1]) / 2
whichsep = "CEsep"

fig, ax = plt.subplots(4, 4, figsize=(16, 15))

for m, model in enumerate(models):
    if "Z" in model:
        Heinter = pd.read_hdf(paths.data / "results.hdf", key = 'intersep_{}_{}_{}_{}'.format(10, 10, "FZ", model))
        COHeinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(11, 10,  "FZ", model))
        COinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(11, 11,  "FZ", model))
        ONeinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(12, 10,  "FZ", model))
    else:
        Heinter = pd.read_hdf(paths.data / "results.hdf", key = 'intersep_{}_{}_{}_{}'.format(10, 10, "F50", model))
        COHeinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(11, 10,  "F50", model))
        COinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(11, 11,  "F50", model))
        ONeinter = pd.read_hdf(paths.data / "results.hdf", key='intersep_{}_{}_{}_{}'.format(12, 10,  "F50", model))
    
    Heavgs = []
    Hecovs = []
    COHeavgs = []
    COHecovs = []
    COavgs = []
    COcovs = []
    ONeavgs = []
    ONecovs = []
    for i in range(num - 1):
        meti = met_bins[i]
        metf = met_bins[i + 1]
    
        Hebin = Heinter.loc[(Heinter.met>=meti)&(Heinter.met<=metf)]
        COHebin = COHeinter.loc[(COHeinter.met>=meti)&(COHeinter.met<=metf)]
        CObin = COinter.loc[(COinter.met>=meti)&(COinter.met<=metf)]
        ONebin = ONeinter.loc[(ONeinter.met>=meti)&(ONeinter.met<=metf)]
        if len(Hebin) != 0:
            Heavgs.append(np.mean(Hebin[whichsep].values))
            Hecovs.append(np.std(Hebin[whichsep].values))
        else:
            Heavgs.append(0.0)
            Hecovs.append(0.0)
    
        if len(COHebin) != 0:
            COHeavgs.append(np.mean(COHebin[whichsep].values))
            COHecovs.append(np.std(COHebin[whichsep].values))
        else:
            COHeavgs.append(0.0)
            COHecovs.append(0.0)
    
        if len(CObin) != 0:
            COavgs.append(np.mean(CObin[whichsep].values))
            COcovs.append(np.std(CObin[whichsep].values))
        else:
            COavgs.append(0.0)
            COcovs.append(0.0)
    
        if len(ONebin) != 0:
            ONeavgs.append(np.mean(ONebin[whichsep].values))
            ONecovs.append(np.std(ONebin[whichsep].values))
        else:
            ONeavgs.append(0.0)
            ONecovs.append(0.0)
    
    Heavgs = np.array(Heavgs)
    Hecovs = np.array(Hecovs)
    COHeavgs = np.array(COHeavgs)
    COHecovs = np.array(COHecovs)
    COavgs = np.array(COavgs)
    COcovs = np.array(COcovs)
    ONeavgs = np.array(ONeavgs)
    ONecovs = np.array(ONecovs)
    
    ax[m, 0].plot(
        np.log10(met_mids[Heavgs>0]), Heavgs[Heavgs>0]/1e3, 
        color='xkcd:tomato red', lw=3, ls='-', label='He + He', 
        drawstyle='steps-mid', rasterized=True
    )
                           
    ax[m, 0].fill_between(
        np.log10(met_mids[Heavgs>0]), (Heavgs[Heavgs>0]-Hecovs[Heavgs>0])/1e3, 
        (Heavgs[Heavgs>0]+Hecovs[Heavgs>0])/1e3, alpha=0.3, 
        color='xkcd:tomato red', zorder=0, step='mid', 
        label='$1\sigma$', rasterized=True
    )
                           
    ax[m, 2].plot(
        np.log10(met_mids[COavgs>0]), COavgs[COavgs>0]/1e3, 
        color='xkcd:pink', lw=3, ls='-', label='CO + CO', 
        drawstyle='steps-mid', rasterized=True
    )
    
    ax[m, 2].fill_between(
        np.log10(met_mids[COavgs>0]), (COavgs[COavgs>0]-COcovs[COavgs>0])/1e3, 
        (COavgs[COavgs>0]+COcovs[COavgs>0])/1e3, 
        alpha=0.3, color='xkcd:pink', zorder=0, step='mid', 
        label='$1\sigma$', rasterized=True
    )
    
    ax[m, 1].plot(
        np.log10(met_mids), COHeavgs/1e3, 
        color='xkcd:blurple', lw=3, ls='-', label='CO + He', 
        drawstyle='steps-mid', rasterized=True
    )
    
    ax[m, 1].fill_between(
        np.log10(met_mids[COHeavgs>0]), (COHeavgs[COHeavgs>0]-COHecovs[COHeavgs>0])/1e3, 
        (COHeavgs[COHeavgs>0]+COHecovs[COHeavgs>0])/1e3, 
        alpha=0.3, color='xkcd:blurple', zorder=0, step='mid', 
        label='$1\sigma$', rasterized=True
    )
    
    ax[m, 3].plot(
        np.log10(met_mids[ONeavgs>0]), ONeavgs[ONeavgs>0]/1e3, 
        color='xkcd:light blue', lw=3, label='ONe + X', 
        drawstyle='steps-mid', rasterized=True
    )
    
    ax[m, 3].fill_between(
        np.log10(met_mids[ONeavgs>0]), (ONeavgs[ONeavgs>0]-ONecovs[ONeavgs>0])/1e3,
        (ONeavgs[ONeavgs>0]+ONecovs[ONeavgs>0])/1e3, 
        alpha=0.3, color='xkcd:light blue', zorder=0, 
        step='mid', label='$1\sigma$', rasterized=True)
    
    for i in range(4):
        ax[m, i].text(0.05, 0.875, model_labels[m], fontsize=18, transform=ax[m, i].transAxes)
        ax[m, i].set_xticks([-3., -2., -1., 0., 1.])
        ax[m, i].tick_params(labelsize=15)
        ax[m, i].set_xlim(np.log10(met_mids[0]), np.log10(met_mids[-1]))
        if m == 3:
            ax[m, i].set_xlabel('Log$_{10}$(Z/Z$_\odot$)', fontsize=18)
        if m == 0:
            ax[m, i].legend(
                loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=2, 
                borderaxespad=0, frameon=False, fontsize=15, markerscale=0.5)
        ax[m, i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[m, i].yaxis.set_minor_locator(AutoMinorLocator())
        
    if whichsep == 'CEsep':
        ax[m, 0].set_ylabel(r'$\overline{\mathit{a}}_{\rm{CE}}$  [10$^3$ R$_\odot$]', fontsize=16)
    else:
        ax[m, 0].set_ylabel(r'$\overline{\mathit{a}}_{\rm{RLO}}$  [10$^3$ R$_\odot$]', fontsize=16)
    
    if "fiducial" in model:
        ax[m, 0].set_yticks(np.arange(0.1, 0.6, 0.1))
        ax[m, 0].set_ylim(0.09, 0.51)
        ax[m, 0].set_yticklabels(['0.10', '0.20', '0.30', '0.40', '0.50'])
        ax[m, 1].set_yticks(np.arange(0.25, 1.5, 0.25))
        ax[m, 1].set_ylim(top=1.28)
        ax[m, 2].set_yticks(np.arange(0.25, 2.75,0.5))
        ax[m, 2].set_ylim(top=2.3)
        ax[m, 3].set_yticks(np.arange(1.0, 3.5, 0.5))
        ax[m, 3].set_yticklabels(['1.00', '1.50', '2.00', '2.50', '3.00'])
        ax[m, 3].set_ylim(top=3.05)
        
    if "alpha25" in model:
        ax[m, 0].set_yticks(np.arange(0.15, 0.55, 0.1))
        ax[m, 0].set_yticklabels(['0.15', '0.25', '0.35', '0.45'])
        ax[m, 1].set_yticks(np.arange(0.5, 3.0, 0.75))
        ax[m, 1].set_ylim(top=2.9)
        ax[m, 2].set_yticks(np.arange(0.75, 4.0, 0.75))
        ax[m, 2].set_ylim(0.15, 3.4)
        ax[m, 3].set_yticks([1., 2., 3., 4.])
        ax[m, 3].set_yticklabels(['1.00', '2.00', '3.00', '4.00'])
        
    if "alpha5" in model:
        ax[m, 0].yaxis.set_minor_locator(AutoMinorLocator(3))
        ax[m, 1].set_yticks([0.0, 0.1, 0.2, 0.3])
        ax[m, 1].set_yticklabels(['0.00', '0.10', '0.20', '0.30'])
        ax[m, 2].set_yticks([0.0, 0.15, 0.3, 0.45])
        ax[m, 2].yaxis.set_minor_locator(AutoMinorLocator(4))
        ax[m, 3].set_yticks(np.arange(0.0, 0.9, 0.2))
        ax[m, 3].set_yticklabels(['0.00', '0.20', '0.40', '0.60', '0.80'])
        
    if "q3" in model:
        ax[m, 0].set_yticks(np.arange(0.04, 0.16, 0.02))
        ax[m, 0].set_ylim(top=0.145)
        ax[m, 1].set_yticks(np.arange(0.04, 0.2, 0.04))
        ax[m, 1].yaxis.set_minor_locator(AutoMinorLocator(5))
        ax[m, 3].set_yticks(np.arange(0.0, 1.5, 0.3))
        ax[m, 3].set_yticklabels(['0.00', '0.30', '0.60', '0.90', '1.20'])
        ax[m, 3].set_ylim(bottom=-0.05)

        
plt.tight_layout()

plt.savefig(paths.figures / "CEsep.pdf", dpi=100)