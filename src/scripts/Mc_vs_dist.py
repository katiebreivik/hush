from matplotlib.ticker import AutoMinorLocator
import legwork.utils as utils
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
import seaborn as sns
import paths

models = ["fiducial", "alpha25", "alpha5", "q3"]
model_Zs = ["fiducial_Z", "alpha25_Z", "alpha5_Z", "q3_Z"]

fig, ax = plt.subplots(4, 4, figsize=(16,16))
levels = [0.05, 0.25, 0.50, 0.75]
label_y = [0.35, 0.49, 0.95, 1.6]
colors = ['#add0ed', '#2b5d87', '#4288c2', '#17334a']
figlabels = ['He + He', 'CO + He', 'CO + CO', 'ONe + X']

ii = 0
for model, model_Z in zip(models, model_Zs):

    resolved_dat_FZ = pd.read_hdf(
       paths.data / "results.hdf", key="resolved_DWDs_{}_{}".format("FZ", model_Z)
    
    )
    resolved_dat_FZ = resolved_dat_FZ.loc[resolved_dat_FZ.resolved_chirp == 1.0]
    resolved_dat_F50 = pd.read_hdf(
        paths.data / "results.hdf", key="resolved_DWDs_{}_{}".format("F50", model)
    )
    resolved_dat_F50 = resolved_dat_F50.loc[resolved_dat_F50.resolved_chirp == 1.0]
    
    Heplot = resolved_dat_FZ.loc[
        (resolved_dat_FZ.kstar_1 == 10) & (resolved_dat_FZ.kstar_2 == 10)
    ]
    COHeplot = resolved_dat_FZ.loc[
        ((resolved_dat_FZ.kstar_1 == 11) & (resolved_dat_FZ.kstar_2 == 10)) |
        ((resolved_dat_FZ.kstar_2 == 11) & (resolved_dat_FZ.kstar_1 == 10))
    ]
    COplot = resolved_dat_FZ.loc[
        (resolved_dat_FZ.kstar_1 == 11) & (resolved_dat_FZ.kstar_2 == 11)
    ]
    ONeplot = resolved_dat_FZ.loc[
        ((resolved_dat_FZ.kstar_1 == 12) & (resolved_dat_FZ.kstar_2.isin([10, 11, 12]))) |
        ((resolved_dat_FZ.kstar_2 == 12) & (resolved_dat_FZ.kstar_1.isin([10, 11, 12])))
    ]
    
    
    Heplot_F50 = resolved_dat_F50.loc[
        (resolved_dat_F50.kstar_1 == 10) & (resolved_dat_F50.kstar_2 == 10)
    ]
    COHeplot_F50 = resolved_dat_F50.loc[
        ((resolved_dat_F50.kstar_1 == 11) & (resolved_dat_F50.kstar_2 == 10)) |
        ((resolved_dat_F50.kstar_2 == 11) & (resolved_dat_F50.kstar_1 == 10))
    ]
    COplot_F50 = resolved_dat_F50.loc[
        (resolved_dat_F50.kstar_1 == 11) & (resolved_dat_F50.kstar_2 == 11)
    ]
    ONeplot_F50 = resolved_dat_F50.loc[
        ((resolved_dat_F50.kstar_1 == 12) & (resolved_dat_F50.kstar_2.isin([10, 11, 12]))) |
        ((resolved_dat_F50.kstar_2 == 12) & (resolved_dat_F50.kstar_1.isin([10, 11, 12])))
    ]
    
    dists = [x.dist_sun.values for x in [Heplot, COHeplot, COplot, ONeplot]]
    dists_F50 = [
        x.dist_sun.values for x in [Heplot_F50, COHeplot_F50, COplot_F50, ONeplot_F50]
    ]
    M_c = [
        utils.chirp_mass(x.mass_1.values * u.M_sun, x.mass_2.values * u.M_sun).value
        for x in [Heplot, COHeplot, COplot, ONeplot]
    ]
    M_c_F50 = [
        utils.chirp_mass(x.mass_1.values * u.M_sun, x.mass_2.values * u.M_sun).value
        for x in [Heplot_F50, COHeplot_F50, COplot_F50, ONeplot_F50]
    ]
    
    for dist, Mc, dist_F50, Mc_F50, jj in zip(dists[:1], M_c[:1], dists_F50[:1], M_c_F50[:1], range(len(dists))[:1]):
        sns.kdeplot(
            x=dist, 
            y=Mc, 
            fill=False, 
            ax=ax[ii, jj], 
            color=colors[0], 
            zorder=3, 
            linewidths=3.5, 
            label='FZ', 
            levels=levels
        )
        sns.kdeplot(
            x=dist_F50, 
            y=Mc_F50, 
            fill=False, 
            ax=ax[ii, jj], 
            color=colors[1], 
            zorder=6,
            linewidths=3.5, 
            linestyles='--', 
            label='F50', 
            levels=levels 
        )
        handles, labels = ax[ii, jj].get_legend_handles_labels(legend_handler_map=None)
        ax[ii, jj].legend(handles=handles,
                      labels=labels,
                      loc=(0, 1.01),
                      prop={'size': 15},
                      ncol=2, 
                      frameon=False)

    ax[ii, 0].set_ylabel('Chirp Mass [M$_\odot$]', fontsize=18)
    for i, name in zip(range(4), figlabels):
        if ii == 1:
            ax[ii, i].set_xlabel(r'Distance [kpc]', fontsize=18)
        ax[ii, i].text(0.05, 0.9, name, fontsize=18, horizontalalignment='left',
                   transform=ax[ii, i].transAxes)
        ax[ii, i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[ii, i].yaxis.set_minor_locator(AutoMinorLocator())
        ax[ii, i].tick_params(labelsize=15)
    
    
    for j in range(4):
        ax[ii, j].set_xlim(0, 30)
    
    if "fiducial" in model:
        ax[ii, 0].set_yticks(np.arange(0.2, 0.42, 0.05))
        ax[ii, 0].set_ylim(0.17, 0.36)
        
        ax[ii, 1].set_yticks([0.25, 0.3, 0.35, 0.4, 0.45, 0.5])
        ax[ii, 1].set_ylim(0.23, 0.525)
        
        ax[ii, 2].set_yticks(np.arange(0.25, 1.05, 0.15))
        ax[ii, 2].set_ylim(0.345, 0.95)
        
        ax[ii, 3].set_yticks(np.arange(0.3, 1.6, 0.3))
        ax[ii, 3].set_yticklabels(['0.30', '0.60', '0.90', '1.20', '1.50'])
        ax[ii, 3].set_ylim(0.175, 1.55)
    
    ii+=1
    print(ii)

plt.tight_layout()
#plt.subplots_adjust(wspace=0.25)

plt.savefig(paths.figures / "Mc_vs_dist.pdf", dpi=100)
