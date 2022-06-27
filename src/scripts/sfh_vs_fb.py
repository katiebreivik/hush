"""
Generate figure 1 of the paper based on galaxy m12i and the Moe+19 metallicity-dependent binary fraction
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from showyourwork.paths import user as Paths

# instantiate the paths
paths = Paths()


def get_FeH_from_Z(Z, Z_sun=0.02):
    """
    Converts from Z to FeH under the assumption that
    all stars have the same abundance as the sun

    INPUTS
    ----------------------
    Z [array]: array of metallicities to convert
    Z_sun [float]: solar metallicity

    RETURNS
    ----------------------
    Z [array]: array of FeH values
    """
    FeH = np.log10(Z) - np.log10(Z_sun)
    return FeH


def get_binfrac_of_Z(Z):
    """
    Calculates the theoretical binary fraction as a
    function of metallicity.

    INPUTS
    ----------------------
    Z [array]: metallicity Z values

    RETURNS
    ----------------------
    binfrac [array]: binary fraction values
    """
    FeH = get_FeH_from_Z(Z)
    FeH_low = FeH[np.where(FeH <= -1.0)]
    FeH_high = FeH[np.where(FeH > -1.0)]
    binfrac_low = -0.0648 * FeH_low + 0.3356
    binfrac_high = -0.1977 * FeH_high + 0.2025
    binfrac = np.append(binfrac_low, binfrac_high)
    return binfrac


rcParams["font.family"] = "serif"
rcParams["font.size"] = 14
rcParams["mathtext.default"] = "regular"

met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
met_arr = np.round(met_arr, 8)
met_arr = np.append(0.0, met_arr)

Z_sun = 0.02

FIRE = pd.read_hdf(paths.data / "FIRE.h5")
fig, ax = plt.subplots()
plt.grid(lw=0.25, which="both")
bins = np.append(met_arr[1:-1] / Z_sun, FIRE.met.max())
bins = np.append(FIRE.met.min(), bins)
bins = np.log10(bins)
ax2 = ax.twinx()
h, bins, _ = ax2.hist(
    np.log10(FIRE.met),
    bins=bins,
    histtype="step",
    lw=2,
    color="xkcd:tomato red",
    label="Latte m12i",
)
ax2.set_yscale("log")
ax2.legend(
    loc="lower left",
    bbox_to_anchor=(0.7, 1.01),
    ncol=4,
    borderaxespad=0,
    frameon=False,
)
ax.scatter(
    np.log10(met_arr[1:] / Z_sun),
    get_binfrac_of_Z(met_arr[1:]),
    color="k",
    s=15,
    zorder=2,
    label="COSMIC Z grid",
)
met_plot = np.linspace(FIRE.met.min() * Z_sun, FIRE.met.max() * Z_sun, 10000)
ax.plot(np.log10(met_plot / Z_sun), get_binfrac_of_Z(met_plot), color="k", label="FZ")
ax.set_xlim(bins[1] - 0.17693008, bins[-2] + 2 * 0.17693008)
ax.legend(
    loc="lower left",
    bbox_to_anchor=(-0.07, 1.01),
    ncol=5,
    borderaxespad=0,
    frameon=False,
)
   
ax.set_zorder(ax2.get_zorder() + 1)
ax.patch.set_visible(False)
ax.set_xlabel("Log$_{10}$(Z/Z$_\odot$)", size=18)
ax.set_ylabel("Binary Fraction", size=18)
ax2.set_ylabel(r"M$_{\rm{stars}}$ per Z bin (M$_\odot$)", size=18)
ax2.set_yticks([1e4, 1e5, 1e6, 1e7])
ax2.set_yticklabels(["7e7", "7e8", "7e9", "7e10"])
ax2.tick_params(labelsize=15)
ax.tick_params(labelsize=15)
ax2.set_ylim(1e4-0.1e4, 1e7+0.1e7)
plt.savefig(paths.figures / "sfh_vs_fb.pdf",  bbox_inches='tight', dpi=100)