import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u


def func(x, a, b, c, d, e):
    return a + b * x + c * x ** 2 + d * x ** 3 + e * x ** 4

model = "fiducial"
colors = ["#add0ed", "#2b5d87", "#4288c2", "#17334a"]
Tobs = 4 * u.yr

power_dat_F50 = pd.read_hdf(
    "../data/results.hdf", key="total_power_DWDs_{}_{}".format("F50", model)
)
popt_F50 = pd.read_hdf("../data/results.hdf",  key="conf_fit_DWDs_{}_{}".format("F50", model))
popt_F50 = popt_F50.values.flatten()

power_dat_FZ = pd.read_hdf(
    "../data/results.hdf", key="total_power_DWDs_{}_{}".format("FZ", model)
)
popt_FZ = pd.read_hdf("../data/results.hdf", key="conf_fit_DWDs_{}_{}".format("FZ", model))
popt_FZ = popt_FZ.values.flatten()

conf_fit_FZ = (
    10
    ** func(
        x=np.log10(np.linspace(1e-4, 1e-1, 100000)),
        a=popt_FZ[0],
        b=popt_FZ[1],
        c=popt_FZ[2],
        d=popt_FZ[3],
        e=popt_FZ[4],
    )
    * Tobs.to(u.s).value
)

conf_fit_F50 = (
    10
    ** func(
        x=np.log10(np.linspace(1e-4, 1e-1, 100000)),
        a=popt_F50[0],
        b=popt_F50[1],
        c=popt_F50[2],
        d=popt_F50[3],
        e=popt_F50[4],
    )
    * Tobs.to(u.s).value
)


fig, (ax1) = plt.subplots(1, 1, figsize=(6, 4.2))
plt.plot(
    power_dat_F50.f_gw[::10],
    power_dat_F50.strain_2[::10] * Tobs.to(u.s).value,
    c=colors[1],
    lw=1,
    alpha=1,
    rasterized=True,
)
plt.plot(
    power_dat_FZ.f_gw[::10],
    power_dat_FZ.strain_2[::10] * Tobs.to(u.s).value,
    c=colors[0],
    lw=1,
    alpha=0.8,
    rasterized=True,
)
plt.plot(
    np.linspace(1e-4, 1e-1, 100000),
    conf_fit_F50,
    c=colors[3],
    ls="--",
    lw=2,
    label=r"F50",
)
plt.plot(
    np.linspace(1e-4, 1e-1, 100000),
    conf_fit_FZ,
    c=colors[2],
    ls="--",
    lw=2,
    label=r"FZ",
)
plt.xscale("log")
plt.yscale("log")

ax1.set_ylabel(r"PSD [Hz$^{-1}$]", size=15)
ax1.set_xlabel(r"f$_{\rm{GW}}$ [Hz]", size=15)
ax1.tick_params(labelsize=12)
ax1.set_yticks([1e-38, 1e-37, 1e-36, 1e-35, 1e-34])
plt.xlim(1e-4, 3e-2)
plt.ylim(1e-38, 5e-34)
plt.legend(prop={"size": 12}, ncol=2, frameon=False, loc=(0, 1))
plt.tight_layout()
plt.savefig("PSD.pdf", dpi=100)
ax1.set_ylabel(r'PSD [Hz$^{-1}$]', size=15)
