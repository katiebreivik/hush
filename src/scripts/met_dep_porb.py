from scipy.stats import gaussian_kde, norm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import interp1d
from dutils import get_Z_from_FeH, get_FeH_from_Z, get_binfrac_of_Z
import scipy.integrate as integrate
from matplotlib import rcParams
import paths


def get_porb_kde(Z, close_logP=4.0, wide_logP=6.0, binfrac_tot_solar=0.53, Z_sun=0.02):
    '''Returns a kde with weights consistent with Fig 19 of Moe+19
    for the orbital period distribution
    
    INPUTS
    ------------------------
    Z [array]: metallicity values
        
    close_logP [float]: divding line beween close and intermediate orbits
    
    wide_logP [float]: dividing line between intermediate and wide orbits
    
    binfrac_tot [float]: integrated total binary fraction at solar metallicity
    
    RETURNS
    ------------------------
    logp [array]: orbital periods drawn from the distribution
    '''
    # import the interpolation package
    from scipy.interpolate import interp1d
    from scipy.stats import norm
    from scipy.integrate import trapz
    
    # fix to values used in Moe+21
    logP_lo_lim=0
    logP_hi_lim=9
    log_P = np.linspace(logP_lo_lim, logP_hi_lim, 10000)
    
    # hard coded from Moe+21
    logP_pdf = norm.pdf(log_P, loc=4.4, scale=2.1)
    
    # set up the wide binary fraction inflection point
    norm_wide = binfrac_tot_solar/trapz(logP_pdf, log_P)
    
    # set up the close binary fraction inflection point
    FeHclose = np.linspace(-3.0, 0.5, 100)
    fclose = -0.0648 * FeHclose + 0.3356
    fclose[FeHclose > -1.0] = -0.1977 * FeHclose[FeHclose > -1.0] + 0.2025

    Zclose = get_Z_from_FeH(FeHclose, Z_sun=Z_sun)
    
    fclose_interp = interp1d(Zclose, fclose)
    
    fclose_Z = fclose_interp(Z)
    norm_close_Z = fclose_Z/trapz(logP_pdf[log_P < close_logP], log_P[log_P < close_logP])
    
    return norm_wide, norm_close_Z
   
def get_logP_dist(nsamp, norm_wide, norm_close):
    prob_wide = norm.pdf(np.linspace(wide_logP, logP_hi_lim, neval), loc=4.4, scale=2.1)*norm_wide
    prob_close = norm.pdf(np.linspace(logP_lo_lim, close_logP, neval), loc=4.4, scale=2.1)*norm_close
    slope = -(prob_close[-1] - prob_wide[0]) / (wide_logP - close_logP)
    prob_intermediate = slope * (np.linspace(close_logP, wide_logP, neval) - close_logP) + prob_close[-1]
    prob_interp_int = interp1d(np.linspace(close_logP, wide_logP, neval), prob_intermediate)

    log_p_success = []
    n_success = 0
    while n_success < nsamp:
        logP_samp = np.random.uniform(logP_lo_lim, logP_hi_lim, nsamp*5)
        logP_prob = np.random.uniform(0, 1, nsamp*5)
        
        logP_samp_lo = logP_samp[logP_samp<close_logP]
        logP_prob_lo = logP_prob[logP_samp<close_logP]
        log_p_success.extend(logP_samp_lo[np.where(logP_prob_lo < norm.pdf(logP_samp_lo, loc=4.4, scale=2.1)*norm_close)])
        
        logP_samp_int = logP_samp[(logP_samp>=close_logP) & (logP_samp<wide_logP)]
        logP_prob_int = logP_prob[(logP_samp>=close_logP) & (logP_samp<wide_logP)]
        log_p_success.extend(logP_samp_int[np.where(logP_prob_int < prob_interp_int(logP_samp_int))])
    
        logP_samp_hi = logP_samp[(logP_samp>=wide_logP)]
        logP_prob_hi = logP_prob[(logP_samp>=wide_logP)]

        log_p_success.extend(logP_samp_hi[np.where(logP_prob_hi < norm.pdf(logP_samp_hi, loc=4.4, scale=2.1)*norm_wide)])

        n_success = len(log_p_success)
    log_p_success = np.array(log_p_success)[np.random.randint(0,n_success,nsamp)]    
    return log_p_success




met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
met_arr = np.round(met_arr, 8)

logP_lo_lim=0
logP_hi_lim=9
close_logP=4.0
wide_logP=6.0
Z_sun=0.02
colors = sns.color_palette("mako", n_colors=len(met_arr))

FeH_list = np.log10(met_arr/Z_sun)
norm_wide_list = []
norm_close_list = []
for Z in met_arr:
    norm_wide, norm_close_Z = get_porb_kde(Z)
    norm_wide_list.append(norm_wide)
    norm_close_list.append(norm_close_Z)
    
    
rcParams["font.family"] = "serif"
rcParams["font.size"] = 14
rcParams["mathtext.default"] = "regular"


logP_lo_lim=0
logP_hi_lim=9
log_P = np.linspace(logP_lo_lim, logP_hi_lim, 10000)

fig, ax = plt.subplots(figsize=(7.5,5))

for norm_wide, norm_close, Z, FeH, color in zip(norm_wide_list, norm_close_list, met_arr, FeH_list, colors):
    neval = 500
    prob_wide = norm.pdf(np.linspace(wide_logP, logP_hi_lim, neval), loc=4.4, scale=2.1)*norm_wide
    prob_close = norm.pdf(np.linspace(logP_lo_lim, close_logP, neval), loc=4.4, scale=2.1)*norm_close
    slope = -(prob_close[-1] - prob_wide[0]) / (wide_logP - close_logP)
    prob_intermediate = slope * (np.linspace(close_logP, wide_logP, neval) - close_logP) + prob_close[-1]
    prob_interp_int = interp1d(np.linspace(close_logP, wide_logP, neval), prob_intermediate)

    x_dat = np.hstack([np.linspace(logP_lo_lim, close_logP, neval),
                       np.linspace(close_logP, wide_logP, neval),
                       np.linspace(wide_logP, logP_hi_lim, neval)])
    y_dat = np.hstack([prob_close, 
                       prob_interp_int(np.linspace(close_logP, wide_logP, neval)), 
                       prob_wide])
    
    if (np.round(FeH, 2) == -2.3) or (np.round(FeH, 2) == 0.18):
        plt.plot(x_dat, y_dat, c=color, label='[Fe/H] = {}'.format(np.round(FeH, 2)))
    else:
        plt.plot(x_dat, y_dat, c=color)

plt.legend(ncol=2, loc=[0, 1.01], frameon=False, prop={'size':16})
plt.tick_params(labelsize=15)
plt.xlabel(r'log$_{10}(P_{\rm{orb}}/day)$', size=18)
plt.ylabel(r'$f_{\rm{logP}}$', size=18)
plt.xlim(0, 8.5)
plt.ylim(0, 0.23)
plt.tight_layout()
plt.savefig(paths.figures / 'met_dep_porb.pdf', facecolor='white', dpi=180)
    
