
import utils as dutil

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.time import Time
import astropy.constants as const
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.optimize import curve_fit
import tqdm
from schwimmbad import MultiPool

from legwork import psd, strain, utils
import legwork.source as source

pd.options.mode.chained_assignment = None


# Specific to Thiele et al. (2021), here are the used metallicity
# array, the associated binary fractions for each Z value, and the ratios 
# of mass in singles to mass in binaries of the Lband with each specific 
# binary fraction as found using COSMIC's independent samplers
# (See Binary_Fraction_Modeling.ipynb for Tutorials). All values were
# rounded to 4 significant digits except metallicity which used 8:

met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
met_arr = np.round(met_arr, 8)
met_arr = np.append(0.0, met_arr)

binfracs = np.array([0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                     0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                     0.3076, 0.2726, 0.2376, 0.2027, 0.1677])

ratios = np.array([0.68, 0.71, 0.74, 0.78, 0.82, 
                   0.86, 0.9, 0.94, 1.05, 1.22, 
                   1.44, 1.7 , 2.05, 2.51, 3.17])

ratio_05 = 0.64

# LEGWORK uses astropy units so we do also for consistency
G = const.G.value  # gravitational constant
c = const.c.value  # speed of light in m s^-1
M_sol = const.M_sun.value  # sun's mass in kg
R_sol = const.R_sun.value  # sun's radius in metres
sec_Myr = u.Myr.to('s')  # seconds in a million years
m_kpc = u.kpc.to('m')  # metres in a kiloparsec
L_sol = const.L_sun.value  # solar luminosity in Watts
Z_sun = 0.02  # solar metallicity
sun = coord.get_sun(Time("2021-04-23T00:00:00", scale='utc'))  # sun coordinates
sun_g = sun.transform_to(coord.Galactocentric)
sun_yGx = sun_g.galcen_distance.to('kpc').value
sun_zGx = sun_g.z.to('kpc').value
M_astro = 7070  # FIRE star particle mass in solar masses


#===================================================================================
# Lband and Evolution Functions:
#===================================================================================

def beta_(pop):
    '''
    Beta constant from page 8 of Peters(1964) used in the evolution 
    of DWDs due to gravitational waves.
    
    INPUTS
    ----------------------
    pop  [pandas dataframe]: DF of population which includes component 
                             masses in solar masses
    
    RETURNS
    ----------------------
    beta [array]: array of beta values
    '''
    m1 = pop.mass_1 * M_sol
    m2 = pop.mass_2 * M_sol
    beta = 64 / 5 * G ** 3 * m1 * m2 * (m1 + m2) / c ** 5
    return beta


def a_of_t(pop, t):
    '''
    Uses Peters(1964) equation (5.9) for circular binaries to find separation.
    as a function of time.

    INPUTS 
    ----------------------
    pop [pandas dataframe]: population subset from COSMIC. 
    t [array]: time at which to find separation. Must be in Myr.

    RETURNS
    ----------------------
    array of separation at time t in solar radii.
    '''
    t = t * sec_Myr
    beta = beta_(pop)
    a_i = pop.sep * R_sol
    a = (a_i ** 4 - 4 * beta * t) ** (1/4)
    return a / R_sol


def porb_of_a(pop, a):
    '''
    Converts semi-major axis "a" to orbital period using Kepler's equations.

    INPUTS
    ----------------------
    pop [pandas dataframe]: population from COSMIC. 
    a [array]:   semi-major axis of systems. Must be in solar radii and an array of 
         the same length as the dateframe pop.

    RETURNS
    t [array]:   orbital period in days.
    '''
    a = a * R_sol
    m1 = pop.mass_1 * M_sol
    m2 = pop.mass_2 * M_sol
    P_sqrd = 4 * np.pi ** 2 * a ** 3 / G / (m1 + m2)
    P = np.sqrt(P_sqrd)
    P = P / 3600 / 24  # converts from seconds to days
    return P


def t_of_a(pop, a):
    '''
    Finds time from SRF at which a binary would have a given separation after
    evolving due to gw radiation. (Re-arrangement of a_of_t(pop, t)).
    
    INPUTS
    ----------------------
    pop [pandas dataframe]: population subset from COSMIC.
    a [array]: separation to find time for. Must be in solar radii.

    RETURNS
    ----------------------
    t [array]: time in Myr where DWD reaches separation "a"
    '''
    beta = beta_(pop)
    a_i = pop.sep * R_sol
    a = a * R_sol
    t = (a_i ** 4 - a ** 4) / 4 / beta
    t = t / sec_Myr
    return t


def t_merge(pop):
    '''
    Uses Peters(1964) equation (5.10) to determine the merger time of a circular
    DWD binary from time of SRF.
    
    INPUTS
    ----------------------
    pop [pandas dataframe]: population subset from COSMIC

    RETURNS
    ----------------------
    t [array]: time in Myr.
    '''
    a_0 = pop.sep * R_sol
    beta = beta_(pop)
    T = a_0 ** 4 / 4 / beta
    T / sec_Myr
    return T


def a_of_RLOF(pop):
    '''
    Finds separation when lower mass WD overflows its
    Roche Lobe. Taken from Eq. 23 in "Binary evolution in a nutshell" 
    by Marc van der Sluys, which is an approximation of a fit
    done of Roche-lobe radius by Eggleton (1983).
    
    INPUTS
    ----------------------
    pop [pandas dataframe]: population subset from COSMIC
    
    RETURNS
    ----------------------
    a [array]: RLO separations of pop
    '''
    m1 = pop.mass_1
    m2 = pop.mass_2
    primary_mass = np.where(m1>m2, m1, m2)
    secondary_mass = np.where(m1>m2, m2, m1)
    secondary_radius = np.where(m1>m2, pop.rad_2, pop.rad_1)
    R2 = secondary_radius
    q = secondary_mass / primary_mass
    num = 0.49 * q ** (2/3)
    denom = 0.6 * q ** (2/3) + np.log(1 + q ** (1/3))
    a = denom * R2 / num
    return a


def random_sphere(R, num):
    '''
    Generates "num" number of random points within a
    sphere of radius R. It picks random x, y, z values
    within a cube and discards it if it's outside the
    sphere.

    INPUTS
    ----------------------
    R [array]: Radius in kpc
    num [int]: number of points to generate

    RETURNS
    ----------------------
    X, Y, Z arrays of length num
    '''
    X = []
    Y = []
    Z = []
    while len(X) < num:
        x = np.random.uniform(-R, R)
        y = np.random.uniform(-R, R)
        z = np.random.uniform(-R, R)
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        if r > R:
            continue
        if r <= R:
            X.append(x)
            Y.append(y)
            Z.append(z)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    return X, Y, Z


def rad_WD(M):
    '''
    Calculates the radius of a WD as a function of mass M in solar masses.
    Taken from Eq. 91 in Hurley et al. (2000), from Eq. 17 in Tout et al. (1997)
    
    INPUTS
    ----------------------
    M [array]: masses of the WDs in solar masses
    
    RETURNS
    ----------------------
    rad[array]: radii of the WDs in solar radii
    '''
    M_ch = 1.44
    R_NS = 1.4e-5*np.ones(len(M))
    A = 0.0115 * np.sqrt((M_ch/M)**(2/3) - (M/M_ch)**(2/3))
    rad = np.max(np.array([R_NS, A]), axis=0)
    return rad 

def evolve(pop_init):
    '''
    Evolve an initial population of binary WD's using
    GW radiation.
    
    INPUTS
    ----------------------
    pop_init [pandas dataframe]: initial population from COSMIC. 
                                 Must include assigned FIRE star
                                 particle age columns.
                                
    RETURNS
    ----------------------
    pop_init [pandas dataframe]: input pop with present-day parameter
                                 columns added with evolution time and
                                 present day separation, orbital period
                                 and GW frequency.
    
    '''
    t_evol = pop_init.age * 1000 - pop_init.tphys
    sep_f = a_of_t(pop_init, t_evol)
    porb_f = porb_of_a(pop_init, sep_f)
    f_gw = 2 / (porb_f * 24 * 3600)
    pop_init['t_evol'] = t_evol
    pop_init['sep_f'] = sep_f
    pop_init['porb_f'] = porb_f
    pop_init['f_gw'] = f_gw
    return pop_init


def position(pop_init):
    '''
    Assigning random microchanges to positions to
    give each system a unique position for identical
    FIRE star particles
    
    INPUTS
    ----------------------
    pop_init [pandas dataframe]: initial population from COSMIC. 
                                 Must include assigned FIRE star
                                 particle columns.
     
    RETURNS
    ----------------------
    pop_init [pandas dataframe]: input pop with columns added for
                                 galactocentric coordinates, and 
                                 Sun-to-DWD distance.                            
    '''
    R_list = pop_init.kern_len.values
    xGx = pop_init.xGx.values.copy()
    yGx = pop_init.yGx.values.copy()
    zGx = pop_init.zGx.values.copy()
    x, y, z = random_sphere(1.0, len(R_list))
    X = xGx + (x * R_list)
    Y = yGx + (y * R_list)
    Z = zGx + (z * R_list)
    pop_init['X'] = X
    pop_init['Y'] = Y
    pop_init['Z'] = Z
    pop_init['dist_sun'] = (X ** 2 + (Y - sun_yGx) ** 2 + (Z - sun_zGx) ** 2) ** (1/2)   
    return pop_init
  
    
def merging_pop(pop_init):
    '''
    Identifies DWD systems which will merge before present day,
    defined as those in which their delay time is less than their
    assigned FIRE star particle age.
    
    INPUTS
    ----------------------
    pop_init [pandas dataframe]:  initial population from COSMIC. 
                                  Must include assigned FIRE star
                                  particle age columns.
                                 
    RETURNS
    ----------------------
    pop_init [pandas dataframe]:  input pop with merged systems 
                                  discarded
    pop_merge [pandas dataframe]: merged population which can be
                                  saved separately 
    '''
    t_m = t_merge(pop_init)
    pop_init['t_delay'] = t_m + pop_init.tphys.values
    pop_merge = pop_init.loc[pop_init.t_delay <= pop_init.age * 1000]
    pop_init = pop_init.loc[pop_init.t_delay >= pop_init.age * 1000]
    return pop_init, pop_merge


def RLOF_pop(pop_init):
    '''
    Identifies DWD systems in which the lower mass WD will overflow
    its Roche Lobe before present day, i.e when the system's RLO time 
    is less than its assigned FIRE star particle age.
    
    INPUTS
    ----------------------
    pop_init [pandas dataframe]:  initial population from COSMIC. 
                                  Must include assigned FIRE star
                                  particle age columns.
                                 
    RETURNS
    ----------------------
    pop_init [pandas dataframe]:  input pop with merged systems 
                                  discarded
    pop_RLOF [pandas dataframe]: RLO population which can be
                                  saved separately 
    '''
    a_RLOF = a_of_RLOF(pop_init)
    t_RLOF = t_of_a(pop_init, a_RLOF)
    pop_init['t_RLOF'] = t_RLOF
    pop_RLOF = pop_init.loc[t_RLOF + pop_init.tphys <= pop_init.age * 1000]
    pop_init = pop_init.loc[t_RLOF + pop_init.tphys >= pop_init.age * 1000]
    return pop_init, pop_RLOF


def filter_population(dat):
    '''
    discards systems which have any of [formation times, delay times, RLOF times] 
    less than their FIRE age. Evolves the remaining systems to present day. Selects
    systems orbiting in the LISA band.
    
    INPUTS
    ----------------------
    dat [list] containing (in order)...
    - pop_init [pandas dataframe]: initial population from COSMIC. 
                                   Must include assigned FIRE star
                                   particle columns.
    - i [int]:          bin number for metallicity bin in [0, 15]
    - label [str]:      label for the DWD type for LISAband file names
    - ratio [float]:    ratio of mass in singles to mass in binaries formed for 
                        metallicity bin i
    - binfrac [float]:  binary fraction, either calculated from model FZ for bin i,
                        or 0.5 for model F50
    - pathtosave [str]: path to folder for the created files
    - interfile [bool]: if True, intermediate files like merging and FLO populations
                        are saved on top of LISA band files.
    
    OUTPUTS:
    ----------------------
    LISA_band [pandas dataframe]: evolved DWDs orbiting in the LISA freq. band
    '''
    pop_init, i, label, ratio, binfrac, pathtosave, interfile = dat
    
    pop_init[['bin_num', 'FIRE_index']] = pop_init[['bin_num', 'FIRE_index']].astype('int64')
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                     met_arr[i+1],
                                                                                     binfrac),
                                                   key='pop_init', format='t', append=True)    
    # Now that we've obtained an initial population, we make data cuts
    # of systems who wouldn't form in time for their FIRE age, or would
    # merge or overflow their Roche Lobe before present day.
    pop_init = pop_init.loc[pop_init.tphys <= pop_init.age * 1000]
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                     met_arr[i+1], 
                                                                                     binfrac), 
                                                    key='pop_age', format='t', append=True)
    
    pop_init, pop_merge = merging_pop(pop_init)
    if interfile == True:
        pop_merge[['bin_num', 'FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                      met_arr[i+1], 
                                                                                      binfrac), 
                                                    key='pop_merge', format='t', append=True)    
    
        pop_init[['bin_num', 'FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                     met_arr[i+1], 
                                                                                     binfrac), 
                                key='pop_nm', format='t', append=True)
    
    pop_merge = pd.DataFrame()
    pop_init, pop_RLOF = RLOF_pop(pop_init)
    
    if interfile == True:
        pop_RLOF[['bin_num','FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                           met_arr[i+1], 
                                                                                           binfrac), 
                                                  key='pop_RLOF', format='t', append=True)
    
        pop_init[['bin_num', 'FIRE_index']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                            met_arr[i+1], 
                                                                                            binfrac), 
                                key='pop_nRLOF', format='t', append=True)
    pop_RLOF = pd.DataFrame()
    
    # We now have a final population which we can evolve
    # using GW radiation
    pop_init = evolve(pop_init)
    
    # Assigning random microchanges to positions to
    # give each system a unique position for identical
    # FIRE star particles
    pop_init = position(pop_init)
    
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index', 'X', 'Y', 'Z']].to_hdf(pathtosave + 'Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                                    met_arr[i+1], 
                                                                                                    binfrac), 
                                                                  key='pop_f', format='t', append=True)    
    if binfrac == 0.5:
        binfrac_write = 0.5
    else:
        binfrac_write = 'variable'
    
    # Assigning weights to population to be used for histograms.
    # This creates an extra columns which states how many times
    # a given system was sampled from the cosmic-pop conv df.
    pop_init = pop_init.join(pop_init.groupby('bin_num')['bin_num'].size(), 
                             on='bin_num', rsuffix='_pw')
    
    # Systems detectable by LISA will be in the frequency band
    # between f_gw's 0.01mHz and 1Hz.
    LISA_band = pop_init.loc[(pop_init.f_gw >= 1e-4)]
    if len(LISA_band) == 0:
        print('No LISA sources for source {} and met {} and binfrac {}'.format(label, met_arr[i+1], binfrac))
        return []
    else:
        pop_init = pd.DataFrame()
        LISA_band = LISA_band.join(LISA_band.groupby('bin_num')['bin_num'].size(), 
                                   on='bin_num', rsuffix='_Lw')
            
        return LISA_band
    
def make_galaxy(dat, verbose=False):
    '''
    Creates populations of DWDs orbiting in the LISA band for a given
    DWD type and metallicity.
    
    INPUTS: 
    dat [list] containing (in order)...
    - pathtodat [str]:  path to COSMIC dat files with BPS DWD populations
    - fire_path [str]:  path to FIRE file with metallicity-dependent SFH data
    - pathtosave [str]: path to folder for the created galaxy files
    - filename [str]:   name of dat file for given DWD type and metallicity bin
    - i [int]:          bin number for metallicity bin in [0, 15]
    - label [str]:      label for the DWD type for LISAband file names
    - ratio [float]:    ratio of mass in singles to mass in binaries formed for 
                        metallicity bin i
    - binfrac [float]:  binary fraction, either calculated from model FZ for bin i,
                        or 0.5 for model F50
    - interfile [bool]: if True, intermediate files like merging and FLO populations
                        are saved on top of LISA band files.
    - nproc:            number of processes to allow if using on compute cluster
    
    OUTPUTS:
    No direct function outputs, but saves the following:
    - HDF file with LISA band systems
    - If interfile is True, HDF file with intermediate populations
    '''
    pathtodat, fire_path, pathtosave, filename, i, label, ratio, binfrac, interfile, model, nproc = dat
    
    if binfrac < 0.5:
        var_label = 'FZ'
    else:
        var_label = 'F50'
    Lkey = 'Lband_{}_{}'.format(var_label, model)
    Rkey = 'rand_seed_{}_{}'.format(var_label, model)
    Lsavefile = 'Lband_{}_{}_{}_{}.hdf'.format(label, var_label, model, i)
    FIREfile = 'FIRE.h5' 
    try:
        pd.read_hdf(pathtosave / Lsavefile, key=Lkey)
        
        return [], [], []
    except:
        FIRE = pd.read_hdf(fire_path / FIREfile).sort_values('met')
    
        rand_seed = np.random.randint(0, 100, 1)
        np.random.seed(rand_seed)
        
        rand_seed = pd.DataFrame(rand_seed)
        rand_seed.to_hdf(pathtosave / Lsavefile, key=Rkey)
    
        # Choose metallicity bin
        met_start = met_arr[i] / Z_sun
        met_end = met_arr[i+1] / Z_sun
        
        # Load DWD data at formation of the second DWD component
        conv = pd.read_hdf(pathtodat / filename, key='conv')
        
        # Limit to 1000 Rsun
        conv = conv.loc[conv.sep < 1000]
        
        # get bin_num indices if needed
        if 'bin_num' not in conv.columns:
            conv.index = conv.index.rename('index')
            conv['bin_num'] = conv.index.values
        
        # overwrite COSMIC radii
        conv['rad_1'] = rad_WD(conv.mass_1.values)
        conv['rad_2'] = rad_WD(conv.mass_2.values)    
        
        # Use ratio to scale to astrophysical pop w/ specific binary frac.
        try:
            mass_binaries = pd.read_hdf(pathtodat/filename, key='mass_stars').iloc[-1]
        except:
            print('m_binaries key')
            mass_binaries = pd.read_hdf(pathtodat/filename, key='mass_binaries').iloc[-1]
            
        mass_total = (1 + ratio) * mass_binaries  # total ZAMS mass of galaxy
        
        # Set up LISAband key to append to:
        final_params = ['bin_num','mass_1','mass_2','kstar_1','kstar_2','sep','met',
                            'tphys','rad_1','rad_2','xGx','yGx','zGx','FIRE_index','f_gw',
                            'dist_sun']
        d0 = pd.DataFrame(columns=final_params)
        d0.to_hdf(pathtosave / Lsavefile, key=Lkey, format='t', append=True)
        
        # Get DWD formatioon efficiency and number of binaries per star particle
        DWD_per_mass = len(conv) / mass_total
        N_astro = DWD_per_mass * M_astro
        
        # Choose FIRE bin based on metallicity:
        FIRE['FIRE_index'] = FIRE.index
        if met_end * Z_sun == met_arr[-1]:
            FIRE_bin = FIRE.loc[FIRE.met >= met_start]
        else:
            FIRE_bin = FIRE.loc[(FIRE.met >= met_start)&(FIRE.met <= met_end)]
        FIRE = []
        
        # We sample by the integer number of systems per star particle,
        # as well as a probabilistic approach for the fractional component
        # of N_astro:
        N_astro_dec = N_astro % 1
        p_DWD = np.random.rand(len(FIRE_bin))
        N_sample_dec = np.zeros(len(FIRE_bin))
        N_sample_dec[p_DWD <= N_astro_dec.values] = 1.0  # assign extra DWD to star particles
        num_sample_dec = int(N_sample_dec.sum())
        if verbose:
            print('we will sample {} stars from the decimal portion'.format(num_sample_dec))
        sample_dec = pd.DataFrame.sample(conv, num_sample_dec, replace=True)
        FIRE_bin_dec = FIRE_bin.loc[N_sample_dec == 1.0]
        
        params_list = ['bin_num', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'porb', 'sep', 
                'met', 'age', 'tphys', 'rad_1', 'rad_2', 'kern_len', 'xGx', 'yGx', 'zGx', 
                       'FIRE_index']
        
        pop_init_dec = pd.concat([sample_dec.reset_index(), FIRE_bin_dec.reset_index()], axis=1)
        sample_dec = pd.DataFrame()
        FIRE_bin_dec = pd.DataFrame()
        
        # get dat list and the population of DWDs orbiting in the LISA band for 
        # systems added from the decimal component of N_astro
        dat = [pop_init_dec[params_list], i, label, ratio, binfrac, pathtosave, interfile]
        LISA_band = filter_population(dat)
    
        if len(LISA_band) > 0:
            LISA_band = LISA_band[final_params]
            LISA_band.to_hdf(pathtosave / Lsavefile, key=Lkey, format='t', append=True)
        
        # now sampling by tthe integer number of systems per star particle:
        N_sample_int = int(N_astro) * len(FIRE_bin)
        if verbose:
            print('we will sample {} stars from the integer portion'.format(N_sample_int))
            print('getting FIRE values')
            
        FIRE_int = pd.DataFrame(np.repeat(FIRE_bin.values, int(N_astro), axis=0))
        FIRE_int.columns = FIRE_bin.columns
        FIRE_bin = pd.DataFrame()
        
        # if the number of populations to be sampled is large, we create galaxies iteratively
        # by looping through.
        Nsamp_split = 5e6
        if N_sample_int < Nsamp_split:
            sample_int = pd.DataFrame.sample(conv, N_sample_int, replace=True)
            pop_init_int = pd.concat([sample_int.reset_index(), 
                                      FIRE_int.reset_index()], axis=1)
            N = len(sample_int)
            sample_int = pd.DataFrame()
            FIRE_int = pd.DataFrame()
            dat = [pop_init_int[params_list], i, label, ratio, binfrac, pathtosave, interfile]
            LISA_band = filter_population(dat)
            
            if len(LISA_band) > 0:
                LISA_band = LISA_band[final_params]
                LISA_band.to_hdf(pathtosave / Lsavefile, key=Lkey, format='t', append=True)
        
        elif N_sample_int > Nsamp_split:
            if verbose:
                print('looping the integer population')
            N = 0
            j = 0
            jlast = int(Nsamp_split)
            dat_filter = []
            while j < N_sample_int:
                if verbose:
                    print('j: ', j)
                    print('jlast: ', jlast)
                    print('sampling {} systems'.format(int(jlast-j)))
                sample_int = pd.DataFrame.sample(conv, int(jlast-j), replace=True)
                N += len(sample_int)
                pop_init_int = pd.concat([sample_int.reset_index(), 
                                          FIRE_int.iloc[j:jlast].reset_index()], axis=1)
                dat_filter.append([pop_init_int[params_list], i, label, ratio, binfrac, pathtosave, interfile])
                j += Nsamp_split
                j = int(j)
                jlast += Nsamp_split
                jlast = int(jlast)
                if jlast > N_sample_int:
                    jlast = N_sample_int
            LISA_band_list = []
            if i < 13:
                for dat in dat_filter:
                    LISA_band_list.append(filter_population(dat))
                for LISA_band in LISA_band_list:
                    LISA_band = LISA_band[final_params]
                    LISA_band.to_hdf(pathtosave / Lsavefile, key=Lkey, format='t', append=True) 
            else:
                with MultiPool(processes=nproc) as pool:
                    LISA_band_list = list(pool.map(filter_population, dat_filter))
            
                for LISA_band in LISA_band_list:
                    LISA_band = LISA_band[final_params]
                    LISA_band.to_hdf(pathtosave / Lsavefile, key=Lkey, format='t', append=True)    
           
        if N != N_sample_int:
            print('loop is incorrect')
    
        FIRE_repeat = pd.DataFrame()
        sample_int = pd.DataFrame()
    
    return DWD_per_mass, label, i, binfrac


def save_full_galaxy(DWD_list, pathtodat, fire_path, pathtosave, interfile, model, nproc):
    '''
    Creates Milky Way-like galaxies of DWDs orbiting in the LISA band,
    over 15 metallicity bins and 4 DWD subtypes.
    
    INPUTS
    ----------------------
    DWD_list [list, str]: list of DWD types in [He_He, CO_He, CO_CO, ONe_X]
    pathtodat [str]:      path to COSMIC dat files with BPS DWD populations
    fire_path [str]:      path to FIRE file with metallicity-dependent SFH data
    pathtoLband [str]:    path to folder for the created galaxy files
    interfile [bool]:     if True, intermediate files like merging and FLO populations
                          are saved on top of LISA band files.
    nproc:                number of processes to allow if using on compute cluster
    
    RETURNS
    ----------------------
    No direct function outputs, but saves 60 HDF files with LISA band systems 
    for each metallicity bin and DWD type at pathtoLband.
    '''
    # Generate array of metallicities:
    
    met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
    met_arr = np.round(met_arr, 8)
    met_arr = np.append(0.0, met_arr)
    
    # Corresponding binary fractions and single-to-binaries mass ratios:
    # (These can be generated by get_binfrac_of_Z(Z) and get_ratios(binfracs)
    # functions)
    binfracs = np.array([0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                         0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                         0.3076, 0.2726, 0.2376, 0.2027, 0.1677])
    
    ratios = np.array([0.68, 0.71, 0.74, 0.78, 0.82, 
                       0.86, 0.9, 0.94, 1.05, 1.22, 
                       1.44, 1.7, 2.05, 2.51, 3.17])
    
    ratio_05 = 0.64

    
    # Run Code:
    # Run through all metallicities for metallicity-dependent
    # binary fraction and binary fraction of 0.5
    
    dat = []
    dat_long = []
    for DWD in DWD_list:
        if DWD == 'He_He':
            kstar1 = '10'
            kstar2 = '10'    
        elif DWD == 'CO_He':
            kstar1 = '11'
            kstar2 = '10'
        elif DWD == 'CO_CO':
            kstar1 = '11'
            kstar2 = '11'
        elif DWD == 'ONe_X':
            kstar1 = '12'
            kstar2 = '10_12'
            
        fnames, label = dutil.getfiles(kstar1=kstar1, kstar2=kstar2, model=model)
        i = 0
        for f, ratio, binfrac in zip(fnames, ratios, binfracs):
            if 'Z' in model:
                dat_list = [pathtodat, fire_path, pathtosave, f, i, label, 0.0, binfrac, interfile, model, nproc]
            else:
                dat_list = [pathtodat, fire_path, pathtosave, f, i, label, ratio_05, 0.5, interfile, model, nproc]

            if i < 13:
                dat.append(dat_list)
            else:
                dat_long.append(dat_list)
            i += 1
    
    with MultiPool(processes=nproc) as pool:
        _ = list(pool.map(make_galaxy, dat))
            
    for d in dat_long:
        _ = make_galaxy(d)
    
    return


def get_interactionsep_and_numLISA(pathtocosmic, pathtoLISA, pathtoresults, model, var, FIREmin=0.00015, FIREmax=13.346, Z_sun=0.02, verbose=False):
    '''
    Creates plot files with the interaction separation of DWDs. In order to
    run this, there must be already-existing LISA band files created with
    make_galaxy() at pathtoLband.
    
    INPUTS
    ----------------------
    pathtocosmic [str]: path to folder containing COSMIC datfiles 
    pathtoLISA [str]: path to folder containing LISA band files
    pathtoresults [str]: path to folder to save results to
    
    RETURNS
    ----------------------
    No direct function outputs, but saves interraction separation data for all
    DWD types and metallicity bins to files contained in pathtosave.
    '''
    def intersep(pathtocosmic, datfile, pathtoresults, interkey, LISA_data, result_file):
        if verbose:
            print('dat file: ' + datfile)
        dat = pd.read_hdf(pathtocosmic / datfile, key='bpp') 
        if 'bin_num' not in dat.columns:
            dat.index = dat.index.rename('index')
            dat['bin_num'] = dat.index.values
        dat = dat[['tphys', 'evol_type', 'sep', 'mass_1', 'mass_2', 'porb', 'bin_num']]
    
        RLOFsep = dat.loc[dat.evol_type==3].groupby('bin_num', as_index=False).first()
        RLOFsep = RLOFsep.loc[RLOFsep.bin_num.isin(LISA_data.bin_num)]
        data_RLOF = LISA_data.loc[LISA_data.bin_num.isin(RLOFsep.bin_num)]
        RLOFsep['weights'] = data_RLOF.bin_num.value_counts().sort_index().values
        
        CEsep = dat.loc[dat.evol_type==7].groupby('bin_num', as_index=False).first()
        CEsep = CEsep.loc[CEsep.bin_num.isin(LISA_data.bin_num)]
        data_CE = LISA_data.loc[LISA_data.bin_num.isin(CEsep.bin_num)]
        CEsep['weights'] = data_CE.bin_num.value_counts().sort_index().values
    
        dat = []
        data_RLOF['CEsep'] = np.repeat(CEsep['sep'], CEsep['weights']).values
        data_RLOF['CEporb'] = np.repeat(CEsep['porb'], CEsep['weights']).values
        data_RLOF['CEmass1'] = np.repeat(CEsep['mass_1'], CEsep['weights']).values
        data_RLOF['CEmass2'] = np.repeat(CEsep['mass_2'], CEsep['weights']).values
        data_RLOF['CEtime'] = np.repeat(CEsep['tphys'], CEsep['weights']).values
        data_RLOF['RLOFsep'] = np.repeat(RLOFsep['sep'], RLOFsep['weights']).values
        data_RLOF['RLOFtime'] = np.repeat(RLOFsep['tphys'], RLOFsep['weights']).values
        
        data_RLOF = data_RLOF[['CEsep', 'CEporb', 'CEmass1', 'CEmass2', 'CEtime', 'RLOFsep', 'RLOFtime', 'met']]    
        print('CEsep results path: ', pathtoresults)
        data_RLOF.to_hdf(pathtoresults / result_file, key=interkey, append=True) 
        return
    

    
    # FZ:
    kstar1_list = ['10', '11', '11', '12']
    kstar2_list = ['10', '10', '11', '10_12']
    k1_ints = [10, 11, 11, 12]
    k2_ints = [10, 10, 11, 10]
    if var:
        var_label = 'FZ'
    else:
        var_label = 'F50'
        
    result_file = 'results_{}_{}.hdf'.format(var_label, model)
    Lkey = 'Lband_{}_{}'.format(var_label, model)
    
    
    
    # get number of LISA sources per metallicity
    num = 30
    met_bins = np.logspace(np.log10(FIREmin), np.log10(FIREmax), num)*Z_sun
    num_LISA = []
    
    for kstar1, kstar2, k1, k2 in zip(kstar1_list, kstar2_list, k1_ints, k2_ints):
        files, label = dutil.getfiles(kstar1, kstar2, model)
        mets = []
        for f, i in zip(files, range(len(files))):
            Lsavefile = 'Lband_{}_{}_{}_{}.hdf'.format(label, var_label, model, i)
            # read in LISA data
            LISA_data = pd.read_hdf(pathtoLISA / Lsavefile, key=Lkey)
    
            if verbose:
                print('i = {}'.format(i))
            binfrac = binfracs[i]
            met = met_arr[i+1]
            
            if len(LISA_data) == 0:
                continue
            else:
                interkey = 'intersep_{}_{}_{}_{}'.format(k1, k2, var_label, model)
                _ = intersep(pathtocosmic, f, pathtoresults, interkey, LISA_data, result_file)
                
                mets.extend(LISA_data.met.values)
        
        mets = np.array(mets)
        nums, bins = np.histogram(mets*Z_sun, bins=met_bins)
        
        num_LISA.append(nums)
    
    
    numLISA_30bins = pd.DataFrame(np.array(num_LISA).T, 
                                  columns=['He', 'COHe', 'CO', 'ONe'])
    if var:
        key = 'numLISA_30bins_{}_{}'.format("FZ", model)
    else:
        key = 'numLISA_30bins_{}_{}'.format("F50", model)
    print('numLISA results path: ', pathtoresults)

    numLISA_30bins.to_hdf(pathtoresults / result_file, key=key)
    
    return
    
def get_resolvedDWDs(pathtoLISA, pathtosave, var, model, window):
    '''
    Creates plot files with the DWDs that have SNR>7, and notes whether or
    not they are chirping. In order to run this, there must be already-existing 
    LISA band files created with make_galaxy() at pathtoLband.
    
    INPUTS
    ----------------------
    pathtoLband [str]: path to folder containing LISA band files
    pathtosave [str]: path to folder to save data to
    var [bool]: if True, uses model FZ. if False, model F50.
    
    RETURNS
    ----------------------
    No direct function outputs, but saves resolved DWD data for all
    DWD types and metallicity bins to files contained in pathtosave.
    '''
    def func(x, a, b, c, d, e):
        return a + b*x + c*x**2 + d*x**3 + e*x**4
    
    def cosmic_confusion_var_fiducial(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('FZ', 'fiducial_Z'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_var_alpha25(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('FZ', 'alpha25_Z'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_var_alpha5(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('FZ', 'alpha5_Z'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_var_q3(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('FZ', 'q3_Z'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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

    def cosmic_confusion_50_fiducial(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('F50', 'fiducial'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) &
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_50_alpha25(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('F50', 'alpha25'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_50_alpha5(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('F50', 'alpha5'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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
    
    def cosmic_confusion_50_q3(f, L, t_obs=4 * u.yr, approximate_R=True, confusion_noise=None):
        power_dat = pd.read_hdf(pathtoLISA / result_file, key='total_power_DWDs_{}_{}'.format('F50', 'q3'))

        power_dat_median = power_dat.rolling(window).median()
        power_dat_median = power_dat_median[window:]

        power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & 
                                                    (power_dat_median.f_gw <= 1.2e-3)]

        popt, pcov = curve_fit(func,
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))

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

    kstar1_list = ['10', '11', '11', '12']
    kstar2_list = ['10', '10', '11', '10_12']
    
    lisa_bins = np.arange(1e-9, 1e-1, 1/(4 * 3.155e7))
    Tobs = 4 * u.yr
    
    if var:
        Rkey = 'resolved_DWDs_{}_{}'.format('FZ', model)
        Pkey = 'total_power_DWDs_{}_{}'.format('FZ', model)
        Ckey = 'conf_fit_DWDs_{}_{}'.format('FZ', model)
        result_file = 'results_{}_{}.hdf'.format('FZ', model)
    else:
        Rkey = 'resolved_DWDs_{}_{}'.format('F50', model)
        Pkey = 'total_power_DWDs_{}_{}'.format('F50', model)
        Ckey = 'conf_fit_DWDs_{}_{}'.format('F50', model)
        result_file = 'results_{}_{}.hdf'.format('F50', model)
    
    if var:
        var_label = 'FZ'
    else:
        var_label = 'F50'
    Lkey = 'Lband_{}_{}'.format(var_label, model)
    
    dat = []
    for kstar1, kstar2 in zip(kstar1_list, kstar2_list):
        for i in range(len(met_arr)-1):
            if kstar1 == '12':
                label='12'
            else:
                label='{}_{}'.format(kstar1, kstar2)
            
            Lsavefile = 'Lband_{}_{}_{}_{}.hdf'.format(label, var_label, model, i)
            if len(dat) == 0:
                Ldat = pd.read_hdf(pathtoLISA / Lsavefile, key=Lkey)
                dat = Ldat[['mass_1', 'mass_2', 'dist_sun', 'f_gw', 'kstar_1', 'kstar_2']]
            else:
                Ldat = pd.read_hdf(pathtoLISA / Lsavefile, key=Lkey)
                dat = dat.append(Ldat[['mass_1', 'mass_2', 'dist_sun', 'f_gw', 'kstar_1', 'kstar_2']])
    
    sources = source.Source(m_1=dat.mass_1.values * u.Msun, 
                            m_2=dat.mass_2.values * u.Msun,  
                            ecc=np.zeros(len(dat.mass_1)), 
                            dist=dat.dist_sun.values * u.kpc, 
                            f_orb=dat.f_gw.values/2 * u.Hz,
                            interpolate_g=True, 
                            interpolate_sc=True, 
                            sc_params={"instrument": "LISA",
                                       "t_obs": Tobs,
                                       "L": 2.5e9 * u.m,
                                       "approximate_R": True,
                                       "confusion_noise": None})
    
    dat['h_0'] = sources.get_h_0_n(harmonics=[2])
    dat['power'] = sources.get_h_0_n(harmonics=[2])**2
    dat['digits'] = np.digitize(dat.f_gw, lisa_bins)
    sources = []    
        
    power = dat.groupby('digits').power.sum()
    power_foreground = np.zeros(len(lisa_bins))
    power_foreground[np.array(power.index.astype(int))] = power
    
    power_dat = pd.DataFrame(np.vstack([lisa_bins, power_foreground]).T, 
                             columns=['f_gw', 'strain_2'])
    
    power_dat_median = power_dat.rolling(window).median()
    power_dat_median = power_dat_median[window:]
    power_dat_downsample = power_dat.sample(int(len(power_dat)/10), replace=False)
    power_dat_downsample = power_dat_downsample.sort_values('f_gw')
    power_dat_downsample.to_hdf(pathtosave / result_file, key=Pkey)
    #power_dat.to_hdf(pathtoLISA / result_file, key=Pkey)
    power_dat = []
    
    dat = dat[['mass_1', 'mass_2', 'dist_sun', 'f_gw', 'kstar_1', 'kstar_2', 'h_0']]
    power_dat_median_fit = power_dat_median.loc[(power_dat_median.strain_2 > 0) & (power_dat_median.f_gw <= 1.2e-3)]
    
    popt, pcov = curve_fit(func, 
                           xdata=np.log10(power_dat_median_fit.f_gw.values),
                           ydata=np.log10(power_dat_median_fit.strain_2.values))
    
    
    if var:
        if model == 'fiducial_Z':
            cosmic_confusion = cosmic_confusion_var_fiducial
        elif model == 'alpha25_Z':
            cosmic_confusion = cosmic_confusion_var_alpha25
        elif model == 'alpha5_Z':
            cosmic_confusion = cosmic_confusion_var_alpha5
        elif model == 'q3_Z':
            cosmic_confusion = cosmic_confusion_var_q3
    else:
        if model == 'fiducial':
            cosmic_confusion = cosmic_confusion_50_fiducial
        elif model == 'alpha25':
            cosmic_confusion = cosmic_confusion_50_alpha25
        elif model == 'alpha5':
            cosmic_confusion = cosmic_confusion_50_alpha5
        elif model == 'q3':
            cosmic_confusion = cosmic_confusion_50_q3
            
    
    n_cut = int(len(dat)/2)
    snrs = []
    chirps = []
    
    for nstart, nstop in zip([0,n_cut], [n_cut, len(dat)]):
    
        sources_conf = source.Source(m_1=dat.mass_1.values[nstart:nstop] * u.Msun, 
                                     m_2=dat.mass_2.values[nstart:nstop] * u.Msun,  
                                     ecc=np.zeros(len(dat.mass_1[nstart:nstop])), 
                                     dist=dat.dist_sun.values[nstart:nstop] * u.kpc, 
                                     f_orb=dat.f_gw.values[nstart:nstop]/2 * u.Hz,
                                     stat_tol = 1/(Tobs.to(u.s).value),
                                     interpolate_g=True, 
                                     interpolate_sc=False, 
                                     sc_params={"instrument": "custom",
                                                "t_obs": Tobs,
                                                "L": 2.5e9 * u.m,
                                                "approximate_R": True,
                                                "confusion_noise": None,
                                                "custom_psd": cosmic_confusion})
        snr = sources_conf.get_snr(t_obs=Tobs, verbose=False, custom_psd=cosmic_confusion)
        snrs.extend(snr)
    dat['snr'] = snrs
    
    dat['chirp'] = utils.fn_dot(utils.chirp_mass(dat.mass_1.values * u.Msun, dat.mass_2.values * u.Msun), dat.f_gw.values/2 * u.Hz, np.zeros(len(dat)), n=2).to(u.s**(-2)).value
    dat = dat.loc[dat.snr > 7]
    dat['resolved_chirp'] = np.zeros(len(dat))
    dat.loc[dat.chirp > 1/((Tobs.to(u.s).value)**2), 'resolved_chirp'] = 1.0
    
    dat.to_hdf(pathtosave / result_file, key=Rkey)    
    sources_conf = []
    print('sources written')
    
    dat = []
    popt_df = pd.DataFrame(popt)
    popt_df.to_hdf(pathtosave / result_file, key=Ckey)
     
    return


def get_formeff(pathtodat, pathtosave, model, var):
    '''
    Creates plot files with formation efficiency of DWDs. In order to
    run this, there must be already-existing LISA band files created with
    make_galaxy() at pathtoLband.

    INPUTS
    ----------------------
    pathtodat [str]: path to folder containing datfiles
    pathtosave [str]: path to folder to save formation efficiency information to

    RETURNS
    ----------------------
    No direct function outputs, but saves formation efficiency data for all
    DWD types and metallicity bins to files contained in pathtosave.
    '''

    def formeff(datfiles, pathtodat, model):
        ratios = np.array([0.68, 0.71, 0.74, 0.78, 0.82,
                           0.86, 0.9, 0.94, 1.05, 1.22,
                           1.44, 1.7, 2.05, 2.51, 3.17])

        ratio_05 = 0.64

        lenconv = []
        masslist = []
        for i in range(15):
            if 'Z' in model:
                # this was covered in the simulations!
                ratio = 0
            else:
                ratio = ratio_05
            mass_binaries = pd.read_hdf(pathtodat / datfiles[i], key='mass_stars').iloc[-1]
            mass = (1 + ratio) * mass_binaries
            conv = pd.read_hdf(pathtodat / datfiles[i], key='conv')
            conv = conv.loc[conv.sep < 1000]
            masslist.append(mass)
            lenconv.append(len(conv))

        lenconv = np.array(lenconv)
        masslist = np.concatenate(np.array(masslist))
        eff = lenconv / masslist
        return eff

    kstar1_list = ['10', '11', '11', '12']
    kstar2_list = ['10', '10', '11', '10_12']
    labels = ['10_10', '11_10', '11_11', '12']

    eff = []
    for kstar1, kstar2, label in tqdm.tqdm(zip(kstar1_list, kstar2_list, labels)):
        files, lab = dutil.getfiles(kstar1=kstar1, kstar2=kstar2, model=model)
        eff.append(formeff(files, pathtodat, model))
        #print('finished {}'.format(label))
    
    if var:
        var_label='FZ'
    else:
        var_label='F50'
    result_file = 'results_{}_{}.hdf'.format(var_label, model)
    DWDeff = pd.DataFrame(np.array(eff).T, columns=['He', 'COHe', 'CO', 'ONe'])
    print('DWDEff results save', pathtosave)
    DWDeff.to_hdf(pathtosave / result_file, key='DWDeff_{}'.format(model))
    
    return
