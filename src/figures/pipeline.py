#=========================================================================
# This script creates LISA DWD galaxies across 15 metallicity bins,
# incorporating the metallicity-dependent binary fraction as
# discussed in Thiele et al. (2021).
#
# Authors: Sarah Thiele & Katelyn Breivik
# Last updated: Oct 14th, 2021
#=========================================================================

import numpy as np
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import argparse
import postproc as pp


DWD_list = ['He_He', 'CO_He', 'CO_CO', 'ONe_X']
dat_path = "../data/cosmic_dat/"
LISA_path = "../data/postproc/"
FIRE_path = "../data/"
results_path = "../data/"
models = ['fiducial', 'fiducial_Z', 'alpha25', 'alpha25_Z', 'alpha5', 'alpha5_Z', 'q3', 'q3_Z']
interfile = False
nproc = 128

for model in models:
    print(model)
    #pp.save_full_galaxy(
    #    DWD_list, dat_path, FIRE_path, LISA_path, interfile, model, nproc
    #)
    #print('Gx done!')

    if model in ['fiducial_Z', 'fiducial']:
        pp.get_formeff(
            pathtodat=dat_path, pathtosave=results_path, model=model
        )
    print('formation efficiency done')
    
    if 'Z' in model:
        pp.get_interactionsep_and_numLISA(
            pathtocosmic=dat_path, pathtoLISA=LISA_path, pathtoresults=results_path, model=model, var=True,     
        )
        print('interaction sep done')
        
        pp.get_resolvedDWDs(
            pathtoLISA=LISA_path, pathtosave=results_path, var=True, model=model, window=1000
        )
        print('resolved FZ done')
    else:
        pp.get_interactionsep_and_numLISA(
            pathtocosmic=dat_path, pathtoLISA=LISA_path, pathtoresults=results_path, model=model, var=False,     
        )
        print('interaction sep done')
    
    
        pp.get_resolvedDWDs(
           pathtoLISA=LISA_path, pathtosave=results_path, var=False, model=model, window=1000
        )
    
        print('resolved done')
    
