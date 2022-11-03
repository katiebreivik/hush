#=========================================================================
# This script creates LISA DWD galaxies across 15 metallicity bins,
# incorporating the metallicity-dependent binary fraction as
# discussed in Thiele et al. (2021).
#
# Authors: Sarah Thiele & Katelyn Breivik
# Last updated: Oct 14th, 2021
#=========================================================================

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import argparse
import postproc as pp
from schwimmbad import MultiPool
import h5py
import os
import paths


DWD_list = ['He_He', 'CO_He', 'CO_CO', 'ONe_X']
dat_path = paths.data / "cosmic_dat"
LISA_path = paths.data / "postproc"
FIRE_path = paths.data
results_path = paths.data
models = ['fiducial', 'fiducial_Z', 'alpha25', 'alpha25_Z', 'alpha5', 'alpha5_Z', 'q3', 'q3_Z']
interfile = False
nproc_gx = 40
nproc_results = len(models)
  

# Get the FIRE metallicity data for the results file
pp.get_FIRE_met_dat(FIRE_path, results_path)

print('FIRE done!')
# Generate Milky Way galaxies for each model
for model in models:
    pp.save_full_galaxy(
        DWD_list, dat_path, FIRE_path, LISA_path, interfile, model, nproc_gx
    )
    print('Gx done!')

# Get the results for each model
def get_results(dat):
    
    dat_path, LISA_path, results_path, var, model, window = dat
    pp.get_formeff(
        pathtodat=dat_path, pathtosave=results_path, model=model, var=var,
    )
    print('formation efficiency done for model: {}'.format(model))
    
    pp.get_interactionsep_and_numLISA(
        pathtocosmic=dat_path, pathtoLISA=LISA_path, pathtoresults=results_path, model=model, var=var,     
    )
    print('interaction sep done')
    
    pp.get_resolvedDWDs(
        pathtoLISA=LISA_path, pathtosave=results_path, var=var, model=model, window=window
    )
    print('resolved sources done')
    
    return []

dat = []
for model in models:
    if 'Z' in model:
        dat.append([dat_path, LISA_path, results_path, True, model, 1000])    
    else:
        dat.append([dat_path, LISA_path, results_path, False, model, 1000])
    
with MultiPool(processes=nproc_results) as pool:
    _ = list(pool.map(get_results, dat))
    
#Save the results into a single file            
for model in models:
    if 'Z' in model:
        var_label = 'FZ'
    else:
        var_label = 'F50'
    result_file = 'results_{}_{}.hdf'.format(var_label, model)
    file = h5py.File(results_path / result_file, "r")
    for key in file.keys():
        df = pd.read_hdf(results_path / result_file, key)
        df.to_hdf(results_path / 'results.hdf', key=key)
        #os.remove(results_path / result_file)
