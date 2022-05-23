# =========================================================================
# This script creates LISA DWD galaxies across 15 metallicity bins,
# incorporating the metallicity-dependent binary fraction as
# discussed in Thiele et al. (2021).
#
# Authors: Sarah Thiele & Katelyn Breivik
# Last updated: Oct 14th, 2021
# =========================================================================

import numpy as np
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import argparse
import postproc as pp
import paths


DWD_list = ["He_He", "CO_He", "CO_CO", "ONe_X"]
dat_path = str(paths.data) + "/"
FIRE_path = str(paths.data) + "/"
models = ["fiducial", "alpha25", "alpha5", "q3"]
interfile = False
nproc = 4

for model in models:
    pp.save_full_galaxy(
        DWD_list, dat_path, FIRE_path, dat_path, interfile, model, nproc
    )
    print("Gx done!")

    if model == "fiducial":
        pp.get_formeff(pathtodat=dat_path, pathtosave=dat_path)
    print("formation efficiency done")

    pp.get_interactionsep_and_numLISA(
        pathtodat=dat_path, pathtosave=dat_path, model=model, var=True
    )
    print("interaction sep FZ done")

    pp.get_interactionsep_and_numLISA(
        pathtodat=dat_path, pathtosave=dat_path, model=model, var=False
    )
    print("interaction sep F50 done")

    pp.get_resolvedDWDs(dat_path, dat_path, var=True, model=model, window=1000)
    print("resolved FZ done")
    pp.get_resolvedDWDs(dat_path, dat_path, var=False, model=model, window=1000)

    print("resolved F50 done")
