import numpy as np

#===================================================================================
# File and Organizing Functions:
#===================================================================================

def getfiles(kstar1, kstar2, model):
    '''
    Get the list of COSMIC dat files used to create 
    DWD populations.
    
    INPUTS
    ----------------------
    kstar1 [10, 11, 12]: type of first WD
    kstar2 [10, 11, 12]: type of second WD
    model: model variation from 'fiducial', 'q3', 'alpha25', 'alpha5'
    
    RETURNS
    ----------------------
    filename_list [list, str]: list of files names
    label [str]: DWD type label
    '''
    met_list = [0.0001, 0.00015029, 0.00022588, 0.00033948, 0.00051021, 
                0.00076681, 0.00115245, 0.00173205, 0.00260314, 0.00391233, 
                0.00587993, 0.0088371, 0.01328149, 0.01996108, 0.03]
    
    filename_list = []
    for ii in range(len(met_list)):
        fname = '{}_{}_{}_{}.h5'.format(model, kstar1, kstar2, ii)
        filename_list.append(fname)
        
    if kstar1 == '12':
        label='12'
    else:
        label='{}_{}'.format(kstar1, kstar2)
        
    return filename_list, label


def Lband_files(kstar1, kstar2, model, var=True):
    '''
    Creates list of LISA band filenames for given DWD type.
    
    INPUTS
    ----------------------
    kstar1 [10, 11, 12]: type of first WD
    kstar2 [10, 11, 12]: type of second WD
    model ['fiducial', 'alpha25', 'alpha5', 'q3']: specifies binary evolution assumptions
    var [bool]: if True, creates keys for FZ, else for F50
    
    RETURNS
    ----------------------
    files [list, str]: list of keys
    '''
    met_list = [0.0001, 0.00015029, 0.00022588, 0.00033948, 0.00051021, 
                0.00076681, 0.00115245, 0.00173205, 0.00260314, 0.00391233, 
                0.00587993, 0.0088371, 0.01328149, 0.01996108, 0.03]
    if var:
        var_list = [0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                    0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                    0.3076, 0.2726, 0.2376, 0.2027, 0.1677]
      
        if kstar1 != '12':
            files = ['Lband_{}_{}_{}_{}_{}.hdf'.format(kstar1, kstar2, int(met*1e5), int(var*1e5), model) for var, met in zip(var_list, met_list)]
        else:
            files = ['Lband_{}_{}_{}_{}.hdf'.format(kstar1, int(met*1e5), int(var*1e5), model) for var, met in zip(var_list, met_list)]
    else:
        if kstar1 != '12':
            files = ['Lband_{}_{}_{}_0.5_{}.hdf'.format(kstar1, kstar2, int(met*1e5), model) for met in met_list]
        else:
            files = ['Lband_{}_{}_0.5_{}.hdf'.format(kstar1, int(met*1e5), model) for met in met_list]
    
    return keys


def get_binfrac_of_Z(Z):
    '''
    Calculates the theoretical binary fraction as a 
    function of metallicity.
    
    INPUTS
    ----------------------
    Z [array]: metallicity Z values
    
    RETURNS
    ----------------------
    binfrac [array]: binary fraction values
    '''
    FeH = get_FeH_from_Z(Z)
    FeH_low = FeH[np.where(FeH<=-1.0)]
    FeH_high = FeH[np.where(FeH>-1.0)]
    binfrac_low = -0.0648 * FeH_low + 0.3356
    binfrac_high = -0.1977 * FeH_high + 0.2025
    binfrac = np.append(binfrac_low, binfrac_high)
    return binfrac

def get_Z_from_FeH(FeH, Z_sun=0.02):
    '''
    Converts from FeH to Z under the assumption that
    all stars have the same abundance as the sun
    
    INPUTS
    ----------------------
    FeH [array]: array of Fe/H values to convert
    Z_sun [float]: solar metallicity
    
    RETURNS
    ----------------------
    Z [array]: array of metallicities
    '''
    Z = 10**(FeH + np.log10(Z_sun))
    return Z

def get_FeH_from_Z(Z, Z_sun=0.02):
    '''
    Converts from Z to FeH under the assumption that
    all stars have the same abundance as the sun
    
    INPUTS
    ----------------------
    Z [array]: array of metallicities to convert
    Z_sun [float]: solar metallicity
    
    RETURNS
    ----------------------
    Z [array]: array of FeH values
    '''
    FeH = np.log10(Z) - np.log10(Z_sun)
    return FeH