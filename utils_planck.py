from datetime import datetime
import time
import os
import pandas as pd
import numpy as np
import healpy as hp
import seaborn as sns
from astropy.io import fits
import matplotlib.pyplot as plt
from healpy.fitsfunc import read_map
import matplotlib.style as style
import camb
from camb import model, initialpower


style.use('tableau-colorblind10')
colors = sns.color_palette("colorblind").as_hex()
#style.use('seaborn-colorblind')
today = datetime.today().strftime('%Y-%b-%d')


# define dictionaries to map gal cover and frequency to indices
GALCOVER_DICT = {
    'GAL020': 0,
    'GAL040': 1,
    'GAL060': 2,
    'GAL070': 3,
    'GAL080': 4,
    'GAL090': 5,
    'GAL097': 6,
    'GAL099': 7
}

FREQ_PS_DICT = {
    'F100': 0,
    'F143': 1,
    'F217': 2,
    'F353': 3,
    'F545': 4,
    'F857': 5,
} 

# define utils functions
def extract_data(filepath, hdu=1):
    '''
    Extracts input HDU from map Planck FITS file (we only need the first one)
    '''
    with fits.open(filepath) as hdul:
        hdul = fits.open(filepath)
        data = hdul[hdu].data                                              # hdul is a list of HDU objects
    
    return data


def clean_skydata(sky_array, 
                  bad_data_threshold=1e12):
    '''
    Checks whether the average of the map is close to 0. Ff not, removes bad data points.
    '''
        
    print(f'Carrying out data quality of sky map:')
    print('-'*80)
    
    init_mean = np.mean(sky_array)
    print(f'Before data clean, the mean of the sky map corresponds to:  {init_mean:.2e}')
    if np.abs(init_mean) > 1:
        print('Skymap contains bad data')
        
        sky_array[np.abs(sky_array) >  bad_data_threshold] = 0
    
        print(f'After data clean, the mean of the sky map corresponds to: {np.mean(sky_array):.2e}')
    
    else:
        print(f'Data is already clean.')
        print(f'Currently, the mean is {np.mean(sky_array):.2e}')
    
    mean_sky = np.mean(sky_array)
    print(f'Data Quality check completed.') 
    print('='*80)
    
    return(sky_array)


def find_spherical_harmonics(sky_i, gp_mask, ps_mask, 
                             lmax=4000, apply_mask=True):
    '''
    Computes the spherical harmonics from a map,
    masking the stokes parameters from the Planck maps.
    '''
    
    start = time.time()
    # rescale the sky map to mKelvin
    sky_i *= 10**6
    
    if apply_mask:
        x_masked = sky_i*gp_mask*ps_mask
    else:
        x_masked = sky_i
    
    # compute spherical harmonics
    a_lm =  hp.sphtfunc.map2alm(x_masked, pol=False, lmax=lmax)             
    
    end = time.time()
    print(f'It took {end-start:.2f} s. to compute a_lm')
    
    return a_lm


def find_power_spectrum(alm_1, alm_2, bl_1=None, bl_2=None,
                        M_ll=1, f_l=1, n_l=0, NSIDE=2048,
                        which_return='both'):
    '''
    Finds the cross-power spectrum (via pseudo cross power spectrum) given the 
    coefficients of the spherical harmonics.
    '''
    
    start = time.time()
    # compute pseudo power spectrum
    print(len(alm_1))
    D_l = hp.sphtfunc.alm2cl(alm_1, alm_2)
    
    if which_return=='pseudo':
        end = time.time()
        print(f'It took {end-start:.2f} s. to compute D_l')
        
        return D_l
    elif which_return=='actual' or which_return=='both':
        # now define all the instrument related biases
        bl_1, bl_2 = bl_1.astype(np.float), bl_2.astype(np.float)
        p_l = hp.sphtfunc.pixwin(NSIDE, lmax=len(bl_1)-1)
                
        # calculate actual power spectrum
        C_l = D_l/(p_l**2 * bl_1*bl_2)
        
        end = time.time()
        print(f'It took {end-start:.2f} s. to compute C_l')
        
        if which_return=='both':
            return (D_l, C_l)
        else:
            return C_l


def clean_power_spectrum(C_l, ind_low=20, ind_max=2000,
                         decimate=True, n_dec=2):
    '''
    Bins the computed power spectrumm removes the edges and decimates it if set to True
    '''
    # define multipole array
    C_l = C_l[ind_low:ind_max]
    l = np.arange(len(C_l))
        
    Cl_planck = l*(l+1)*C_l/(2*np.pi)

    # decimate each array
    if decimate:
        n_dec = 2
        Cl_planck, l_planck = Cl_planck[::n_dec], l[::n_dec]
    
    return (l_planck, Cl_planck)

def peak_finder(x,y):
    '''
    Finds the value of x for which y(x) is max (peaks)
    '''
    ind = np.where(y == np.max(y))
    x_peak = x[ind][0]
    
    return x_peak

def model_cl(l, a_0, a_2, l_p):
    '''
    Second order model with offset and no linear term
    '''
    return a_2*(l-l_p)**2 + a_0

def boltzmann_cl(omega_bary, omega_cdm, omega_k, H_0=67.5):
    '''
    Finds the cross power spectrum for a Boltzmann code given
    input density parameters
    '''
    # Initialize parameters of CAMB model
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H_0, ombh2=omega_bary, omch2=omega_cdm, omk=omega_k,  mnu=0.06,  tau=0.06)
    #pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
    pars.set_for_lmax(2500, lens_potential_accuracy=0)

    #calculate results for these parameters
    results = camb.get_results(pars)
    Cl_allpol = results.get_cmb_power_spectra(pars, CMB_unit='muK', spectra={'total'})['total']     # the C_l returned is already binned as l*(l+1)+C_l/(2pi)

    Cl_boltz = Cl_allpol[:,0]                                                          # we only need the TT/II polarization (i.e. the temperature)
    l_boltz = np.arange(Cl_boltz.shape[0])
    
    return l_boltz, Cl_boltz