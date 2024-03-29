import numpy as np
from astropy.io import fits
from healpy.fitsfunc import read_map
import time

APO = 5
NSIDE = 2048
FREQ = 143
GAL_COVER = 80
HDU = 1 

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


# define useful functions
def extract_data(filepath, hdu=1):
    '''
    Extracts input HDU from map Planck FITS file (we only need the first one)
    Document function once fully tested
    '''
    with fits.open(filepath) as hdul:
        hdul = fits.open(filepath)
        data = hdul[hdu].data                                              # hdul is a list of HDU objects
    
    return data

def fits2npy(filename, 
            field, hdu=1, gp=False, beamfunc=False):
    '''
    Extracts useful data from FITS file and saves it as .npy file
    Document function once fully tested
    '''
    
    print(f'Extracting data from {filename}:')
    print('-'*75)

    short_name = filename.split('.fits')[0]

    if beamfunc:
        array = extract_data('data/'+filename)
    else:
        array = read_map('data/'+filename, field=field, hdu=HDU)

    array = np.array(array.astype(np.float64))

    if gp:
        output_name = 'data_test/'+short_name+f'_fsky={GAL_COVER/100}.npy'

    else:
        output_name = 'data_test/'+short_name+'.npy'

    np.save(output_name, array)
    print(f"Data was extracted and saved into {output_name} succesfully")
    print('-'*70)

    return None

start = time.time()
# define filenames as downloaded from planck release: https://pla.esac.esa.int/#home
filepath_hm1 = 'HFI_SkyMap_143_2048_R3.01_halfmission-1.fits'         # map for half mission 1 143 GHZ
filepath_hm2 = 'HFI_SkyMap_143_2048_R3.01_halfmission-2.fits'         # map for half mission 2 143 GHZ
filepath_gp = f'HFI_Mask_GalPlane-apo{APO}_2048_R2.00.fits'           # galactic plane mask for input apodization length
filepath_ps = 'HFI_Mask_PointSrc_2048_R2.00.fits'                     # mask point source
filepath_beamwin1 = 'Bl_T_R3.01_fullsky_143hm1x143hm1.fits'           # beam transfer function for hm1
filepath_beamwin2 = 'Bl_T_R3.01_fullsky_143hm1x143hm2.fits'           # beam transfer function for hm2
         
start1 = time.time()
print(f'Started Planck data extraction:')
print('='*80)
# extract data by using read_map funcitn from healpy
skies2npy = False
if skies2npy:
    sky_hm1, sky_hm2 = fits2npy(filepath_hm1, field=0, hdu=HDU), fits2npy(filepath_hm2, field=0, hdu=HDU)

field_gp = GALCOVER_DICT[f'GAL0{GAL_COVER}']
gp2npy = False
if gp2npy:
    mask_gp = fits2npy(filepath_gp, field=field_gp, hdu=HDU, gp=True)

field_ps = FREQ_PS_DICT['F143']
ps2npy = False
if ps2npy:
    mask_ps = fits2npy(filepath_ps, field=field_ps, hdu=HDU)

# use the extract_data function to find the beam window function
beamfunc2npy = False
if beamfunc2npy:
    beam_hm1, beam_hm2 = fits2npy(filepath_beamwin1, field=None, beamfunc=True), fits2npy(filepath_beamwin2, field=None, beamfunc=True)

end = time.time()
if (skies2npy or gp2npy or ps2npy or beamfunc2npy):
    print(f'='*80)
    print(f'File(s) were loaded and saved succesfully in {end-start:.2f} s.')
else:
    print('No files were selected for conversion from .fits to .npy')