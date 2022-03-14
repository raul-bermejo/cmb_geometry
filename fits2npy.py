import numpy as np
from astropy.io import fits

# define useful functions
def extract_data(filepath, file_type,
                 hdu=1):
    '''
    Extracts input HDU from map Planck FITS file (we only need the first one)
    Document function once fully tested
    '''
    
    with fits.open(filepath) as hdul:
        data = hdul[hdu].data                    # hdul is a list of HDU objects

        hdul.info()                              # prints metadata about FITS file (HDUs)
        print('='*90)
    
    return data


# define filenames as downloaded from planck release: https://pla.esac.esa.int/#home
filenames = ['HFI_SkyMap_143_2048_R3.01_halfmission-1.fits',        # map for half mission 1 143 GHZ
             'HFI_SkyMap_143_2048_R3.01_halfmission-2.fits',        # map for half mission 2 143 GHZ
            f'HFI_Mask_GalPlane-apo{APO}_2048_R2.00.fits',          # galactic plane mask for input apodization length
             'HFI_Mask_PointSrc_2048_R2.00.fits',                   # mask point source
             'Bl_T_R3.01_fullsky_143hm1x143hm1.fits',               # beam transfer function for hm1
             'Bl_T_R3.01_fullsky_143hm1x143hm2.fits',               # beam transfer function for hm2
            ]
