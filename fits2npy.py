import numpy as np
from astropy.io import fits
import time

APO = 5
NSIDE = 2048
FREQ = 143
COVER = 60

# define useful functions
def extract_data(filename, file_kind,
                 hdu=1):
    '''
    Extracts useful data from FITS file and saves it as .npy file
    Document function once fully tested
    '''
    
    print(f'Extracting data from {filename}:')
    print(f'='*80)
    with fits.open('data/'+filename) as hdul:
        data_table = hdul[hdu].data                    # hdul is a list of HDU objects
        short_name = filename.split('.fits')[0]

        if file_kind == 'map':
            array = data_table['I_STOKES']

        elif file_kind == 'mask_point_source':
            array = data_table[f'F{FREQ}']

        elif file_kind == 'mask_galactic_plane':
            array = data_table[f'GAL0{COVER}']
            short_name += f'_GAL0{COVER}_'
        
        else:
            raise SyntaxError("Invalid file_kind, must be either: ['map', 'mask_point_source', 'mask_galactic_plane' ]")
        
        np.save('data/'+short_name+'.npy', array)

        print(f"Data was extracted and saved into {'data/'+short_name+'.npy'} succesfully")
        print('-'*90)
    
    return None


# define filenames as downloaded from planck release: https://pla.esac.esa.int/#home
filenames = ['HFI_SkyMap_143_2048_R3.01_halfmission-1.fits',        # map for half mission 1 143 GHZ
             'HFI_SkyMap_143_2048_R3.01_halfmission-2.fits',        # map for half mission 2 143 GHZ
            f'HFI_Mask_GalPlane-apo{APO}_2048_R2.00.fits',          # galactic plane mask for input apodization length
             'HFI_Mask_PointSrc_2048_R2.00.fits',                   # mask point source
             'Bl_T_R3.01_fullsky_143hm1x143hm1.fits',               # beam transfer function for hm1
             'Bl_T_R3.01_fullsky_143hm1x143hm2.fits',               # beam transfer function for hm2
            ]

start = time.time()

# Call function for each of the files
extract_data(filenames[0], file_kind='map')
extract_data(filenames[1], file_kind='map')

extract_data(filenames[2], file_kind='mask_galactic_plane')

extract_data(filenames[3], file_kind='mask_point_source')