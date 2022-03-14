import numpy as np
from astropy.io import fits
from healpy.fitsfunc import read_map
import time

APO = 5
NSIDE = 2048
FREQ = 143
COVER = 60

# define useful functions
def masks2npy(filename, file_kind,
              hdu=1):
    '''
    Extracts useful data from FITS file and saves it as .npy file
    Document function once fully tested
    '''
    
    print(f'Extracting data from {filename}:')
    print(f'='*130)
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
        
        elif file_kind == 'beam_mask':
            array = data_table
        
        else:
            raise SyntaxError("Invalid file_kind, must be either: ['map', 'mask_point_source', 'mask_galactic_plane', beam ]")
        
        array = np.array(array.astype(np.float))
        np.save('data/'+short_name+'.npy', array)

        print(f"Data was extracted and saved into {'data/'+short_name+'.npy'} succesfully")
        print('-'*100)
    
    return None

def maps2npy(filename, hdu = 1, field=0):

    print(f'Extracting data from {filename}:')
    print(f'='*130)

    array = read_map('data/'+filename, field=field, hdu=hdu, h=False)
    short_name = filename.split('.fits')[0]

    np.save('data/'+short_name+'.npy', array)

    print(f"Data was extracted and saved into {'data/'+short_name+'.npy'} succesfully")
    print('-'*100)

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
#maps2npy(filenames[0])
#maps2npy(filenames[1])

masks2npy(filenames[2], file_kind='mask_galactic_plane')
masks2npy(filenames[3], file_kind='mask_point_source')

masks2npy(filenames[4], file_kind='beam_mask')
masks2npy(filenames[5], file_kind='beam_mask')

end = time.time()
print(f'='*130)
print(f'All files were loaded and saved succesfully in {end-start:.2f}')