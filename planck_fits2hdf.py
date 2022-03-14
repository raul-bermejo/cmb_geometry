import time
from fits2hdf.io.fitsio import read_fits
from fits2hdf.io.hdfio import export_hdf

# Define each filepath
filenames = ['HFI_SkyMap_143_2048_R3.01_halfmission-1.fits',        # map for half mission 1 143 GHZ
             'HFI_SkyMap_143_2048_R3.01_halfmission-2.fits',        # map for half mission 2 143 GHZ
             'HFI_Mask_GalPlane-apo5_2048_R2.00.fits',              # galactic plane mask for input apodization
             'HFI_Mask_PointSrc_2048_R2.00.fits',                   # mask point source
             'Bl_T_R3.01_fullsky_143hm1x143hm1.fits',               # beam transfer function for hm1
             'Bl_T_R3.01_fullsky_143hm1x143hm2.fits',               # beam transfer function for hm2
            ]

start = time.time()
print('The conversion process has started:')
print('='*80)
# export each FITS file in the list to hdf
for file in filenames:
    r = read_fits('data/'+file)
    name = file.split('.fits')[0]
    export_hdf(r, 'data/'+name+'.hdf')
    print(f'{file} was succesfully converted.')


end = time.time()
print(f'All the Planck FITS files were succesfully converted to HDF format. The conversion of all files took {end-start:.2f}.')