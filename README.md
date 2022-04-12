# cmb_geometry

Anisotropy analysis of Cosmic Microwave Background (CMB) from Planck's Release 3 to estimate the density curvature parameter.
The statistical analysis consists of computing the _cross power spectrum_ from two Planck sky maps and associated masks. 

## Data 

The Planck satellite was launched by the European Space Agency (ESA) in May 2019 with the goal to measure the temperature fluctuations of the CMB across the sky. The Planck temperature sky maps used for this study are publicly available at the [ESA Planck mission website](http://pla.esac.esa.int/pla/#maps). The following files were downloaded:

> 'HFI_SkyMap_143_2048_R3.01_halfmission-1.fits' 
> 
> 'HFI_SkyMap_143_2048_R3.01_halfmission-2.fits'
> 
> 'HFI_Mask_GalPlane-apo5_2048_R2.00.fits'
> 
> 'HFI_Mask_PointSrc_2048_R2.00.fits'
> 
> 'Bl_T_R3.01_fullsky_143hm1x143hm1.fits'
> 
> 'Bl_T_R3.01_fullsky_143hm1x143hm2.fits' 
> 

As seen above, two half-mission single-frequency maps of 143 GHz were used and a few masks were applied to reduce signal-to-noise ratio.


## Usage

After downloading the files above, move into `data` directory. Then, run the `fits2npy.py` to convert the data files from `.fits` to `.npy` for convenience. This will create .npy files with the same name and save into the `data` directory.
Finally, run `planck_analysis.ipynb` to produce the anisotropy results. Plots will be produced in the `plots` directory.

## Contributions
 
 Pull requests are welcome. Please ensure that the pipeline `planck_analysis.ipynb` without errors.
 
## Acknowledgements

This project was funded by the [Laby Research Scholars Program](https://physics.unimelb.edu.au/study/Scholarships/laby-research-scholars-program).
