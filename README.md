# README for UVIS_IR_catalogs.py

# USAGE
./UVIS_IR_Catalogs.py [-h] par

par is the ParID and can include any text (Par123, par123, field123,
123, etc) as it only reads the numbers


export PATH="${PATH}:$HOME/[pathtoscripts]/wisp-iruvis-cats/"
export PYTHONPATH=$PYTHONPATH:$HOME/[pathtoscripts]/wisp-iruvis-cats/

# 1. need to convolve images where bad pixels are removed
# 2. need to run SE in dual image mode for all filters
#       - F110W: for detection, we remove all problem pixels from RMS maps
#           (replace all pixels > 0.01 or something)
#       - for analysis, leave bad pixels in
#   do not convolve rms maps

# REQUIREMENTS
numpy
scipy
astropy
pyraf

These are all in requirements.txt so if necessary you can do:
    pip install -r requirements.txt


# NOTES:
- It assumes the IR and UVIS aligned data are in ParXXX/DATA/UVIS/IRtoUVIS. 
  If your data are somewhere else, you can change datadir in line 20.

- The code uses ConfigParser for SExtractor parameters and the thresholds
  for image convolution. This allows you to set parameters simulataneously 
  for all filters (in parameters.cfg). Parameters in this file will override 
  those in the config.sex

- It will convolve each image to the psf of the reddest filter, so the psfs 
  must be present for the calculation of the convolution kernel.

- I've only carefully checked the convolution threshold for the following 
  filters:
     F110W --> F160W
     F475X --> F160W
     F600LP --> F160W
     F606W --> F160W
  You may want to check these yourself. The convolution is not great around 
  bright point sources, where there is a lot of fringing. Marc and I have 
  discussed this at length and we think this is because I'm using pure TinyTim 
  psfs. I do not smooth the psfs because I wasn't sure how this would affect 
  the measurements.

- Because the IR has been drizzled to the UVIS pixscale, sextracted objects can 
  sometimes look strangely pixelated in the segmentation map. I've found that 
  using a Gaussian filter helps. The default in the code is gauss_5.0_9x9.conv,
  but this can be changed in parameters.cfg.

# STEPS:
1. Finds _sci.fits images for all filters available in the Par. 
   Uses IR and UVIS images drizzled onto the UVIS pixscale and aligned in
   pixel space for use with SExtractor's dual-image mode

2. Cleans all images.
   Replaces edges and chip gaps with random noise generated from the 
   distribution of noise in the image
   (cleanedges_general.py)

3. Replaces 0's in the _rms.fits images with NaNs. 
   Bad pixels are incorrectly identified with value=0 in the rms images
   output by the pipeline. SExtractor expects bad pixels in an RMS map 
   to have very large rms. Source detection is better if these bad pixels
   are replaced with NaNs.

4. Convolves images to the psf of the reddest filter.
   IRAF's psfmatch first computes the psf matching function. (This is 
   where the convolution threshold is used.) Astropy's convolve_fft then
   convolves the higher resolution image to the lower res psf. 
   Outputs: 
        ker_XXXtoYYY.fits - the computed kernel (for convolving FXXXW to FYYYW)
        XXX_convtoYYY.fits - the convolved image
   (convolvebypsf.py)

5. Runs SExtractor on the detection image.
   The reddest available filter is used as the detection image.
   
6. Runs SE on all remaining images in dual-image mode with the detection image.

# SEXTRACTOR:
- GAIN is set to the EXPTIME of each image
- MAG_ZEROPOINT is found based on the observation date
- A segmentation map is created for the detection image
- Catalogs are output as fits files and include ISO, AUTO, and APER 
  (radius = 0.2, 0.3, 0.4, 0.5, 0.6 arcsec) photometry. 
  Change aperture diamteres in parameters.cfg
- backphoto_think is currently set to a square 1.5" on a side
- Filter for detection is currently set to gauss_5.0_9x9.conv


