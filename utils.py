#! /usr/bin/env python
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


def match_cats(RA, Dec, refRA, refDec):
    """Match a catalog of RA's and Dec's to a reference catalog.

       Return the indices of the reference catalog that mach each
       source in the input catalog, and the on-sky separation
       between each source's closest match in arcsec.
    """
    # create SkyCoord objects to use with the matching
    SCcat = SkyCoord(ra=RA, dec=Dec, frame='icrs', unit=(u.deg,u.deg))
    SCrefcat = SkyCoord(ra=refRA, dec=refDec, frame='icrs', unit=(u.deg,u.deg))

    # idx    - indices of matched sources in reference cat
    # sep2d  - on-sky angular separation between closest match
    # dist3d - 3D distance between closest matches
    idx,sep2d,dist3d = match_coordinates_sky(SCcat, SCrefcat)

    return (idx, sep2d)


def region_wcs(filename, cat, color='green'):
    """Create a ds9 region file with elliptical apertures in WCS coords. """
    f = open(filename, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('fk5\n')
    print cat['X_WORLD'].shape
    for i in range(cat['X_WORLD'].shape[0]):
        f.write('ellipse(%f,%f,%f,%f,%f) # color=%s text={%i}\n' %
                (cat['X_WORLD'][i], cat['Y_WORLD'][i], cat['A_WORLD'][i]*2.,
                 cat['B_WORLD'][i]*2., cat['THETA_WORLD'][i]*(-1.), color,
                 cat['NUMBER'][i]))
    f.close()


def plot_scatter(ax, xx, yy):
    """ """
    ax.scatter(xx, yy, marker='o', edgecolor='none', alpha=0.5, color='b')
    ax.plot([14,26], [1,1], 'k', lw=2)
    faintlim = np.max(xx)+0.5
    ax.set_xlim(17, faintlim)


def check_conv_phot(origcat, newcat):
    """ """
    # read in original cat
    o = np.genfromtxt(origcat, dtype=[('num',float), ('x_world',float),
                                      ('y_world',float), ('a_world',float),
                                      ('b_world',float), ('theta_world',float),
                                      ('mag',float), ('emag',float)], 
                               usecols=(1,7,8,9,10,11,12,13))

    # read in new convolved cat
    f = fits.open(newcat)
    d = f[1].data
    f.close()
    mag = d['mag_auto']
    emag = d['magerr_auto']
    a = d['a_world']
    b = d['b_world']
    theta = d['theta_world']
    
    # match catalogs
    idx,separc = match_cats(d['x_world'],d['y_world'],o['x_world'],o['y_world'])
    match = (separc.value*3600. <= 0.1)
    # matches in original catalog are [idx[match]]
    # matches in new catalogs are [match]

    fig = plt.figure()
    gs = GridSpec(2,6)
    ax1 = fig.add_subplot(gs[0,:-3]) 
    ax2 = fig.add_subplot(gs[0,-3:]) 
    ax3 = fig.add_subplot(gs[1,:2]) 
    ax4 = fig.add_subplot(gs[1,2:4]) 
    ax5 = fig.add_subplot(gs[1,4:]) 

    plot_scatter(ax1, mag[match], o['mag'][idx[match]]/mag[match])    
    plot_scatter(ax2, mag[match], o['emag'][idx[match]]/emag[match])    
    plot_scatter(ax3, mag[match], o['a_world'][idx[match]]/a[match])    
    plot_scatter(ax4, mag[match], o['b_world'][idx[match]]/b[match])    
    plot_scatter(ax5, mag[match], o['theta_world'][idx[match]]/theta[match])    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(os.path.dirname(newcat), 'check_conv_phot.pdf'))
                    

def calc_errors(zeropoint, exptime, bckthick):
    """ """
    # convert zeropoint, image, and rms map from e-/s to e-
    zpt = zeropoint + 2.5*np.log10(exptime)
    im = im * exptime
    rms = rms * exptime
    


# PROBLEM WITH ERRORS IN IR IMAGES  - CALCULATE MYSELF?
# FIX SE parameters for getting all WISP objects and correct DEBLENDING
# FINISH CODE TO CHECK PHOTOMETRY AND PUT AS AN OPTION

