#! /usr/bin/env python
import os
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


def region_wcs(filename, ra, dec, a, b, theta, objid, color='green'):
    """Create a ds9 region file with elliptical apertures in WCS coords. """
    f = open(filename, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('fk5\n')
    print ra.shape
    for i in range(ra.shape[0]):
        f.write('ellipse(%f,%f,%f,%f,%f) # color=%s text={%i}\n' %
                (ra[i], dec[i], a[i]*3.5, b[i]*3.5, theta[i]*(-1.), 
                 color, objid[i]))
    f.close()


def plot_scatter(ax, xx, yy):
    """ """
    ax.scatter(xx, yy, marker='o', edgecolor='none', alpha=0.5, color='b')
    ax.plot([14,29], [1,1], 'k', lw=2)
    faintlim = np.max(xx)+0.5
    ax.set_xlim(17, faintlim)


def check_conv_phot(origcat, newcat, filt):
    """ """
    # read in original cat
    o = np.genfromtxt(origcat, dtype=[('num',float), ('x_world',float),
                                      ('y_world',float), ('a_world',float),
                                      ('b_world',float), ('theta_world',float),
                                      ('mag',float), ('emag',float)], 
                               usecols=(1,7,8,9,10,11,12,13))

    # read in new convolved cat
    d = fits.getdata(newcat)
    mag = d['mag_auto_%s' % filt]
    emag = d['magerr_auto_%s' % filt]
    a = d['a_world']
    b = d['b_world']
    t = d['theta_world']
    
    # match old and new catalog
    idx,separc = match_cats(d['x_world'],d['y_world'],o['x_world'],o['y_world'])
    match = (separc.value*3600. <= 0.1)
    # matches in original catalog are [idx[match]]
    # matches in new catalogs are [match]

    fig = plt.figure(figsize=(8.5,11))
    gs = GridSpec(5,6)
    ax1 = fig.add_subplot(gs[:2,:-3]) 
    ax2 = fig.add_subplot(gs[:2,-3:]) 
    ax3 = fig.add_subplot(gs[2,:2]) 
    ax4 = fig.add_subplot(gs[2,2:4]) 
    ax5 = fig.add_subplot(gs[2,4:]) 
    ax6 = fig.add_subplot(gs[3:,1:-1])
    
    ax1.set_title(filt)

    # plot ratios
    if filt == 'F110W':
        w = np.where((o['emag'][idx[match]] != 0) & 
                     (o['num'][idx[match]] < 2000))
    elif filt == 'F160W':
        w = np.where((o['emag'][idx[match]] != 0) & 
                     ((o['num'][idx[match]] >= 2000) |
                      (o['num'][idx[match]] <= 1000)))
    plot_scatter(ax1, mag[match][w], mag[match][w]/o['mag'][idx[match]][w])    
    plot_scatter(ax2, mag[match][w], emag[match][w]/o['emag'][idx[match]][w])
    plot_scatter(ax3, mag[match][w], a[match][w]/o['a_world'][idx[match]][w])
    plot_scatter(ax4, mag[match][w], b[match][w]/o['b_world'][idx[match]][w])
    plot_scatter(ax5, mag[match][w], t[match][w]/o['theta_world'][idx[match]][w])
    st1 = np.std(mag[match][w]/o['mag'][idx[match]][w])
    st2 = np.std(emag[match][w]/o['emag'][idx[match]][w])
    st5 = np.std(t[match][w]/o['theta_world'][idx[match]][w])
    print st1   

    # compare the errors on the convolved images:
    #   SE errors with those from calc_error
    c = fits.getdata(os.path.join(os.path.dirname(newcat), '%s_cat.fits'%filt))
    id2,sep2 = match_cats(d['x_world'],d['y_world'],c['x_world'],c['y_world'])
    mat2 = (sep2.value*3600. <= 0.1)
    calc_eiso = d['magerr_iso_%s' % filt][mat2]
    calc_eauto = d['magerr_auto_%s' % filt][mat2]
    calc_eaper = d['magerr_aper_%s' % filt][mat2]
    calc_eaper = calc_eaper[:,0]
    SE_eiso = c['magerr_iso'][id2[mat2]]
    SE_eauto = c['magerr_auto'][id2[mat2]]
    SE_eaper = c['magerr_aper'][id2[mat2]]
    SE_eaper = SE_eaper[:,0]
    # find all sources with good errors in SE catalog
    wi = np.where(SE_eiso < 1000.)
    wa = np.where(SE_eauto < 1000.)
    w2 = np.where(SE_eaper < 1000.)

    ax6.plot([0,10], [0,10], 'k', lw=2)
    ax6.scatter(SE_eiso[wi], calc_eiso[wi], color='r', edgecolor='none',
                alpha=0.7, label='ISO')
    ax6.scatter(SE_eauto[wa], calc_eauto[wa], color='b', edgecolor='none',
                alpha=0.7, label='AUTO')
    ax6.scatter(SE_eaper[w2], calc_eaper[w2], color='g', edgecolor='none',
                alpha=0.7, label='APER')
    ax6.legend(scatterpoints=1, loc=2)
    # standard deviation for plot limits
    st6 = np.std(SE_eauto[wa])

    for ax in [ax1,ax2,ax3,ax4,ax5]:
        ax.set_xlabel(r'new $m_{\mathrm{AUTO}}$')
        ax.set_xlim(17,29)
    ax1.set_ylabel(r'$m_{\mathrm{AUTO}}$ new / old')
    ax2.set_ylabel(r'$\sigma_{m,\mathrm{AUTO}}$ new / old')
    ax3.set_ylabel(r'$a$ [wcs], new / old')
    ax4.set_ylabel(r'$b$ [wcs], new / old')
    ax5.set_ylabel(r'$\theta$, new / old')
    ax6.set_xlabel(r'$\sigma_{mag}$ old')
    ax6.set_ylabel(r'$\sigma_{mag}$ new')
    
    # plot limits
    y2 = np.min([1+st2*2, 10])
    ax2.set_ylim(0, y2)
    ax5.set_ylim(1-st5*3, 1+st5*3)
    if filt == 'F110W':
        ax1.set_ylim(1-st1*10, 1+st1*10)
        ax6.set_xlim(0, 1.5)
        ax6.set_ylim(0, 1.5)
    elif filt == 'F160W':
        ax1.set_ylim(1-0.5*st1, 1+0.5*st1)
        ax6.set_xlim(0, 6)
        ax6.set_ylim(0, 6)

    plt.tight_layout()
    plt.savefig(os.path.join(os.path.dirname(newcat), 'check_conv_phot_%s.pdf'%filt))
    

def aperture(xc, yc, x, y, cat, objid, radius=0.0, ap='circle'):
    """Identify pixels inside a circular or elliptical aperture"""
    # distance of every pixel from star's position
    if ap == 'circle':
        r2 = (x - xc)**2 + (y - yc)**2
        # find all pixels inside a circular aperture
        #mask = r2 <= (radius + 0.5)**2
        #mask = np.where(r2 <= (radius-0.5)**2)
        mask = np.where(r2 <= (radius)**2)

    elif ap == 'ellipse':
        cxx = cat['CXX_IMAGE'][objid]
        cxy = cat['CXY_IMAGE'][objid]
        cyy = cat['CYY_IMAGE'][objid]
        r2 = cxx*(x - xc)**2 + cxy*(x - xc)*(y - yc) + cyy*(y - yc)**2
        # find all pixels inside an elliptical aperture
        # pad apertuers a bit to match SE errors
        mask = np.where(r2 <= (3.5)**2)
        #mask = np.where(r2 <= (3.5+1.5)**2)

    return mask


def calc_ap_error(xc, yc, objidx, cat, x, y, img, rms, seg, phot, skysz, objid=None, ap=None):
    """
       SE flux uncertainties:
         sqrt(A * sig_i^2 + F/g)
         A = area of aperture
         sig = average sky value
         F = flux
         g = gain
       
    """
    # get sky errors
    # using semi-global approach so should be same for all ap types
    objpix = np.where(seg == objid)
    imgbox = img[yc-skysz:yc+skysz, xc-skysz:xc+skysz]
    segbox = seg[yc-skysz:yc+skysz, xc-skysz:xc+skysz]
    rmsbox = rms[yc-skysz:yc+skysz, xc-skysz:xc+skysz]
    # consider pixels that are 
    #   1. not associated with an object
    #   2. on-image
    #   3. not bad pixels
    goodsky = np.where((segbox == 0) & (imgbox != 0) & (rmsbox != 1.e10))
    nsky = imgbox[goodsky].shape[0]
    if nsky > 0:
        sky = imgbox[goodsky]
        # get sigma of sky 
        skysig = np.std(sky)
        skyerr = skysig / nsky
    else:
        skyerr = 0

    # pixels corresponding to aperture
    if phot == 'ISO':
        appix = np.where(seg == objid)
    elif phot == 'AUTO':
        appix = aperture(xc, yc, x, y, cat, objidx, ap='ellipse')
    elif phot == 'APER':
        appix = aperture(xc, yc, x, y, cat, objidx, radius=ap)

    # rms added in quadrature
    # this assumes they are not correlated, but rms map has already
    # been corrected for correlated pixels
    aprms = rms[appix]
    # remove bad pixels from error consideration
    goodrms = aprms[aprms != 1.e10]
    nrms = goodrms.shape[0]
    # this is the dominant source of error
    rmserr = np.sqrt(np.sum(goodrms**2))

    return rmserr, nrms*skyerr


def calc_errors(cat, imgfile, rmsfile, segfile, expt, phot='AUTO', skysz=256):
    """ """
    objids = cat['NUMBER']
    xx = cat['X_IMAGE']
    yy = cat['Y_IMAGE']
    flux = cat['FLUX_%s'%phot]
    eflux = cat['FLUXERR_%s'%phot]
    magerr = cat['MAGERR_%s'%phot]
    img = fits.getdata(imgfile)# * expt
    rms = fits.getdata(rmsfile)# * expt
    seg = fits.getdata(segfile)

    # get x and y vectors of image
    y,x = np.ogrid[:img.shape[0], :img.shape[1]]

    # circular aperture sizes
    aps = np.array([10, 15, 20, 25, 30]) / 2.

    err = eflux.copy()

    if phot != 'APER':
        for i in range(eflux.shape[0]):
            # only need to recalculate for objects that are affected
#            if eflux[i] < 1000.:
            objid = objids[i]
            ef = calc_ap_error(xx[i], yy[i], i, cat, x, y, img, rms, seg, 
                               phot, skysz, objid=objid)
            # ef[0] = sumerr
            # ef[1] = npix*skyerr
            err[i] = np.sqrt(ef[0]**2 + ef[1]**2 + flux[i]/expt)
    
    if phot == 'APER':
        for j in range(eflux.shape[1]):
            for i in range(eflux.shape[0]):
                objid = objids[i]
#                if eflux[i,j] < 1000.:
                ef = calc_ap_error(xx[i], yy[i], i, cat, x, y, img, rms,
                                   seg, phot, skysz, ap=aps[j], objid=objid)
                err[i,j] = np.sqrt(ef[0]**2 + ef[1]**2 + flux[i,j]/expt)

#    emag = err / (flux*expt) * 2.5 / np.log(10)
    emag = err / flux * 2.5 / np.log(10)

    return err, emag

