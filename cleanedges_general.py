#! /usr/bin/env python
import argparse
import os
import numpy as np
import astropy.io.fits as fits
from astropy.convolution import Gaussian2DKernel,convolve
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def clean_image(image, output, check=False, showGauss=False):
    '''Cleans edges of an image. Also works for images with chip gaps.
       Replaces edges and chip gaps with random noise generated 
       from the distribution of noise in the image.

       Based on cleanedges_general.pro written by Marc Rafelski
    '''
    drz,hdr = fits.getdata(image, header=True)

    # replace all data with 1
    cln = np.copy(drz)
    cln[drz != 0] = 1

    # create Gaussian kernel
    gauss = Gaussian2DKernel(stddev=3.8)#, x_size=30, y_size=30)
    filt = convolve(cln, gauss, boundary=None, normalize_kernel=True)
    if check:
        fits.writeto('test.fits', filt, header=hdr, clobber=True)

    # bin data
    dat = drz[drz > 0]
    hist,bins = np.histogram(dat, bins=np.arange(0,0.2,0.001))
    rbins = 0.5 * (bins[1:] + bins[:-1])
    binsize = bins[1] - bins[0]
    nbins = hist.shape[0]

    # mirror distribution
    peak = np.argmax(hist)
    right = hist[peak:]
    left = right[::-1]
    tot = np.append([left], [right])
    ntot = tot.shape[0]
    lbins = rbins*-1.
    bincens = np.append([lbins[::-1]], [rbins])

    # fit a gaussian
    gaussfunc = lambda x,a,mu,sig: a * np.exp(-(x-mu)**2 / (2.*sig**2))
    p0 = [np.max(tot), bincens[ntot/2.],(np.max(bincens)-np.min(bincens))/4.]
    popt,pcov = curve_fit(gaussfunc, bincens, tot, p0=p0)
    mu = popt[1] if popt[1] >= 0. else 0.
    sigma = np.abs(popt[2])
    if showGauss:
        plt.bar(rbins, hist, align='center', width=binsize, alpha=0.3, 
                color='b', linewidth=0)
        plt.bar(bincens, tot, align='center', width=binsize, alpha=0.3, 
                color='k', linewidth=0)
        plt.plot(bincens, gaussfunc(bincens, *popt), 'k', linewidth=2)
        plt.show()

    # replace edges with random noise in the image
    bad = np.where((filt != 0) & (filt != 1))

    # random numbers with mean of 0 and sigma of 1*sig
    noise = np.random.normal(0, 1, drz[bad].shape[0])

    drz[bad] = mu + noise * sigma
    fits.writeto(output, drz, hdr, clobber=True)


def main():
    parser = argparse.ArgumentParser(description=
        'Clean edges and chip gap (if applicable) of an image')
    parser.add_argument('image', type=str, help='Filename of image to clean')
    parser.add_argument('--out', 
        help='Output name for cleaned file ([image]_cln.fits)')
    parser.add_argument('--check', action='store_true',
        help='Check filtered image of 1s and 0s')
    args = parser.parse_args()
    image = args.image
    output = args.out if args.out else os.path.splitext(image)[0]+'_cln.fits'
    check = args.check

    clean_image(image, output, check)

if __name__ == '__main__':
    main()
