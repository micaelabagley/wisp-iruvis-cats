#! /usr/bin/env python
import os
import argparse
from astropy.io import fits
from astropy.convolution import convolve_fft
from pyraf.iraf import imarith, imcopy, imdelete, psfmatch

def zeropadfits(smfits, bigfits, padfits):
    '''Pads smfits with zeros to match size of bigfits.
       Result is padfits, centered as was smfits.
       Assumes smfits & bigfits are squares w/ odd # of pixels across.
    '''
    NY,NX = fits.getheader(bigfits)['NAXIS2'],fits.getheader(bigfits)['NAXIS1']
    ny,nx = fits.getheader(smfits)['NAXIS2'],fits.getheader(smfits)['NAXIS1']
    print "\nPadding 'smfits' at %ix%i to match 'bigfits' at %ix%i\n" % \
            (nx,ny,NX,NY)
    center = (NY+1)/2
    border = ny / 2
    lo = center - border
    hi = center + border
    croprange = '[%d:%d,%d:%d]' % (lo,hi,lo,hi)

    imarith(bigfits, '*', 0, padfits)
    imcopy(smfits, padfits+croprange)


def convolve(lopsf, hipsf, outker, threshold, highresimg, outname):
    '''Match the higher-res psf to the lower-res psf and 
    ### note for optical and UV, a threshold of 0.14 was found to be ideal.
    ### for two IR ones, a threshold of 0.03 seemed to work better. 
    '''
    # compute the psf matching function and save it to outker
    psfmatch(hipsf, lopsf, hipsf, outker, convolution='psf',
             background='none', threshold=threshold)

    # read in kernel
    k = fits.getdata(outker)
    # read in high res image to be convolved
    im,hdr = fits.getdata(highresimg, header=True)

    # convolve highresimg with the psf matching function in outker 
    # to produce output image
    convimg = convolve_fft(im, k)

    fits.writeto(outname, convimg, header=hdr, clobber=True)


def main():
    parser = argparse.ArgumentParser(description=
        'Convolve a higher resolution image to a lower resolution')
    parser.add_argument('lopsf', type=str, help='Filename of lower res psf')
    parser.add_argument('hipsf', type=str, help='Filename of higher res psf')
    parser.add_argument('outker', type=str,
        help='Name of kernel to be created by psfmatch for convolution')
    parser.add_argument('threshold', type=float, 
        help='Low freq threshold in frac of total input image spectrum power for filtering option "replace"')
    parser.add_argument('highresimg', type=str,
        help='Filename of image to be convolved')
    parser.add_argument('outname', type=str,
        help='Filename of output, convolved image')
    args = parser.parse_args()
    lopsf = args.lopsf
    hipsf = args.hipsf
    outker = args.outker
    threshold = args.threshold
    highresimg = args.highresimg
    outname = args.outname

    padpsf = 'padpsf'

    if os.path.exists(padpsf+'.fits'):
        imdelete(padpsf)
    if os.path.exists(outker+'.fits'):
        imdelete(outker)

    # check that psfs are the same size
    lnx,lny = fits.getheader(lopsf)['NAXIS1'], fits.getheader(lopsf)['NAXIS2']
    hnx,hny = fits.getheader(hipsf)['NAXIS1'], fits.getheader(hipsf)['NAXIS2']
    # the psfs should be square
    if (lnx != lny) | (hnx != hny):
        print '\npsfs not square.'
        exit()
    if (lnx != hnx):
        # find which is smaller
        print '\npsfs not equal size. Need to pad one with zeros.'
        print 'Not implemented yet. ...'
        exit()
    else:
        convolve(lopsf, hipsf, outker, threshold, highresimg, outname) 


if __name__ == '__main__':
    main()
