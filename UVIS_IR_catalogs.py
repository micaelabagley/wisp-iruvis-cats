#! /usr/bin/env python
import argparse
import os
import subprocess
import re
import ConfigParser
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table
from datetime import datetime
from cleanedges_general import clean_image
from convolvebypsf import convolve


class SECatalogs():

    def __init__(self, par):
        self.par = par
        self.datadir = os.path.join('..', self.par, 'DATA/UVIS/IRtoUVIS')
    
        # check that top-level directory exists for output
        topdir = '../UVIScatalogs'
        if os.path.isdir(topdir) is False:
            os.mkdir(topdir)
        # create output directory
        self.outdir = os.path.join(topdir, self.par)
        try:
            os.mkdir(self.outdir)
        except OSError:
            print '\n%s already exists.' % self.outdir
            print '\tPlease remove before continuing.\n'
            exit()

        # config file for convolution thresholds & SExtractor parameters
        self.Config = ConfigParser.ConfigParser()
        self.Config.read('parameters.cfg')
        # get thresholds for convolution
        options = self.Config.options('thresholds')
        self.thresh = {}
        for option in options:
            self.thresh[option.upper()] = self.Config.get('thresholds', option)
            
        # get list of images, filters, rms maps, exptimes and 
        # zero points for this field
        self.par_info = self.setup()
        # reddest filter to be used for detection
        self.detect_filt = self.detection_image()


    def get_zp(self, image):
        '''Get zero point for WFC3 filters based on observation date'''
        # obs date for new photometry
        photdate = '2012-03-06'
        HSTdate = datetime.strptime(photdate, '%Y-%m-%d')
        # get DATE-OBS from header
        obsdate = fits.getheader(image)['DATE-OBS']
        date = datetime.strptime(obsdate, '%Y-%m-%d')

        zp = {}
        if date.date() >= HSTdate.date():
            # new zero points
            zp['F110W'] = 26.8223
            zp['F140W'] = 26.4524
            zp['F160W'] = 25.9463
            zp['F475X'] = 26.1579
            zp['F600LP'] = 25.8746
            zp['F606W'] = 26.0691
            zp['F814W'] = 25.0985
        if date.date() < HSTdate.date():
            # old zero points
            zp['F110W'] = 26.83
            zp['F140W'] = 26.46
            zp['F160W'] = 25.96
            zp['F475X'] = 26.15
            zp['F600LP'] = 25.85
            zp['F606W'] = 26.08
            zp['F814W'] = 25.09

        return zp

    
    def replace_zeros(self, image):
        '''Replace all 0s in an RMS map with NaNs for proper use
           with Sextractor
        '''
        im,hdr = fits.getdata(image, header=True)
        rms = np.where(im != 0, im, np.nan)
        tmp = os.path.splitext(image)[0] + '_nan.fits'
        output = os.path.join(self.outdir, os.path.basename(tmp))
        fits.writeto(output, rms, header=hdr, clobber=True)
        return output


    def detection_image(self):
        '''Determine reddest filter to use for detection'''
        ir_images = glob(os.path.join(self.datadir, 'F1*_UVIS_sci.fits'))
        ir_images.sort()
        return fits.getheader(ir_images[-1])['FILTER']


    def setup(self):
        '''Return lists of all images and filters used in this Par.
           We will use the unrotated images for use with a single psf
           Image filenames:
             ParXXX/DATA/UVIS/IRtoUVIS/FYYY_UVIS_sci.fits
             ParXXX/DATA/UVIS/IRtoUVIS/FYYY_UVIS_rms.fits
        '''
        images = glob(os.path.join(self.datadir, 'F*_UVIS_sci.fits'))

        # dictionary of zero points
        zps = self.get_zp(images[0])

        # build table
        t = Table(data=None, names=['filt', 'image', 'rms', 'exptime', 'zp'],
                  dtype=['S10', 'S60', 'S60', float, float])
        for image in images:
            filt = fits.getheader(image)['FILTER']

            # clean image for convolution
            print 'Cleaning %s' % os.path.basename(image)
            tmp = os.path.splitext(image)[0] + '_cln.fits'
            image_cln = os.path.join(self.outdir, os.path.basename(tmp))
#            clean_image(image, image_cln)

            # replace zeros with NaNs in rms images
            rms_0 = image.split('_sci.fits')[0] + '_rms.fits'
            rms_nan = self.replace_zeros(rms_0)

            exptime = fits.getheader(image)['EXPTIME']
            zp = zps[filt]

            t.add_row([filt, image_cln, rms_nan, exptime, zp])
        return t


    def convolve_images(self):
        '''Convolve each image to the psf of the detection image.

           IRAF's psfmatch first computes the psf matching function
           Astropy's convolve_fft convolves the higher resolution image
           to the lower resolution psf.
        '''
        # lower resolution filter for convolution is always detect_filt
        lopsf = '%s_psf.fits' % self.detect_filt
        
        # get filter number for filenames
        check = re.search('\d+', self.detect_filt)
        df = check.group(0)

        for i in range(len(self.par_info['filt'])):
            if self.par_info['filt'][i] != self.detect_filt:
                # get filter number for filenames
                check = re.search('\d+', self.par_info['filt'][i])
                f = check.group(0)

                # setup parameters for image convolution
                # filename of higher res psf
                hipsf = '%s_psf.fits' % self.par_info['filt'][i]

                # name of kernel to be created by psfmatch for convolution
                outker = os.path.join(self.outdir, 'ker_%sto%s.fits'%(f, df))

                # low freq threshold in frac of total input image 
                # spectrum power for filtering option "replace"
                threshold = float(self.thresh[self.par_info['filt'][i]])

                # filename of image to be convolved
                highresimg = self.par_info['image'][i]
                # also have to convolve the rms maps
                highresrms = self.par_info['rms'][i]

                # filename of output, convolved image
                outname = os.path.join(self.outdir, '%s_convto%s.fits' % 
                                        (self.par_info['filt'][i], df))

                # convolve image
#                convolve(lopsf, hipsf, outker, threshold, highresimg, outname)
                

    def single_SE(self, updates):
        '''Run SE in single image mode'''
        wdet = np.where(self.par_info['filt'] == self.detect_filt)
        image = self.par_info['image'][wdet][0]
        cat = '%s_cat.fits' % self.par_info['filt'][wdet][0]
        cat = os.path.join(self.outdir, cat)
        seg = '%s_seg.fits' % self.par_info['filt'][wdet][0]
        seg = os.path.join(self.outdir, seg)
    
        # get parameters for this filter
        options = self.Config.options(self.detect_filt)
        params = {}
        for option in options:
            params[option] = self.Config.get(self.detect_filt, option)
        # override any parameters or add new ones?
        params.update(updates)

        # set parameters specific to this image
        params['-gain'] = '%.1f' % self.par_info['exptime'][wdet]
        params['-mag_zeropoint'] = '%.4f' % self.par_info['zp'][wdet]
        params['-weight_image'] = self.par_info['rms'][wdet][0]
        params['-catalog_name'] = cat
        params['-checkimage_name'] = seg
        params['-checkimage_type'] = 'SEGMENTATION'

        # set up SE parameters
        args = ['sex', image]
        for key,value in params.iteritems():
            args.append(key)
            args.append(value)
        # run SE
#        subprocess.check_call(args) 
        print image
        print args
        print

 
    def dual_SE(self, updates):
        '''Run SE in dual image mode'''
        # get filter number for filenames
        check = re.search('\d+', self.detect_filt)
        df = check.group(0)

        wdet = np.where(self.par_info['filt'] == self.detect_filt)

        # run SE on all images except detection image
        wphot = np.where(self.par_info['filt'] != self.detect_filt)
        filts = self.par_info['filt'][wphot]
        images = []
        for i in range(filts.shape[0]):
            images.append(os.path.join(self.outdir, '%s_convto%s.fits' % 
                                (self.par_info['filt'][wphot][i], df)))

        #images = self.par_info['image'][wphot]
        rmss = self.par_info['rms'][wphot]
        exptimes = self.par_info['exptime'][wphot]
        zps = self.par_info['zp'][wphot]

        for i in range(len(images)):
            image = images[i]
            rms = rmss[i]
            cat = '%s_cat.fits' % filts[i]
            cat = os.path.join(self.outdir, cat)
            options = self.Config.options(filts[i])
            params = {}
            for option in options:
                params[option] = self.Config.get(filts[i], option)
            # detection parameters should be from detection image section
            for option in ['-detect_minarea', '-detect_thresh', '-filter',\
                           '-filter_name', '-deblend_nthresh', 
                           '-deblend_mincont', '-back_filtersize']:
                params[option] = self.Config.get(self.detect_filt, option)
            # update rms maps for dual image mode
            params['-weight_image'] ='%s,%s'%(self.par_info['rms'][wdet][0],rms)
            params.update(updates)

            # set parameters specific to this image
            params['-gain'] = '%.1f'%exptimes[i]
            params['-mag_zeropoint'] = '%.4f'%zps[i]
            params['-catalog_name'] = cat

            # set up SE parameters
            args = ['sex', self.par_info['image'][wdet][0], image]
            for key,value in params.iteritems():
                args.append(key)
                args.append(value)
#            subprocess.check_call(args)
            print image
            print args
            print


    def __str__(self):
        print '\n%s:\n' % self.par
        
        ret_str = ''
        for i in range(len(self.par_info['filt'])):
            ret_str += '%s\n   Image:   %s\n   Rms:     %s\n   Exptime: %.1f\n   Zp:      %.4f\n' % \
                (self.par_info['filt'][i], 
                 os.path.basename(self.par_info['image'][i]),
                 os.path.basename(self.par_info['rms'][i]),
                 self.par_info['exptime'][i], self.par_info['zp'][i])
        return ret_str


    def __repr__(self):
        return str(self)


def main():
    parser = argparse.ArgumentParser(description='' )
    parser.add_argument('par', type=str, help='')
    args = parser.parse_args()
    
    # user can enter 'Par123', '123' 'field123', etc. 
    # so use just the par number 
    check = re.search('\d+', args.par)
    par = 'Par%s' % check.group(0)   

    # clean images
    Cat = SECatalogs(par)
    print Cat

    # convolve all images and rms maps to detection image psf 
    Cat.convolve_images()

    # run SExtractor on detection image
    # update any SE parameter values?
    updates = {}
    Cat.single_SE(updates)

    # run SExtractor on all images in dual image mode with the
    # reddest filter used for detection
    updates = {}
    Cat.dual_SE(updates)
    


if __name__ == '__main__':
    main()
