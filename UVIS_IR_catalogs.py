#! /usr/bin/env python
import timing
import argparse
import os
import inspect
import subprocess
import re
import ConfigParser
from glob import glob
import numpy as np
from astropy.io import fits
from astropy.table import Table, Column
from datetime import datetime

from cleanedges import clean_image
from convolvebypsf import convolve
from utils import *

class SECatalogs():

    def __init__(self, par, datadir, topdir, SE_only, no_conv):
        self.par = par
        self.datadir = datadir
        self.imdir = os.path.join(datadir, self.par, 'DATA/UVIS/IRtoUVIS')
        self.scriptdir = os.path.dirname(inspect.getfile(clean_image))

        if topdir is None:
            self.topdir = os.path.join(datadir, 'UVIScatalogs')
        else:
            self.topdir = topdir

        # check that top-level directory exists for output
        if os.path.isdir(self.topdir) is False:
            os.mkdir(self.topdir)
        # create output directory
        self.outdir = os.path.join(self.topdir, self.par)
        try:
            os.mkdir(self.outdir)
        except OSError:
            if SE_only is False:
                print '\n%s already exists.' % self.outdir
                print '\tPlease remove before continuing.\n'
                exit()

        if (SE_only is True) & (no_conv is False):
            # check that cleaned and convolved images alredy exist
            conv_images = glob(os.path.join(self.outdir, '*convto*.fits'))
            if len(conv_images) == 0:
                print '\nImages have not been convolved. Ignoring SE_only ' +\
                      'and beginning full process.\n'
                SE_only = False

        # reddest filter to be used for convolution
        self.reddest_filt = self.reddest_filter()
        # detection filter is always F110W
        self.detect_filt = 'F110W'
        self.detect_rms = None

        # config file for convolution thresholds & SExtractor parameters
        self.Config = ConfigParser.ConfigParser()
        self.Config.read(os.path.join(self.scriptdir,'parameters.cfg'))
        # get thresholds for convolution
        options = self.Config.options('thresholds')
        self.thresh = {}
        for option in options:
            self.thresh[option.upper()] = self.Config.get('thresholds', option)
            
        # get list of images, filters, rms maps, exptimes and 
        # zero points for this field
        self.par_info = self.setup(SE_only, no_conv)


    def get_zp(self, image):
        """Get zero point for WFC3 filters based on observation date"""
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

    
    def fix_rms_map(self, image, output, value=1.e10, rmstype='analysis', whtim=None):
        """Fix an RMS map for either detection or analysis with SExtractor.
           
           rmstype: either 'analysis' or 'detection'

           Detection:
               De-weight bad pixels in WHT map (replace <100 with 0.01 in RMS)

           Analysis: 
               Replace all 0s with large number (1e10) for proper use 
               with SExtractor
               (SExtractor v5 does not work properly with NaNs)
        """
        # read in rms map
        im,hdr = fits.getdata(image, header=True)
        if rmstype == 'analysis':
            # rms map for analysis: replace 0's with large number
            rms = np.where(im != 0, im, value)
        elif rmstype == 'detection':
            # rms map for detection: de-weight bad pixels in wht map
            if whtim is None:
                print '\nNo weight image provided\n'
                exit()
            wht = fits.getdata(whtim)
            rms = np.where(wht >= 100., im, value)
        else:
            print "\nUnknown rmstype: %s. ('analysis' or 'detection')\n"%rmstype
            exit()
        fits.writeto(output, rms, header=hdr, clobber=True)


    def reddest_filter(self):
        """Determine reddest filter to use for image convolution"""
        ir_images = glob(os.path.join(self.imdir, 'F1*_UVIS_sci.fits'))
        ir_images.sort()
        return fits.getheader(ir_images[-1])['FILTER']


    def setup(self, SE_only, no_conv):
        """Return lists of all images and filters used in this Par.
           We will use the unrotated images for use with a single psf
           Image filenames:
             ParXXX/DATA/UVIS/IRtoUVIS/FYYY_UVIS_sci.fits
             ParXXX/DATA/UVIS/IRtoUVIS/FYYY_UVIS_rms.fits
        """
        images = glob(os.path.join(self.imdir, 'F*_UVIS_sci.fits'))

        # dictionary of zero points
        zps = self.get_zp(images[0])

        # build table
        t = Table(data=None, 
                  names=['filt','image','convim','rms','wht','exptime','zp'],
                  dtype=['S10', 'S60', 'S60', 'S60', 'S60', float, float])
        for image in images:
            filt = fits.getheader(image)['FILTER']

            # weight map
            wht = image.split('_sci.fits')[0] + '_wht.fits'

            # clean image for convolution
            tmp = os.path.splitext(image)[0] + '_cln.fits'
            image_cln = os.path.join(self.outdir, os.path.basename(tmp))
            if SE_only is False:
                print 'Cleaning %s' % os.path.basename(image)
                if no_conv:
                    clean_image(image, image_cln, cln_by_wht=False, whtim=wht)
                else:
                    clean_image(image, image_cln, cln_by_wht=True, whtim=wht)
            
            # names of convolved images
            if filt == self.reddest_filt:
                convim = image_cln
            else:
                check = re.search('\d+', self.reddest_filt)
                rf = check.group(0)
                convim = os.path.join(self.outdir,'%s_convto%s.fits'%(filt,rf))

            # replace zeros with 1.e9 in rms analysis maps
            rms0 = image.split('_sci.fits')[0] + '_rms.fits'
            tmp = os.path.splitext(rms0)[0] + '_analysis.fits'
            rms_analysis = os.path.join(self.outdir, os.path.basename(tmp))
            self.fix_rms_map(rms0, rms_analysis, value=1.e10,rmstype='analysis')

            # for detection image, create detection RMS map as well
            if filt == self.detect_filt:
                tmp2 = os.path.splitext(rms0)[0] + '_detection.fits'
                rms_detect = os.path.join(self.outdir, os.path.basename(tmp2))
                self.fix_rms_map(rms0, rms_detect, value=0.01, 
                                 rmstype='detection', whtim=wht)
            
            exptime = fits.getheader(image)['EXPTIME']
            zp = zps[filt]

            t.add_row([filt, image_cln, convim, rms_analysis, wht, exptime, zp])
        # set detection RMS map
        self.detect_rms = rms_detect
        return t


    def convolve_images(self, compute_kernel, SE_only):
        """Convolve each image to the psf of the reddest image.

           IRAF's psfmatch first computes the psf matching function
           Astropy's convolve_fft convolves the higher resolution image
           to the lower resolution psf.
        """
        # lower resolution filter for convolution is always reddest_filt
        lopsf = os.path.join(self.scriptdir, '%s_psf.fits' % self.reddest_filt)
        
        # get filter number for filenames
        check = re.search('\d+', self.reddest_filt)
        rf = check.group(0)

        for i in range(len(self.par_info['filt'])):
            if self.par_info['filt'][i] != self.reddest_filt:
                # get filter number for filenames
                check = re.search('\d+', self.par_info['filt'][i])
                f = check.group(0)

                # filename of output, convolved image
                outname = self.par_info['convim'][i]

                if SE_only:
                    continue

                # setup parameters for image convolution
                # filename of higher res psf
                hipsf = os.path.join(self.scriptdir, 
                                     '%s_psf.fits' % self.par_info['filt'][i])

                if compute_kernel is False:
                    # name of kernel to be used for convolution
                    outker = os.path.join(self.scriptdir, 
                                          'ker_%sto%s.fits'%(f, rf))
                else:
                    # name of kernel to be created by psfmatch for convolution
                    outker = os.path.join(self.outdir, 
                                          'ker_%sto%s.fits'%(f, rf))

                # low freq threshold in frac of total input image 
                # spectrum power for filtering option "replace"
                threshold = float(self.thresh[self.par_info['filt'][i]])

                # filename of image to be convolved
                highresimg = self.par_info['image'][i]
                # also have to convolve the rms maps
                highresrms = self.par_info['rms'][i]

                # convolve image
                print 'Convolving %s to %s' % (os.path.basename(highresimg), 
                                               self.reddest_filt)
                convolve(lopsf, hipsf, outker, threshold, highresimg, outname,
                         compute_kernel)

 
    def run_SE(self, updates, no_conv):
        """Run SE in dual image mode"""
        filts = self.par_info['filt']

        # get list of images, all of which will be run in dual image mode
        if no_conv:
            images = self.par_info['image']
        else:
            images = self.par_info['convim']

        # detection image
        wdet = np.where(self.par_info['filt'] == self.detect_filt)
        detim = images[wdet][0]

        for i in range(len(images)):
            image = images[i]

            rms = self.par_info['rms'][i]
            cat = os.path.join(self.outdir, '%s_cat.fits' % filts[i])   
            options = self.Config.options(filts[i])
            params = {}
            for option in options:
                params[option] = self.Config.get(filts[i], option)
            # detection parameters should be from detection image section
            for option in ['-detect_minarea', '-detect_thresh', '-filter',\
                           '-filter_name', '-deblend_nthresh', \
                           '-deblend_mincont', '-back_filtersize', \
                           '-thresh_type', '-back_size']:
                params[option] = self.Config.get(self.detect_filt, option)

            # update rms maps for dual image mode
            params['-weight_image'] = '%s,%s' % (self.detect_rms, rms)
            
            params.update(updates)
            
            # fix filenames of SExtractor required files
            for key in ['-c', '-parameters_name', '-filter_name', \
                        '-starnnw_name']:
                params[key] = os.path.join(self.scriptdir, params[key])

            # set parameters specific to this image
            params['-gain'] = '%.1f' % self.par_info['exptime'][i]
            params['-mag_zeropoint'] = '%.4f' % self.par_info['zp'][i]
            params['-catalog_name'] = cat
            # segmentation map, filtered and background images
            if i == wdet[0][0]:
                seg = os.path.join(self.outdir, '%s_seg.fits' % filts[i])   
                bck = os.path.join(self.outdir, '%s_bck.fits' % filts[i])   
                fil = os.path.join(self.outdir, '%s_fil.fits' % filts[i])   
                aps = os.path.join(self.outdir, '%s_aps.fits' % filts[i])
                checkstr = 'SEGMENTATION,BACKGROUND,FILTERED,APERTURES'
                params['-checkimage_type'] = checkstr
                params['-checkimage_name'] = '%s,%s,%s,%s' % (seg,bck,fil,aps)

            # set up SE parameters
            args = [os.path.join(self.scriptdir, 'sex'), detim, image]
            for key,value in params.iteritems():
                args.append(key)
                args.append(value)
            subprocess.check_call(args)


    def join_cats(self, no_conv, check_phot):
        """ """
        if no_conv:
            images = self.par_info['image']
        else:
            images = self.par_info['convim']

        # read in F110W catalog (all UVIS fields have F110W)
        f110cat = os.path.join(self.outdir, 'F110W_cat.fits')
        f = fits.open(f110cat)
        uvisdata = f[1].data
        f.close()
        # get all other catalogs
        cats = glob(os.path.join(self.outdir, 'F*cat.fits'))
        cats = [x for x in cats if x != f110cat]

        # read in original catalog
        c = os.path.join(self.datadir,self.par,'DATA/DIRECT_GRISM/fin_F110.cat')
        refdata = np.genfromtxt(c)
        idx,separc = match_cats(uvisdata['x_world'], uvisdata['y_world'], 
                                refdata[:,7], refdata[:,8])
        match = (separc.value*3600. <= 0.2)
        # d['x_world'][match]
        # dref[:,7][idx[match]]
        t = Table(uvisdata[match])  
     
        print
        for phot in ['ISO', 'AUTO', 'APER']: # 'PETRO'
            # fix errors
            w110 = np.where(self.par_info['filt'] == 'F110W')
            print 'fix errors for F110W, %s' % phot
            eflux,emag = calc_errors(t, images[w110][0], 
                                self.par_info['rms'][w110][0],
                                os.path.join(self.outdir,'F110W_seg.fits'), 
                                self.par_info['exptime'][w110][0], phot=phot)

            #eflux=eflux.reshape(t['FLUXERR_%s'%phot].shape)
            t['FLUXERR_%s'%phot] = eflux
            t['MAGERR_%s'%phot] = emag

            # rename F110 photometry columns
            t.rename_column('FLUX_%s'%phot, 'FLUX_%s_F110W'%phot)
            t.rename_column('FLUXERR_%s'%phot, 'FLUXERR_%s_F110W'%phot)
            t.rename_column('MAG_%s'%phot, 'MAG_%s_F110W'%phot)
            t.rename_column('MAGERR_%s'%phot, 'MAGERR_%s_F110W'%phot)
            
        # add column of indices from fin_F110.cat 
        t.add_column(Column(data=refdata[:,1][idx[match]], name='WISP_NUMBER'),
                     index=0)
        
        # add in photometry from other catalogs
        for i, cat in enumerate(cats):
            i += 1
            print cat
            f = fits.getdata(cat)
            d = f[match]
            filtstr = os.path.basename(cat).split('_')[0]
            for j,phot in enumerate(['ISO', 'AUTO', 'APER']):
                index = (i*12 + 2) + (j * 4)

                # fix errors - UVIS errors are fine
                if filtstr == 'F160W':
                    w160 = np.where(self.par_info['filt'] == 'F160W')
                    print 'fix errors for %s, %s' % (filtstr, phot)
                    eflux,emag = calc_errors(d, images[w160][0],
                                    self.par_info['rms'][w160][0],
                                    os.path.join(self.outdir,'F110W_seg.fits'), 
                                    self.par_info['exptime'][w160][0], 
                                    phot=phot)
                    d['FLUXERR_%s'%phot] = eflux
                    d['MAGERR_%s'%phot] = emag

                t.add_columns([Column(data=d['FLUX_%s'%phot], 
                                  name='FLUX_%s_%s'%(phot,filtstr)),
                               Column(data=d['FLUXERR_%s'%phot],
                                  name='FLUXERR_%s_%s'%(phot,filtstr)),
                               Column(data=d['MAG_%s'%phot],
                                  name='MAG_%s_%s'%(phot,filtstr)),
                               Column(data=d['MAGERR_%s'%phot], 
                                  name='MAGERR_%s_%s'%(phot,filtstr))], 
                               indexes=[index, index, index, index])

        # sort by WISP number
        t.sort(['WISP_NUMBER'])
        output = os.path.join(self.outdir, '%s_cat.fits'%self.par)
        t.write(output, format='fits')
        """
        hdu0 = fits.PrimaryHDU()
        hdu1 = fits.BinTableHDU(np.array(t))
        hdulist = fits.HDUList([hdu0, hdu1])
        output = os.path.join(self.outdir, '%s_cat.fits'%self.par)
        hdulist.writeto(output, clobber=True)
        """

        # make region file
        region_wcs(os.path.splitext(output)[0]+'.reg', 
                   t['X_WORLD'], t['Y_WORLD'], t['A_WORLD'], t['B_WORLD'],
                   t['THETA_WORLD'], t['NUMBER'])
        # make region file of original F110W catalog
        region_wcs(os.path.join(self.outdir, 'F110W_orig.reg'), 
                   refdata[:,7], refdata[:,8], refdata[:,9], refdata[:,10],
                   refdata[:,11], refdata[:,1], color='red')

        if check_phot is True:
            print c
            print f110cat
            check_conv_phot(c, output, 'F110W')
            c = os.path.join(self.datadir,self.par,'DATA/DIRECT_GRISM/fin_F160.cat')
            print c
            check_conv_phot(c, output, 'F160W')


    def __str__(self):
        print '\n%s:\n' % self.par
        
        ret_str = ''
        ret_str += 'Detection filter:   %s\nConvolution filter: %s\n\n' % \
                    (self.detect_filt, self.reddest_filt)
        for i in range(len(self.par_info['filt'])):
            ret_str += '%s\n   Images:   %s, %s\n   Rms:     %s\n   Exptime: %.1f\n   Zp:      %.4f\n' % \
                (self.par_info['filt'][i], 
                 os.path.basename(self.par_info['image'][i]),
                 os.path.basename(self.par_info['convim'][i]),
                 os.path.basename(self.par_info['rms'][i]),
                 self.par_info['exptime'][i], self.par_info['zp'][i])
        return ret_str


    def __repr__(self):
        return str(self)


def main():
    parser = argparse.ArgumentParser(description='' )
    parser.add_argument('par', type=str, help='')
    parser.add_argument('--datadir', type=str, default='./',
        help='Path to WISP v5.0 reductions. Default is working directory.')
    parser.add_argument('--topdir', type=str, default='UVIScatalogs',
        help='Desired path for output. Default is UVIScatalogs.')
    parser.add_argument('--compute_kernel', action='store_true',
        help='Compute new kernels for convolutions?')
    parser.add_argument('--SE_only', action='store_true',
        help='Run SExtractor only. Assumes images already cleaned & convolved')
    parser.add_argument('--no_conv', action='store_true',
        help='Run SExtractor on unconvolved images? By default, SE runs on ' +\
             'convolved images')
    parser.add_argument('--check_phot', action='store_true',
        help='Compare convolved phot with original')
    args = parser.parse_args()
    
    # determine par number as user can enter 'Par123', '123' 'field123', etc.  
    check = re.search('\d+', args.par)
    par = 'Par%s' % check.group(0)   

    # clean images
    Cat = SECatalogs(par, args.datadir, args.topdir, args.SE_only, args.no_conv)

    # convolve all images to reddest image psf 
    if args.no_conv is False:
        Cat.convolve_images(args.compute_kernel, args.SE_only)
    print Cat

    # run SExtractor on all image in dual image mode
    #   F110W is used for detection, with the detection RMS map
    #   The analysis image (including F110W) is run with analysis RMS
    # update any SE parameter values?
    updates = {}
    Cat.run_SE(updates, args.no_conv)

    # match to original F110W catalog and join all photometry into output
    Cat.join_cats(args.no_conv, args.check_phot)
    


if __name__ == '__main__':
    main()
