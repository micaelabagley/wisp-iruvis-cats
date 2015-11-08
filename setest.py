#! /usr/bin/env python
import subprocess
import ConfigParser
from astropy.table import Table


#image = '../UVIScatalogs/Par340_1/F160W_UVIS_sci_cln.fits'
image = '../UVIScatalogs/Par333/F160W_UVIS_sci_cln.fits'

# get parameters for filter
Config = ConfigParser.ConfigParser()
Config.read('parameters.cfg')
options = Config.options('F160W')
params = {}
for option in options:
    params[option] = Config.get('F160W', option)

# alter detection thresholds
#detthresh = ['4', '3', '2', '1.5', '1.25', '1', '0.75', '0.5', '0.4', '0.3', '0.2']
detthresh = ['0.2']
# colors for region files
c = ['yellow', 'green', 'red', 'cyan', 'magenta', 'blue', 'orange', 'purple', 'white', 'red', 'green']

for j in range(len(detthresh)):
    det = detthresh[j]
    cat = '../UVIScatalogs/tests/F160W_dt_%s_cat.fits'%det
    seg = '../UVIScatalogs/tests/F160W_dt_%s_seg.fits'%det
 
    params['-filter_name'] = 'gauss_6.0_27x27.conv'#
    params['-detect_minarea'] = '30'
    params['-detect_thresh'] = det   
    params['-catalog_name'] = cat
    params['-checkimage_type'] = 'SEGMENTATION'
    params['-checkimage_name'] = seg
    params['-weight_image'] = '../UVIScatalogs/Par333/F160W_UVIS_rms.fits'
    
    # set up SE parameters
    args = ['sex', image]
    for key,value in params.iteritems():
        args.append(key)
        args.append(value)
    subprocess.check_call(args)

    # write ds9 region files
    rcat = Table.read(cat, format='fits')

    f = open('../UVIScatalogs/tests/dt_%s.reg'%det, 'w')
    f.write('global color=green dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')
    for i in range(rcat['NUMBER'].shape[0]):
        f.write('circle(%f,%f,10) # color=%s \n' % 
                (rcat['X_IMAGE'][i], rcat['Y_IMAGE'][i], c[j]))
    f.close()


