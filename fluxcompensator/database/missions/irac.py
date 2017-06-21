import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../'

from ...psf import GaussianPSF, FilePSF, FunctionPSF
from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on IRAC detector with Filter objects, PSF objects including resolution and zero-magnitude fluxes.
'''

IRAC1_FILTER = Filter(name = 'IRAC1', filter_file = ROOT + 'missions/filter/IRAC1.txt', waf_0 = 3.550, alpha = 1., beta = 0.)
IRAC2_FILTER = Filter(name = 'IRAC2', filter_file = ROOT + 'missions/filter/IRAC2.txt', waf_0 = 4.493, alpha = 1., beta = 0.)
IRAC3_FILTER = Filter(name = 'IRAC3', filter_file = ROOT + 'missions/filter/IRAC3.txt', waf_0 = 5.731, alpha = 1., beta = 0.)
IRAC4_FILTER = Filter(name = 'IRAC4', filter_file = ROOT + 'missions/filter/IRAC4.txt', waf_0 = 7.872, alpha = 1., beta = 0.)
# filter_file from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/calibrationfiles/spectralresponse/, Quijada et al. (2004)
# waf_0 from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/18/
# alpha, beta from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/18/, Reach et al. (2005), Hora et al. (2008)

IRAC1_FILTER_PLOT = OneFilter(name = 'IRAC1', filter_file = ROOT + 'missions/filter/IRAC1.txt', beta = 0)
IRAC2_FILTER_PLOT = OneFilter(name = 'IRAC2', filter_file = ROOT + 'missions/filter/IRAC2.txt', beta = 0)
IRAC3_FILTER_PLOT = OneFilter(name = 'IRAC3', filter_file = ROOT + 'missions/filter/IRAC3.txt', beta = 0)
IRAC4_FILTER_PLOT = OneFilter(name = 'IRAC4', filter_file = ROOT + 'missions/filter/IRAC4.txt', beta = 0)

IRAC1_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/irac_ch1_flight.fits', oversampled = 4)
IRAC2_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/irac_ch2_flight.fits', oversampled = 4)
IRAC3_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/irac_ch3_flight.fits', oversampled = 4)
IRAC4_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/irac_ch4_flight.fits', oversampled = 4)
# psf_file from http://dirty.as.arizona.edu/~kgordon/mips/conv_psfs/conv_psfs.html
# oversampled from psf_file

IRAC1_PSF_RESOLUTION = 1.221
IRAC2_PSF_RESOLUTION = 1.213
IRAC3_PSF_RESOLUTION = 1.222
IRAC4_PSF_RESOLUTION = 1.220
# from psf_file

IRAC1_ZERO = 280.9
IRAC2_ZERO = 179.7
IRAC3_ZERO = 115.0
# zero-point magnitude (Jy) from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/, Reach et al. (2005)
IRAC4_ZERO = 64.9
# zero-point magnitude (Jy) from http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/17/