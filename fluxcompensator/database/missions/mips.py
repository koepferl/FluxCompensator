import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../'

from ...psf import GaussianPSF, FilePSF, FunctionPSF
from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on MIPS detector with Filter objects, PSF objects including resolution and zero-magnitude fluxes.
'''

MIPS1_FILTER = Filter(name = 'MIPS1', filter_file = ROOT + 'missions/filter/MIPS1.txt', waf_0 = 23.68, alpha = -2., beta = -1.)
MIPS2_FILTER = Filter(name = 'MIPS2', filter_file = ROOT + 'missions/filter/MIPS2.txt', waf_0 = 71.42, alpha = -2., beta = -1.)
MIPS3_FILTER = Filter(name = 'MIPS3', filter_file = ROOT + 'missions/filter/MIPS3.txt', waf_0 = 155.9, alpha = -2., beta = -1.)
# filter_file from http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/calibrationfiles/spectralresponse/
# waf_0 from http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/49/
# alpha from http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/51/
# beta = educated guess

MIPS1_FILTER_PLOT = OneFilter(name = 'MIPS1', filter_file = ROOT + 'missions/filter/MIPS1.txt', beta = -1)
MIPS2_FILTER_PLOT = OneFilter(name = 'MIPS2', filter_file = ROOT + 'missions/filter/MIPS2.txt', beta = -1)
MIPS3_FILTER_PLOT = OneFilter(name = 'MIPS3', filter_file = ROOT + 'missions/filter/MIPS3.txt', beta = -1)

MIPS1_PSF_10K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_10K.fits'  , oversampled = 5)
MIPS1_PSF_25K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_25K.fits'  , oversampled = 5)
MIPS1_PSF_50K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_50K.fits'  , oversampled = 5)
MIPS1_PSF_75K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_75K.fits'  , oversampled = 5)
MIPS1_PSF_100K  = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_100K.fits' , oversampled = 5)
MIPS1_PSF_3000K = FilePSF(psf_file = ROOT + 'missions/PSF/mips_24_3000K.fits', oversampled = 5)

MIPS2_PSF_10K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_10K.fits'  , oversampled = 5)
MIPS2_PSF_25K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_25K.fits'  , oversampled = 5)
MIPS2_PSF_50K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_50K.fits'  , oversampled = 5)
MIPS2_PSF_75K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_75K.fits'  , oversampled = 5)
MIPS2_PSF_100K  = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_100K.fits' , oversampled = 5)
MIPS2_PSF_3000K = FilePSF(psf_file = ROOT + 'missions/PSF/mips_70_3000K.fits', oversampled = 5)

MIPS3_PSF_10K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_10K.fits'  , oversampled = 5)
MIPS3_PSF_25K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_25K.fits'  , oversampled = 5)
MIPS3_PSF_50K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_50K.fits'  , oversampled = 5)
MIPS3_PSF_75K   = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_75K.fits'  , oversampled = 5)
MIPS3_PSF_100K  = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_100K.fits' , oversampled = 5)
MIPS3_PSF_3000K = FilePSF(psf_file = ROOT + 'missions/PSF/mips_160_3000K.fits', oversampled = 5)
# psf_file from http://dirty.as.arizona.edu/~kgordon/mips/conv_psfs/conv_psfs.html
# oversampled from psf_file

MIPS1_PSF_RESOLUTION = 2.49
MIPS2_PSF_RESOLUTION = 9.85
MIPS3_PSF_RESOLUTION = 16.00
# from psf_file

MIPS1_ZERO = 7.17
MIPS2_ZERO = 0.778
MIPS3_ZERO = 0.159
# zero-point magnitude (Jy) from http://irsa.ipac.caltech.edu/data/SPITZER/docs/mips/mipsinstrumenthandbook/49/
