import os
ROOT = os.path.dirname(os.path.abspath(__file__))+ '/../'

from ...psf import GaussianPSF, FilePSF, FunctionPSF
from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on PACS detector with Filter objects, PSF objects including resolution and zero-magnitude fluxes.
'''

PACS1_FILTER = Filter(name = 'PACS1', filter_file = ROOT + 'missions/filter/PACS1.txt', waf_0 = 70. , alpha = 1., beta = -1.)
PACS2_FILTER = Filter(name = 'PACS2', filter_file = ROOT + 'missions/filter/PACS2.txt', waf_0 = 100., alpha = 1., beta = -1.)
PACS3_FILTER = Filter(name = 'PACS3', filter_file = ROOT + 'missions/filter/PACS3.txt', waf_0 = 160., alpha = 1., beta = -1.)
# filter_file from https://nhscsci.ipac.caltech.edu/sc/index.php/Pacs/FilterCurves
# waf_0 from http://herschel.esac.esa.int/Docs/PACS/html/ch03s02.html
# alpha from http://herschel.esac.esa.int/Docs/PACS/html/ch03s03.html
# beta = educated guess

PACS1_FILTER_PLOT = OneFilter(name = 'PACS1', filter_file = ROOT + 'missions/filter/PACS1.txt', beta = -1.)
PACS2_FILTER_PLOT = OneFilter(name = 'PACS2', filter_file = ROOT + 'missions/filter/PACS2.txt', beta = -1.)
PACS3_FILTER_PLOT = OneFilter(name = 'PACS3', filter_file = ROOT + 'missions/filter/PACS3.txt', beta = -1.)

PACS1_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/PSF_blue_slope-1_small.fits' , oversampled = 10)
PACS2_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/PSF_green_slope-1_small.fits', oversampled = 10)
PACS3_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/PSF_red_slope-1_small.fits'  , oversampled = 10)
# psf_file from http://dirty.as.arizona.edu/~kgordon/mips/conv_psfs/conv_psfs.html
# oversampled from psf_file

PACS1_PSF_RESOLUTION = 4.000
PACS2_PSF_RESOLUTION = 4.000
PACS3_PSF_RESOLUTION = 4.000
# from psf_file

PACS1_ZERO = 0.78
PACS2_ZERO = 0.38
PACS3_ZERO = 0.14
# zero-point magnitude (Jy) from http://svo2.cab.inta-csic.es/theory/fps/index.php

