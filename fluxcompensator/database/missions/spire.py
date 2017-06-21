import os
ROOT = os.path.dirname(os.path.abspath(__file__))+ '/../'

from ...psf import GaussianPSF, FilePSF, FunctionPSF
from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on SPIRE detector with Filter objects, PSF objects including resolution and zero-magnitude fluxes.
'''

SPIRE1_FILTER = Filter(name = 'SPIRE1', filter_file = ROOT + 'missions/filter/SPIRE1.txt', waf_0 = 250., alpha = 1., beta = -1.)
SPIRE2_FILTER = Filter(name = 'SPIRE2', filter_file = ROOT + 'missions/filter/SPIRE2.txt', waf_0 = 350., alpha = 1., beta = -1.)
SPIRE3_FILTER = Filter(name = 'SPIRE3', filter_file = ROOT + 'missions/filter/SPIRE3.txt', waf_0 = 500., alpha = 1., beta = -1.)
# filter_file from Tom
# waf_0 from http://herschel.esac.esa.int/Docs/SPIRE/html/
# alpha from http://herschel.esac.esa.int/hcss-doc-11.0/load/spire_drg/html/ch05s07.html
# beta = educated guess

SPIRE1_FILTER_PLOT = OneFilter(name = 'SPIRE1', filter_file = ROOT + 'missions/filter/SPIRE1.txt', beta = -1)
SPIRE2_FILTER_PLOT = OneFilter(name = 'SPIRE2', filter_file = ROOT + 'missions/filter/SPIRE2.txt', beta = -1)
SPIRE3_FILTER_PLOT = OneFilter(name = 'SPIRE3', filter_file = ROOT + 'missions/filter/SPIRE3.txt', beta = -1)

SPIRE1_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/theoretical_spire_beam_model_psw_V0_2.fits', oversampled = 10)
SPIRE2_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/theoretical_spire_beam_model_pmw_V0_2.fits', oversampled = 10)
SPIRE3_PSF = FilePSF(psf_file = ROOT + 'missions/PSF/theoretical_spire_beam_model_plw_V0_2.fits', oversampled = 10)
# psf_file from http://dirty.as.arizona.edu/~kgordon/mips/conv_psfs/conv_psfs.html
# oversampled from psf_file

SPIRE1_PSF_RESOLUTION = 6.000
SPIRE2_PSF_RESOLUTION = 6.000
SPIRE3_PSF_RESOLUTION = 6.000
# from psf_file

SPIRE1_ZERO = 0.06
SPIRE2_ZERO = 0.03
SPIRE3_ZERO = 0.01
# zero-point magnitude (Jy) from http://svo2.cab.inta-csic.es/theory/fps/index.php

