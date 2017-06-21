import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../'

from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on 2MASS survey with Filter objects, zero-magnitude fluxes and telescope diameter.
'''

diameter = 130.
# diameter from http://www.ipac.caltech.edu/2mass/releases/first/doc/sec3_1a.html

J_2MASS_FILTER = Filter(name = 'J_2MASS', filter_file = ROOT + 'missions/filter/J.txt', waf_0 = 1.235, alpha = 1., beta = 0.)
H_2MASS_FILTER = Filter(name = 'H_2MASS', filter_file = ROOT + 'missions/filter/H.txt', waf_0 = 1.662, alpha = 1., beta = 0.)
K_2MASS_FILTER = Filter(name = 'K_2MASS', filter_file = ROOT + 'missions/filter/K.txt', waf_0 = 2.159, alpha = 1., beta = 0.)
# filter_file from http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html, Cohen et al. (2003)
# waf_0 from http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html
# alpha, beta = educated guess

J_2MASS_FILTER_PLOT = OneFilter(name = 'J_2MASS', filter_file = ROOT + 'missions/filter/J.txt', beta = 0)
H_2MASS_FILTER_PLOT = OneFilter(name = 'H_2MASS', filter_file = ROOT + 'missions/filter/H.txt', beta = 0)
K_2MASS_FILTER_PLOT = OneFilter(name = 'K_2MASS', filter_file = ROOT + 'missions/filter/K.txt', beta = 0)

J_2MASS_ZERO = 1594.
H_2MASS_ZERO = 1024.
K_2MASS_ZERO = 666.7
# zero-point magnitude (Jy) from http://www.ipac.caltech.edu/2mass/releases/allsky/doc/sec6_4a.html, Cohen et al. (2003)
