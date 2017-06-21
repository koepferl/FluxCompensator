import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../'

from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on WISE survey with Filter objects, zero-magnitude fluxes and telescope diameter.
'''

diameter = 40.
# diameter from http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec3_2.html

WISE1_FILTER = Filter(name = 'WISE1', filter_file = ROOT + 'missions/filter/WISE1.txt', waf_0 = 3.3526 , alpha = 2., beta = 0.)
WISE2_FILTER = Filter(name = 'WISE2', filter_file = ROOT + 'missions/filter/WISE2.txt', waf_0 = 4.6028 , alpha = 2., beta = 0.)
WISE3_FILTER = Filter(name = 'WISE3', filter_file = ROOT + 'missions/filter/WISE3.txt', waf_0 = 11.5608, alpha = 2., beta = 0.)
WISE4_FILTER = Filter(name = 'WISE4', filter_file = ROOT + 'missions/filter/WISE4.txt', waf_0 = 22.0883, alpha = 2., beta = 0.)
# filter_file, waf_0 from http://wise2.ipac.caltech.edu/docs/release/prelim/expsup/sec4_3g.html, Wright et al. (2010)
# alpha from http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
# beta = educated guess

WISE1_FILTER_PLOT = OneFilter(name = 'WISE1', filter_file = ROOT + 'missions/filter/WISE1.txt', beta = 0)
WISE2_FILTER_PLOT = OneFilter(name = 'WISE2', filter_file = ROOT + 'missions/filter/WISE2.txt', beta = 0)
WISE3_FILTER_PLOT = OneFilter(name = 'WISE3', filter_file = ROOT + 'missions/filter/WISE3.txt', beta = 0)
WISE4_FILTER_PLOT = OneFilter(name = 'WISE4', filter_file = ROOT + 'missions/filter/WISE4.txt', beta = 0)

WISE1_ZERO = 390.540
WISE2_ZERO = 171.787
WISE3_ZERO = 31.674
WISE4_ZERO = 8.363
# zero-point magnitude (Jy) from Jarrett et al. (2011)