import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/../'

from ...filter import Filter
from ...utils.one_filter import OneFilter

'''
Information on IRAS detector with Filter objects and zero-magnitude fluxes.
'''

diameter = 57. 
# diameter from http://irsa.ipac.caltech.edu/IRASdocs/iras_mission.html

IRAS1_FILTER = Filter(name = 'IRAS1', filter_file = ROOT + 'missions/filter/IRAS1.txt', waf_0 = 12. , alpha = 1., beta = -1.)
IRAS2_FILTER = Filter(name = 'IRAS2', filter_file = ROOT + 'missions/filter/IRAS2.txt', waf_0 = 25. , alpha = 1., beta = -1.)
IRAS3_FILTER = Filter(name = 'IRAS3', filter_file = ROOT + 'missions/filter/IRAS3.txt', waf_0 = 60. , alpha = 1., beta = -1.)
IRAS4_FILTER = Filter(name = 'IRAS4', filter_file = ROOT + 'missions/filter/IRAS4.txt', waf_0 = 100., alpha = 1., beta = -1.)
# filter_file, waf_0 from http://irsa.ipac.caltech.edu/IRASdocs/exp.sup/ch2/tabC5.html
# alpha from http://lambda.gsfc.nasa.gov/product/iras/docs/exp.sup/ch6/C3.html
# beta = educated guess


IRAS1_FILTER_PLOT = OneFilter(name = 'IRAS1', filter_file = ROOT + 'missions/filter/IRAS1.txt', beta = -1)
IRAS2_FILTER_PLOT = OneFilter(name = 'IRAS2', filter_file = ROOT + 'missions/filter/IRAS2.txt', beta = -1)
IRAS3_FILTER_PLOT = OneFilter(name = 'IRAS3', filter_file = ROOT + 'missions/filter/IRAS3.txt', beta = -1)
IRAS4_FILTER_PLOT = OneFilter(name = 'IRAS4', filter_file = ROOT + 'missions/filter/IRAS4.txt', beta = -1)

IRAS1_ZERO = 30.88
IRAS2_ZERO = 7.26
IRAS3_ZERO = 1.11
IRAS4_ZERO = 0.39
# zero-point magnitudes (Jy) from http://svo2.cab.inta-csic.es/theory/fps/index.php