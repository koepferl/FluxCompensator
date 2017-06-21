import pytest

import numpy as np
import os

from ..utils.one_filter import OneFilter

def test_run_default():
    path = os.path.join(os.path.dirname(__file__), 'IRAC4.txt')
    o = OneFilter(name='IRAC4_default', filter_file=path, beta=None)
    o.onefilter_plot(line='-r', style='linear', normalized=None, unit=None)
    o.onefilter_save_plot(dpi=None)

def test_run_style():
    path = os.path.join(os.path.dirname(__file__), 'IRAC4.txt')
    
    for style in ['linear', 'loglog', 'logx', 'logy']:
    
        # style
        o = OneFilter(name='IRAC4_style_' + style, filter_file=path, beta=None)
        o.onefilter_plot(line='-r', style=style, normalized=None, unit=None)
        o.onefilter_save_plot(dpi=None)
    
        # dic
        o = OneFilter(name='IRAC4_style_dic_' + style, filter_file=path, beta=None)
        o.onefilter_plot(line={'color' : 'r', 'linestyle' : '--', 'linewidth' : 2}, style=style, normalized=None, unit=None)
        o.onefilter_save_plot(dpi=None)

def test_run_normalize():
    path = os.path.join(os.path.dirname(__file__), 'IRAC4.txt')
    o = OneFilter(name='IRAC4_normalize', filter_file=path, beta=None)
    o.onefilter_plot(line='-r', style='linear', normalized=True, unit=None)
    o.onefilter_save_plot(dpi=None)

def test_run_unit():
    path = os.path.join(os.path.dirname(__file__), 'IRAC4.txt')
    path_MIPS = os.path.join(os.path.dirname(__file__), 'MIPS1.txt')
    
    # beta = None
    with pytest.raises(Exception) as exc:
        o = OneFilter(name='IRAC4_unit_ext', filter_file=path, beta=None)
        o.onefilter_plot(line='-r', style='linear', normalized=None, unit='unit energy')
    
    # unit energy
    o = OneFilter(name='IRAC4_unit_ener_beta0', filter_file=path, beta=0) # ok
    o.onefilter_plot(line='-r', style='linear', normalized=None, unit='unit energy')
    o.onefilter_save_plot(dpi=None)
    
    o = OneFilter(name='MIPS1_unit_ener_beta-1', filter_file=path_MIPS, beta=-1)
    o.onefilter_plot(line='-r', style='linear', normalized=None, unit='unit energy')
    o.onefilter_save_plot(dpi=None)
        
    # unit energy
    o = OneFilter(name='IRAC4_unit_elec_beta0', filter_file=path, beta=0)
    o.onefilter_plot(line='-r', style='linear', normalized=None, unit='unit electron')
    o.onefilter_save_plot(dpi=None)
    
    o = OneFilter(name='MIPS1_unit_elec_beta-1', filter_file=path_MIPS, beta=-1) # ok
    o.onefilter_plot(line='-r', style='linear', normalized=None, unit='unit electron')
    o.onefilter_save_plot(dpi=None)
