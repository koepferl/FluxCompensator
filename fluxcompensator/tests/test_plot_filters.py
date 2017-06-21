import pytest

import numpy as np
import os

from ..utils.plot_filters import PlotFilters

def test_run_default():
    pf = PlotFilters(style='linear', normalized=None, unit=None)
    pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)

def test_run_unit():
    pf = PlotFilters(style='linear', normalized=None, unit='unit energy')
    pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.collect_filters(filter_database='MIPS1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    
    pf = PlotFilters(style='linear', normalized=None, unit='unit electron')
    pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.collect_filters(filter_database='MIPS1_FILTER_PLOT', name=None, filter_file=None, beta=None)

def test_run_construct():
    pf = PlotFilters(style='linear', normalized=None, unit=None)
    pf.collect_filters(filter_database=None, name='IRAC1_constract', filter_file=os.path.join(os.path.dirname(__file__), 'IRAC1.txt'), beta=0)

def test_run_get_axis():
    
    for style in ['linear', 'loglog', 'logx', 'logy']:
        pf = PlotFilters(style=style, normalized=None, unit='unit energy')
        pf.collect_filters(filter_database='IRAC1_FILTER_PLOT', name=None, filter_file=None, beta=None)
        pf.get_axis()
        
def test_run_plot_filters_plot_name():
    pf = PlotFilters(style='loglog', normalized=True, unit='unit energy')
    i1 = pf.collect_filters(filter_database='IRAC1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i2 = pf.collect_filters(filter_database='IRAC2_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i3 = pf.collect_filters(filter_database='IRAC3_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i4 = pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)
    m1 = pf.collect_filters(filter_database='MIPS1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.get_axis()
    pf.plot(data_filter=i1, line='r:')
    pf.plot(data_filter=i2, line='r-.')
    pf.plot(data_filter=i3, line='r--')
    pf.plot(data_filter=i4, line='r-')
    pf.plot(data_filter=m1, line='k-')
    pf.save_plot(name='plot_filters', dpi=100)
    
def test_run_plot_filters_plot_name2():
    pf = PlotFilters(style='loglog', normalized=True, unit='unit electron')
    i1 = pf.collect_filters(filter_database='IRAC1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i2 = pf.collect_filters(filter_database='IRAC2_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i3 = pf.collect_filters(filter_database='IRAC3_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i4 = pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)
    m1 = pf.collect_filters(filter_database='MIPS1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.get_axis()
    pf.plot(data_filter=i1, line='r:')
    pf.plot(data_filter=i2, line='r-.')
    pf.plot(data_filter=i3, line='r--')
    pf.plot(data_filter=i4, line='r-')
    pf.plot(data_filter=m1, line='k-')
    pf.save_plot(name='plot_filters2', dpi=100)

    
def test_run_plot_filters_no_plot_name():
    pf = PlotFilters(style='loglog', normalized=True, unit='unit electron')
    i1 = pf.collect_filters(filter_database='IRAC1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i2 = pf.collect_filters(filter_database='IRAC2_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i3 = pf.collect_filters(filter_database='IRAC3_FILTER_PLOT', name=None, filter_file=None, beta=None)
    i4 = pf.collect_filters(filter_database='IRAC4_FILTER_PLOT', name=None, filter_file=None, beta=None)
    m1 = pf.collect_filters(filter_database='MIPS1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.get_axis()
    pf.plot(data_filter=i1, line='r:')
    pf.plot(data_filter=i2, line='r-.')
    pf.plot(data_filter=i3, line='r--')
    pf.plot(data_filter=i4, line='r-')
    pf.plot(data_filter=m1, line='k-')
    pf.save_plot(name=None, dpi=100)
    
def test_run_plot_filters_no_unit():
    pf = PlotFilters(style='loglog', normalized=None, unit=None)
    i1 = pf.collect_filters(filter_database='IRAC1_FILTER_PLOT', name=None, filter_file=None, beta=None)
    pf.get_axis()