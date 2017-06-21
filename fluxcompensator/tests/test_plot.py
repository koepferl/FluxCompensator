import pytest

from ..cube import SyntheticCube
from hyperion.model.image import Image
from hyperion.util.constants import c, pc, au
import numpy as np
import os

from ..utils.plot import MakePlots


class TestPlot(object):

    def setup_method(self, method):
        np.random.seed(3)
        
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random((10,20,30))
        
        self.image = Image(nu=nu, val=val, units='Jy')
        self.image.distance = 3 * pc
        self.image.x_min = - 10 * au
        self.image.x_max =  10 * au
        self.image.y_min = - 10 * au
        self.image.y_max =  10 * au
        self.cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        print self.cube.filter

    def test_plot_exceptions(self):
        
        with pytest.raises(Exception) as exc:
            a = MakePlots(input_array=self.cube, wav_interest=None, prefix=None, name=None, multi_cut=None, single_cut=None, set_cut=None, dpi=None)
        
        with pytest.raises(Exception) as exc:
            a = MakePlots(input_array=self.cube, wav_interest=None, prefix=None, name=None, multi_cut=True, single_cut=None, set_cut=None, dpi=None)
        
        with pytest.raises(Exception) as exc:
            a = MakePlots(input_array=self.cube, wav_interest=None, prefix='name_prefix_exeptions', name=None, multi_cut=True, single_cut=None, set_cut=None, dpi=None)
        with pytest.raises(Exception) as exc:
            sed = self.cube.get_rough_sed()
            a = MakePlots(input_array=sed, wav_interest=None, prefix='name_prefix_exeptions', name=None, multi_cut=True, single_cut=None, set_cut=None, dpi=None)
                    
    
    def test_multi_cut(self):
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image_filtered = self.cube.convolve_filter(filter_input)
        
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix='name_prefix_multi', name=None, multi_cut=True, single_cut=None, set_cut=None, dpi=None)
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix=None, name='name_no_prefix_multi', multi_cut=True, single_cut=None, set_cut=None, dpi=None)

    def test_single_cut(self):
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image_filtered = self.cube.convolve_filter(filter_input)
        
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix='name_prefix_single', name=None, multi_cut=None, single_cut=90., set_cut=None, dpi=None)
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix=None, name='name_no_prefix_single', multi_cut=None, single_cut=90., set_cut=None, dpi=None)

    def test_histogram(self):
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image_filtered = self.cube.convolve_filter(filter_input)
        
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix='name_prefix1', name=None, multi_cut=None, single_cut=90., set_cut=None, dpi=None)
        a.histogram_cuts()
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix='name_prefix2', name=None, multi_cut=None, single_cut=90., set_cut=None, dpi=None)
        a.histogram_cuts(dpi=100)
        
        a = MakePlots(input_array=image_filtered, wav_interest=None, prefix=None, name='name_no_prefix3', multi_cut=None, single_cut=90., set_cut=None, dpi=None)
        a.histogram_cuts()
        
        