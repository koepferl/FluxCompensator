import pytest

from ..cube import SyntheticCube
from ..sed import SyntheticSED
from hyperion.model.sed import SED
from hyperion.model.image import Image
from hyperion.util.constants import c, pc, au
import numpy as np
import os


class TestSyntheticFlux(object):

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
        
        
    
    ################
    # init
    ################
                
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_test_init_read(self, units):
    
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
        flux = cube.get_total_val(wav_1=10., wav_2=30.)      
        
    def test_prop(self):
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        flux = cube.get_total_val(wav_1=10., wav_2=30.)     
        
        assert flux.pixel == (None, None)
        assert flux.shape == '()' 
          
                
    
    
    ################
    # extinction
    ################
                
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_extinction_run(self, units):
    
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
        flux = cube.get_total_val(wav_1=10., wav_2=30.)
        
        ext = flux.extinction(0)
        assert np.allclose(flux.val, ext.val) 
     

        ext1 = flux.extinction(1)
        assert flux.val.shape == ext1.val.shape
        #assert not np.allclose(flux.val, ext1.val)  
    
    
        path = os.path.join(os.path.dirname(__file__), 'extinction_law.txt')
        ext1i = flux.extinction(1, input_opacities=path)
        assert np.array_equal(ext1.val, ext1i.val)  

        path2 = os.path.join(os.path.dirname(__file__), 'extinction_law_reversed.txt')
        ext1r = flux.extinction(1, input_opacities=path2)
        assert np.array_equal(ext1i.val, ext1r.val)  
                
    
    ################
    # magnitudes
    ################
                
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_magnitudes(self, units):
    
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
        flux = cube.get_total_val(wav_1=10., wav_2=30.)
        
        mag = flux.magnitudes(zero_point=10.)
        
        
    ####################
    # plot_sed_filter
    ####################
    
    def test_plot_sed_filter(self):
        np.random.seed(3)
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random(30)
        
        sed = SED(nu=nu, val=val, units='Jy')
        sed.distance = 3 * pc
        
        sed_sim = SyntheticSED(input_array=sed, name='test_sed')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        filtered = sed_sim.convolve_filter(filter_input, plot_rebin=True)
        
        a = filtered.plot_sed_filter(wav_sed=sed_sim.wav, val_sed=sed_sim.val, ymin=10. ** (-10), my_own_filter=None, dpi=None)
           

    def test_plot_sed_filter_collapsed(self):
        np.random.seed(3)
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random(30)
        
        sed = SED(nu=nu, val=val, units='Jy')
        sed.distance = 3 * pc
        
        sed_sim = SyntheticSED(input_array=sed, name='test_sed')
        
        collapsed = sed_sim.get_total_val(wav_1=1., wav_2=10.)
        a = collapsed.plot_sed_filter(wav_sed=sed_sim.wav, val_sed=sed_sim.val, ymin=10. ** (-10), my_own_filter=None, dpi=None)