import pytest

from ..cube import SyntheticCube
from hyperion.model.image import Image
from hyperion.util.constants import c, pc, au
import numpy as np
import os

from ..utils.pseudo import Pseudo


class TestSyntheticCube(object):

    def setup_method(self, method):
        np.random.seed(3)
        
        self.wav = np.logspace(-1,3,30)
        nu = c / (self.wav * 1e-4)
        self.val = np.random.random((10,20,30))
        
        self.image = Image(nu=nu, val=self.val, units='Jy')
        self.image.distance = 3 * pc
        self.image.x_min = - 10 * au
        self.image.x_max =  10 * au
        self.image.y_min = - 10 * au
        self.image.y_max =  10 * au
    
    ################
    # init
    ################
                
    def test_run_origin(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        new_cube = Pseudo(wav=self.wav, val=cube.val, origin=cube)
        
        assert np.allclose(cube.wav, new_cube.wav)   
        assert np.allclose(cube.val, new_cube.val)   
        
    def test_run_exception(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        with pytest.raises(Exception) as exc:
            new_cube = Pseudo(wav=self.wav, val=cube.val)
            
    def test_run_no_origin(self):
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        new_cube = Pseudo(wav=self.wav, val=cube.val, origin=None, units=2, distance=3, x_min=4, x_max=5, y_min=6, y_max=7, lon_min=8, lon_max=9, lat_min=10, lat_max=11, pix_area_sr=12, ap_min=13, ap_max=14)
        
        assert np.allclose(cube.wav, new_cube.wav)   
        assert np.allclose(cube.val, new_cube.val)
    
    def test_run_sed_and_flux(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        flux = cube.get_total_val(wav_1=10., wav_2=30.)
        
        with pytest.raises(Exception) as exc:
            new_flux_ori = Pseudo(wav=flux.wav, val=flux.val, origin=flux)

        new_flux_ori = Pseudo(wav=flux.wav, val=flux.val, origin=flux, ap_min=13, ap_max=14)
        assert np.allclose(flux.wav, new_flux_ori.wav)   
        assert np.allclose(flux.val, new_flux_ori.val)   
        

        new_flux_no_ori = Pseudo(wav=flux.wav, val=flux.val, origin=None, units=2, distance=3, x_min=4, x_max=5, y_min=6, y_max=7, lon_min=8, lon_max=9, lat_min=10, lat_max=11, pix_area_sr=12, ap_min=13, ap_max=14)
        
        assert np.allclose(flux.wav, new_flux_no_ori.wav)   
        assert np.allclose(flux.val, new_flux_no_ori.val)   
        
