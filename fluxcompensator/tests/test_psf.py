import pytest

import numpy as np
import os

class TestPSFs(object):

    def setup_method(self, method):
        np.random.seed(3)
    
        self.array_3D = np.random.random((10,20,30))
        self.array_2D = np.random.random((10,11))
            
    def test_convolve_astropy(self):
        from fluxcompensator.psf import convolve_astropy
        
        with pytest.raises(Exception) as exc:
            convolve_astropy(array=self.array_3D, psf=self.array_2D, mode='typo')

    def test_FilePSF(self):
        
        path = os.path.join(os.path.dirname(__file__), 'mips_24_100K.fits')
        
        from fluxcompensator.psf import FilePSF
        file_psf = FilePSF(psf_file=path)
        with pytest.raises(Exception) as exc:
            file_psf.convolve(array=self.array_3D)     
        
    def test_FunctionPSF(self):
        
        # Function PSF
        def my_psf(X,Y,wavelength):
            # resolution in rad/pixel
            resolution = self.sim_image.FOV[0] / self.sim_image.distance / self.sim_image.pixel[0]

            # telescope diameter in cm
            D_telescope = 350.

            # standard deviation in pixel
            sigma = 0.44 * wavelength / D_telescope / resolution

            Z = np.exp(-(X**2 + Y**2)/(2 * sigma**2))

            return Z
            
        
        from fluxcompensator.psf import FunctionPSF
        psf_object = FunctionPSF(psf_function=my_psf, width=32)
        
        with pytest.raises(Exception) as exc:
            psf_object.convolve(wav=20., array=self.array_3D)     
        
