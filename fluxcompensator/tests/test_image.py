import pytest

from ..cube import SyntheticCube
from ..image import SyntheticImage
from hyperion.model.image import Image
from hyperion.util.constants import c, pc, au
import numpy as np
import os

class TestSyntheticImage(object):

    def setup_method(self, method):
        np.random.seed(3)
    
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random((10,20,30))
        
        self.output = Image(nu=nu, val=val, units='Jy')
        self.output.distance = 3 * pc
        self.output.x_min = - 10 * au
        self.output.x_max =  10 * au
        self.output.y_min = - 10 * au
        self.output.y_max =  10 * au
        
        cube = SyntheticCube(input_array=self.output, name='test_image')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        self.sim_image = cube.convolve_filter(filter_input)
    
    ################
    # init
    ################
    
    def test_run(self):
        assert self.sim_image.val.ndim == 2
        assert self.sim_image.wav.ndim == 0
        assert self.sim_image.wav == np.array(70.)
        
        assert self.sim_image.shape == '(x, y)'
        assert self.sim_image.pixel != (None, None)
        
        
    def test_SyntheticImage(self):
        assert isinstance(self.sim_image, SyntheticImage) is True
        
    def test_wav_spacing(self):
        assert self.sim_image.spacing_wav is None
                
    ################
    # extinction
    ################
                
    def test_extinction_run(self):
        
        print self.sim_image.val.shape
        
        ext = self.sim_image.extinction(0)
        assert np.array_equal(self.sim_image.val, ext.val)  
        
        ext1 = self.sim_image.extinction(1)
        assert not np.array_equal(ext.val, ext1.val)  
                

    def test_extinction_run_input(self):
        
        ext1 = self.sim_image.extinction(1)
        
        path = os.path.join(os.path.dirname(__file__), 'extinction_law.txt')
        ext1i = self.sim_image.extinction(1, input_opacities=path)
        assert np.array_equal(ext1.val, ext1i.val)  

        path2 = os.path.join(os.path.dirname(__file__), 'extinction_law_reversed.txt')
        ext1r = self.sim_image.extinction(1, input_opacities=path2)
        assert np.array_equal(ext1i.val, ext1r.val)  
    
    ####################
    # change_resolution
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_change_resolution(self, units):
        
        cube = SyntheticCube(input_array=self.output, name='test_image', unit_out=units)
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image = cube.convolve_filter(filter_input)
        
        
        zoom_in = image.change_resolution(new_resolution=image.resolution['arcsec'] * 0.4)
        
        assert (zoom_in.shape[0] >= image.shape[0]) 
        assert (zoom_in.shape[1] >= image.shape[1]) 

        zoom_out = image.change_resolution(new_resolution=image.resolution['arcsec'] * 2.5, grid_plot=True)

        assert (zoom_out.shape[0] <= image.shape[0]) 
        assert (zoom_out.shape[1] <= image.shape[1]) 
    
    ####################
    # central_pixel
    ####################
    def test_central_pixel(self):
        
        center = self.sim_image.central_pixel(dx=0.2, dy=0.3)
        
        assert np.sum(self.sim_image.val) == np.sum(center.val)
    
    
    
    ####################
    # convolve_psf
    ####################
    def test_psf_Gauss(self):
        
        # Gaussian PSF
        from fluxcompensator.psf import GaussianPSF
        psf_object = GaussianPSF(diameter=350.)
        psf = self.sim_image.convolve_psf(psf_object)
    
    def test_psf_Database(self):
        
        # Database PSF
        import fluxcompensator.database.missions as PSFs
        psf_object = getattr(PSFs, 'PACS1_PSF')
        psf = self.sim_image.convolve_psf(psf_object)
    
    def test_psf_File(self):
        
        path = os.path.join(os.path.dirname(__file__), '../database/missions/PSF/irac_ch1_flight.fits')
        
        # File PSF
        from fluxcompensator.psf import FilePSF
        psf_object = FilePSF(psf_file=path, oversampled=4)
        psf = self.sim_image.convolve_psf(psf_object)
    
    def test_psf_Function(self):
        
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
        psf = self.sim_image.convolve_psf(psf_object)
    

    ################
    # add_noise
    ################
                
    def test_add_noise_read(self):
        
        noise_0 = self.sim_image.add_noise(mu_noise=0., sigma_noise=0.)
        noise_1 = self.sim_image.add_noise(mu_noise=0., sigma_noise=0.1, diagnostics=True)
        
        assert np.allclose(self.sim_image.val, noise_0.val)  
        assert not np.allclose(self.sim_image.val, noise_1.val)

    ####################
    # get_total_val
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_get_total_val(self, units):
        
        cube = SyntheticCube(input_array=self.output, name='test_cube', unit_out=units)
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image = cube.convolve_filter(filter_input)
        
        tot1 = image.get_total_val()
        assert tot1.shape == '()'
        
    ################
    # plot_image
    ################
                
    def test_plot_image_read(self):
        
        self.sim_image.plot_image(name='init', set_cut=(1e-13, 9.4e-12),
                        single_cut=None, multi_cut=None, dpi=300)
        
        with pytest.raises(Exception) as exc:
            self.sim_image.plot_image(prefix=None, name=None, set_cut=(1e-13, 9.4e-12),
                            single_cut=None, multi_cut=None, dpi=300)
        assert exc.value.args[0] == "If prefix name is not given, you need to give the a name to enable the default naming chain."
        #
                        
        self.sim_image.plot_image(prefix='prefix_name_image', set_cut=(1e-13, 9.4e-12),
                        single_cut=None, multi_cut=None, dpi=300)
        

        with pytest.raises(Exception) as exc:
            self.sim_image.plot_image(prefix='prefix_name_multiTrue_image', set_cut=(1e-13, 9.4e-12),
                            single_cut=None, multi_cut=True, dpi=300)
        assert exc.value.args[0] == "If prefix naming is enabled only one plotting option can be chosen."
        
        with pytest.raises(Exception) as exc:
            self.sim_image.plot_image(prefix='prefix_name_multiNone_image', set_cut=(1e-13, 9.4e-12),
                            single_cut=True, multi_cut=None, dpi=300)
        assert exc.value.args[0] == "If prefix naming is enabled only one plotting option can be chosen."        
        
        
    ###############
    # add_to_obs
    ###############
    
    def test_add_to_observations_exceptions(self):
        
        # default
        with pytest.raises(Exception) as exc:
            self.sim_image.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                           name='output_fits_file', position_pix=None, position_world=None, zero_edges=None)

        # zero_edges & wrong resolution
        with pytest.raises(Exception) as exc:
            self.sim_image.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                                          name='output_fits_file', position_pix=(3000,2500), position_world=None, zero_edges=True)

        # postion_world & wrong resolution
        with pytest.raises(Exception) as exc:
            self.sim_image.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                                          name='output_fits_file', position_pix=None, position_world=(0.2, 0.2), zero_edges=True)

        
        # position_pix & wrong resolution
        with pytest.raises(Exception) as exc:
             self.sim_image.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                                               name='output_fits_file', position_pix=(3000,2500), position_world=None, zero_edges=None)

        # position_pix & wrong resolution
        with pytest.raises(Exception) as exc:
             self.sim_image.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                                               name='output_fits_file', position_pix=(3000,2500), position_world=None, zero_edges=None)
                                               
    def test_add_to_observations_resolution(self):
        
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random((20,20,30))
        
        self.output = Image(nu=nu, val=val, units='Jy')
        self.output.distance = 3 * pc
        self.output.x_min = - 10 * au
        self.output.x_max =  10 * au
        self.output.y_min = - 10 * au
        self.output.y_max =  10 * au
        
        cube = SyntheticCube(input_array=self.output, name='test_image')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        image = cube.convolve_filter(filter_input)
    
        res = image.change_resolution(new_resolution=3.2)
        res.add_to_observation(fits_file=os.path.join(os.path.dirname(__file__), 'pacs70.fits'), 
                               name='output_fits_file', position_pix=(3000,2500), position_world=None, zero_edges=None)
    
