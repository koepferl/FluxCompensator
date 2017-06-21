import pytest

from ..cube import SyntheticCube
from hyperion.model.image import Image
from hyperion.util.constants import c, pc, au
import numpy as np
import os


class TestSyntheticCube(object):

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
                
    def test_init_read(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_init_units(self, units):
        
        # units works?
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
        
    
    ################
    # extinction
    ################
                
    def test_extinction_run(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
    
        ext = cube.extinction(0)
        assert np.allclose(cube.val, ext.val)  

        ext1 = cube.extinction(1)
        assert not np.allclose(cube.val, ext1.val)  
        

    def test_extinction_run_input(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
    
        ext = cube.extinction(1)
        assert not np.allclose(cube.val, ext.val)  
        
        path = os.path.join(os.path.dirname(__file__), 'extinction_law.txt')
        ext1 = cube.extinction(1, input_opacities=path)
        assert np.allclose(ext.val, ext1.val)  

        path2 = os.path.join(os.path.dirname(__file__), 'extinction_law_reversed.txt')
        ext1r = cube.extinction(1, input_opacities=path2)
        assert np.allclose(ext1.val, ext1r.val)  

    
    ####################
    # change_resolution
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_change_resolution(self, units):#, tmpdir):
        
        #os.chdir(tmpdir.strpath)
        
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
    
        zoom_in = cube.change_resolution(new_resolution=cube.resolution['arcsec'] * 0.4)
        
        assert (zoom_in.shape[0] >= cube.shape[0]) 
        assert (zoom_in.shape[1] >= cube.shape[1]) 

        zoom_out = cube.change_resolution(new_resolution=cube.resolution['arcsec'] * 2.5)

        assert (zoom_out.shape[0] <= cube.shape[0]) 
        assert (zoom_out.shape[1] <= cube.shape[1]) 
        
    def test_change_resolution_plot(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out='Jy')
    
        zoom_in = cube.change_resolution(new_resolution=cube.resolution['arcsec'] * 0.4, grid_plot='grid_plot')    
    
    ####################
    # central_pixel
    ####################
    #@pytest.mark.skip(reason="fix later")
    def test_central_pixel(self):
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        center = cube.central_pixel(dx=0.2, dy=0.3)
        
        assert np.sum(cube.val) == np.sum(center.val)
    
    
    
    ####################
    # convolve_psf
    ####################
    def test_psf_Gauss(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        # Gaussian PSF
        from sim_obs.psf import GaussianPSF
        psf_object = GaussianPSF(diameter=350.)
        psf = cube.convolve_psf(psf_object)
    
    def test_psf_Database(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        # Database PSF
        import sim_obs.database.missions as PSFs
        psf_object = getattr(PSFs, 'PACS1_PSF')
        psf = cube.convolve_psf(psf_object)
    
    def test_psf_File(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        path = os.path.join(os.path.dirname(__file__), '../database/missions/PSF/irac_ch1_flight.fits')
        
        # File PSF
        from sim_obs.psf import FilePSF
        psf_object = FilePSF(psf_file=path, oversampled=4)
        psf = cube.convolve_psf(psf_object)
    
    def test_psf_Function(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        # Function PSF
        def my_psf(X,Y,wavelength):
            # resolution in rad/pixel
            resolution = cube.FOV[0] / cube.distance / cube.pixel[0]

            # telescope diameter in cm
            D_telescope = 350.

            # standard deviation in pixel
            sigma = 0.44 * wavelength / D_telescope / resolution

            Z = np.exp(-(X**2 + Y**2)/(2 * sigma**2))

            return Z
        
        from sim_obs.psf import FunctionPSF
        psf_object = FunctionPSF(psf_function=my_psf, width=32)
        psf = cube.convolve_psf(psf_object)
    


    ###################
    # convolve_filter
    ###################
                
    def test_convolve_filter_read(self):
        
        # noise works?
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        
        import sim_obs.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        filtered = cube.convolve_filter(filter_input)
        
        assert len(filtered.val.shape) == 2
        
        filtered_plot = cube.convolve_filter(filter_input, plot_rebin=True, plot_rebin_dpi=100)

    ################
    # add_noise
    ################
                
    def test_add_noise_read(self):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube')
    
        noise_0 = cube.add_noise(mu_noise=0., sigma_noise=0.)
        noise_1 = cube.add_noise(mu_noise=0., sigma_noise=0.1, diagnostics=True)
        
        assert np.allclose(cube.val, noise_0.val)  
        assert not np.allclose(cube.val, noise_1.val)


    ####################
    # get_rough_sed
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_get_rough_sed(self, units):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
        r_sed = cube.get_rough_sed()
        
        assert r_sed.val.shape == r_sed.wav.shape
        
    ####################
    # get_total_val
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_get_total_val(self, units):
        
        cube = SyntheticCube(input_array=self.image, name='test_cube', unit_out=units)
    
        tot1 = cube.get_total_val(wav_1=10., wav_2=30.)
        assert tot1.shape == '()'
        
        tot2 = cube.get_total_val(wav_1=30., wav_2=10.)
        assert tot2.shape == '()'
        
        assert np.allclose(tot1.val, tot2.val)
        
    ################
    # plot_image
    ################
                
    def test_plot_image_read(self):
    
        cube = SyntheticCube(input_array=self.image, name='test_cube')
        cube.plot_image(name='init', wav_interest=60., set_cut=(1e-13, 9.4e-12),
                        single_cut=None, multi_cut=None, dpi=300)
        
        with pytest.raises(Exception) as exc:
            cube.plot_image(prefix=None, name=None, wav_interest=60., set_cut=(1e-13, 9.4e-12),
                            single_cut=None, multi_cut=None, dpi=300)
        assert exc.value.args[0] == "If prefix name is not given, you need to give the a name to enable the default naming chain."
        #
                        
        cube.plot_image(prefix='prefix_name_cube', wav_interest=60., set_cut=(1e-13, 9.4e-12),
                        single_cut=None, multi_cut=None, dpi=300)
        

        with pytest.raises(Exception) as exc:
            cube.plot_image(prefix='prefix_name_multiTrue_cube', wav_interest=60., set_cut=(1e-13, 9.4e-12),
                            single_cut=None, multi_cut=True, dpi=300)
        assert exc.value.args[0] == "If prefix naming is enabled only one plotting option can be chosen."
        
        with pytest.raises(Exception) as exc:
            cube.plot_image(prefix='prefix_name_multiNone_cube', wav_interest=60., set_cut=(1e-13, 9.4e-12),
                            
                            single_cut=True, multi_cut=None, dpi=300)
        assert exc.value.args[0] == "If prefix naming is enabled only one plotting option can be chosen."