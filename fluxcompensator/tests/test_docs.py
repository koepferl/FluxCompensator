import pytest
import os

class TestDOCS(object):
    
    def setup_method(self, method):
        
        import numpy as np
        
        from hyperion.model import ModelOutput
        from hyperion.util.constants import kpc
        from fluxcompensator.cube import *
        
        # read in from HYPERION
        m = ModelOutput(os.path.join(os.path.dirname(__file__), 'hyperion_output.rtout'))
        array = m.get_image(group=0, inclination=0, distance=10*kpc, units='ergs/cm^2/s')
                            
        # initial FluxCompensator array
        self.FC_object = SyntheticCube(input_array=array, unit_out='ergs/cm^2/s', name='test_cube')
    
    
    def test_docs_cube(self):
        
        # plot c.val at 60 microns
        self.FC_object.plot_image(name='init', wav_interest=60., set_cut=(1e-14, 9.4e-13),
                     single_cut=None, multi_cut=None, dpi=300)

    def test_docs_sed(self):
       
       import numpy as np
    
       from hyperion.model import ModelOutput
       from hyperion.util.constants import pc
       from fluxcompensator.sed import *
    
       # read in from HYPERION
       m = ModelOutput(os.path.join(os.path.dirname(__file__), 'B5_class2_45.rtout'))
       array = m.get_sed(group=0, inclination=0, distance=300*pc,
                         units='ergs/cm^2/s')
       
       # initial FluxCompensator array
       s = SyntheticSED(input_array=array, unit_out='ergs/cm^2/s', name='test_sed')
    
    
    def test_docs_extinction(self):
       
       # dered with provided extinction law
       ext = self.FC_object.extinction(A_v=20.)
       
       # plot ext.val (3D) at 60 microns
       ext.plot_image(name='ext', wav_interest=60., set_cut=(1e-14, 9.4e-13),
                      single_cut=None, multi_cut=None, dpi=300)
                      
    def test_docs_resolution(self):
       
       # change resolution
       zoom = self.FC_object.change_resolution(new_resolution=6., grid_plot=True)
       
       # plot zoom.val (3D) at 60 microns
       zoom.plot_image(name='zoom', wav_interest=60., set_cut=(1e-13, 9.4e-12),
                       single_cut=None, multi_cut=None, dpi=300)
    
    def test_docs_psf_file(self):
        
        from fluxcompensator.psf import FilePSF
        
        # create PSF from my own file
        file_psf = FilePSF(psf_file=os.path.join(os.path.dirname(__file__), 'mips_24_100K.fits'), oversampled=5)
        

    def test_docs_psf_file_database(self):
        
        from fluxcompensator.psf import FilePSF
        
        # create PSF from my own file
        file_psf = FilePSF(psf_file=os.path.join(os.path.dirname(__file__), 'mips_24_100K.fits'), oversampled=5)
        import fluxcompensator.database.missions as PSFs

        # call object from the psf database
        psf_object = getattr(PSFs, 'PACS1_PSF')
        
        # convolve with PSF
        psf = self.FC_object.convolve_psf(psf_object)
        
        # psf_object is FilePSF
        # plot psf.val (3D) at 60 microns
        psf.plot_image(name='psf_file', wav_interest=60., set_cut=(1e-14, 9.4e-13),
                       single_cut=None, multi_cut=None, dpi=300)
                       
    def test_docs_psf_function(self):
        
        import numpy as np
        
        def my_psf(X,Y,wavelength):
            # resolution in rad/pixel
            resolution = self.FC_object.FOV[0] / self.FC_object.distance / self.FC_object.pixel[0]

            # telescope diameter in cm
            D_telescope = 350.

            # standard deviation in pixel
            sigma = 0.44 * wavelength / D_telescope / resolution

            Z = np.exp(-(X**2 + Y**2)/(2 * sigma**2))

            return Z
            
        from fluxcompensator.psf import FunctionPSF

        # create PSF from a Function
        psf_object = FunctionPSF(psf_function=my_psf, width=32)
        
        # convolve with PSF
        psf = self.FC_object.convolve_psf(psf_object)

        # psf_object is FunctionPSF
        # plot psf.val (3D) at 60 microns
        psf.plot_image(name='psf_func', wav_interest=60., set_cut=(1e-14, 9.4e-13),
                       single_cut=None, multi_cut=None, dpi=300)

    def test_docs_psf_gaussian(self):
        
        from fluxcompensator.psf import GaussianPSF

        # create Gaussian PSF
        psf_object = GaussianPSF(diameter=350.)
        
        # convolve with PSF
        psf = self.FC_object.convolve_psf(psf_object)
        
        # psf_object is GaussianPSF
        # plot psf.val (3D) at 60 microns
        psf.plot_image(name='psf_gauss', wav_interest=60., set_cut=(1e-14,
                       9.4e-13), single_cut=None, multi_cut=None, dpi=300)
        
    def test_docs_filter_database(self):
        
        import fluxcompensator.database.missions as filters
    
        # call object from the filter database
        filter_input = getattr(filters, 'PACS1_FILTER')
        
        # convolve with filter object
        filtered = self.FC_object.convolve_filter(filter_input, plot_rebin=None,
                                             plot_rebin_dpi=None)
        
        # plot filtered.val (3D)
        filtered.plot_image(name='filter', set_cut=(1e-14, 9.4e-13),
                            single_cut=None, multi_cut=None, dpi=300)
    
    def test_docs_filter_database_plot_rebin(self):
        
        import fluxcompensator.database.missions as filters
    
        # call object from the filter database
        filter_input = getattr(filters, 'PACS1_FILTER')
        
        # convolve with filter object
        filtered = self.FC_object.convolve_filter(filter_input, plot_rebin=True, plot_rebin_dpi=300)    
                                             
    def test_docs_filter_read(self):
        
        import fluxcompensator.database.missions as filters
    
        # call object from the filter database
        filter_input_database = getattr(filters, 'IRAC1_FILTER')
        
        from fluxcompensator.filter import Filter
    
        # create own filter object
        filter_input = Filter(name='my_IRAC1', filter_file='IRAC1.txt',
                              waf_0=3.550, alpha=1, beta=0)
                              
        assert filter_input_database.waf_0 == filter_input.waf_0

    def test_docs_noise(self):
        
        # add noise
        noise = self.FC_object.add_noise(mu_noise=0, sigma_noise=1e-13, seed=2, diagnostics=None)
        
        # plot noise.val at 60 microns
        noise.plot_image(name='noise', wav_interest=60., set_cut=(1e-14, 9.4e-13),
                         single_cut=None, multi_cut=None, dpi=300)
                         
    def test_docs_example(self):
        import numpy as np

        from hyperion.model import ModelOutput
        from hyperion.util.constants import pc

        from fluxcompensator.cube import *
        from fluxcompensator.psf import *
        from fluxcompensator.utils.resolution import *

        # read in from HYPERION
        m = ModelOutput(os.path.join(os.path.dirname(__file__), 'hyperion_output.rtout'))
        array = m.get_image(group=0, inclination=0, distance=300*pc,
                            units='ergs/cm^2/s')

        # initial FluxCompensator array
        c = SyntheticCube(input_array=array, unit_out='ergs/cm^2/s',
                          name='test_cube')

        # dered with provided extinction law
        ext = c.extinction(A_v=20.)

        # change resolution to 10-times of the initial
        zoom = ext.change_resolution(new_resolution=10*ext.resolution['arcsec'],
                                     grid_plot=True)

        import fluxcompensator.database.missions as PSFs

        # call object from the psf database
        psf_object = getattr(PSFs, 'PACS1_PSF')

        # convolve with PSF
        psf = zoom.convolve_psf(psf_object)

        import fluxcompensator.database.missions as filters

        # call object from the filter database
        filter_input = getattr(filters, 'PACS1_FILTER')

        # convolve with filter
        filtered = psf.convolve_filter(filter_input, plot_rebin=None,
                                       plot_rebin_dpi=None)

        # add noise
        noise = filtered.add_noise(mu_noise=0, sigma_noise=5e-15, diagnostics=None)

    def test_docs_output_rough(self):
        
        # collapse 3D cube to rough SED
        rough = self.FC_object.get_rough_sed()

    def test_docs_output_total(self):
        
        # collapse FC_object.val (3D or 1D) within 30 and 60 microns
        FC_object_tot = self.FC_object.get_total_val(wav_1=30., wav_2=60.)
        
        # convert to magnitudes
        FC_object_mag = FC_object_tot.magnitudes(zero_point=64.9)
        
    def test_docs_output_image_plot(self):
        
        # plot FC_object.val (3D) at 60 microns with default naming
        self.FC_object.plot_image(name='plot', wav_interest=60., set_cut=(1e-14,
                             9.4e-13), single_cut=80., multi_cut=True, dpi=300)
        
        # plot FC_object.val (3D) at 60 microns with prefix
        self.FC_object.plot_image(prefix='prefix1', wav_interest=60., multi_cut=True,
                             dpi=300)
        self.FC_object.plot_image(prefix='prefix2', wav_interest=60., set_cut=(1e-14,
                             9.4e-13), dpi=300)
        self.FC_object.plot_image(prefix='prefix3', wav_interest=60., single_cut=80.,
                             dpi=300)
                             
    def test_docs_output_sed_plot(self):
        
        # collapse 3D cube to rough SED
        s = self.FC_object.get_rough_sed()
        print s.val

        import fluxcompensator.database.missions as filters

        # call object from the filter database
        filter_input = getattr(filters, 'PACS1_FILTER')

        # convolve with filter object
        filtered = self.FC_object.convolve_filter(filter_input, plot_rebin=None,
                                     plot_rebin_dpi=None)

        # collapse filtered.val
        f = filtered.get_total_val()

        # plot c_tot.val with c_sed.val in one plot
        f.plot_sed_filter(wav_sed=s.wav, val_sed=s.val, ymin=1.e-5, dpi=300)

    def test_docs_tutorial_multi_filter(self):
        
        import numpy as np
        
        # collapse 3D cube to rough SED
        s = self.FC_object.get_rough_sed()

        import fluxcompensator.database.missions as filters

        # empty arrays for storage
        val_array = np.array(())
        wav_array = np.array(())
        filter_array = np.array(())

        for loop_filter in ['J_2MASS', 'H_2MASS', 'K_2MASS', 'IRAS1', 'IRAS2',
                            'IRAS3', 'IRAS4','IRAC1', 'IRAC2', 'IRAC3', 'IRAC4',
                            'MIPS1', 'MIPS2', 'MIPS3','WISE1', 'WISE2', 'WISE3',
                            'WISE4', 'PACS1', 'PACS2', 'PACS3', 'SPIRE1', 'SPIRE2',
                            'SPIRE3']:

            # call object from the filter database
            filter_input = getattr(filters, loop_filter + '_FILTER')

            # convolve with filter object
            filtered = self.FC_object.convolve_filter(filter_input, plot_rebin=None,
                                         plot_rebin_dpi=None)

            # collapse FC_object.val
            f = filtered.get_total_val()

            # store f.val, f.wav and filter_input in arrays
            # for plot_sed_multi_filter()
            val_array = np.append(val_array, f.val)
            wav_array = np.append(wav_array, f.wav)
            filter_array = np.append(filter_array, f.filter['name'] + '_FILTER_PLOT')

        # plot all filters in loop with f.val and s.val
        s.plot_sed_multi_filter(multi_filter_val=val_array,
                                multi_filter_wav=wav_array, names=filter_array,
                                ymin=1e-5, filter_label_size=None, dpi=300)

    def test_docs_database(self):
        
        # telescopes
        from fluxcompensator.database.missions import *
        assert wise.diameter == 40.

        from fluxcompensator.database.missions import wise
        assert wise.diameter == 40.
        
        # PSFs
        from fluxcompensator.database.missions import *
        assert IRAC4_PSF_RESOLUTION == 1.22
        assert IRAC4_ZERO == 64.9

        from fluxcompensator.database.missions import irac
        assert irac.IRAC4_PSF_RESOLUTION == 1.22
        assert irac.IRAC4_ZERO == 64.9

        from fluxcompensator.database.missions import spitzer
        assert spitzer.IRAC4_PSF_RESOLUTION == 1.22
        assert spitzer.IRAC4_ZERO == 64.9

        # Access objects from database
        
        import fluxcompensator.database.missions as filters
        filter_object = getattr(filters, 'IRAC4_FILTER')

        import fluxcompensator.database.missions as PSFs
        psf_object = getattr(PSFs, 'IRAC4_PSF')

        import fluxcompensator.database.missions as Plots
        plot_filter = getattr(Plots, 'IRAC4_FILTER_PLOT')
        
        # Access attributes from objects in database
        
        from fluxcompensator.database.missions import *

        # name of Filter object
        assert IRAC4_FILTER.name == 'IRAC4'

        # central wavelength in microns
        assert IRAC4_FILTER.waf_0 == 7.872

        # power law exponent alpha
        assert IRAC4_FILTER.alpha == 1.0

        # power law exponent beta
        assert IRAC4_FILTER.beta == 0.0