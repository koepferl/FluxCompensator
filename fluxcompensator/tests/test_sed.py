import pytest

from ..sed import SyntheticSED
from hyperion.model.sed import SED
from hyperion.util.constants import c, pc, au
import numpy as np
import os


class TestSyntheticSED(object):

    def setup_method(self, method):
        np.random.seed(3)
        nu = c / (np.logspace(-1,3,30) * 1e-4)
        val = np.random.random(30)
        
        self.sed = SED(nu=nu, val=val, units='Jy')
        self.sed.distance = 3 * pc
    
    ################
    # init
    ################
                
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2'])
    def test_init_read(self,units):
        
        if units in ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy']:
            sed = SyntheticSED(input_array=self.sed, name='test_sed', unit_out=units)            
            
        else:
            with pytest.raises(Exception) as exc:
                sed = SyntheticSED(input_array=self.sed, name='test_sed', unit_out=units)
            assert exc.value.args[0] == "WARNING: Input or Output units needs to differ from MJy/sr if Input_array is not SyntheticCube or SyntheticSED or HyperionCube"
        
        
    
    ################
    # extinction
    ################
                
    def test_extinction_run(self):
        
        sed = SyntheticSED(input_array=self.sed, name='test_cube')
    
        ext = sed.extinction(0)
        assert np.allclose(sed.val, ext.val)  
    
        ext1 = sed.extinction(1)
        assert not np.allclose(sed.val, ext1.val)  
        
        assert ext.shape == '(wav)'
        assert ext.pixel == (None, None)                

    def test_extinction_run_input(self):
        
        sed = SyntheticSED(input_array=self.sed, name='test_cube')
        ext1 = sed.extinction(1)
        
        path = os.path.join(os.path.dirname(__file__), 'extinction_law.txt')
        ext1i = sed.extinction(1, input_opacities=path)
        assert np.array_equal(ext1.val, ext1i.val)  

        path2 = os.path.join(os.path.dirname(__file__), 'extinction_law_reversed.txt')
        ext1r = sed.extinction(1, input_opacities=path2)
        
    
    ###################
    # convolve_filter
    ###################
                
    def test_convolve_filter_read(self):
        
        sed = SyntheticSED(input_array=self.sed, name='test_sed')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        filtered = sed.convolve_filter(filter_input, plot_rebin=True)
        
        assert filtered.val.ndim == 0
        assert filtered.wav == 70.
        assert filtered.shape == '()'
        assert filtered.pixel == (None, None)
        
    
    ####################
    # get_total_val
    ####################
    
    @pytest.mark.parametrize('units', ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy'])
    def test_get_total_val(self, units):
        
        sed = SyntheticSED(input_array=self.sed, name='test_sed', unit_out=units)
    
        tot1 = sed.get_total_val(wav_1=10., wav_2=30.)
        assert tot1.shape == '()'
        
        tot2 = sed.get_total_val(wav_1=30., wav_2=10.)
        assert tot2.shape == '()'
        
        assert np.allclose(tot1.val, tot2.val)
        
    
    ####################
    # plot_sed_multi_filter
    ####################
    
    def test_plot_sed_multi_filter(self):
        
        sed = SyntheticSED(input_array=self.sed, name='test_sed')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'PACS1_FILTER')
        filtered = sed.convolve_filter(filter_input, plot_rebin=True)

        a = sed.plot_sed_multi_filter(multi_filter_val=filtered.val, multi_filter_wav=filtered.wav, names=['PACS1_FILTER_PLOT'], filter_label_size=True, ymin=10. ** (-5), dpi=None)
        
    def test_plot_sed_multi_filter2(self):
        
        sed = SyntheticSED(input_array=self.sed, name='test_sed')
        
        import fluxcompensator.database.missions as filters
        filter_input = getattr(filters, 'MIPS1_FILTER')
        filtered = sed.convolve_filter(filter_input, plot_rebin=True)

        a = sed.plot_sed_multi_filter(multi_filter_val=filtered.val, multi_filter_wav=filtered.wav, names=['MIPS1_FILTER_PLOT'], filter_label_size=None, ymin=10. ** (-5), dpi=None)
        
