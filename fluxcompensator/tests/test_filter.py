import pytest

import numpy as np
import os

from fluxcompensator.filter import Filter


class TestFilter(object):
      
    def setup_method(self, method):
        np.random.seed(3)
        
        path = os.path.join(os.path.dirname(__file__), 'MIPS1.txt')
        self.filter = Filter(name='my_filter', filter_file=path, waf_0=24., alpha=1, beta=0)
                
    def test_rebin_incl_hidden(self):
        
        with pytest.raises(Exception) as exc:
            self.filter.rebin(wav=np.tile(np.arange(3), 3), val=np.tile(np.arange(3), 3), inclination=None)
        
    def test_rebin_incl_clip(self):
        
        with pytest.raises(Exception) as exc:
            self.filter.rebin(wav=np.tile(np.arange(50), 3), val=10 * np.tile(np.arange(50), 3), inclination='0:50')
            
    def test_rebin_flux(self):
        wav = np.arange(10,50)
        flux = np.random.random(len(wav)) * 10
        self.filter.rebin(wav=wav, val=flux, inclination=None)

    def test_rebin_boundaries(self):
        wav = np.arange(1,10)
        flux = np.random.random(len(wav)) * 10
        with pytest.raises(Exception) as exc:
            self.filter.rebin(wav=wav, val=flux, inclination=None)

    def test_plot(self):
        wav = np.arange(10,50)
        flux = np.random.random(len(wav)) * 10
        self.filter.rebin(wav=wav, val=flux, inclination=None)
        self.filter.plot()
        
        