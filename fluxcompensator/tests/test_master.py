import pytest
import os
import numpy as np
from hyperion.util.constants import pc, kpc, au
from hyperion.model import ModelOutput

from ..interface import Interface2FITS
from ..database.compact_pipeline import PACS1

class TestInterface2FITS(object):

    def setup_method(self, method):
        path = os.path.join(os.path.dirname(__file__), 'hyperion_output.rtout')
        m = ModelOutput(path)
        self.image = m.get_image(group=0, inclination=0, distance=8.5*kpc, units='ergs/cm^2/s')
        

    def test_run_interface(self):
        path = os.path.join(os.path.dirname(__file__), 'pacs70.fits')
        
        synobs = Interface2FITS(obs=path, 
                model=self.image, compact_pipeline=PACS1, exposure=10, A_v=20)
        print np.shape(synobs.val)
        synobs.save2fits('test_py_interface_fits')
        synobs.add2observation('test_py_interface_2obs', position_pix=(3000,2500))