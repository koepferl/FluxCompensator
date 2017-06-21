import pytest
import numpy as np
from hyperion.util.constants import pc, kpc, au

from ..utils.units import ConvertUnits

class TestUnits(object):

    def setup_method(self, method):
        
        self.arr_0D = np.array(1)
        self.arr_1D = np.arange(3)
        self.arr_1D_long = np.arange(4)
        self.arr_2D = np.ones((3,3))
        self.arr_3D = np.ones((3,4,4))


    def test_units_exceptions(self):
        
        with pytest.raises(Exception) as exc:
            ConvertUnits(wav=self.arr_1D, val=self.arr_0D)
        
        with pytest.raises(Exception) as exc:
            ConvertUnits(wav=self.arr_1D, val=self.arr_1D_long)

        with pytest.raises(Exception) as exc:
            ConvertUnits(wav=self.arr_1D, val=self.arr_2D)

        un = ConvertUnits(wav=self.arr_1D, val=self.arr_3D)
        un.get_unit(in_units='Jy', out_units='mJy')

    def test_units_string_exceptions(self):
        
        un = ConvertUnits(wav=self.arr_1D, val=self.arr_1D)
        
        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='Jy', out_units='typo')

        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='typo', out_units='Jy')
            
    def test_units_MJysr_exceptions(self):
        un = ConvertUnits(wav=self.arr_1D, val=self.arr_1D)
        un.get_unit(in_units='MJy/sr', out_units='Jy', distance = 3 * kpc, pixel=(10, 10), FOV=(3 * pc, 3 * pc))

        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='MJy/sr', out_units='Jy')

        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='MJy/sr', out_units='Jy', distance = 3 * kpc, pixel=(None, None), FOV=(3 * pc, 3 * pc))

        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='MJy/sr', out_units='Jy', distance = 3 * kpc, pixel=(None, None), FOV=(4 * pc, 3 * pc))

        with pytest.raises(Exception) as exc:
            un.get_unit(in_units='MJy/sr', out_units='Jy', distance = 3 * kpc, pixel=(10, 11), FOV=(3 * pc, 3 * pc))

    def test_units_mag(self):
        un = ConvertUnits(wav=self.arr_0D, val=self.arr_0D)
        un.get_unit(in_units='mag', out_units='Jy', zero_point=2.)


#def test_run_properties():
#    
#    
#    # flux
#    properties(wav=arr_0D, val=arr_0D)
#    
#    # sed
#    properties(wav=arr_1D, val=arr_1D)
#    
#    # image
#    properties(wav=arr_0D, val=arr_2D)
#    
#    # cube
#    properties(wav=arr_1D, val=arr_3D)
#    
#
#def test_run_grid_units():
#    assert grid_units(x= 2 * au)['grid_unit_name'] == 'au'
#    assert grid_units(x= 2000 * au)['grid_unit_name'] == 'pc'
#    assert grid_units(x= 100 * pc)['grid_unit_name'] == 'pc'
#    assert grid_units(x= 200 * pc)['grid_unit_name'] == 'kpc'
#
#def test_get_slices_exceptions():
#    
#    arr_0D = np.array(1)
#    arr_2D = np.ones((3,3))
#    
#    # image
#    with pytest.raises(Exception) as exc:
#        get_slices(wav=arr_0D, val=arr_2D, wav_1=2., wav_2=10.)
#    
#    # flux
#    with pytest.raises(Exception) as exc:
#        get_slices(wav=arr_0D, val=arr_0D, wav_1=2., wav_2=10.)
#
#def test_average_collapse_exceptions():
#    
#    arr_0D = np.array(1)
#    arr_2D = np.ones((3,3))
#    
#    # image
#    with pytest.raises(Exception) as exc:
#        average_collapse(val=arr_2D)
#    
#    # flux
#    with pytest.raises(Exception) as exc:
#        average_collapse(val=arr_0D)