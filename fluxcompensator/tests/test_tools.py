import pytest
import numpy as np
from hyperion.util.constants import pc, kpc, au

from ..utils.tools import properties, grid_units, get_slices, average_collapse

def test_run_properties():
    
    arr_0D = np.array(1)
    arr_1D = np.arange(3)
    arr_2D = np.ones((3,3))
    arr_3D = np.ones((3,3,3))
    
    # flux
    properties(wav=arr_0D, val=arr_0D)
    
    # sed
    properties(wav=arr_1D, val=arr_1D)
    
    # image
    properties(wav=arr_0D, val=arr_2D)
    
    # cube
    properties(wav=arr_1D, val=arr_3D)
    

def test_run_grid_units():
    assert grid_units(x= 2 * au)['grid_unit_name'] == 'au'
    assert grid_units(x= 2000 * au)['grid_unit_name'] == 'pc'
    assert grid_units(x= 100 * pc)['grid_unit_name'] == 'pc'
    assert grid_units(x= 200 * pc)['grid_unit_name'] == 'kpc'

def test_get_slices_exceptions():
    
    arr_0D = np.array(1)
    arr_2D = np.ones((3,3))
    
    # image
    with pytest.raises(Exception) as exc:
        get_slices(wav=arr_0D, val=arr_2D, wav_1=2., wav_2=10.)
    
    # flux
    with pytest.raises(Exception) as exc:
        get_slices(wav=arr_0D, val=arr_0D, wav_1=2., wav_2=10.)

def test_average_collapse_exceptions():
    
    arr_0D = np.array(1)
    arr_2D = np.ones((3,3))
    
    # image
    with pytest.raises(Exception) as exc:
        average_collapse(val=arr_2D)
    
    # flux
    with pytest.raises(Exception) as exc:
        average_collapse(val=arr_0D)