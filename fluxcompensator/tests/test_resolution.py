import pytest
import numpy as np

from ..utils.resolution import ConservingZoom, central

def test_ConservingZoom():
    
    arrays = [(5,5), (4,4), (4,5), (5,4)]
    for i in arrays:
        array = np.ones(i)
        regrid = ConservingZoom(array=array, initial_resolution=1., new_resolution=1.5)
        zoomout = regrid.zoom()
        
        regrid = ConservingZoom(array=array, initial_resolution=1., new_resolution=0.75)
        zoomin = regrid.zoom()
        
        # total sum is conserved?
        assert np.sum(array) == np.sum(zoomout)
        assert np.sum(array) == np.sum(zoomin)       

def test_zoom_grid_else():
    
    array = np.ones((5,5))
    regrid = ConservingZoom(array=array, initial_resolution=1., new_resolution=1.5)
    zoomout = regrid.zoom()
    regrid.zoom_grid(name='test_zoom_grid_else_1', dpi=100, zoom_in=1.)
    regrid.zoom_grid(name='test_zoom_grid_else_2', dpi=100, zoom_in=2.)
        
def test_central_exceptions():
    
    val = np.ones((10))
    
    with pytest.raises(Exception) as exc:
        central(array=val, dx=0.5, dy=0.5)
    
    val2 = np.ones((10,10))
    
    with pytest.raises(Exception) as exc:
        central(array=val2, dx=0.6, dy=0.5)

def test_central_zero_shift():
    
    val = np.ones((10,10))
    
    cent_x = central(array=val, dx=0., dy=0.5)
    cent_y = central(array=val, dx=0.5, dy=0.)
        
    assert np.sum(val) == np.sum(cent_x)
    assert np.sum(val) == np.sum(cent_y) 
    