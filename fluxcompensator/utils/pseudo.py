import numpy as np

class Pseudo(object):

    '''
    Mimics the properties of an ModelOutput array.
    
    Parameters
    ----------
    wav : numpy.ndarray
        * SyntheticCube like: The 1D wavelength of val cube in microns.
        * SyntheticImage like: The 0D wavelength of val image in microns.
        * SyntheticSED like: The 1D wavelengths of val vector in microns.
        * SyntheticFlux like: The 0D wavelength of val scalar in microns.
    
    val : numpy.ndarray
        * SyntheticCube like: The 3D array cube with shape (x, y, wav).
        * SyntheticImage like: The 2D array image with shape (x, y).
        * SyntheticSED like: The 1D array vector with shape like wav.
        * SyntheticFlux like: The 0D array scalar with shape like wav.

    units : str
        Unit of val array at input. Valid options are:

            * ``'ergs/cm^2/s'``
            * ``'ergs/cm^2/s/Hz'``
            * ``'Jy'``
            * ``'mJy'``
            * ``'MJy/sr'`` (not for SyntheticSED)
        
    distance : str
        Distance to the observed object in cm.
        
    x_min : float 
        Physical offset from axis origin in FOV in cm.
        
    x_max : float 
        Physical offset from axis origin in FOV in cm.
    
    y_min : float
        Physical offset from axis origin in FOV in cm.
    
    y_max : float 
        Physical offset from axis origin in FOV in cm.
    
    lon_min : float
        Minimal longitudinal angle.
        
    lon_max : float
        Maximal longitudinal angle.
        
    lat_min : float
        Minimal latitudinal angle.
        
    lat_max : float
        Maximal latitudinal angle.
    
    pix_area_sr : float
        Pixel area per sr.
        
    ap_min : float
        Minimal aperture of the telescope. Default is ``None``.
        
    ap_max : float
        Maximal aperture of the telescope. Default is ``None``.
        
        
    Returns
    -------
    
    cube : Pseudo like ModelOutput 3D
        3D val array with ModelOutput properties.
    sed : Pseudo like ModelOutput 1D
        1D val array with ModelOutput properties.
    '''

    def __init__(self, wav, val, origin=None, units=None, distance=None, x_min=None, x_max=None, y_min=None, y_max=None, lon_min=None, lon_max=None, lat_min=None, lat_max=None, pix_area_sr=None, ap_min=None, ap_max=None):

        self.wav = np.array(wav)
        self.val = np.array(val)
        
        if origin is not None:
            self.units = origin.units
            self.distance = origin.distance
            self.x_min = origin.x_min
            self.x_max = origin.x_max
            self.y_min = origin.y_min
            self.y_max = origin.y_max
            self.lon_min = origin.lon_min
            self.lon_max = origin.lon_max
            self.lat_min = origin.lat_min
            self.lat_max = origin.lat_max
            self.pix_area_sr = origin.pix_area_sr
            
            from ..cube import SyntheticCube
            from ..image import SyntheticImage
            
            if not isinstance(origin, SyntheticImage) and not isinstance(origin, SyntheticCube):
                if ap_min is None or ap_max is None:
                    raise Exception('ap_min and ap_max need to be defined.')
                else:
                    self.ap_min = ap_min
                    self.ap_max = ap_max
            else:
                self.ap_min = None
                self.ap_max = None
                

        else:
            if units is not None and distance is not None and x_min is not None and x_max is not None and y_min is not None and y_max is not None and lon_min is not None and lon_max is not None and lat_min is not None and lat_max is not None and pix_area_sr is not None:
                self.units = units
                self.distance = distance
                self.x_min = x_min
                self.x_max = x_max
                self.y_min = y_min
                self.y_max = y_max
                self.lon_min = lon_min
                self.lon_max = lon_max
                self.lat_min = lat_min
                self.lat_max = lat_max
                self.pix_area_sr = pix_area_sr
                self.ap_min = ap_min
                self.ap_max = ap_max
            else:
                raise Exception('Either origin must be defined or units, ... , pix_area_sr must be defined.')
