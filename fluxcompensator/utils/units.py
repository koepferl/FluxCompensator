import numpy as np
from hyperion.util.constants import c


class ConvertUnits(object):

    '''
    Tool to convert units of val arrays.
    
    Parameters
    ----------
    
    wav : numpy.ndarray
        The wavelength of the val image in microns.
    
    val : numpy.ndarray
        Flux array from 3D to 0D.
    
    
    '''

    def __init__(self, wav, val):
        self.val = np.array(val)
        self.wav = np.array(wav)
        self.nu = c / (self.wav * 10. ** (-4))  # Hz

        # if dimension = 3 than the wav needs to be in the same dimension as the last enterie of val
        if np.shape(self.val) != ():
            self.dimension = len(np.shape(self.val))
        else:
            self.dimension = None

        if np.shape(self.wav) != ():
            self.lenght_wav = len(self.wav)
        else:
            self.lenght_wav = None

        if self.dimension is None:
            self.orientation = '()'
            if self.lenght_wav is not None:
                raise Exception('Warning: for scalar val input wav need to be scalar')

        if self.dimension == 1.:
            self.orientation = '(wav)'
            if self.lenght_wav != len(self.val):
                raise Exception('Warning: for vector val input wav need to be vector of same length')

        if self.dimension == 2.:
            self.orientation = '(x,y)'
            if self.lenght_wav is not None:
                raise Exception('Warning: for grid val inputs wav need to be a scalar')

        if self.dimension == 3:
            if np.shape(self.wav) == np.shape(self.val[0, 0, :]):
                self.orientation = '(x, y, wav)'

            elif np.shape(self.wav) == np.shape(self.val[:, 0, 0]):
                self.orientation = '(wav, x, y)'
                self.val = self.val.swapaxes(1, 2).swapaxes(0, 2)
                self.orientation = 'initial (wav, x, y) now (x, y, wav)'

    def get_unit(self, in_units, out_units, distance=None, FOV=None, pixel=None, input_resolution=None, zero_point=None):
        '''
        Convert val units from cgs to Jy like units.
        
        Parameters
        ----------
        
        in_units : str
            Unit of input val array. Valid options are:

                * ``'ergs/cm^2/s'``
                * ``'ergs/cm^2/s/Hz'``
                * ``'Jy'``
                * ``'mJy'``
                * ``'MJy/sr'``
                * ``'mag'``
                
        out_units : str
            Unit of output val array. Valid options are:

                * ``'ergs/cm^2/s'``
                * ``'ergs/cm^2/s/Hz'``
                * ``'Jy'``
                * ``'mJy'``
                * ``'MJy/sr'``
                * ``'mag'``
                
        distance : float
            Only needed if in_units or out_units is ``'MJy/sr'``.
            Distance of observed object in cm. 

        FOV : tuple
            Only needed if in_units or out_units is ``'MJy/sr'``.
            FOV (x,y) of observed object in cm. 

        pixel : tuple
            Only needed if in_units or out_units is ``'MJy/sr'``.
            Pixel in (x,y) direction.

        input_resolution : float
            Only needed if in_units or out_units is ``'MJy/sr'``.
            Force resolution to certain val in arcsec / pixel if any of 
            FOV, pixel and distance are not known.
            
        zero_point : float
            Zero-magnitude flux is need, if one wants to convert to ``mag``. 
            ``zero_point`` needs to be given in ``Jy``. Default is ``None``.
            
        
        Returns
        -------
        val : numpy.ndarray
            Flux array from 3D to 0D in units of out_units.
        '''

        # are the units defined
        if in_units not in ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2', 'mag']:
            raise Exception('This in_units is not defined use only "ergs/cm^2/s", "ergs/cm^2/s/Hz", "Jy", "mJy", "MJy/sr", "Jy/arcsec^2", "mag".')

        if out_units not in ['ergs/cm^2/s', 'ergs/cm^2/s/Hz', 'Jy', 'mJy', 'MJy/sr', 'Jy/arcsec^2', 'mag']:
            raise Exception('This out_units is not defined use only "ergs/cm^2/s", "ergs/cm^2/s/Hz", "Jy", "mJy", "MJy/sr", "Jy/arcsec^2", "mag".')

        # criteria for MJy/sr
        self.resolution = {}
        if in_units == 'MJy/sr' or out_units == 'MJy/sr' or in_units == 'Jy/arcsec^2' or out_units == 'Jy/arcsec^2':
            if input_resolution is None:
                # no distance, physical size and pixel
                if distance is None or FOV is None or pixel is None:
                    raise Exception('WARNING: define either distance, F0V and pixel or give input_resolution val in arcsec/pixel!')
                # no pixel or physical size
                if pixel[0] is None and FOV[0] is not None:
                    raise Exception('WARNING: number of pixel are not defined!')
                # no square
                if (FOV[1] is not None and FOV[1] != FOV[0]) or (pixel[1] is not None and pixel[1] != pixel[0]):
                    raise Exception('WARNING: Can not use MJy/sr => rectengular FOV')

                else:
                    self.resolution['rad'] = FOV[0] / pixel[0] / distance                   # resolution in rad
                    self.resolution['arcsec'] = np.degrees(self.resolution['rad']) * 3600   # resolution in arcsec

            if input_resolution is not None:
                self.resolution['arcsec'] = input_resolution
                self.resolution['rad'] = np.radians(self.resolution['arcsec'] / 3600.)  # resolution arcsec > rad

        if in_units == out_units:
            out = self.val

        else:
            # convert all to Jy
            if in_units == 'ergs/cm^2/s':
                out_Jy = self.val / self.nu * 10. ** (23)
            if in_units == 'ergs/cm^2/s/Hz':
                out_Jy = self.val * 10. ** (23)
            if in_units == 'Jy':
                out_Jy = self.val
            if in_units == 'mJy':
                out_Jy = self.val * 10. ** (-3)
            if in_units == 'MJy/sr':
                out_Jy = self.val * 10. ** (6) * self.resolution['rad'] ** 2
            if in_units == 'Jy/arcsec^2':
                out_Jy = self.val * self.resolution['arcsec'] ** 2
            if in_units == 'mag':
                out_Jy = zero_point * 10 ** (-self.val / 2.5)

            # convert to out_units
            if out_units == 'ergs/cm^2/s':
                out = out_Jy * 10. ** (-23) * self.nu
            if out_units == 'ergs/cm^2/s/Hz':
                out = out_Jy * 10. ** (-23)
            if out_units == 'Jy':
                out = out_Jy
            if out_units == 'mJy':
                out = out_Jy * 10. ** (3)
            if out_units == 'MJy/sr':
                out = out_Jy * 10. ** (-6) / self.resolution['rad'] ** 2
            if out_units == 'Jy/arcsec^2':
                out = out_Jy / self.resolution['arcsec'] ** 2
            if out_units == 'mag':
                out = 2.5 * np.log10(zero_point / out_Jy)

        # bring back to original orientation
        if self.orientation == 'initial (wav, x, y) now (x, y, wav)':
            out = out.swapaxes(0, 2).swapaxes(1, 2)
            self.orientation = 'back to initial (wav, x, y)'
        
        return out
