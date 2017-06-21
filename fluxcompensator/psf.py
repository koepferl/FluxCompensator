from astropy import log as logger
from astropy.io import fits
import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.convolution import convolve_fft

def convolve_astropy(array, psf, mode='nearest'):
    # This does the mode='nearest' option
    if mode == 'nearest': 
        astropy_mode = 'edge'
    else:
        raise Exception('this function only works for mode "nearest" and moves it to "edge"')

    # Get shape of PSF
    ny, nx = psf.shape

    # Pad the image
    array = np.pad(array, nx, mode=astropy_mode)

    # Do the convolution
    array = convolve_fft(array, psf, normalize_kernel=True)

    # Crop the image
    array = array[nx:-nx, nx:-nx]

    return array


class GaussianPSF(object):

    '''
    Convolves 2D val array with Gaussian PSF. 
      
    Parameters
    ----------
    
    diameter : float
        Diameter of the telescope in cm. 
            
    '''

    def __init__(self, diameter):
        self.diameter = diameter

    def convolve(self, wav, array, resolution):
        '''
        Actually convolves 2D val array with Gaussian PSF with 
        ``sigma_psf = 0.44 * wav * 10.**(-4) / diameter / resolution['rad']``.
        
        Parameters
        ----------
        
        wav : float
            Current val of wavelength which is passed by SyntheticCube
            or SyntheticImage. 
        
        array : numpy.ndarray
            Flux array of SyntheticCube or SyntheticImage.
            
        resolution : dict={'rad', 'arcsec'}
            Resolution of val slice passed by the SyntheticCube
            or SyntheticImage. Here 'rad' needs to be defined.
        
                
        Returns
        -------
    
        val : numpy.ndarray
            2D val arrays convolved with PSF.
        
        '''

        sigma_psf = 0.44 * wav * 10. ** (-4) / self.diameter / resolution['rad']  # pixel
        #print wav, sigma_psf

        logger.debug(str('%9.4f' % wav) + ' microns ; sigma = ' + str('%9.4f' % sigma_psf) + ' pixel')

        val = gaussian_filter(array, sigma_psf, mode='nearest')

        return val


class FilePSF(object):

    '''
    Convolves val with PSF from input file.
    
    Parameters
    ----------
            
    psf_file : str
        Location of the file e. g. psf_file.fits.
        WARNING: If file from a detector/pipeline (e.g. IRAC4_PSF) is called, than the all slices of the cube, even outside the corresponding filter limits are convolved with that PSF or different filter could have been used before.
                    
    oversampled : int
        Times of oversampeling in PSF. Default is ``None``.
    '''

    def __init__(self, psf_file, oversampled=None):
        contents = fits.open(psf_file)
        self.psf = contents[0].data

        if oversampled is not None:  # importaint for Spitzer
            psf_pfr = self.psf[::oversampled, ::oversampled]
            self.psf = psf_pfr / np.sum(psf_pfr)

    def convolve(self, array, resolution=None, wav=None):
        '''
        Actually convolves 2D val array with PSF from input file.
            
        Parameters
        ----------
            
        array : numpy.ndarray
            Flux array of SyntheticCube or SyntheticImage.
        
        resolution : ``None``
            Just here to keep convolve_PSF in FluxCompensator working. 
            
        wav : ``None``
            Just here to keep convolve_PSF in FluxCompensator working. 
                
                
        Returns
        -------
        
        val : numpy.ndarray
                2D val arrays convolved with PSF.
            
        '''

        if array.ndim != 2:
            raise Exception('WARNING: FilePSF only works when 2D arrays are past.')

        val = convolve_astropy(array, self.psf, mode='nearest')

        return val


class FunctionPSF(object):

    '''
    Convolves 2D val with defined PSF function.
    
    Parameters
    ----------
            
    psf_function : method
        Method ``psf_function(X, Y, wavelength)`` depends on:
        the pixel coordinates in x and y direction and the
        wavelength of the current slice.
        
    width : int
        Ratio of pixel area of val slice to PSF function slice.
    '''

    def __init__(self, psf_function, width):
        self.width = width
        self.psf_function = psf_function

    def convolve(self, wav, array, resolution=None):
        '''
        Actually convolves 2D val array with defined PSF.
            
        Parameters
        ----------
        
        wav : float
            Wavelength of val array in microns. 
            
        array : numpy.ndarray
            Flux array of SyntheticCube or SyntheticImage.
        
        
        resolution : ``None``
            Just here to keep convolve_PSF in FluxCompensator working.
        
                
        Returns
        -------
        
        val : numpy.ndarray
                2D val arrays convolved with PSF.

        '''

        if array.ndim != 2:
            raise Exception('WARNING: FileFunction only works when 2D arrays are past.')

        # grid of PSF
        end = int(self.width / 2.)
        x = np.arange(- end, end + 1)
        y = np.arange(- end, end + 1)

        X, Y = np.meshgrid(x, y)

        wavelength = wav * 10. ** (-4)
        function_grid = self.psf_function(X, Y, wavelength)

        # normalizing
        norm = np.sum(function_grid)
        function_grid_norm = function_grid / norm

        val = convolve_astropy(array, function_grid_norm, mode='nearest')

        return val
