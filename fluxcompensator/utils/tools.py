from astropy import log as logger
import numpy as np
from hyperion.util.constants import au, pc, kpc


def properties(wav, val):
    '''
    Tool to determine the properties of wav and val passed in the
    FluxCompensator.


    Parameters
    ----------

    wav : numpy.ndarray
        The wavelengths of val array in microns.

    val : numpy.ndarray
        The val array of the FluxCompensator object.


    Returns
    -------

    properties : dict
        Dictionary with relevant properties of the FluxCompensator object.

        ``properties{spacing_wav, pixel, shape}``

             * ``spacing_wav`` : float width between wav entries.
             * ``pixel`` : tuple of (x, y) entries.
             * ``shape`` : shape of cube, image, sed or flux.
    '''
    wav = np.array(wav)
    val = np.array(val)

    # spacing of wav
    if wav.ndim != 0:
        spacing_wav = np.log10(wav[0] / wav[-1]) / (len(wav) - 1)
    else:
        spacing_wav = None

    # pixel
    if val.ndim in (0, 1):
        pixel = (None, None)
    if val.ndim in (2, 3):
        pixel = (val.shape[0], val.shape[1])

    # shape
    if val.ndim == 0:
        shape = '()'
    if val.ndim == 1:
        shape = '(wav)'
    if val.ndim == 2:
        shape = '(x, y)'
    if val.ndim == 3:
        shape = '(x, y, wav)'

    return {'spacing_wav': spacing_wav, 'pixel': pixel, 'shape': shape}


def grid_units(x):
    #input in cgs
    '''
    Tool to determine the astronomical units of x.
    (e. g. grid units of the axis in the FOV).


    Parameters
    ----------

    x : float
        Length in cm.


    Returns
    -------

    gird_units : dict
        Dictionary with information of x.

        ``grid_units{grid_unit, grid_unit_name}``

             * ``grid_unit`` : Estimated unit of x in cm.
             * ``grid_unit_name`` : Name of new unit.

                 Defined possibilities are:
                     * ``'au'``
                     * ``'pc'``
                     * ``'kpc'``
    '''

    if abs(x) <= 1e3 * au:
        grid_unit = au
        grid_unit_name = 'au'
    if abs(x) > 1e3 * au and abs(x) <= 1e2 * pc:
        grid_unit = pc
        grid_unit_name = 'pc'
    if abs(x) > 1e2 * pc:
        grid_unit = kpc
        grid_unit_name = 'kpc'

    return {'grid_unit': grid_unit, 'grid_unit_name': grid_unit_name}


def get_slices(wav, val, wav_1, wav_2):
    '''
    Tool to extract slices of the val cube in SyntheticCube.


    Parameters
    ----------

    wav : numpy.ndarray
        The wavelengths of val array in microns.

    val : numpy.ndarray
        The val array of the FluxCompensator object.

    wav_1 : float
        Wavelength boundary in microns.

    wav_2 : float
        Wavelength boundary in microns.


    Returns
    -------

    slices : dict
        Dictionary relevant for the FluxCompensator object.

        ``slices{wav_short, val_short, filter_index}``

             * ``wav_short`` : Vector like wav but all entries outside
                               the filter boundaries are erased.
             * ``val_short`` : Shape like wav_short for val.
             * ``filter_index`` : Indices of original wav within the filter.
    '''

    wav = np.array(wav)
    val = np.array(val)
    spacing_wav = properties(wav, val)['spacing_wav']

    if val.ndim != 3 and val.ndim != 1:
        raise Exception('Warning: only 3D and 1D arrays can be used to get slices')

    # swapp if wav_1 < wav_2
    if wav_1 < wav_2:
        wav_1, wav_2 = wav_2, wav_1
        
    # get_slice spacing_wav
    cond = (10. ** (np.log10(wav) - spacing_wav / 2.) <= wav_1) & (10. ** (np.log10(wav) + spacing_wav / 2.) >= wav_2)
    filter_index = cond.nonzero()[0]
    
    if val.ndim == 1: 
        wav_short = wav[cond]
        val_short = val[cond]
    
    if val.ndim == 3:
        wav_short = wav[cond]
        val_short = val[:,:,cond]

    return {'wav_short': wav_short, 'val_short': val_short, 'filter_index': filter_index}


def average_collapse(val):
    '''
    Tool to collapse slices of cube or vector averaged into one direction.


    Parameters
    ----------

    val : numpy.ndarray
        The val array of FluxCompensator object in 3D or 1D.


    Returns
    -------

    val : numpy.ndarray
        The 2D val array used in SyntheticImage.


    val : numpy.ndarray
        The 0D val array used in SyntheticFlux.

    '''

    val = np.array(val)

    if val.ndim != 3 and val.ndim != 1:
        raise Exception('Warning only 3D cubes or 1D vectors can be used')

    if val.ndim == 3:
        val = np.sum(val, axis=2) / len(val[0, 0, :])
    if val.ndim == 1:
        val = np.sum(val) / len(val)

    return val


def central_wav(wav):
    '''
    Tool to estimate the central wavelength of an averaged collapsed val.


    Parameters
    ----------

    wav : numpy.ndarray
        The wavelength vector of val array in microns.

    Returns
    -------

    central_wav : float
        Central wavelength in microns.
    '''

    central_wav = 10 ** (0.5 * np.log10(wav[0] * wav[-1]))
    
    return central_wav


def where_is_1D(n, o):
    '''
    Localizes fragment of pixel when grid resolution is changed. See also
    ConservingZoom. Computes the fraction of the old pixel fragment in the
    new resolution pixel.

    Parameters
    ----------

    n : tuple
        Boundaries of pixels new in arcsec.

            ``n = (x0, x1, y0, y1)``


    o : tuple
        Boundaries of old pixels in arcsec.

            ``o = (x0, x1, y0, y1)``

    Returns
    -------

    frac : float
        Fragment of old pixel in new resolution pixel.

    '''
    d_o = o[1] - o[0]
    d_n = n[1] - n[0]

    # debugging comments
    # logger.debug('-'*70)
    #logger.debug('fluxcompensator.utils.tools where_is_1D')
    # logger.debug('-'*70)

    #logger.debug('locaton of n with respect to o')
    #logger.debug('l = left, r = right, m = middle, out = outside')
    # logger.debug('-'*30)

    if n[0] > o[1] or n[1] < o[0]:
        out = True
        frac = 0
        #logger.debug('out: ' + str(frac))
    else:
        out = None

    if n[0] >= o[0] and n[1] >= o[1] and out is not True:
        frac = (o[1] - n[0]) / d_o
        #logger.debug('l  : ' + str(frac))

    elif n[0] <= o[0] and n[1] <= o[1] and out is not True:
        frac = (n[1] - o[0]) / d_o
        #logger.debug('r  : ' + str(frac))

    elif n[0] >= o[0] and n[1] <= o[1] and out is not True:
        frac = d_n / d_o
        #logger.debug('lr : ' + str(frac))

    elif n[0] <= o[0] and n[1] >= o[1] and out is not True:
        frac = 1
        #logger.debug('m  : ' + str(frac))

    return frac