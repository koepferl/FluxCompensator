from copy import deepcopy
import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/'

from astropy import log as logger
from astropy.io import fits
import numpy as np
from numpy.random import normal

from .psf import GaussianPSF, FilePSF, FunctionPSF
from .filter import Filter

from .utils.plot import MakePlots
from .utils.resolution import ConservingZoom, central
from .utils.tools import properties, grid_units, get_slices, average_collapse, central_wav
from .utils.units import ConvertUnits

# submitting PhD thesis today :)

class SyntheticCube(object):

    '''
    SyntheticCube is part the FluxCompensator. It converts 
    input_arrays (e. g. HYPERION ModelOutput) to "realistic"
    synthetic observations (e.g. accounts for PSF, filters & noise). 
    It contains attributes like ModelOutput (see Notes).
    If input_array is already a SyntheticCube object, the attributes are
    passed. If input_array is not a SyntheticCube object, SyntheticCube
    specific attributes are defined and then passed. 
    
    Parameters
    ----------
    input_array : SyntheticCube, ModelOutput, optional
        input_array also reads arrays with ModelOutput like properties.
        
    unit_out : str, optional
        The output units for SyntheticCube val. Valid options are:

            * ``'ergs/cm^2/s'``
            * ``'ergs/cm^2/s/Hz'``
            * ``'Jy'``
            * ``'mJy'``
            * ``'MJy/sr'``
        
        The default is ``'ergs/cm^2/s'``.
                
    name : str
        The name of the FluxCompensator object until another
        input_array is called. The default is ``None``.
        
    
    Attributes
    ----------
    
    wav : numpy.ndarray
        The wavelengths of val cube slices in microns.
    
    val : numpy.ndarray
        The 3D cube with shape (x, y, wav).
    
    units : str
        Current units of the val cube.
        
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
        
     
    Notes
    -----

    unit_in : str
        Unit of val in input_array. Valid options are:
        
            * ``'ergs/cm^2/s'``
            * ``'ergs/cm^2/s/Hz'``
            * ``'Jy'``
            * ``'mJy'``
            * ``'MJy/sr'``
              
    grid_unit : float
        Physical unit of FOV axis in cm. Valid options are:
        
            * ``au`` in cm
            * ``pc`` in cm
            * ``kpc`` in cm
        
    grid_unit_name
        Astronomical unit of FOV axis. Valid options are:
        
            * ``'au'``
            * ``'pc'``
            * ``'kpc'``
            
    FOV : tuple
        Tuple ``FOV(x,y)`` of Field of View pixel entries.
        
            * pixel in x direction: ``FOV[0]``
            * pixel in y direction: ``FOV[1]``
                
    name : str
        The name of the FluxCompensator object until another  
        input_array is called. The default is ``None``.
    
    stage : str
        Gives current operation stage of SyntheticCube.
        E. g. ``'SyntheticCube: convolve_filter'``
        
    log : list
        List of strings of the previous and current stages.
        
    filter : dict
        Dictionary ``filter = {name, waf_0, waf_min, waf_max}`` 
        of the applied filter. 
        
            * name of filter:     ``filter['name']``
            * central wavelength: ``filter['waf_0']``
            * minimal wavelength:  ``filter['waf_min']``
            * maximal wavelength: ``filter['waf_max']``

    
    Returns
    -------
    
    cube : SyntheticCube
        3D val array with SyntheticCube properties.
        
    image : SyntheticImage
        2D val array with SyntheticImage properties.
        
    sed : SyntheticSED
        1D val array (collapsed rough SED) with SyntheticSED properties.
        
    flux : SyntheticFlux
        0D val array (scalar) with SyntheticFlux properties.
    '''

    def __init__(self, input_array, unit_out='ergs/cm^2/s', name=None):

        # Hyperion ModelOutput attributes
        #if input_array.val.ndim == 3:
        self.val = np.array(deepcopy(input_array.val))
        #else:
        #    raise Exception('input_array does not have the right dimensions. numpy array of (x, y, wav) is required.')

        self.wav = np.array(deepcopy(input_array.wav))
        self.units = input_array.units
        self.distance = input_array.distance
        self.x_max = input_array.x_max
        self.x_min = input_array.x_min
        self.y_max = input_array.y_max
        self.y_min = input_array.y_min
        self.lon_min = input_array.lon_min
        self.lon_max = input_array.lon_max
        self.lat_min = input_array.lat_min
        self.lat_max = input_array.lat_max
        self.pix_area_sr = input_array.pix_area_sr

        ##################
        # new attributes #
        ##################

        if isinstance(input_array, SyntheticCube):
            # attributes with are passed, since input_array is SyntheticCube

            # physical values
            self.unit_in = input_array.unit_in
            self.unit_out = input_array.unit_out

            self.grid_unit = input_array.grid_unit
            self.grid_unit_name = input_array.grid_unit_name

            # properties of cube
            self.FOV = deepcopy(input_array.FOV)

            # name
            self.name = input_array.name
            self.stage = input_array.stage
            self.log = deepcopy(input_array.log)

            # filter
            self.filter = deepcopy(input_array.filter)

        else:  # attributes are defined, since input_array is NOT SyntheticCube

            # physical values
            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = grid_units(self.x_max - self.x_min)['grid_unit']
            self.grid_unit_name = grid_units(self.x_max - self.x_min)['grid_unit_name']

            self.FOV = (self.x_max - self.x_min, self.y_max - self.y_min)

            # name
            self.name = name
            self.stage = 'SyntheticCube: initial'
            self.log = [self.stage]

            # filter
            self.filter = {'name': None, 'waf_0': None, 'waf_min': None, 'waf_max': None}

            # convert into val units into unit_out
            s = ConvertUnits(wav=self.wav, val=self.val)
            self.val = s.get_unit(in_units=self.unit_in, out_units=self.unit_out, input_resolution=self.resolution['arcsec'])

            self.units = self.unit_out

    def extinction(self, A_v, input_opacities=None):
        '''
        Accounts for reddening.
        
        Parameters
        ----------
        
        A_v : Value of the visible extinction.
                         
        input_opacities : ``None``, str
            If ``None`` standard extinction law is used. 
            Otherwise a e. g. input_opacities.txt file can be passed 
            as a str to read an opacity file with column #1 wav in microns 
            and column #2 in cm^2/g. 
            Default is ``None``.

                
        Returns
        -------
        
        cube : SyntheticCube    
        '''

        stage = 'SyntheticCube: extinction'

        # read own extinction law
        if input_opacities is None:
            t = np.loadtxt(ROOT + 'database/extinction/extinction_law.txt')

        else:
            t = np.loadtxt(input_opacities)
            
        wav_ext = t[:, 0]
        k_lam = t[:, 1]
        
        # wav_ext monotonically increasing
        if wav_ext[0] > wav_ext[1]:
            wav_ext = wav_ext[::-1]
            k_lam = k_lam[::-1]
        
        k_v = np.interp(0.550, wav_ext, k_lam)

        # interpolate to get A_int for a certain wavelength
        k = np.interp(self.wav, wav_ext, k_lam)
        A_int_lam = A_v * (k / k_v)

        # apply extinction law
        val_ext = np.zeros(shape=np.shape(self.val))
        val_ext[:,:,:len(self.wav)] = self.val[:,:,:len(self.wav)] * 10**(-0.4 * A_int_lam[:len(self.wav)])

        # return SimulateCube
        c = SyntheticCube(self)
        c.val = val_ext
        c.stage = stage
        c.log.append(c.stage)

        return c

    def change_resolution(self, new_resolution, grid_plot=None):
        '''
        Changes the resolution of every slice of the val cube.
        
        Parameters
        ----------
        
        new_resolution : Resolution which the val array should get in ``arcsec/pixel.``
                         
        grid_plot : ``None``, ``True``
            If ``True`` old and new resolution is visualized in a plot.
            Default is ``None``.
                
        Returns
        -------
        
        cube : SyntheticCube    
        '''

        stage = 'SyntheticCube: change_resolution'

        # debugging comment
        logger.debug('-' * 70)
        logger.debug(stage)
        logger.debug('-' * 70)
        logger.debug('total value before zoom : ' + str('%1.4e' % np.sum(self.val)) + ' ' + str(self.units))

        # match resolution of psf and val slice
        f = ConservingZoom(array=self.val, initial_resolution=self.resolution['arcsec'], new_resolution=new_resolution)
        zoomed_val = f.zoom()

        # average after changing resolution for MJy/sr
        if self.units == 'MJy/sr' or self.units == 'Jy/arcsec^2':
            # size of new pixel in units of old pixel
            size = new_resolution ** 2 / self.resolution['arcsec'] ** 2
            zoomed_val = zoomed_val / size

        if grid_plot is not None:
            f.zoom_grid(self.name)

        # debugging comment
        logger.debug('total value after zoom  : ' + str('%1.4e' % np.sum(zoomed_val)) + ' ' + str(self.units))

        # return SimulateCube
        c = SyntheticCube(self)
        c.val = zoomed_val
        c.stage = stage
        c.log.append(c.stage)
        c.FOV = (f.len_nx / f.len_nrx * self.FOV[0], f.len_ny / f.len_nry * self.FOV[1])

        return c
        
        
    def central_pixel(self, dx, dy):
        '''
        Move array right and up to create a central pixel.
        
        Returns
        -------
        
        cube : SyntheticCube   
        '''
        
        stage = 'SyntheticCube: central_pixel'
        
        ce = central(array=self.val, dx=dx, dy=dy)
        
        len_x_old = float(self.pixel[0])
        len_x_new = float(len(ce[:,0]))
        len_y_old = float(self.pixel[1])
        len_y_new = float(len(ce[0,:]))
        old_FOV = self.FOV
        new_FOV = (len_x_new / len_x_old * old_FOV[0], len_y_new / len_y_old * old_FOV[1])        
        
        # return SimulateCube
        c = SyntheticCube(self)
        c.val = ce
        c.stage = stage
        c.log.append(c.stage)
        c.FOV = new_FOV
        
        return c
    

    def convolve_psf(self, psf):
        '''
        Convolves every slice of the val cube with a PSF of choice.
        
        Parameters
        ----------
        
        psf : GaussianPSF, FilePSF, database, FunctionPSF
        
            * GaussianPSF(self, diameter): Convolves val with Gaussian PSF.
              
            * FilePSF(self, psf_file, condensed): Reads PSF from input file.
                  
            * database: PSF object defined in FluxCompensator database. 

            * FunctionPSF(self, psf_function, width): Convolves val with calculated PSF.
            
        
        Returns
        -------
        
        cube : SyntheticCube       
        '''

        stage = 'SyntheticCube: convolve_PSF'

        # debugging comments
        if isinstance(psf, GaussianPSF):
            logger.debug('-' * 70)
            logger.debug(stage + 'with GaussianPSF')
            logger.debug('-' * 70)

        # convolve value with classes GaussianPSF, FilePSF and FunctionPSF
        val = self.val.copy()
        for i in range(len(self.wav)):
            val[:, :, i] = psf.convolve(wav = self.wav[i], array = self.val[:,:, i], resolution = self.resolution)

        # return SimulateCube
        c = SyntheticCube(self)
        c.val = val
        c.stage = stage
        c.log.append(c.stage)

        return c

    def convolve_filter(self, filter_input, plot_rebin=None, plot_rebin_dpi=None):
        '''
        Convolves slice within filter limits into a 2D image.
        
        Parameters
        ----------
        
        filter_input : object
        
            * database : if filter ``name`` from FluxCompensator database is used.            
            * Filter : if own filter is used.
            
        plot_rebin : ``True``, ``None``
            Switch to plot the rebined filter and the original filter in one plot.
            
        plot_rebin_dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the value savefig.dpi 
            in the matplotlibrc file.
        
        Returns
        -------
        
        image : SyntheticImage
        '''

        stage = 'SyntheticCube: convolve_filter'

        # debugging comment
        logger.debug('-' * 70)
        logger.debug(stage)
        logger.debug('-' * 70)

        weight = filter_input.rebin(self.wav, self.val)

        # returns weight{'wav_short' 'val_short' 'Response_new' 'filter_index' 'wavf_0' 'waf_min' 'waf_max' 'filter_name'}
        wav_short = weight['wav_short']
        val_short = weight['val_short']
        filter_index = weight['filter_index']
        Response_new = weight['Response_new']
        waf_0 = weight['waf_0']
        waf_min = weight['waf_min']
        waf_max = weight['waf_max']
        filter_name = weight['filter_name']

        if plot_rebin is not None:
            plot = filter_input.plot(val_name=self.name, dpi=plot_rebin_dpi)

        # weight val_short with rebined response
        val = val_short.copy()
        val[:, :, :len(wav_short)] = val_short[:,:, :len(wav_short)] * Response_new[:len(wav_short)]

        # collapse remaining cube into 2D
        val_2D = np.sum(val, axis=2)

        # return SyntheticImage
        from .image import SyntheticImage
        i = SyntheticImage(self)
        i.log.append(stage)
        i.stage = 'SyntheticImage: initial'
        i.log.append(i.stage)
        i.filter = {'name': filter_name, 'waf_0': waf_0, 'waf_min': waf_min, 'waf_max': waf_max}
        i.val = val_2D
        i.wav = np.array(waf_0)

        return i

    def add_noise(self, mu_noise, sigma_noise, seed=None, diagnostics=None):
        '''
        Adds normal distributed noise to every slice in the val cube 
        of SyntheticCube.
        
        Parameters
        ----------
        
        mu_noise : float
            Mean of the normal distribution.
            Good choice: mu_noise = 0.
            
        sigma_noise : float
            Standard deviation of the normal distribution. Good choice around:
                * ``'ergs/cm^2/s'``     : sigma_noise = 10.**(-13)
                * ``'ergs/cm^2/s/Hz'``  : sigma_noise = 10.**(-26)
                * ``'Jy'``              : sigma_noise = 10.**(-3)
                * ``'mJy'``             : sigma_noise = 10.**(-1)
                * ``'MJy/sr'``          : sigma_noise = 10.**(-10)
                
        seed : float, ``None``
            When float seed fixes the random numbers to a certain
            sequence in order to create reproducible results.
            Default is ``None``.
        
        diagnostics : truetype
            When ``True`` noise array is stored in a fits file.
        
        Returns
        -------
        
        cube : SyntheticCube       
        '''

        stage = 'SyntheticCube: add_noise'

        # add different noise with same mu and sigma to 3D cube
        val = self.val.copy()
        for i in range(len(self.wav)):
            if sigma_noise != 0. and sigma_noise != 0:
                if seed is not None:
                    np.random.seed(seed=seed)
                noise = normal(mu_noise, sigma_noise, self.pixel)
            if sigma_noise == 0. or sigma_noise == 0:
                noise = np.zeros(self.pixel)
            val[:, :, i] = self.val[:,:, i] + noise        
        
        if diagnostics is True:
            fits.writeto(self.name + '_process-output_SC-noise.fits', noise, clobber=True)

        # return SyntheticCube
        c = SyntheticCube(self)
        c.val = val
        c.stage = stage
        c.log.append(c.stage)

        return c

    def get_rough_sed(self):
        '''
        Collapses the current val cube into 1D array (SED).
        
        Returns
        -------
        
        sed : SyntheticSED
        '''

        stage = 'SyntheticCube: get_rough_sed'

        # for MJy/sr convert first, add and then convert back
        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=self.val)
            self.val = s.get_unit(in_units=self.units, out_units='Jy', input_resolution=self.resolution['arcsec'])

        # collapse every slice to one scalar value
        rough_sed = np.sum(np.sum(self.val.copy(), axis=1), axis=0)

        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=rough_sed)
            rough_sed = s.get_unit(in_units='Jy', out_units=self.unit_out, input_resolution=self.resolution['arcsec'] * self.pixel[0])

        # return SyntheticSED
        from .sed import SyntheticSED
        s = SyntheticSED(self)
        s.log.append(stage)
        s.stage = 'SyntheticSED: initial'
        s.log.append(s.stage)
        s.val = rough_sed

        return s

    def get_total_val(self, wav_1, wav_2):
        '''
        Collapses the val of SyntheticCube within the boundaries wav_1
        and wav_2 into a 0D value val. 
        
        WARNING: This tool cannot replace convolve_filter! 
                 But it can be used to produce rough estimates 
                 in-between the processes.
                 
        Parameters
        ----------
        
        wav_1, wav_2 : float
            Boundaries in microns.
            
        
        Returns
        -------
        
        val : SyntheticFlux
        '''
        stage = 'SyntheticCube: get_total_val'

        # slices within boundaries are extracted, averaged collapsed to an 2D image and finally collpased to a single scalar value
        # for MJy/sr convert first, add and then convert back
        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=self.val)
            val = s.get_unit(in_units=self.units, out_units='Jy', input_resolution=self.resolution['arcsec'])
        else: val = self.val

        c = get_slices(wav=self.wav, val=val, wav_1=wav_1, wav_2=wav_2)
        i = average_collapse(val=c['val_short'])
        f_total = np.sum(i)

        # real limits within collaps
        wav_max = 10 ** (np.log10(self.wav[c['filter_index'][0]]) + self.spacing_wav / 2.)
        wav_min = 10 ** (np.log10(self.wav[c['filter_index'][-1]]) - self.spacing_wav / 2.)
        wav_total = central_wav(wav=[wav_min, wav_max])

        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=wav_total, val=f_total)
            f_total = s.get_unit(in_units='Jy', out_units=self.unit_out, input_resolution=self.resolution['arcsec'] * self.pixel[0])

        # return SyntheticFlux
        from .flux import SyntheticFlux
        f = SyntheticFlux(self)
        f.log.append(stage)
        f.stage = 'SyntheticFlux: initial'
        f.log.append(f.stage)
        f.wav = np.array(wav_total)
        f.val = np.array(f_total)
        f.filter = {'name': 'val_tot', 'waf_0': wav_total, 'waf_min': wav_min, 'waf_max': wav_max}

        return f

    def plot_image(self, wav_interest, prefix=None, name=None, multi_cut=None, single_cut=None, set_cut=None, dpi=None):
        '''
        Plots a certain slice close the wav_interest. 
        The wavelength interval of the chosen slice labels the plot.
        
        
        Parameters
        ----------
    
        wav_interest : float, ``None``
        
            * float : wavelength close to slice in microns. 
            * ``None`` : Only if input_array is SyntheticImage like
            
        prefix : str
            Name of the image. Default naming chain is switched off.
        
        name : str
            Name of image within the default naming chain to distinguish the
            plot files. E. g. 'PSF_gaussian'
        
        mulit_cut : ``True``, ``None``
        
            * ``True`` : plots chosen image slice at cuts of [100, 99, 95, 90]%.
            * ``None`` : no mulit-plot is returned.
    
            Default is ``None``.  
        
        single_cut : float, ``None``
        
            * float : cut level for single plot of image slice between 0 and 100.
            * ``None`` : no single plot is returned.
        
        set_cut : tuple, ``None``
        
            * tuple : set_cut(v_min, v_max)
                      Minimal and maximal physical value of val in the colorbars.
            * ``None`` : no plot with minimal and maximal cut is returned.
        
            Default is ``None``.
        
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the value valig.dpi 
            in the matplotlibrc file.

        
        Returns
        -------
        
        cube : SyntheticCube
        '''

        stage = 'SyntheticCube: plot_image'

        if prefix is None and name is None:
            raise Exception('If prefix name is not given, you need to give the a name to enable the default naming chain.')
        
        if prefix is not None:
            if multi_cut is True and (single_cut is not None or set_cut is not None):
                raise Exception('If prefix naming is enabled only one plotting option can be chosen.')
            elif multi_cut is None and (single_cut is not None and set_cut is not None):
                raise Exception('If prefix naming is enabled only one plotting option can be chosen.')
        
        plot = MakePlots(prefix=prefix, name=name, input_array=SyntheticCube(self), wav_interest=wav_interest, multi_cut=multi_cut, single_cut=single_cut, set_cut=set_cut, dpi=dpi)

        # return SyntheticCube
        c = SyntheticCube(self)
        c.stage = stage
        c.log.append(c.stage)

        return c

    @property
    def spacing_wav(self):
        '''
        The property spacing_wav estimates the width of the logarithmic
        spaced wav entries.
        '''

        if self.wav.ndim != 0:
            spacing_wav = np.log10(self.wav[0] / self.wav[-1]) / (len(self.wav) - 1)
        else:
            spacing_wav = None
        return spacing_wav

    @property
    def pixel(self):
        '''
        The property pixel is a tuple which resembles the current pixel in a
        value val. ``pixel(x,y)`` are calls as follows:
         
        ``x = pixel[0]``
        ``y = pixel[1]``
        '''
        if self.val.ndim in (0, 1):
            pixel = (None, None)
        if self.val.ndim in (2, 3):
            pixel = (self.val.shape[0], self.val.shape[1])
        return pixel

    @property
    def shape(self):
        '''
        The property shape is a string, which resembles the current shape of
        the value val. 
        
        scalar: ``'()'`` 
        1D:     ``'(wav)'`` 
        2D:     ``'(x, y)'`` 
        3D:     ``'(x, y , wav)'`` 
        '''
        if self.val.ndim == 0:
            shape = '()'
        if self.val.ndim == 1:
            shape = '(wav)'
        if self.val.ndim == 2:
            shape = '(x, y)'
        if self.val.ndim == 3:
            shape = '(x, y, wav)'
        return shape

    @property
    def resolution(self):
        '''
        The property resolution tells you the current resolution. If we are already 
        in the SED or flux everything is considered as one large pixel.

            resolution in arcsec per pixel : ``resolution['arcsec']``
            resolution in rad per pixel    : ``resolution['rad']``
        
        '''
        resolution = {}
        if self.pixel[0] is None:
            resolution['rad'] = self.FOV[0] / 1. / self.distance
        else:
            resolution['rad'] = self.FOV[0] / self.pixel[0] / self.distance
        resolution['arcsec'] = np.degrees(resolution['rad']) * 3600
        return resolution
