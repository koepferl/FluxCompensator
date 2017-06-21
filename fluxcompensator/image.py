from copy import deepcopy
import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/'

import numpy as np
from numpy.random import normal
from astropy import log as logger
from astropy.io import fits
from astropy.wcs import WCS

from .psf import GaussianPSF, FilePSF, FunctionPSF

from .utils.plot import MakePlots
from .utils.resolution import ConservingZoom, central
from .utils.tools import properties, grid_units, get_slices, average_collapse, central_wav
from .utils.units import ConvertUnits


class SyntheticImage(object):

    '''
    SyntheticImage is part the FluxCompensator. It converts 
    input_arrays (e. g. HYPERION ModelOutput in 2D) to "realistic" synthetic observations 
    (e. g. by accounting for PSF and noise). 
    It contains attributes like ModelOutput (see Notes).
    If input_array is already a SyntheticImage object, the attributes are
    passed. If input_array is not a SyntheticImage object, SyntheticImage
    specific attributes are defined and then passed. 
    
    
    Parameters
    ----------
    input_array : SyntheticImage, ModelOutput, optional
        input_array also reads arrays with ModelOutput like properties.
        
    unit_out : str, optional
        The output units for SyntheticImage val. Valid options are:

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
        The wavelength of the val image in microns.
    
    val : numpy.ndarray
        The 2D image with shape (x, y).
        
    units : str
        Current units of the val image.
        
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
        Tuple ``FOV(x,y)`` of Field of View pixel entries:
        
            * pixel in x direction: ``FOV[0]``
            * pixel in y direction: ``FOV[1]``
    
    name : str
        The name of the FluxCompensator object until another
        input_array is called. The default is ``None``.
    
    stage : str
        Gives current operation stage of SyntheticImage.
        E. g. ``'SyntheticImage: convolve_PSF'``
    
    log : list
        List of strings of the previous and current stages.
    
    filter : dict
        Dictionary filter = ``{name, waf_0, waf_min, waf_max}`` 
        of the applied filter: 
        
            * name of filter:     ``filter['name']``
            * central wavelength: ``filter['waf_0']``
            * minimal wavelength:  ``filter['waf_min']``
            * maximal wavelength: ``filter['waf_max']``
    
    
    Returns
    -------
    
    image : SyntheticImage
        2D val array with SyntheticImage properties.
                
    flux : SyntheticFlux
        0D val array (scalar) with SyntheticFlux properties.
    '''

    def __init__(self, input_array, unit_out='ergs/cm^2/s', name=None):

        # Hyperion ModelOutput attributes
        #print input_array.val.ndim, input_array.val.shape[2]
        #if input_array.val.ndim == 3 and input_array.val.shape[2] == 1:
            #self.val = np.array(deepcopy(input_array.val[:,:,0]))
        #if input_array.val.ndim == 2:
        self.val = np.array(deepcopy(input_array.val))
        #else:
        #    raise Exception('input_array does not have the right dimensions. numpy array of (x, y) or (x, y, 1) is required.')

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
        from .cube import SyntheticCube
        if isinstance(input_array, SyntheticImage) or isinstance(input_array, SyntheticCube):
        # attributes with are passed, since input_array is SyntheticCube or SyntheticImage

            # physical values
            self.unit_in = input_array.unit_in
            self.unit_out = input_array.unit_out

            self.grid_unit = input_array.grid_unit
            self.grid_unit_name = input_array.grid_unit_name

            # properties of image
            self.FOV = deepcopy(input_array.FOV)

            # name
            self.name = input_array.name
            self.stage = input_array.stage
            self.log = deepcopy(input_array.log)

            # filter
            self.filter = deepcopy(input_array.filter)

        else:  # attributes are defined, since input_array is NOT SyntheticCube or Image
            # physical values
            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = grid_units(self.x_max - self.x_min)['grid_unit']
            self.grid_unit_name = grid_units(self.x_max - self.x_min)['grid_unit_name']

            self.FOV = (self.x_max - self.x_min, self.y_max - self.y_min)

            # name
            self.name = name
            self.stage = 'SyntheticImage:  initial'
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
        
        image : SyntheticImage    
        '''

        stage = 'SyntheticImage: extinction'

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
        val_ext = self.val * 10 ** (-0.4 * A_int_lam)
        
        # return SimulateImage
        i = SyntheticImage(self)
        i.val = val_ext
        i.stage = stage
        i.log.append(i.stage)

        return i

    def change_resolution(self, new_resolution, grid_plot=None):
        '''
        Changes the resolution of val image.
        
        Parameters
        ----------
        
        new_resolution : Resolution which the val array should get in
                         ``arcsec/pixel.``
                         
        grid_plot : ``None``, ``True``
            If ``True`` old and new resolution is visualized in a plot.
            Default is ``None``.
                
        Returns
        -------
        
        image : SyntheticImage   
        '''

        stage = 'SyntheticImage: change_resolution'

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
        logger.debug('total val after zoom  : ' + str('%1.4e' % np.sum(zoomed_val)) + ' ' + str(self.units))

        # return SimulateCube
        i = SyntheticImage(self)
        i.val = zoomed_val
        i.stage = stage
        i.log.append(i.stage)
        i.FOV = (f.len_nx / f.len_nrx * self.FOV[0], f.len_ny / f.len_nry * self.FOV[1])

        return i


    def central_pixel(self, dx, dy):
        '''
        Move array right and up to create a central pixel.
        
        Returns
        -------
        
        image : SyntheticImage   
        '''
        
        stage = 'SyntheticImage: central_pixel'
        
        # match resolution of psf and val slice
        ce = central(array=self.val, dx=dx, dy=dy)
        
        len_x_old = float(self.pixel[0])
        len_x_new = float(len(ce[:,0]))
        len_y_old = float(self.pixel[1])
        len_y_new = float(len(ce[0,:]))
        old_FOV = self.FOV
        new_FOV = (len_x_new / len_x_old * old_FOV[0], len_y_new / len_y_old * old_FOV[1])        
        
        # return SimulateCube
        i = SyntheticImage(self)
        i.val = ce
        i.stage = stage
        i.log.append(i.stage)
        i.FOV = new_FOV

        return i


    def convolve_psf(self, psf):
        '''
        Convolves the val image with a PSF of choice.
        
        Parameters
        ----------
        
        psf : GaussianPSF, FilePSF, database, FunctionPSF
        
            * GaussianPSF(self, diameter): Convolves val with Gaussian PSF.
              
            * FilePSF(self, psf_file, condensed) : Reads PSF from input file.
                  
            * database : object
                If PSF ``name_PSF`` from FluxCompensator database is used. 
              
            * FunctionPSF(self, psf_function, width): Convolves defined PSF.
                    2D val image of SyntheticImage.val convolved with PSF.
        
        Returns
        -------
        
        image : SyntheticImage       
        '''

        stage = 'SyntheticImage: convolve_PSF'

        # debugging comments
        if isinstance(psf, GaussianPSF):
            logger.debug('-' * 70)
            logger.debug(stage + 'with GaussianPSF')
            logger.debug('-' * 70)

        # convolve val with classes GaussianPSF, FilePSF and FunctionPSF
        val = psf.convolve(wav=self.wav, array=self.val, resolution=self.resolution)

        # return SyntheticImage
        i = SyntheticImage(self)
        i.stage = stage
        i.log.append(i.stage)
        i.val = np.array(val)

        return i

    def add_noise(self, mu_noise, sigma_noise, seed=None, diagnostics=None):
        '''
        Adds normal distributed noise to the val image of SyntheticImage.
        
        Parameters
        ----------
        
        mu_noise : float
            Mean of the normal distribution.
            Good choice: mu_noise = 0.
            
        sigma_noise : float
            Standard deviation of the normal distribution.
            
            Good choice arround:
                * ``'ergs/cm^2/s'``     : sigma_noise = 10.**(-13)
                * ``'ergs/cm^2/s/Hz'``  : sigma_noise = 10.**(-26)
                * ``'Jy'``              : sigma_noise = 10.**(-3)
                * ``'mJy'``             : sigma_noise = 10.**(-1)
                * ``'MJy/sr'``          : sigma_noise = 10.**(-10)
        
        seed : float, ``None``
            When float seed fixes the random numbers to a certain sequence in order to create reproducible results.
            Default is ``None``.
        
        diagnostics : truetype
            When ``True`` noise array is stored in a fits file.
        
        
        Returns
        -------
        
        image : SyntheticImage      
        '''

        stage = 'SyntheticImage: add_noise'
        
        if sigma_noise != 0. and sigma_noise != 0:
            if seed is not None:
                np.random.seed(seed=seed)
            noise = normal(mu_noise, sigma_noise, self.pixel)
        if sigma_noise == 0. or sigma_noise == 0:
            noise = np.zeros(self.pixel)

        # Get noise.fits file
        if diagnostics is True:
            fits.writeto(self.name + '_' + 'process-output_SI-noise.fits', noise, clobber=True)

        # add noise if val is already collapsed (x, y)
        val = self.val.copy() + noise

        # return SyntheticImage
        i = SyntheticImage(self)
        i.stage = stage
        i.log.append(i.stage)
        i.val = np.array(val)

        return i

    def get_total_val(self):
        '''
        Collapses the val image of SyntheticImage into a 0D val array. 
        
        Returns
        -------
        
        flux : SyntheticFlux
        '''

        stage = 'SyntheticImage: get_total_val'

        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=self.val)
            val = s.get_unit(in_units=self.units, out_units='Jy', input_resolution=self.resolution['arcsec'])
        else: val = self.val

        # collapse 2D image to a single scalar val
        total_val = np.sum(val)

        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=total_val)
            total_val = s.get_unit(in_units='Jy', out_units=self.unit_out, input_resolution=self.resolution['arcsec'] * self.pixel[0])

        # return SyntheticFlux
        from .flux import SyntheticFlux
        f = SyntheticFlux(self)
        f.log.append(stage)
        f.stage = 'SyntheticFlux: initial'
        f.log.append(f.stage)
        f.val = np.array(total_val)

        return f

    def plot_image(self, prefix=None, name=None, multi_cut=None, single_cut=None, set_cut=None, dpi=None):
        '''
        Plots the val image of SyntheticImage. The wavelength interval
        around the central wavelength labels the plot.
        
        
        Parameters
        ----------
        
        prefix : str
            Name of the image. Default naming chain is switched off.
        
        name : str
            Name of image within the default naming chain to distinguish the
            plot files. E. g. 'PSF_gaussian'
        
        mulit_cut : ``True``, ``None``
        
            * ``True`` : plots chosen image slice at cuts of [100, 99, 95, 90]%.
            * ``None`` : no mulit-plot is returned.
    
            Default is ``None``.  
        
        single_cut : float [0,100], ``None``
        
            * float : cut level for single plot of image slice.
            * ``None`` : no single plot is returned.
        
        set_cut : tuple, ``None``
        
            * tuple : set_cut(v_min, v_max)
                      Minimal and maximal physical val presented in the colorbars.
            * ``None`` : no plot with minimal and maximal cut is returned.
        
            Default is ``None``.
        
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
        
        
        Returns
        -------
        
        image : SyntheticImage
        '''

        stage = 'SyntheticImage: plot_image'

        if prefix is None and name is None:
            raise Exception('If prefix name is not given, you need to give the a name to enable the default naming chain.')
        
        if prefix is not None:
            if multi_cut is True and (single_cut is not None or set_cut is not None):
                raise Exception('If prefix naming is enabled only one plotting option can be chosen.')
            elif multi_cut is None and (single_cut is not None and set_cut is not None):
                raise Exception('If prefix naming is enabled only one plotting option can be chosen.')

        plot = MakePlots(prefix=prefix, name=name, input_array=SyntheticImage(self), multi_cut=multi_cut, single_cut=single_cut, set_cut=set_cut, dpi=dpi)

        # return SyntheticImage
        i = SyntheticImage(self)
        i.stage = stage
        i.log.append(i.stage)

        return i

    def add_to_observation(self, fits_file, name, position_pix=None, position_world=None, zero_edges=None):
        '''
        Blends the modeled realistic synthetic observation to a real observation in a fits file.
        
        Parameters
        ----------
        
        fits_file : str
            fits_file of the observation.
        
        name : str
            Name of the output fits file.
        
        position_pix : list, ``None``
            Center position of the model in observation pixel coordinates.
            Default is ``None``.  
        
        position_world : list, ``None``
            Center position of the model in observation world coordinates.
            Default is ``None``.  
        
        zero_edges : ``True``, ``None``
            If ``True`` edges of model are normalized to zero.
            Default is ``None``.  
        
        
        Returns
        -------
        
        image : SyntheticImage
        
        '''

        stage = 'SyntheticImage: add_to_observation'

        # world coordinates from fits_file
        w = WCS(fits_file)

        if position_world is None and position_pix is None:
            raise Exception('WARNING: Position of model center needs to be given either in world or pixel coordinates.')

        if position_pix is not None:
            pos = position_pix
            p_x_pos, p_y_pos = pos[0], pos[1]
        else:
            pos = position_world
            p_x_pos, p_y_pos = w.wcs_world2pix(pos[0], pos[1], 1)

        # center position in pixel and adjust position in current grid
        x_round = np.round(p_x_pos, 0)
        x_int = int(p_x_pos)
        y_round = np.round(p_y_pos, 0)
        y_int = int(p_y_pos)

        # even or odd
        if len(self.val[0]) % 2 == 0 and len(self.val[1]) % 2 == 0:
            pos = np.array([x_round, y_round])
        else:
            if x_int == int(x_round):
                if y_int == int(y_round):
                    pos = np.array([x_round + 0.5, y_round + 0.5])
                else:
                    pos = np.array([x_round + 0.5, y_round - 0.5])
            else:
                if y_int == int(y_round):
                    pos = np.array([x_round - 0.5, y_round + 0.5])
                else:
                    pos = np.array([x_round - 0.5, y_round - 0.5])

        # limits of model in observation
        start_x = pos[0] - len(self.val[0]) / 2.
        stop_x = pos[0] + len(self.val[0]) / 2.
        start_y = pos[1] - len(self.val[1]) / 2.
        stop_y = pos[1] + len(self.val[1]) / 2.

        # normalized that edges are zero
        if zero_edges is True:
            model = self.val.copy() - np.min(self.val)
        else:
            model = self.val.copy()

        # open fits_file
        hdulist = fits.open(fits_file)
        hdu = hdulist[0]
        header = hdu.header
        
        if np.allclose(np.abs(header['CDELT1'] * 3600), self.resolution['arcsec']) is not True:
            raise Exception('WARNING: make sure that resolution of observation and model are the same! E. g. change resolution of FC_object first.')

        image = hdu.data

        # add model to observation
        image[start_y:stop_y, start_x:stop_x] = image[start_y:stop_y, start_x:stop_x] + model

        # store to name.fits file
        fits.writeto(name + '.fits', image, clobber=True)

        # return SyntheticImage
        i = SyntheticImage(self)
        i.stage = stage
        i.log.append(i.stage)

        return i

    def add_field_stars(self, extinction_map, database=None, star_file=None, seed=None, ISMextinction=None):
        '''
        Adds field stars to synthetic image.
        
        Parameters
        ----------
        
        extinction_map : object
            Created with ``fluxcompensator.utils.fieldstars.extract_extinction_map``.
        
        database : dict, ``None``
            Dictionary sets the parameters for field stars loaded for the respective
            band from the built-in database.
            
            dict = {'number':200, 'distance_range':[3*kpc, 50*kpc], 'ground': 0.02}
        
            The dictionary is structured as follows:
        
            * ``'number'``         : int in [0,288]
            * ``'distance_range'`` : list
        
              Distance lower and upper limit in units of cm
        
            * ``'ground'``         : str, float
        
              Distribution of stars before (``'foreground'``) or behind (``'background'``) the synthetic object. 
              When ``'ground'`` is a ``float`` in the limits of [0,1] then this is the fraction of foreground stars.
        
            Default is ``None``.  
        
        star_file : str, ``None``
            To load individual file with field stars in the format of (distance[pc], mag[band]).
            Default is ``None``.  
        
        seed : int, ``None``
            To create reproducible results for the positions of field stars.
            Default is ``None``.  
        
        ISMextinction : float, ``None``
            Optical extinction A_V along the line of sight in units mag/kpc.
            Default is ``None``.  
        
        
        Returns
        -------
        
        image : SyntheticImage
        
        '''

        stage = 'SyntheticImage: add_field_stars'

        # make sure resolution and PSF was not applied before
        if 'SyntheticImage: convolve_PSF' in self.log or 'SyntheticCube: convolve_PSF' in self.log \
                                                      or 'SyntheticImage: change_resolution' in self.log \
                                                      or 'SyntheticCube: change_resolution' in self.log:
            raise Exception('WARNING: Adding field stars should happen before changing resolution or convolution with PSF.')
        
        # make sure that filter was applied before
        if 'SyntheticCube: convolve_filter' not in self.log:
            raise Exception('WARNING: Image must be convolved with the transmission of a detector.')

        if extinction_map.shape != self.val.shape:
            raise Exception('WARNING: Extinction map and val of SyntheticImage do not have the same dimension.')
        
        # load file or give parameters to read from database
        if database is None and star_file is None:
            raise Exception('WARNING: Either database or star_file need to be different from None.')
        
        # read from database
        if database is not None:
            from utils.fieldstars import get_stars_from_database
            mag, star_distance = get_stars_from_database(band=self.filter['name'], number=database['number'], 
                                                          distance_range=database['distance_range'], ground=database['ground'], 
                                                          object_distance=self.distance, seed=seed)

        # read field star data from file, distance in pc, mag in magnitudes
        pc = 3.08568025e+18
        if database is None and star_file is not None:
            print 'CAUTION: only stars in the same band as the image should be loaded.'
            print 'CAUTION: units of distance is in pc, stellar photometry in mag.'
            f = np.loadtxt(star_file)
            star_distance = f[:,0] * pc # pc>cm
            mag = f[:,1]
        
        # ensure that random numbers are the same every time
        if seed is not None:
            np.random.seed(seed)
        
        #print star_distance - self.distance
        # extinction from extinction map for objects
        x = np.random.uniform(0, self.pixel[0], len(mag)).astype('int')
        y = np.random.uniform(0, self.pixel[1], len(mag)).astype('int')
        A_obj = extinction_map[x,y]
        
        # convert to from A_v to A_filter
        print 'CAUTION: Extinction law from Kim et al. is used.'
        wav_ext, k_lam = np.loadtxt(ROOT + 'database/extinction/extinction_law.txt', unpack=True)
        k_v = np.interp(0.550, wav_ext, k_lam)
        k = np.interp(self.wav, wav_ext, k_lam)
        A_filter = A_obj * (k / k_v)
        
        MAG_ext = np.where([star_distance[i] >= self.distance for i in range(len(star_distance))], mag + A_filter, mag)
        if ISMextinction is not None:
            ISM_extinction_filter = ISMextinction * (k / k_v)
            MAG_ext_ISM = np.where([star_distance[i] >= self.distance for i in range(len(star_distance))], MAG_ext + ISM_extinction_filter * star_distance/(1e3 * pc), MAG_ext)
            #print mag[10], MAG_ext[10], MAG_ext_ISM[10]
            MAG_ext = MAG_ext_ISM
        
        # zero-point
        import database.missions as filters
        zero_point = getattr(filters, self.filter['name'] + '_ZERO')
        wav_1D = np.ones(np.shape(MAG_ext))*self.wav

        # converting mag to flux
        flux = ConvertUnits(wav=wav_1D, val=MAG_ext)
        if self.units == 'MJy/sr' or self.units == 'Jy/arcsec^2':
            starflux = flux.get_unit(in_units='mag', out_units=self.units, zero_point=zero_point, input_resolution=self.resolution['arcsec'])
        else:
            starflux = flux.get_unit(in_units='mag', out_units=self.units, zero_point=zero_point)

        # position of star on image
        add_stellar_flux = self.val.copy()
        for i in range(len(starflux)):
            add_stellar_flux[x[i],y[i]] = add_stellar_flux[x[i],y[i]] + starflux[i]

        # return SyntheticImage
        i = SyntheticImage(self)
        i.val = add_stellar_flux
        i.stage = stage
        i.log.append(i.stage)

        return i

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
        val slice. ``pixel(x,y)`` are calls as follows:
         
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
        the val array. 
        
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
        in the SED or val dimension everything is considered as one large pixel.

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
