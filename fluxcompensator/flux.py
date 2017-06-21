from copy import deepcopy
import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/'

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from .utils.one_filter import OneFilter
from .utils.plot import MakePlots
from .utils.plot_filters import PlotFilters
from .utils.tools import properties, grid_units, get_slices, average_collapse, central_wav
from .utils.units import ConvertUnits


class SyntheticFlux(object):

    '''
    SyntheticFlux is part of FluxCompensator. It converts 
    input_arrays (e. g. HYPERION ModelOutput) to "realistic" synthetic observations. 
    It contains attributes like ModelOutput (see Notes).
    If input_array is already a SyntheticFlux object, the attributes are
    passed. If input_array is not a SyntheticFlux object, SyntheticFlux
    specific attributes are defined and then passed. 
    
    
    Parameters
    ----------
    
    input_array : SyntheticFlux, ModelOutput (get_sed, one line), optional
        input_array also reads arrays with ModelOutput like properties.
        
    unit_out : str, optional
        The output units for SyntheticFlux val. Valid options are:

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
        The wavelength of scalar val array in microns.
    
    val : numpy.ndarray
        The scalar val with array shape like wav.
        
    units : str
        Current units of the vector val.
        
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
        Minimal aperture.
    
    ap_max : float
        Maximal aperture.
    
    
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
        Gives current operation stage of SyntheticFlux.
        E. g. ``'SyntheticFlux: plot_sed_filter'``
        
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
    
    flux : SyntheticFlux
        0D val array (scalar) with SyntheticFlux Properties.
    '''

    def __init__(self, input_array, unit_out='ergs/cm^2/s', name=None):

        # Hyperion ModelOutput attributes (SED or image)
        self.val = np.array(deepcopy(input_array.val))
        self.wav = np.array(deepcopy(input_array.wav))
        self.units = input_array.units
        self.distance = input_array.distance

        # Hyperion Image
        try:
            self.x_max = input_array.x_max
            self.x_max = input_array.x_max
            self.x_min = input_array.x_min
            self.y_max = input_array.y_max
            self.y_min = input_array.y_min
            self.lon_min = input_array.lon_min
            self.lon_max = input_array.lon_max
            self.lat_min = input_array.lat_min
            self.lat_max = input_array.lat_max
            self.pix_area_sr = input_array.pix_area_sr
            # switch
            self.hyperion_cube = True

        except AttributeError:
            self.x_max = None
            self.x_max = None
            self.x_min = None
            self.y_max = None
            self.y_min = None
            self.lon_min = None
            self.lon_max = None
            self.lat_min = None
            self.lat_max = None
            self.pix_area_sr = None
            # switch
            self.hyperion_cube = None

        # Hyperion SED
        try:
            self.ap_min = input_array.ap_min
            self.ap_max = input_array.ap_max
            # switch
            self.hyperion_sed = True

        except AttributeError:
            self.ap_min = None
            self.ap_max = None
            # switch
            self.hyperion_sed = None

        ##################
        # new attributes #
        ##################

        from .cube import SyntheticCube
        from .image import SyntheticImage
        from .sed import SyntheticSED

        if isinstance(input_array, SyntheticFlux) or isinstance(input_array, SyntheticSED) or isinstance(input_array, SyntheticImage) or isinstance(input_array, SyntheticCube):
            # attributes with are passed, since input_array is SyntheticCube, SyntheticImage, SyntheticImage, SyntheticFlux

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

        elif not isinstance(input_array, SyntheticFlux) or not isinstance(input_array, SyntheticSED) or not isinstance(input_array, SyntheticImage) or not isinstance(input_array, SyntheticCube) and self.hyperion_sed is True:
            # attributes are defined, since input_array is NOT SyntheticCube, SyntheticImage, SyntheticImage, SyntheticFlux or HyperionOutput sed

            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = None
            self.grid_unit_name = None

            self.FOV = None

            # name
            self.name = name
            self.stage = 'SyntheticFlux:  initial'
            self.log = [self.stage]

            # filter
            self.filter = {'name': None, 'waf_0': None, 'waf_min': None, 'waf_max': None}

            # convert into val units into unit_out
            if self.unit_in == 'MJy/sr' or self.unit_out == 'MJy/sr'  or self.unit_in == 'Jy/arcsec^2' or self.unit_out == 'Jy/arcsec^2':
                raise Exception('WARNING: Input or Output units needs to differ from MJy/sr or Jy/arcsec^2 if Input_array is not SyntheticCube or SyntheticSED or HyperionCube')

            s = ConvertUnits(wav=self.wav, val=self.val)
            self.val = s.get_unit(in_units=self.unit_in, out_units=self.unit_out)
            self.units = self.unit_out

        else:   # attributes are defined, since input_array is NOT FluxCompensator object
            # physical values
            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = grid_units(self.x_max - self.x_min)['grid_unit']
            self.grid_unit_name = grid_units(self.x_max - self.x_min)['grid_unit_name']

            self.FOV = (self.x_max - self.x_min, self.y_max - self.y_min)

            # name
            self.name = name
            self.stage = 'SyntheticFlux:  initial'
            self.log = np.array(self.stage)

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
        
        flux : SyntheticFlux    
        '''

        stage = 'SyntheticFlux: extinction'

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
        
        # return SimulateFlux
        f = SyntheticFlux(self)
        f.val = val_ext
        f.stage = stage
        f.log.append(f.stage)
        
        return f
        
    def magnitudes(self, zero_point):
        '''
        Converts f.val to magnitudes.
        
        Parameters
        ----------
        
        zero_point : zero-magnitude flux in Jy.
            The zero-magnitude flux can be accessed from the database 
            if the filter objects are defined therein.         
                
        Returns
        -------
        
        flux : SyntheticFlux    
        '''

        stage = 'SyntheticFlux: magnitude'

        mag = ConvertUnits(wav=self.wav, val=self.val)
        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            val_mag = mag.get_unit(in_units=self.units, out_units='mag', zero_point=zero_point, input_resolution=self.resolution['arcsec'])
        else:
            val_mag = mag.get_unit(in_units=self.units, out_units='mag', zero_point=zero_point)

        # return SimulateFlux
        f = SyntheticFlux(self)
        f.val = val_mag
        f.stage = stage
        f.log.append(f.stage)

        return f

    def plot_sed_filter(self, wav_sed, val_sed, ymin=10. ** (-10), my_own_filter=None, dpi=None):
        '''
        Reads SED like val from SyntheticCube (from get_rough_sed) or SyntheticSED and plots 
        it with the current passed val of SyntheticFlux in a log-log diagram. 
        That way the quality of the filtered val can be checked. 
        If the used filter is part of the filter database the filter 
        curve will be plotted. Otherwise just the boundaries are given.
        
        Parameters
        ----------
        
        wav_sed : np.ndarray
            1D vector of wavelength from ModelOutput (get_sed) or
            SyntheticCube.rough_sed
            
        val_sed : np.ndarray
            1D vector with val entries from ModelOutput (get_sed) or
            SyntheticCube.rough_sed with the same dimensions like wav_sed.
            
        ymin : float
            Minimal vertical limit for SED plot. 
            Default is 10.**(-10) of the maximum val.
            
        my_own_filter  :  ``None``, OneFilter 
            OneFilter object which is passed from outside to 
            plot the filter curve.

        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
            
            
        Returns
        -------
        
        flux : SyntheticFlux
        '''

        stage = 'SyntheticFlux: plot_sed_filter'

        # plot sed and val within boundaries of sed
        fig = plt.figure()
        font = {'size': 8}

        mpl.rc('font', **font)
        fig.suptitle('output SED with filtered val', fontsize=12.)

        x_max = wav_sed[0]
        x_min = wav_sed[-1]
        y_max = max(val_sed) * 10
        y_min = max(val_sed) * ymin

        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        plt.subplots_adjust(wspace=0, hspace=0)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])

        #ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel('val unit ' + '[' + str(self.unit_out) + ']')
        plt.loglog(wav_sed, val_sed, 'o-r', label='SED', linewidth=2)

        # val from get_total_val
        if self.filter['name'] == 'val_tot':
            # filter
            text = 'val get_total_val'
            plt.loglog(self.filter['waf_0'], self.val, 'ok', label=text)
            ax2 = plt.subplot(gs[1])
            ax2.set_xlim([x_min, x_max])
            ax2.set_ylim([10. ** (-4), 10.])

            plt.loglog([self.filter['waf_min'], self.filter['waf_min'], self.filter['waf_max'], self.filter['waf_max']], [10. ** (-4), 1., 1., 10. ** (-4)], '--k', label='get_total_val boundaries')

        # no filter or get_total_val previousely
        elif self.filter['name'] is None:
            raise Exception('No filter or get_total_val applied previously.')

        # val from filter
        else:
            text = 'filtered val by ' + str(self.filter['name'])

            import fluxcompensator.database.missions as of

            # filter not defined in one_filter
            try:
                the_filter = getattr(of, self.filter['name'] + '_FILTER_PLOT')
            except AttributeError:

                # pass one_filter object from outside
                the_filter = my_own_filter
                if the_filter is None:
                    raise Exception('The filter you are using is not a permanent object of OneFilter. Pass the OneFilter object with my_own_filer.')

            # filter defined OneFilter
            if isinstance(the_filter, OneFilter):
                x = PlotFilters(style='loglog', normalized=True, unit='unit energy')
                plot_filter = x.collect_filters(self.filter['name'] + '_FILTER_PLOT')
                plt.loglog(self.filter['waf_0'], self.val, 'ok', label=text)

                ax2 = plt.subplot(gs[1])
                ax2.set_xlim([x_min, x_max])
                ax2.set_ylim([10. ** (-4), 10.])
                ax2.set_ylabel('filter response [unit energy]', color='k')
                x.plot(plot_filter, line='-k')
                ax2.set_xlabel(r'$\lambda$ [$\mu$m]')
            else:
                raise Exception('my_own_filter is not a member of OneFilter.')

        # legend merged
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, loc=3)

        fig.savefig(str(self.name) + '_' + 'process-output_SF-' + self.filter['name'] + '.png', dpi=dpi)

        # return SyntheticFlux
        f = SyntheticFlux(self)
        f.stage = stage
        f.log.append(f.stage)

        return f

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
        val slice. ``pixel(x,y)`` are calles as follows:
         
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
