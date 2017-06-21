from copy import deepcopy
import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/'

from astropy import log as logger
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from .filter import Filter

from .utils.one_filter import OneFilter
from .utils.plot import MakePlots
from .utils.plot_filters import PlotFilters
from .utils.tools import properties, grid_units, get_slices, average_collapse, central_wav
from .utils.units import ConvertUnits


class SyntheticSED(object):

    '''
    SyntheticSED is part of the FluxCompensator. It converts 
    input_arrays (e. g. HYPERION ModelOutput) to "realistic" synthetic observations. 
    It contains attributes like ModelOutput (see Notes).
    If input_array is already a SyntheticSED object, the attributes are
    passed. If input_array is not a SyntheticSED object, SyntheticSED
    specific attributes are defined and then passed. 
    
    
    Parameters
    ----------
    
    input_array : SyntheticSED, ModelOutput (get_sed), optional
        input_array also reads arrays with ModelOutput like properties.        
    
    unit_out : str, optional
        The output units for SyntheticSED val. Valid options are:

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
        The wavelengths of the val vector entries in microns.
    
    val : numpy.ndarray
        The 1D vector with shape like wav.
        
    units : str
        Current units of the val vector.
        
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
        Gives current operation stage of SyntheticSED.
        E. g. ``'SyntheticSED: convolve_filter'``
        
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
    
    sed : SyntheticSED
        1D val array (collapsed rough SED) with SyntheticSED properties.
        
    flux : SyntheticFlux
        0D val array (scalar) with SyntheticFlux properties.
    '''

    def __init__(self, input_array, unit_out='ergs/cm^2/s', name=None):

        # Hyperion ModelOutput attributes (image and sed)
        if input_array.val.ndim == 2 and input_array.val[:, 0].ndim == 1:
            self.val = np.array(input_array.val[0, :])
        else:
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

        if isinstance(input_array, SyntheticSED) or isinstance(input_array, SyntheticCube):
            # attributes with are passed, since input_array is SyntheticSED or SyntheticCube

            # physical values
            self.unit_in = input_array.unit_in
            self.unit_out = input_array.unit_out

            if self.x_max is not None:
                self.grid_unit = grid_units(self.x_max - self.x_min)['grid_unit']
                self.grid_unit_name = grid_units(self.x_max - self.x_min)['grid_unit_name']
            else:
                self.grid_unit = None
                self.grid_unit_name = None

            # properties of cube
            self.FOV = deepcopy(input_array.FOV)

            # name
            self.name = input_array.name
            self.stage = input_array.stage
            self.log = deepcopy(input_array.log)

            # filter
            self.filter = deepcopy(input_array.filter)

        elif not isinstance(input_array, SyntheticSED) and not isinstance(input_array, SyntheticCube) and self.hyperion_sed is True:
            # attributes with are defined, since input_array is NOT SyntheticSED or SyntheticCube but HyperionSED
            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = None
            self.grid_unit_name = None

            self.FOV = None

            # name
            self.name = name
            self.stage = 'SyntheticSED:  initial'
            self.log = [self.stage]

            # filter
            self.filter = {'name': None, 'waf_0': None, 'waf_min': None, 'waf_max': None}

            # convert into val units into unit_out
            if self.unit_in == 'MJy/sr' or self.unit_out == 'MJy/sr' or self.unit_in == 'Jy/arcsec^2' or self.unit_out == 'Jy/arcsec^2':
                raise Exception('WARNING: Input or Output units needs to differ from MJy/sr if Input_array is not SyntheticCube or SyntheticSED or HyperionCube')

            s = ConvertUnits(wav=self.wav, val=self.val)
            self.val = s.get_unit(in_units=self.unit_in, out_units=self.unit_out)

            self.units = self.unit_out

        else:   # attributes with are defined, since input_array is NOT SyntheticSED or SyntheticCube or HyperionOutput
            # physical values
            self.unit_in = input_array.units
            self.unit_out = unit_out

            self.grid_unit = grid_units(self.x_max)['grid_unit']
            self.grid_unit_name = grid_units(self.x_max)['grid_unit_name']

            self.FOV = (self.x_max - self.x_min, self.y_max - self.y_min)

            # name
            self.name = name
            self.stage = 'SyntheticSED:  initial'
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
        
        sed : SyntheticSED    
        '''

        stage = 'SyntheticSED: extinction'

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
        val_ext[:len(self.wav)] = self.val[:len(self.wav)] * 10 ** (-0.4 * A_int_lam[:len(self.wav)])

        # return SyntheticSED
        s = SyntheticSED(self)
        s.val = val_ext
        s.stage = stage
        s.log.append(s.stage)

        return s

    def convolve_filter(self, filter_input, plot_rebin=None, plot_rebin_dpi=None):
        '''
        Convolves vector val entries within filter limits in a 0D val array.
        
        Parameters
        ----------
        
        filter_input : object
        
            * database : if filter ``name`` from FluxCompensator database is used.            
            * Filter : if own filter is used.
            
        plot_rebin : ``True``, ``None``
            Switch to plot the rebined filter and the original filter in one plot.
            
        plot_rebin_dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
        

        Returns
        -------
        
        flux : SyntheticFlux
        '''

        stage = 'SyntheticSED: convolve_filter'

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
        val[:len(wav_short)] = val_short[:len(wav_short)] * Response_new[:len(wav_short)]

        # collapse remaining vector to val scalar
        val_tot = np.sum(val)

        # return SyntheticFlux
        from .flux import SyntheticFlux
        f = SyntheticFlux(self)
        f.log.append(stage)
        f.stage = 'SyntheticFlux: initial'
        f.log.append(f.stage)
        f.val = val_tot
        f.wav = np.array(waf_0)
        f.filter = {'name': filter_name, 'waf_0': waf_0, 'waf_min': waf_min, 'waf_max': waf_max}

        return f

    def get_total_val(self, wav_1, wav_2):
        '''
        Collapses the val entries in the vector within the boundaries wav_1
        and wav_2 into a 0D val array. 
        
        WARNING: This tool cannot replace convolve_filter! 
                 But it can be used to produce rough estimates 
                 in-between the processes.


        Parameters
        ----------
        
        wav_1, wav_2 : float
            Boundaries in microns.
            
        
        Returns
        -------
        
        flux : SyntheticFlux
        '''

        stage = 'SyntheticSED: get_total_val'

        # for MJy/sr convert first, add and then convert back
        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=self.wav, val=self.val)
            self.val = s.get_unit(in_units=self.units, out_units='Jy', input_resolution=self.resolution['arcsec'])

        # slices within boundaries are extracted, averaged collapsed to a single scalar val
        vec = get_slices(wav=self.wav, val=self.val, wav_1=wav_1, wav_2=wav_2)
        f_total = average_collapse(val=vec['val_short'])

        # real limits within collapse
        wav_max = 10 ** (np.log10(self.wav[vec['filter_index'][0]]) + self.spacing_wav / 2.)
        wav_min = 10 ** (np.log10(self.wav[vec['filter_index'][-1]]) - self.spacing_wav / 2.)
        wav_total = central_wav(wav=[wav_min, wav_max])

        # for MJy/sr convert first, add and then convert back
        if self.unit_out == 'MJy/sr' or self.unit_out == 'Jy/arcsec^2':
            s = ConvertUnits(wav=wav_total, val=f_total)
            f_total = s.get_unit(in_units='Jy', out_units=self.unit_out, input_resolution=self.resolution['arcsec'] * self.pixel[0])

        # return SyntheticFlux
        from .flux import SyntheticFlux
        f = SyntheticFlux(self)
        f.log.append(stage)
        f.stage = 'SyntheticFlux: initial'
        f.log.append(f.stage)
        f.val = np.array(f_total)
        f.wav = np.array(wav_total)
        f.filter = {'name': 'val_tot', 'waf_0': wav_total, 'waf_min': wav_min, 'waf_max': wav_max}

        return f

    def plot_sed_multi_filter(self, multi_filter_val, multi_filter_wav, names, filter_label_size=None, ymin=10. ** (-5), dpi=None):
        '''
        Reads in array of filtered val and plots it with the current passed val of  
        SyntheticSED in a log-log diagram. 
        That way the quality of multiple filtered val can be checked, if filters are
        in database.
        If the used filters are not part of the filter database, define a OneFilter object.
        
        
        Parameters
        ----------
        
        multi_filter_val : np.ndarray
            1D vector with val entries from several filters.
                    
        multi_filter_wav : np.ndarray
            1D vector of central wavelengths from several filters.

        names : np.ndarray
            Original name of several filters, if found in database.
            
        filter_label_size : ``None``, ``True``
            Switch wether to set print labels above the filters.
            
                * ``None``: No labels plotted.
                * ``True``: Plots names[index] as label.
            
        ymin : float
            Minimal vertical limit for SED plot. 
            Default is 10.**(-10) of the maximum val.
            
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
            
            
        Returns
        -------
        
        flux : SyntheticFlux
        '''

        stage = 'SyntheticSED: plot_sed_multi_filter'

        # plot sed and val within boundaries of sed
        fig = plt.figure()
        font = {'size': 8}

        mpl.rc('font', **font)
        suptitle = fig.suptitle('output SED with filtered val', fontsize=12.)

        x_max = self.wav[0]
        x_min = self.wav[-1]
        y_max = 10 * max(self.val)
        y_min = max(self.val) * ymin

        # subplot with SED and collapsed fluval
        gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
        ax = plt.subplot(gs[0])
        plt.subplots_adjust(wspace=0, hspace=0)
        ax.set_xlim([x_min, x_max])
        ax.set_ylim([y_min, y_max])

        #ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel('val ' + '[' + str(self.unit_out) + ']')
        plt.loglog(self.wav, self.val, 'o-k', label='rough SED', linewidth=2)

        # sort c.fluval.wav, filt after c.wav
        # get index of c.wav in right order
        # [::-1] largest cwav comes first
        sort_order = np.argsort(multi_filter_wav)
        multi_filter_val = np.take(multi_filter_val, sort_order)[::-1]
        multi_filter_wav = np.take(multi_filter_wav, sort_order)[::-1]
        names = np.take(names, sort_order)[::-1]

        # color
        f = np.linspace(0, 200, len(names))

        for i in range(len(names)):
            current_waf_0 = multi_filter_wav[i]
            current_val = multi_filter_val[i]
            current_name = names[i]

            color = plt.cm.RdYlBu(int(f[i]))
            #color = plt.cm.autumn(int(f[i]))

            # filter
            text = 'filtered val by ' + str(current_name)
            x = PlotFilters(style='loglog', normalized=True, unit='unit energy')
            plot_filter = x.collect_filters(current_name)
            plt.loglog(current_waf_0, current_val, 'o', color=color)

        # subplot with filter responses
        ax2 = plt.subplot(gs[1])
        ax2.set_xlim([x_min, x_max])
        ax2.set_ylim([10. ** (-4), 10.])
        ax2.set_ylabel('filter response [unit energy]', color='k')

        for i in range(len(names)):
            current_waf_0 = multi_filter_wav[i]
            current_val = multi_filter_val[i]
            current_name = names[i]

            color = plt.cm.RdYlBu(int(f[i]))
            #color = plt.cm.autumn(int(f[i]))
            fancy_line = {'color': color, 'linestyle': '-', 'linewidth': 2}

            x = PlotFilters(style='loglog', normalized=True, unit='unit energy')
            plot_filter = x.collect_filters(current_name)
            x.plot(plot_filter, line=fancy_line)

            # text above filters
            if filter_label_size is not None:
                plt.text(current_waf_0, 3, current_name, horizontalalignment='center', verticalalignment='center', fontsize=filter_label_size)

        ax2.set_xlabel(r'$\lambda$ [$\mu$m]')

        # legend merged
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()

        leg = ax.legend(h1 + h2, l1 + l2, loc=(1.03, -0.5))
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize

        fig.savefig(str(self.name) + '_' + 'process-output_SS-multi-filter.png', dpi=dpi, bbox_inches='tight')

        # return SyntheticSED
        s = SyntheticSED(self)
        s.log.append(stage)
        s.stage = 'SyntheticSED: initial'
        s.log.append(s.stage)

        return s

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
