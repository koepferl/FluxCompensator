import numpy as np
from astropy import log as logger
from astropy.io import fits
import matplotlib.pyplot as plt
from hyperion.util.constants import c
from hyperion.util.integrate import integrate_subset
#val = integrate_subset(x, y, xmin, xmax)

from .utils.tools import properties, get_slices


class Filter(object):

    '''
    Convolves slice of SyntheticCube val within filter limits in a 2D val 
    image or entries of SyntheticSED val into a scalar val.
    
    Parameters
    ----------
    
    name : str
        Name of filter, default is ``None``.
        
    filter_file : str
        Location of e. g. filter_file.txt file with filter function.
        column #1: wave_filter in microns
        column #2: Response_filter
    
        Default is ``None``.
    
    waf_0 : float
        Central wavelength of filter in microns. Default is ``None``.
          
    alpha : float
        Power of nu, when ``nu**(alpha) * F_nu = const``. 
        Default is ``None``.
    
    beta : float
        Power of nu, when ``R_input = R * nu**(beta)``.
        If unit of R_input is per unit energy : ``beta = -1``.
        If unit of R_input is per unit photon : ``beta = 0``.
        Default is ``None``.
    '''

    def __init__(self, name, filter_file, waf_0, alpha, beta):
        self.name = name
        self.filter_file = filter_file
        self.waf_0 = waf_0
        self.alpha = alpha
        self.beta = beta

    def rebin(self, wav, val, inclination=None):
        '''
        Rebins the filter response into bins which are passed by wav and val.
        
        Parameters
        ----------
        
        wav : numpy.ndarray
            The wavelengths of val array slices / entries in microns.
    
        val : numpy.ndarray
            The 3D cube with shape (x, y, wav) or 1D sed with shape like wav.
            
        inclination : str, ``None``
            Gives index range of monotonically decreasing wavelength values for 
            arrays with more than one inclination.
            
            Use indices from 100 to (and excluding) 200 set ``inclination = '100:200'``.

        
        Returns
        -------
        
        weight : dict
            Dictionary relevant for FluxCompensator.convolve_filter.
            
            ``weight{wav_short, val_short, Response_new, filter_index, waf_0, waf_min, waf_max}``
             
                 * ``wav_short`` : Vector like wav but all entries outside
                                   the filter boundaries are erased.
                 * ``val_short`` : Shape like wav_short for val.
                 * ``Response_new`` : New response in bins of wav_short.
                 * ``filter_index`` : Indices of original wav within filter.
                 * ``waf_0`` : Central wavelength of filter.
                 * ``waf_min`` : Lower boundary of filter.
                 * ``waf_max`` : Upper boundary of filter.
                 * ``filter_name`` : Name of filter. E. g. IRAC4_FILTER
        '''
        
        # deal with different inclinations in one file
        if val.ndim == 1 and inclination is not None:
            bound = inclination.split(':')
            wav = wav[int(bound[0]):int(bound[1])]
            val = val[int(bound[0]):int(bound[1])]

        if val.ndim == 1 and (len(wav[wav[:-1]<wav[1:]]) < len(wav)-1):
            raise Exception('More than one inclination, use parameter "inclination" the index range of your preferred segment, where the wavelength decreases monotonically.')
        

        self.wav = np.array(wav).astype(float)
        self.val = np.array(val).astype(float)

        nuf_0 = c / (self.waf_0 * 10. ** (-4))

        # read in filter properties
        f = np.loadtxt(self.filter_file)
        wave_filter = f[:, 0]
        R_filter = f[:, 1]

        if wave_filter[0] < wave_filter[1]:
            wave_filter = wave_filter[::-1]
            R_filter = R_filter[::-1]

        self.wave_filter = np.array(wave_filter)  # (350)
        self.R_filter = np.array(R_filter)  # (350)

        nu_filter = c / (wave_filter * 10. ** (-4))
        self.nu_filter = np.array(nu_filter)

        # wavelength bins in filter
        self.waf_max = wave_filter[0]
        self.waf_min = wave_filter[-1]
        nuf_min = c / (self.waf_max * 10. ** (-4))
        nuf_max = c / (self.waf_min * 10. ** (-4))

        self.nu = c / (self.wav * 10. ** (-4))
        self.nu = np.array(self.nu)

        # Getting the proper filter definition
        integrad = R_filter / nu_filter ** (1. + self.alpha + self.beta)
        integral = integrate_subset(nu_filter, integrad, nuf_min, nuf_max)
        R_proper = (1. / nuf_0 ** self.alpha) * (R_filter / nu_filter ** (1. + self.beta) / integral)

        # wavelength bins in model
        spacing_wav = properties(wav=self.wav, val=self.val)['spacing_wav']

        # eliminate wavelength outside filter
        cube_short = get_slices(wav=self.wav, val=self.val, wav_1=self.waf_min, wav_2=self.waf_max)
        wav_short = cube_short['wav_short']
        val_short = cube_short['val_short']
        filter_index = cube_short['filter_index']
        
        if len(wav_short) == 0:
            raise Exception('Wavelength range of radiative transfer model lies outside of filter boundaries.')

        # weight by filter response
        Response_new = np.array(())
        self.plot_Response = np.array(())
        self.plot_lam = np.array(())
        for i in range(len(wav_short)):
            int_max = c / (10. ** (-4) * 10. ** (np.log10(wav_short[i]) - spacing_wav / 2.))
            int_min = c / (10. ** (-4) * 10. ** (np.log10(wav_short[i]) + spacing_wav / 2.))
            if int_max > nuf_max:
                int_max = nuf_max
            if int_min < nuf_min:
                int_min = nuf_min
            R_int = integrate_subset(nu_filter, R_proper, int_min, int_max)
            Response_new = np.append(Response_new, R_int)

            # for the plot
            self.plot_Response = np.append(self.plot_Response, R_int)
            self.plot_Response = np.append(self.plot_Response, R_int)
            wmax = 10. ** (4) * c / int_min
            self.plot_lam = np.append(self.plot_lam, wmax)
            wmin = 10. ** (4) * c / int_max
            self.plot_lam = np.append(self.plot_lam, wmin)

        # debugging comment
        logger.debug('-' * 30)
        logger.debug('rebined response for ' + self.name)
        logger.debug('-' * 30)
        for i in range(len(wav_short)):
            logger.debug(str('%8.4f' % wav_short[i]) + ', ' + str('%2.5f' % Response_new[i]))
        logger.debug('-' * 30)
        logger.debug('sum of rebinded response: ' + str(sum(Response_new)))
        logger.debug('sum of original response: ' + str(sum(R_filter)))
        logger.debug('max of original response: ' + str(max(R_filter)))

        return {'wav_short': wav_short, 'val_short': val_short, 'Response_new': Response_new, 'filter_index': filter_index, 'waf_0': self.waf_0, 'waf_min': self.waf_min, 'waf_max': self.waf_max, 'filter_name': self.name}

    def plot(self, val_name=None, dpi=None):
        '''
        Plots the original filter curves with the new rebined response.
        
        Parameters
        ----------
        
        val_name : str
            Name of the filter which will appear in the plot file's name.
            It should be equal with name in SyntheticCube, SyntheticImage,
            SyntheticSED or SyntheticFlux.
        
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
        '''

        # visualize filter approximation by monochromatic filter
        fig = plt.figure()
        titel = 'filter approximation ' + self.name
        fig.suptitle(titel, fontsize=16.)
        ax = fig.add_subplot(1, 1, 1)

        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(r'Response')

        plt.plot(self.plot_lam, self.plot_Response, 'o-r')
        plt.plot(self.wave_filter, self.R_filter, '-k')

        # Set view limits
        ax.set_xlim(self.waf_max, self.waf_min)
        ax.set_ylim(min(self.R_filter), 1.)

        if val_name is None:
            val_name = 'no_val_name_defined'

        fig.savefig(val_name + '_process-output_FR-' + self.name + '-rebin.png', dpi=dpi)
