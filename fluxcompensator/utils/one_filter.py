import numpy as np
import matplotlib.pyplot as plt
from hyperion.util.constants import c


class OneFilter(object):

    '''
    Tool to plot one filter of choise.
        
    
    Parameters
    ----------
    
    name : str
        Name of the filter.
    
    filter_file : str
        Location of file e. g. filter_file.txt where wavelength (column 1) and response
        (column 2) of the filter is stored.
        
    beta : float
        Only important if multiple filters with PlotFilters are plotted. 
        Default is ``None``.
        Power of nu, when ``R_input = R * nu**(beta)``.
        If unit of R_input is per unit energy : ``beta = -1``.
        If unit of R_input is per unit photon : ``beta = 0``.
    
    Returns
    -------

    wave_filter : numpy.ndarray
        Wavelength vector of filter with descending wavelength values 
        in microns.
        
    R_filter : numpy.ndarray
        Response of the filter in input units. beta could give hints 
        of the units.
        
    Rmax : float
        Maximum of the vertical component of the filter curve.
        
    nu_Rmax : float
        Frequency at maximum of the vertical component of the filter curve
        in Hz.
    
    '''

    def __init__(self, name, filter_file, beta=None):
        self.name = name
        self.beta = beta

        # read in filter properties
        f = np.loadtxt(filter_file)
        wave_filter = f[:, 0]
        R_filter = f[:, 1]

        if wave_filter[0] < wave_filter[1]:
            wave_filter = wave_filter[::-1]
            R_filter = R_filter[::-1]

        self.wave_filter = np.array(wave_filter)  # (350)
        self.R_filter = np.array(R_filter)  # (350)
        self.Rmax = np.max(self.R_filter)
        self.nu_Rmax = c / (self.wave_filter[self.R_filter == self.Rmax] * 10. ** (-4))
        

    def onefilter_plot(self, line='-r', style='linear', normalized=None, unit=None):
        '''
        Only important if multiple filters with PlotFilters are plotted.             
        
        Parameters
        ----------
        
        line : str, dict
            Linestyle of the filter curve. Valid are all python commands like:
            
                * ``'-r'`` : red with solid line
                * ``'--g'`` : green with dashed line
                * ``'.-b'`` : blue with dot dashed line 
                * dict : line{'color' : color, 'linestyle' : '-', 'linewidth' : 2}
                
            Default is ``'-r'``.
        
        style : str
            Style of the plotting axis. Valid options are:
            
                * ``'linear'`` : red with solid line
                * ``'loglog'`` : green with dashed line
                * ``'logx'`` : blue with dot dashed line 
                * ``'logy'`` : blue with dot dashed line 

            Default is ``'linear'``.
            
        normalized : ``None``, ``True``
                
                * ``None`` : Plot original scale.
                * ``True`` : Plot normalized curve with unity.
            
            Default is ``None``.
            
        unit : str, ``None``
            
            Valid options are: 
                
                * ``'unit electrons'``
                * ``'unit energy'``
            
            Default is ``None``, which might plot filters in different units.
        
        Returns
        -------

        plot : plt command
            Plotting command to pass to PlotFilters.
        '''

        self.style = style
        self.normalized = normalized
        self.unit = unit
        self.nu = c / (self.wave_filter * 10. ** (-4))
        self.line = line

        if unit is not None:
            if self.beta is None:
                raise Exception('Warning: in onefilter_plot, beta needs to be defined')

            if unit == 'unit energy' and self.beta == 0.:
                self.R_filter_unit = self.R_filter / self.nu

            if unit == 'unit energy' and self.beta == -1.:
                self.R_filter_unit = self.R_filter

            if unit == 'unit electron' and self.beta == 0.:
                self.R_filter_unit = self.R_filter

            if unit == 'unit electron' and self.beta == -1.:
                self.R_filter_unit = self.R_filter * self.nu
        else:
            self.normalized = True
            self.R_filter_unit = self.R_filter 

        if self.normalized is not None:
            self.Rmax = max(self.R_filter_unit)
            self.R_filter_unit = self.R_filter_unit / self.Rmax
        
        self.Rmax = np.max(self.R_filter_unit)

        if type(line) != dict:
            if self.style == 'logx':
                pass_back_plot = plt.semilogx(self.wave_filter, self.R_filter_unit, line, linewidth=2, label=self.name)
            if self.style == 'logy':
                pass_back_plot = plt.semilogy(self.wave_filter, self.R_filter_unit, line, linewidth=2, label=self.name)
            if self.style == 'loglog':
                pass_back_plot = plt.loglog(self.wave_filter, self.R_filter_unit, line, linewidth=2, label=self.name)
            if self.style == 'linear':
                pass_back_plot = plt.plot(self.wave_filter, self.R_filter_unit, line, linewidth=2, label=self.name)

        if type(line) == dict:
            if self.style == 'logx':
                pass_back_plot = plt.semilogx(self.wave_filter, self.R_filter_unit, label=self.name, **line)
            if self.style == 'logy':
                pass_back_plot = plt.semilogy(self.wave_filter, self.R_filter_unit, label=self.name, **line)
            if self.style == 'loglog':
                pass_back_plot = plt.loglog(self.wave_filter, self.R_filter_unit, label=self.name, **line)
            if self.style == 'linear':
                pass_back_plot = plt.plot(self.wave_filter, self.R_filter_unit, label=self.name, **line)

        return pass_back_plot

    def onefilter_save_plot(self, dpi=None):
        '''
        Only important if one filter wants to be plotted directly with 
        original scale.             
        
        Parameters
        ----------
        
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
        
        
        Returns
        -------

        plot : e. g. plot.png
            Plot of filter with chosen style.
        '''

        # visualize filter approximation by monochromatic filter
        self.fig = plt.figure()
        titel = 'filter curve ' + self.name
        self.fig.suptitle(titel, fontsize=16.)
        ax = self.fig.add_subplot(1, 1, 1)

        ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        ax.set_ylabel(r'Response')

        if type(self.line) != dict:
            if self.style == 'logx':
                one_plot = plt.semilogx(self.wave_filter, self.R_filter_unit, self.line, linewidth=2, label=self.name)
            if self.style == 'logy':
                one_plot = plt.semilogy(self.wave_filter, self.R_filter_unit, self.line, linewidth=2, label=self.name)
            if self.style == 'loglog':
                one_plot = plt.loglog(self.wave_filter, self.R_filter_unit, self.line, linewidth=2, label=self.name)
            if self.style == 'linear':
                one_plot = plt.plot(self.wave_filter, self.R_filter_unit, self.line, linewidth=2, label=self.name)

        if type(self.line) == dict:
            if self.style == 'logx':
                one_plot = plt.semilogx(self.wave_filter, self.R_filter_unit, label=self.name, **self.line)
            if self.style == 'logy':
                one_plot = plt.semilogy(self.wave_filter, self.R_filter_unit, label=self.name, **self.line)
            if self.style == 'loglog':
                one_plot = plt.loglog(self.wave_filter, self.R_filter_unit, label=self.name, **self.line)
            if self.style == 'linear':
                one_plot = plt.plot(self.wave_filter, self.R_filter_unit, label=self.name, **self.line)

        # Set view limits
        ax.set_xlim(np.min(self.wave_filter), np.max(self.wave_filter))
        if self.style == 'logy' or self.style == 'loglog':
            ax.set_ylim(10. ** (-4), 1)
        else:
            ax.set_ylim(0, self.Rmax)
        ax.legend(loc='best')

        self.fig.savefig('filter_' + self.name + '.png', dpi=dpi)
