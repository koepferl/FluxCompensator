import numpy as np
import matplotlib.pyplot as plt

from .one_filter import OneFilter


class PlotFilters(object):

    '''
    Tool to plot multiple filter of choice.
    
        
    Parameters
    ----------
    
    style : str
        Style of the plotting axis. Valid options are:
        
            * ``'linear'`` : red with solid line
            * ``'loglog'`` : green with dashed line
            * ``'logx'`` : blue with dot dashed line 
            * ``'logy'`` : blue with dot dashed line 
        
    normalized : ``None``, ``True``
            
            * ``None`` : Plot original scale.
            * ``True`` : Plot normalized curve with unity.
                
    unit : str, ``None``
        
        Valid options are: 
            
            * ``'unit electrons'``
            * ``'unit energy'``
            * ``None`` : Which might plot filters in different units.
    '''

    def __init__(self, style, normalized, unit):
        self.style = style
        self.normalized = normalized
        self.unit = unit
        self.R_filter_lim = []
        self.wave_filter_lim = []

    def collect_filters(self, filter_database, name=None, filter_file=None, beta=None):
        '''
        Collect filters of choice.
        
    
        Parameters
        ----------

        filter_database : str, ``None``
            
            * ``None`` : Choose own filter, therefore name, filter_file and
                         beta need to be defined.
            * str : Choose filter ``name_PLOT`` from database like in FluxCompensator.
            
        name : str
            Name of the own filter. Default is ``None``.
            
            
        filter_file : str, ``None``
            Location of own file e. g. filter_file.txt where wavelength 
            (column 1) and response (column 2) of the filter is stored.         
            Default is ``None``. Than filter_database is used.
        
        beta : float, ``None``
            Default is ``None``. Then filter_database is used.
            Power of nu, when ``R_input = R * nu**(beta)``.
            If unit of R_input is per unit energy : ``beta = -1``.
            If unit of R_input is per unit photon : ``beta = 0``.
    
        Returns
        -------

        wave_filter_lim : numpy.ndarray
            Vector of the wavelength limits of all filters red in so far.
        
        R_filter_lim : numpy.ndarray
            Vector of the response maxima red in so far.
            
        data_filter : object
            Returns database like object, even if own filter was red in.
        
        '''
        # use database or own filter
        import fluxcompensator.database.missions as of
        if filter_database is not None and isinstance(getattr(of, filter_database), OneFilter):
            data_filter = getattr(of, filter_database)
        else:
            data_filter = OneFilter(name=name, filter_file=filter_file, beta=beta)

        # store relevant valures for plotting limits
        if self.unit is not None:
            self.normalize = True
                
            if self.unit == 'unit energy' and data_filter.beta == 0.:
                R_max = data_filter.Rmax / data_filter.nu_Rmax
            
            if self.unit == 'unit energy' and data_filter.beta == -1.:
                R_max = data_filter.Rmax
            
            if self.unit == 'unit electron' and data_filter.beta == 0.:
                R_max = data_filter.Rmax
            
            if self.unit == 'unit electron' and data_filter.beta == -1.:
                R_max = data_filter.Rmax * data_filter.nu_Rmax
            
        if self.unit is None:
            R_max = data_filter.Rmax
            
        self.R_filter_lim.append(np.max(R_max))
        self.wave_filter_lim.append(data_filter.wave_filter[0])
        self.wave_filter_lim.append(data_filter.wave_filter[-1])

        return data_filter

    def get_axis(self):
        '''
        Sets up plotting frame and defines axis.
        '''

        self.fig = plt.figure()
        if self.unit is not None:
            titel = 'filter curves in ' + self.unit
        else:
            titel = 'filter curves'

        self.fig.suptitle(titel, fontsize=16.)
        self.ax = self.fig.add_subplot(1, 1, 1)

        # labels
        self.ax.set_xlabel(r'$\lambda$ [$\mu$m]')
        self.ax.set_ylabel(r'Response')

        # limits
        w = np.sort(self.wave_filter_lim)
        R = np.max(self.R_filter_lim)

        self.w = [w[0], w[-1]]
        self.interval = [round(self.w[0], 1), round(self.w[-1], 1)]

        self.ax.set_xlim(self.w)
        
        if self.unit is not None:
            if self.style == 'logy' or self.style == 'loglog':
                self.ax.set_ylim(10. ** (-4), 1)
            else:
                self.ax.set_ylim(0, R)

    def plot(self, data_filter, line):
        '''
        Plots filter curves in empty frame.
        
    
        Parameters
        ----------

        data_filter : object
            Database like object returned from collect_filters.
            
        line : str, dict
            Linestyle of the filter curve. Valid are all python commands like:
            
                * ``'-r'`` : red with solid line
                * ``'--g'`` : green with dashed line
                * ``'.-b'`` : blue with dot dashed line 
                * dict : line{'color' : 'r', 'linestyle' : '-', 'linewidth' : 2}
    
        '''

        plot = data_filter.onefilter_plot(line=line, style=self.style, normalized=self.normalized, unit=self.unit)

    def save_plot(self, name=None, dpi=None):
        '''
        Saves the plots of the filter response curves.
        
    
        Parameters
        ----------

        name : str
            Name of the filter collection e. g. Spitzer.
            
        dpi  :  ``None``, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.
        '''

        self.ax.legend(loc='best')

        if name is None:
            self.fig.savefig('filter_' + str(self.interval[0]) + '_' + str(self.interval[1]) + '.png', dpi=dpi)
            plt.close()
        if name is not None:
            self.fig.savefig('filter_' + name + '_' + str(self.interval[0]) + '_' + str(self.interval[1]) + '.png', dpi=dpi)
            plt.close()
