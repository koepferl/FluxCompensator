import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from hyperion.util.integrate import integrate_subset
from scipy import stats

from .units import ConvertUnits


class MakePlots(object):

    '''
    Plots slices of the val cube from SyntheticCube a the closest slice
    to wav_interest or val images of SyntheticImage. The boundaries of
    the slice or image is ploted as well. The plot has 4 different cut 
    levels at 100, 99, 95 and 90 percent.
    
    
    Parameters
    ----------
    
    input_array : SyntheticCube, SyntheticImage, optional
        input_array also reads arrays with SyntheticCube and
        SyntheticImage properties.
        
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
                  Minimal and Maximal physical value presented in the colorbars.
        * ``None`` : no plot with minimal and maximal cut is returned.
        
        Default is ``None``.
        
    dpi  :  ``None``, scalar > 0 
        The resolution in dots per inch. 
        ``None`` is default and will use the value savefig.dpi 
        in the matplotlibrc file.
        
    
    Returns
    -------
    
    cube : SyntheticCube
        3D val array with SyntheticCube properties.
        
    image : SyntheticImage
        2D val array with SyntheticImage properties.    
    '''

    def __init__(self, input_array, wav_interest=None, prefix=None, name=None, multi_cut=None, single_cut=None, set_cut=None, dpi=None):
        
        if multi_cut == None and single_cut == None and set_cut == None:
            raise Exception('At least one plotting routine (multi_cut, single_cut or set_cut == None) has to be chosen.')

        self.prefix = prefix

        if self.prefix is None and name is None:
            raise Exception('If prefix name is not given, you need to give the a name to enable the default naming chain.')

        if input_array.val.ndim in (2, 3):
            # input_array properties
            self.input_name = name

            self.name = input_array.name
            self.unit_out = input_array.unit_out

            self.val = input_array.val
            self.wav = input_array.wav
            self.wav_interest = wav_interest
            self.filter = input_array.filter
            #print self.filter

            self.grid_unit = input_array.grid_unit
            self.grid_unit_name = input_array.grid_unit_name

            # measure of the image
            self.FOV = input_array.FOV
            self.x_min = input_array.x_min
            self.x_max = input_array.x_max
            self.y_min = input_array.y_min
            self.y_max = input_array.y_max

            self.pixel = input_array.pixel
            self.pixel_2D = self.pixel[0] * self.pixel[1]
            #print self.val.shape

            # condition to find slice close to wav_interest
            if input_array.val.ndim == 3:
                if self.wav_interest is None and self.filter['waf_min'] is None:
                    raise Exception('WARNING: wav_interest or waf_0 need to be defined if 3D cube is pasted.')
                
                find_minimum = np.abs(self.wav - self.wav_interest)
                num = np.arange(len(self.wav))
                
                index = num[find_minimum == np.min(find_minimum)][0]
                wav_min = 10. ** (np.log10(self.wav[index]) - input_array.spacing_wav / 2.)
                wav_max = 10. ** (np.log10(self.wav[index]) + input_array.spacing_wav / 2.)
                self.val_2D = self.val[:, :, index]
                self.wav_real = (round(wav_min, 2), round(wav_max, 2))                

            # condition for image
            if input_array.val.ndim == 2:
                self.val_2D = self.val.copy()
                self.wav_real = (round(self.filter['waf_min'], 2), round(self.filter['waf_max'], 2))

        else:
            raise Exception('WARNING: MakePlots only can use SyntheticCube or SyntheticImage.')

        # creat cut levels
        self.val_sort = np.sort(self.val_2D.ravel())
        self.xx = np.linspace(0, len(self.val_sort), len(self.val_sort))

        # statstics
        self.median = stats.scoreatpercentile(self.val_sort, 50)
        self.min_0 = stats.scoreatpercentile(self.val_sort, 0)
        self.min_5 = stats.scoreatpercentile(self.val_sort, 5)
        self.max_95 = stats.scoreatpercentile(self.val_sort, 95)
        self.max_100 = stats.scoreatpercentile(self.val_sort, 100)

        # grid of X, Y plot
        x = np.linspace(self.x_min / self.grid_unit, self.x_max / self.grid_unit, self.pixel[0])
        y = np.linspace(self.y_min / self.grid_unit, self.y_max / self.grid_unit, self.pixel[1])
        X, Y = np.meshgrid(y,x)
        label = 'Flux [' + self.unit_out + ']'

        # titel of plot
        titel = self.name + ' ' + str(self.wav_real) + ' micron'

        # ploting multiple plot
        if multi_cut is not None:
            fig2 = plt.figure()
            fig2.suptitle(titel, fontsize=10.)
            fig2.subplots_adjust(hspace=0.3, wspace=0.3)

            font = {'size': 6}
            mpl.rc('font', **font)

            a = np.array([100, 99, 95, 90])
            b = np.array([1, 2])
            c = np.array([1, 1])

            for l in range(len(a)):
                k = l + 1
                ax = fig2.add_subplot(2, 2, k)
                title = str(int(a[l])) + ' %'
                plt.title(title)
                self.percentage = a[l]
                self.percent = self.percentage / 2.

                lower_cut = (100 - a[l]) / 2.
                upper_cut = lower_cut + a[l]
                self.min = stats.scoreatpercentile(self.val_sort, lower_cut)
                self.max = stats.scoreatpercentile(self.val_sort, upper_cut)
                vmin = self.min
                vmax = self.max

                self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
                c = ax.pcolormesh(X, Y, self.val_2D, cmap=plt.cm.gist_heat, norm=self.norm)

                ax.set_xlim(x[0], x[-1])
                ax.set_ylim(y[0], y[-1])
                ax.set_xlabel('x [' + self.grid_unit_name + ']')
                ax.set_ylabel('y [' + self.grid_unit_name + ']')
                cb = fig2.colorbar(c)
                cb.set_label(label)

            if self.prefix is None:
                self.plot_name = self.name + '_image_' + self.input_name + '_multi_cut_' + str(self.wav_real[0]) + '_' + str(self.wav_real[1]) + '.png'
            if self.prefix is not None:
                self.plot_name = self.prefix + '.png'

            fig2.savefig(self.plot_name, bbox_inches='tight', dpi=dpi)

        # single plot for certain cut if cut is not None
        if single_cut is not None:
            fig3 = plt.figure()
            fig3.suptitle(titel, fontsize=10.)
            fig3.subplots_adjust(hspace=0.3, wspace=0.3)

            font = {'size': 6}
            mpl.rc('font', **font)

            a = np.array([single_cut])

            ax = fig3.add_subplot(1, 1, 1)
            title = str(int(a[0])) + ' %'
            plt.title(title)

            lower_cut = (100 - single_cut) / 2.
            upper_cut = lower_cut + single_cut
            self.min = stats.scoreatpercentile(self.val_sort, lower_cut)
            self.max = stats.scoreatpercentile(self.val_sort, upper_cut)
            vmin = self.min
            vmax = self.max

            self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            c = ax.pcolormesh(X, Y, self.val_2D, cmap=plt.cm.gist_heat, norm=self.norm)

            ax.set_xlim(x[0], x[-1])
            ax.set_ylim(y[0], y[-1])
            ax.set_xlabel('x [' + self.grid_unit_name + ']')
            ax.set_ylabel('y [' + self.grid_unit_name + ']')
            cb = fig3.colorbar(c)
            cb.set_label(label)

            if self.prefix is None:
                self.plot_name = self.name + '_image_' + self.input_name + '_single_cut_' + str(single_cut) + '%_' + str(self.wav_real[0]) + '_' + str(self.wav_real[1]) + '.png'

            if self.prefix is not None:
                self.plot_name = self.prefix + '.png'

            fig3.savefig(self.plot_name, bbox_inches='tight', dpi=dpi)

        # single plot for certain cut if cut is not None
        if set_cut is not None:
            fig4 = plt.figure()
            fig4.suptitle(titel, fontsize=10.)
            fig4.subplots_adjust(hspace=0.3, wspace=0.3)

            font = {'size': 6}
            mpl.rc('font', **font)

            # min max
            vmin = set_cut[0]
            vmax = set_cut[1]

            ax = fig4.add_subplot(1, 1, 1)
            title = '[' + str("%0.2e" % vmin) + ', ' + str("%0.2e" % vmax) + ']'
            title2 = 'flux values within ' + title + ' ' + self.unit_out
            plt.title(title2)

            self.norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
            #print X.shape, Y.shape, self.val_2D.shape
            c = ax.pcolormesh(X, Y, self.val_2D, cmap=plt.cm.gist_heat, norm=self.norm)

            ax.set_xlim(x[0], x[-1])
            ax.set_ylim(y[0], y[-1])
            ax.set_xlabel('x [' + self.grid_unit_name + ']')
            ax.set_ylabel('y [' + self.grid_unit_name + ']')
            cb = fig4.colorbar(c)
            cb.set_label(label)

            if self.prefix is None:
                self.plot_name = self.name + '_image_' + self.input_name + '_set_cut_' + str("%0.2e" % vmin) + '_' + str("%0.2e" % vmax) + '_' + str(self.wav_real[0]) + '_' + str(self.wav_real[1]) + '.png'

            if self.prefix is not None:
                self.plot_name = self.prefix + '.png'

            fig4.savefig(self.plot_name, bbox_inches='tight', dpi=dpi)

    def histogram_cuts(self, dpi=None):
        '''
        DS9 like histograms of the cuts can be checked here.
        
        
        Parameters
        ----------
        
        dpi  :  None, scalar > 0 
            The resolution in dots per inch. 
            ``None`` is default and will use the val savefig.dpi 
            in the matplotlibrc file.        
        '''

        # Histograms of cut levels
        fig = plt.figure()
        fig.suptitle(self.name, fontsize=10.)
        fig.subplots_adjust(hspace=0.3)

        ax1 = fig.add_subplot(2, 1, 1)
        plt.semilogy(self.val_sort[::-1], self.xx, 'b-')
        plt.semilogy([self.median, self.median], [self.xx[0] + 1., self.xx[-1]], 'r-')
        # min max
        plt.semilogy([self.min_0, self.min_0], [self.xx[0] + 1., self.xx[-1]], 'g-')
        plt.semilogy([self.min_5, self.min_5], [self.xx[0] + 1., self.xx[-1]], 'y-')
        plt.semilogy([self.max_95, self.max_95], [self.xx[0] + 1., self.xx[-1]], 'y-')
        plt.semilogy([self.max_100, self.max_100], [self.xx[0] + 1., self.xx[-1]], 'g-')
        ax1.set_xlabel('val distribution')
        ax1.set_ylabel('Number of pixels')

        ax2 = fig.add_subplot(2, 1, 2)
        plt.plot(self.val_sort)
        ax2.set_xlabel('Number of pixels')
        ax2.set_ylabel('val distribution')
        ax2.set_xlim(self.xx[0], self.xx[-1])
        ax2.set_ylim(self.val_sort[0], self.val_sort[-1])

        if self.prefix is None:
            self.plot_name = self.name + '_image_' + self.input_name + '_' + str(self.wav_real[0]) + '_' + str(self.wav_real[1]) + '_histogram.png'
        if self.prefix is not None:
            self.plot_name = self.prefix + '_histogram.png'

        fig.savefig(self.plot_name, bbox_inches='tight', dpi=dpi)
