from astropy import log as logger
from astropy.io import fits
import numpy as np
import matplotlib.pylab as plt
from pylab import *
from datetime import datetime

from .tools import where_is_1D

class ConservingZoom(object):

    '''
    Changes resolution of arrays under conservation of the total sum.
    E. g. total val conservation.

    Parameters
    ----------

    array : numpy.ndarray
        2D array, where resolution should be changed.

    initial_resolution : float
        Initial resolution of the entries in the array.
        Units are ``arcsec/pixel``.

    new_resolution : float
        Wanted resolution of the entries in the array.
        Units are ``arcsec/pixel``.
    '''

    def __init__(self, array, initial_resolution, new_resolution):

        # if val not 3D currently change it to 2D
        if len(np.shape(array)) == 2:
            self.fake_cube = True
            array = np.array([array]).swapaxes(0, 2).swapaxes(0, 1)
        else:
            self.fake_cube = None

        self.val = array
        self.len_wav = len(self.val[0, 0, :])
        self.ro = initial_resolution
        self.rn = new_resolution

        # number of pixels
        self.len_ox = len(self.val[:, 0, 0])
        self.len_oy = len(self.val[0, :, 0])
        # real new lenght
        self.len_nrx = self.len_ox * self.ro / self.rn
        self.len_nry = self.len_oy * self.ro / self.rn
        # rounded up new lenght
        self.len_nx = int(round(self.len_nrx + 0.4999))
        self.len_ny = int(round(self.len_nry + 0.4999))

        # boundaries of pixels
        self.pix_nx = np.arange(0, (self.len_nx + 1), 1)[1:] * self.rn
        self.pix_ny = np.arange(0, (self.len_ny + 1), 1)[1:] * self.rn

        self.pix_ox = np.arange(0, (self.len_ox + 1), 1)[1:] * self.ro
        self.pix_oy = np.arange(0, (self.len_oy + 1), 1)[1:] * self.ro

        self.diff_x = self.pix_nx[-1] - self.pix_ox[-1]
        self.diff_y = self.pix_ny[-1] - self.pix_oy[-1]

        self.pix_nx = self.pix_nx - self.diff_x / 2.
        self.pix_ny = self.pix_ny - self.diff_y / 2.

        # size of one pixel
        self.size_o = self.ro ** 2
        self.size_n = self.rn ** 2

    def zoom(self):
        '''
        Actual change of resolution or zooming. Is an O(n^2) operation.

        Returns
        ---------

        array : numpy.ndarray
            3D array with changed resolution.
            Note: now number of x and y entireties may have changed.

        '''
        
        #print np.sum(self.val)

        startTime = datetime.now()

        # x direction
        x = np.zeros(shape=(self.len_nx, self.len_ox))
        for nx in range(self.len_nx):
            if nx == 0:
                n = (0, self.pix_nx[nx])
            if nx > 0:
                n = (self.pix_nx[nx - 1], self.pix_nx[nx])

            for ox in range(self.len_ox):
                if ox == 0:
                    o = (0, self.pix_ox[ox])
                if ox > 0:
                    o = (self.pix_ox[ox - 1], self.pix_ox[ox])

                frac = where_is_1D(n, o)
                x[nx, ox] = frac

        # stores enteries which are not zero on diagonal matrix
        start_x = []
        stop_x = []
        for i in range(self.len_nx):
            for j in range(self.len_ox):
                if (j > 0 and x[i][j] != 0. and x[i][j - 1] == 0.) or (j == 0 and x[i][j] != 0.):
                    start_x.append(j)
                if (j == (self.len_ox - 1) and x[i][self.len_ox - 1] != 0.) or (j < (self.len_ox - 1) and x[i][j] != 0. and x[i][j + 1] == 0):
                    stop_x.append(j)

        # y direction
        y = np.zeros(shape=(self.len_ny, self.len_oy))
        for ny in range(self.len_ny):
            if ny == 0:
                n = (0, self.pix_ny[ny])
            if ny > 0:
                n = (self.pix_ny[ny - 1], self.pix_ny[ny])

            for oy in range(self.len_oy):
                if oy == 0:
                    o = (0, self.pix_oy[oy])
                if oy > 0:
                    o = (self.pix_oy[oy - 1], self.pix_oy[oy])

                frac = where_is_1D(n, o)
                y[ny, oy] = frac

        # stores enteries which are not zero on diagonal matrix
        start_y = []
        stop_y = []

        for i in range(self.len_ny):
            for j in range(self.len_oy):
                if (j > 0 and y[i][j] != 0. and y[i][j - 1] == 0.) or (j == 0 and y[i][j] != 0.):
                    start_y.append(j)
                if (j == (self.len_oy - 1) and y[i][self.len_oy - 1] != 0.) or (j < (self.len_oy - 1) and y[i][j] != 0. and y[i][j + 1] == 0):
                    stop_y.append(j)

        nn = np.zeros(shape=(self.len_nx, self.len_ny, self.len_wav))

        for i in range(self.len_nx):
            for j in range(self.len_ny):
                ox_a = start_x[i]
                ox_e = stop_x[i]
                ox = np.arange(ox_a, ox_e + 1, 1)

                oy_a = start_y[j]
                oy_e = stop_y[j]
                oy = np.arange(oy_a, oy_e + 1, 1)

                #new = 0.
                for k in ox:
                    for l in oy:
                        for w in range(self.len_wav):
                            nn[i][j][w] = nn[i][j][w] + x[i][k] * y[j][l] * self.val[k][l][w]

        # debugging statments
        # logger.debug('-'*30)
        # logger.debug('ConservingZoom')
        # logger.debug('-'*30)
        logger.debug('total val initial array : ' + str('%4.4e' % np.sum(self.val)))
        logger.debug('total val new array     : ' + str('%4.4e' % np.sum(nn)))
        #logger.debug('compilation time         : ' + str((datetime.now()-startTime)))

        if self.fake_cube is True:
            nn = nn[:, :, 0]
        
        #print np.sum(nn)

        return nn

    def zoom_grid(self, name, dpi=None, zoom_in=1.):
        '''
        Visualization of old and new grid.

        Parameters
        ----------

        name : str
            Name of the current project.
            E.g. use name from the FluxCompensator run to match with the file names.

        dpi  :  ``None``, scalar > 0
            The resolution in dots per inch.
            ``None`` is default and will use the val savefig.dpi
            in the matplotlibrc file.

        zoom_in  :  float <= 1.
            If there are too many pixels it is helpful,
            to zoom in at the origin of the grid. The float
            ``zoom_in`` is the fraction of the initial grid size in arcsec.
            Default is ``1.``.


        '''
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

        nx = np.sort(np.append(self.pix_nx, self.pix_nx[0] - self.rn))
        ny = np.sort(np.append(self.pix_ny, self.pix_ny[0] - self.rn))
        ox = np.sort(np.append(self.pix_ox, self.pix_ox[0] - self.ro))
        oy = np.sort(np.append(self.pix_oy, self.pix_oy[0] - self.ro))

        # new figure
        gca().add_patch(Rectangle((nx[0], ny[0]), nx[-1] + abs(self.diff_x / 2.), ny[-1] + abs(self.diff_y / 2.), facecolor='#FFCC00', label='new', edgecolor='k', lw=3))

        # old figure
        gca().add_patch(Rectangle((ox[0], oy[0]), ox[-1], oy[-1], facecolor='#990000', label='old', edgecolor='w', lw=1))

        for i in range(len(nx)):
            plt.plot(np.ones(shape=np.shape(ny)) * nx[i], ny, 'k-', linewidth=3)
        for j in range(len(ny)):
            plt.plot(nx, np.ones(shape=np.shape(nx)) * ny[j], 'k-', linewidth=3)
        for i in range(len(ox)):
            plt.plot(np.ones(shape=np.shape(oy)) * ox[i], oy, 'w-', linewidth=1)
        for j in range(len(oy)):
            plt.plot(ox, np.ones(shape=np.shape(ox)) * oy[j], 'w-', linewidth=1)

        ax.set_xlabel('x (arcsec)')
        ax.set_ylabel('y (arcsec)')
        legend = ax.legend(loc='upper left', bbox_to_anchor=(0.965, 1.13))
        legend.get_frame().set_facecolor('#808080')
        
        if self.len_nx >= self.len_ox: 
            if self.len_nx > 5: zoom_in = 5. / self.len_nx
        if self.len_nx < self.len_ox: 
            if self.len_ox > 5: zoom_in = 5. / self.len_ox
        
        if zoom_in == 1.:
            ax.set_xlim((self.pix_nx[0] - self.rn - 1), (self.pix_nx[-1] + 1))
            ax.set_ylim((self.pix_ny[0] - self.rn - 1), (self.pix_ny[-1] + 1))
        else:
            ax.set_xlim((self.pix_nx[0] - self.rn - 1), self.pix_ox[-1] * zoom_in)
            ax.set_ylim((self.pix_ny[0] - self.rn - 1), self.pix_oy[-1] * zoom_in)

        # old in pixel
        aox = ax.twinx()
        aox.set_ylabel('x (pixel)', color='#990000')
        if zoom_in == 1.:
            aox.set_ylim(-(1 + oy[0] - ny[0]) / self.ro, ((ny[-1] - oy[-1] + 1) / self.ro + self.len_oy))
        else:
            aox.set_ylim(-(1 + oy[0] - ny[0]) / self.ro, self.len_oy * zoom_in)

        # new in pixel
        anx = ax.twiny()
        anx.set_xlabel('x (pixel)', color='#FFCC00')
        if zoom_in == 1.:
            anx.set_xlim(-1 / self.rn, (1 / self.rn + self.len_nx))
        else:
            anx.set_xlim(-1 / self.rn, (ox[0] - nx[0]) / self.rn + self.len_ox * zoom_in * self.ro / self.rn)

        fig.savefig(name + '_process-output_zoom_grid_' + str(self.ro) + '_to_' + str(self.rn) + '_arcsec_pixel.png', dpi=dpi)



def central(array, dx=0.5, dy=0.5):
    #print np.sum(array)
    
    # 1D arrays
    if len(np.shape(array)) < 2:
        raise Exception('central_pixel does not work for 1D arrays.')
        
    # 2D arrays change to 3D arrays
    if len(np.shape(array)) == 2:
        fake_cube = True
        array = np.array([array]).swapaxes(0, 2).swapaxes(0, 1)
    else:
        fake_cube = None
    
    # shift parameter space
    if abs(dx) > 0.5 or abs(dy) > 0.5:
        raise Exception('Shift parameter dx or dy out of range, allowed are values between -0.5 and 0.5. Default is 0.5, which means a pixel shift of 0.5 to the right and down.')
    
    
    # shift types
    shift_type={}
    if dx > 0.: shift_type['x'] = 'x_pos'
    if dx == 0.: shift_type['x'] = 'x_zero'
    if dx < 0.: shift_type['x'] = 'x_neg'
    if dy > 0.: shift_type['y'] = 'y_pos'
    if dy == 0.: shift_type['y'] = 'y_zero'
    if dy < 0.: shift_type['y'] = 'y_neg'
    
    dx = abs(dx)
    dy = abs(dy)
    
    shift = {}
    shift['x_pos'] = [1.-dx,dx]
    shift['x_neg'] = [dx,1.-dx]
    shift['y_pos'] = [1.-dy,dy]
    shift['y_neg'] = [dy,1.-dy]
    
    # move array in x direction
    new_x = np.zeros(shape=(len(array[:,0,0])+1, len(array[0,:,0]), len(array[0,0,:])))
    
    if shift_type['x'] != 'x_zero':
        
        new_x[0,:,:] = array[0,:,:] * shift[shift_type['x']][1]
        new_x[-1,:,:] = array[-1,:,:] * shift[shift_type['x']][0]
        
        for nx in range(len(new_x[:,0,:])-2):
            new_x[nx+1,:,:] = array[nx,:,:] * shift[shift_type['x']][0] + array[nx+1,:,:]* shift[shift_type['x']][1]
    
    if shift_type['x'] == 'x_zero':
        new_x[:-1,:,:] = array
    
    #print np.sum(new_x)    
    
    # move array in y direction
    new_y = np.zeros(shape=(len(new_x[:,0,0]), len(new_x[0,:,0])+1,len(new_x[0,0,:])))
        
    if shift_type['y'] != 'y_zero':
        new_y[:,0] = new_x[:,0,:] * shift[shift_type['y']][1]
        new_y[:,-1] = new_x[:,-1,:] * shift[shift_type['y']][0]
        
        for ny in range(len(new_y[0,:,:])-2):
            new_y[:,ny+1,:] = new_x[:,ny,:] * shift[shift_type['y']][0] + new_x[:,ny+1,:] * shift[shift_type['y']][1]
        
    if shift_type['y'] == 'y_zero':
        new_y[:,:-1,:] = new_x
    
    new_array = new_y
    #print np.sum(new_array)
    
    if fake_cube is True:
        new_array = new_array[:, :, 0]
        
    return new_array
