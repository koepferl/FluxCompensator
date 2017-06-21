# Master plan 
 
import numpy as np 
 
from astropy.io import fits 
from astropy.wcs import WCS 
 
from fluxcompensator.cube import * 
from fluxcompensator.image import * 
from fluxcompensator.filter import * 
from fluxcompensator.psf import * 
from fluxcompensator.utils.resolution import * 
 
class Interface2FITS(object): 
    ''' 
    Reads information from real observation and radiative transfer model and creates  
    a synthetic observation of the model with the specifications of the observations. 
     
    Parameters 
    ---------- 
     
    obs : str 
        Stores the observation of a certain telescope (e.g 'SF-region_IRAC4.fits') 
        The file needs to be in fits format and WCS.  
     
    model : ModelOutput.get_image 
        Radiative transfer model from HYPERION. 
     
    compact_pipeline : CompactPipeline object 
        Stores filter and PSF information.   
     
    exposure : float 
        Exposure time of telescope. 
         
    A_v : float 
        Optical extinction coefficient.  
        Default is ``A_v=0``, which disables the extinction. 
        Hence intrinsic flux is observed. 
         
     
    Returns 
    ------- 
     
    image : SyntheticImage 
     
    ''' 
 
    def __init__(self, obs, model, compact_pipeline, exposure, A_v=0): 
 
        # information from real observation 
        self.obs_file = obs 
        self.wcs_obs = WCS(self.obs_file) 
        hdulist = fits.open(self.obs_file) 
        hdu = hdulist[0] 
        self.image_obs = hdu.data 
        self.bunit = hdu.header['BUNIT'] 
        self.obs_res = abs(hdu.header['CDELT1'] * 3600.) 
 
        if abs(hdu.header['CDELT1']) != abs(hdu.header['CDELT2']): 
            print 'WARNING: resolution in x direction is different than resolution in y direction.' 
            print '         This code will only use squared grids and the model is rescaled with ' 
            print '         res_x in both directions.' 
 
 
        # information from radiative transfer model 
        self.model = model 
 
        # information from chosen compact_pipeline 
        self.filter_input = compact_pipeline.filter_object 
        self.psf_original = np.array([compact_pipeline.PSF_object.psf]).swapaxes(0, 2).swapaxes(0, 1) 
        self.psf_original_res = compact_pipeline.PSF_resolution 
 
        # information about extinction 
        self.A_v = A_v 
 
 
        ########## 
        # simobs 
        ########## 
 
        # initial FluxCompensator array 
        c = SyntheticCube(input_array=self.model, unit_out=self.bunit, name='interface') 
 
        # dered with provided extinction law 
        ext = c.extinction(A_v=self.A_v) 
 
        # convolve with filter 
        filtered = ext.convolve_filter(self.filter_input) 
 
        if self.psf_original_res < filtered.resolution['arcsec']: 
            # change resolution of PSF 
            psf_original_zoom = ConservingZoom(array=self.psf_original, 
                                               initial_resolution=self.psf_original_res, 
                                               new_resolution=filtered.resolution['arcsec']) 
            psf_zoom = psf_original_zoom.zoom()[:,:,0] 
            fits.writeto('psf_zoom.fits', psf_zoom, clobber=True) 
 
            psf_object = FilePSF(psf_file='psf_zoom.fits') 
 
        if self.psf_original_res > filtered.resolution['arcsec']: 
            # change resolution of image 
            zoom_before_psf = filtered.change_resolution(new_resolution=self.psf_original_res) 
 
            psf_object = compact_pipeline.PSF_object 
 
 
        # convolve with PSF 
        psf = filtered.convolve_psf(psf_object) 
 
        # change resolution to obs.res 
        if psf.resolution['arcsec'] > self.obs_res: 
            raise Exception('Resolution from observation is smaller than resolution of PSF. Use different PSF to convolve.') 
        zoom = psf.change_resolution(new_resolution=self.obs_res) 
 
        # information about noise from exposure 
        mean = np.mean(zoom.val) 
        self.sigma_exp = exposure#0.1 * mean/sqrt(exposure) 
        print 'WARNING: sigma is curretly set to exposure value. Will be fixed in later versions.' 
 
 
        # add noise 
        noise = zoom.add_noise(mu_noise=0, sigma_noise=self.sigma_exp) 
 
 
        ############### 
        # header 
        ############### 
 
        # Initialize WCS information 
        wcs_model = WCS(naxis=2) 
 
        # Set the x pixel scale (in deg/pix) 
        scale_x = self.obs_res / 3600. 
        wcs_model.wcs.cdelt = [-scale_x, scale_x] 
 
        # Set the coordinate system 
        wcs_model.wcs.ctype = ['GLON-CAR', 'GLAT-CAR'] 
 
        # Use the center of the image as projection center 
        if noise.val.shape[1]%2 == 0: 
            center_y = noise.val.shape[1] / 2. + 0.5 
        else: 
            center_y = noise.val.shape[1] / 2. 
 
        if noise.val.shape[0]%2 == 0: 
            center_x = noise.val.shape[0] / 2. + 0.5 
        else: 
            center_x = noise.val.shape[0] / 2. 
        wcs_model.wcs.crpix = [center_y, center_x] 
 
        # And produce a FITS header 
        self.header_model = wcs_model.to_header() 
 
 
        ######### 
        # debug 
        ######### 
 
        #print np.shape(c.val) 
        #print np.shape(ext.val) 
        #print np.shape(filtered.val) 
        #print '_'*12 
        #print self.psf_original_res 
        #print filtered.resolution['arcsec'] 
        #print np.shape(self.psf_original) 
        #if self.psf_original_res < filtered.resolution['arcsec']: 
        #    print np.shape(psf_zoom) 
        #print np.shape(psf_object.psf) 
        #print '_'*12 
        #print np.shape(psf.val), psf.resolution 
        #print np.shape(zoom.val), zoom.resolution 
        #print np.shape(noise.val) 
 
 
        ################## 
        # Simulated image 
        ################## 
 
        self.val = np.array(deepcopy(noise.val)) 
        self.wav = np.array(deepcopy(noise.wav)) 
        self.units = noise.units 
        self.distance = noise.distance 
        self.x_max = noise.x_max 
        self.x_min = noise.x_min 
        self.y_max = noise.y_max 
        self.y_min = noise.y_min 
        self.lon_min = noise.lon_min 
        self.lon_max = noise.lon_max 
        self.lat_min = noise.lat_min 
        self.lat_max = noise.lat_max 
        self.pix_area_sr = noise.pix_area_sr 
        self.res = noise.resolution 
 
        self.unit_in = noise.unit_in 
        self.unit_out = noise.unit_out 
        self.grid_unit = noise.grid_unit 
        self.grid_unit_name = noise.grid_unit_name 
        self.FOV = deepcopy(noise.FOV) 
        self.name = noise.name 
        self.stage = noise.stage 
        self.log = deepcopy(noise.log) 
        self.filter = deepcopy(noise.filter) 
 
 
    def save2fits(self, name): 
        ''' 
        Saves synthetic observation of radiative transfer model in a fits file. 
        Updated information is stored in a header of WCS format. 
         
        Parameters 
        ---------- 
         
        name : str 
            Name of the output fits file. 
         
         
        Returns 
        ------- 
         
        name.fits : Synthetic observation from radiative transfer model 
         
        ''' 
 
        # store model2simobs.val in name.fits 
        fits.writeto(name + '.fits', self.val, header=self.header_model, clobber=True) 
 
 
    def add2observation(self, name, position_pix=None, position_world=None, zero_edges=None): 
        ''' 
        Blends the modeled realistic observation to the real observation. 
         
        Parameters 
        ---------- 
         
        name : str 
            Name of the output fits file. 
         
        position_pix : list, ``None`` 
            Center position of the model in observation pixel coordinates. 
            Default is ``None``.   
         
        position_world : list, ``None`` 
            Center position of the model in observation world coordinates. 
            Default is ``None``.   
         
        zero_edges : ``True``, ``None`` 
            If ``True``, edges of model are normalized to zero. 
            Default is ``None``.   
         
         
        Returns 
        ------- 
         
        name.fits : Combination of synthetic observation from radiative transfer model and real observation 
         
        ''' 
 
        if position_world is None and position_pix is None: 
            raise Exception('WARNING: Position of model center needs to be given either in world or pixel coordinates.') 
 
        if position_pix is not None: 
            pos = position_pix 
            p_x_pos, p_y_pos = pos[0], pos[1] 
        else: 
            pos = position_world 
            p_x_pos, p_y_pos = self.wcs_obs.wcs_world2pix(pos[0], pos[1], 1) 
 
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
        hdulist = fits.open(self.obs_file) 
        hdu = hdulist[0]
 
        file_header = hdu.header 
 
        print 'resolution: real obs                       ', np.abs(file_header['CDELT1'] * 3600)
        print 'resolution: realistic synthetic observation', self.res['arcsec'] 
 
        if np.allclose(np.abs(file_header['CDELT1'] * 3600), self.res['arcsec']) is not True: 
            raise Exception('WARNING: make sure that resolution of observation and model are the same! E. g. change resolution of FC_object first.') 
 
        image = hdu.data 
 
        # add model to observation 
        image[start_y:stop_y, start_x:stop_x] = image[start_y:stop_y, start_x:stop_x] + model 
 
        ## create header 
        file_header['SYNOBS'] = 'FluxCompensator realistic synthetic observation'
        file_header['SYNOBS_X'] = (x_round, 'center pixel coordinates')
        file_header['SYNOBS_Y'] = (y_round, 'center pixel coordinates')
        #file_header['SYNOBSlog'] = str(self.log)
        #new_header = self.header_model.copy() 
        #new_header['CRPIX1'] = hdu.header['CRPIX1'] 
        #new_header['CRPIX2'] = hdu.header['CRPIX2'] 
        #new_header['CTYPE1'] = hdu.header['CTYPE1'] 
        #new_header['CTYPE2'] = hdu.header['CTYPE2'] 
        #new_header['CRVAL1'] = hdu.header['CRVAL1'] 
        #new_header['CRVAL2'] = hdu.header['CRVAL2'] 
 
        # store to name.fits file 
        fits.writeto(name + '.fits', image, header=file_header, clobber=True) 
 