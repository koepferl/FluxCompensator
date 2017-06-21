import numpy as np
import os
ROOT = os.path.dirname(os.path.abspath(__file__)) + '/'

def set_up_extinction_map():
    '''
    Text description of how to extract a synthetic extinction map from a radiative transfer calculation.
    '''
    
    print 'Set up a Hyperion model, like the module that produced your ideal'
    print 'image. We will make a 1cm image that with a temperature of T=10K. '
    print 'We only care about the dust at that wavelength, this is why initial'
    print '(temperature iterations), imaging (scattering) and '
    print 'raytracing_sources are switched off or set to zero.'
    print ''
    print 'To set this up only change the following things of your copied script:'
    print 'image.set_wavelength_range(1,9999., 10000.)'
    print 'm.set_minimum_temperature(10.)'
    print 'm.set_n_initial_iterations(0)'
    print 'm.set_raytracing(True)'
    print 'm.set_n_photons(imaging=0, raytracing_sources=0, raytracing_dust=1e6)'
    print ''
    print 'To prevent overwriting change the name of the *.rtin and *.rtout file.'
    print ''
    

def extract_extinction_map(SO_cm, SO, dust):
    '''
    Extracting extinction map from radiative transfer calculation.
    
    Parameters
    ----------
    
    SO_cm : SyntheticImage 
        FluxCompensator object of cm observation
    
    SO : SyntheticImage 
        FluxCompensator object of synthetic observation where fieldstars should be added to.
    
    dust : str
        Path and name of dust file.
    
    
    Returns
    -------
    
    A_v : numpy.ndarray
        Optical extinction map.
    
    '''
    
    # S0_cm needs to be at 1cm
    w0 = SO_cm.wav[0]
    if w0 > 10000. or w0 < 9999.: 
        raise Exception('WARNING: Dust extinction image is not at 1cm.')
    
    # check if resolution is the same
    if SO.resolution['rad'] != SO_cm.resolution['rad']:
        raise Exception('WARNING: Extinction image at 1cm and image from pipeline do not have the same resolution.')
    
    # check if units are correct
    if SO_cm.units != 'ergs/cm^2/s/Hz':
        raise Exception('WARNING: Units of SO_cm need to be ergs/cm^2/s/Hz.')
        
    if isinstance(dust,str):
        # load dust properties from hyperion
        from hyperion.dust import SphericalDust
        d = SphericalDust(dust)
        kappa = d.optical_properties.kappa
        chi = d.optical_properties.chi
        wav = d.optical_properties.wav
    
    else:
        # load dust_kappa, dust_chi & dust_wav 
        # from tuple dust={kappa, chi, wav}
        kappa = dust['kappa']
        chi =   dust['chi']
        wav =   dust['wav']
        
        if wav[0] < wav[1]: 
            kappa = kappa[::-1]
            chi = chi[::-1]
            wav = wav[::-1]
    
    # extrapolate kappa at 1cm and chi at all wavelengths of SO
    kappa_cm = np.interp(w0, wav[::-1], kappa[::-1]) # cm^2/g
    chi_wav = np.interp(SO.wav, wav[::-1], chi[::-1]) # cm^2/g
    
    # constants in cgs
    T = 10.
    pi = 3.141592653589793
    h = 6.626068e-27
    k = 1.3806503e-16
    c = 29979245800.0
    
    # Plank function and Jnu at 1cm with T=10K
    nu = c/(w0 * 1e-4)
    Bnu = 2 * h* nu**3/ c**2 * (np.exp(h*nu/k/T)-1)**(-1)    
    Jnu = Bnu * kappa_cm
    
    # surface density * chi = tau
    # surface density = surface brightnes / Jnu
    sr = SO_cm.resolution['rad']**2
    tau = SO_cm.val / Jnu / sr * chi_wav
    
    # extinction at 1 cm
    A_lam = 1.086 * tau
    
    # convert to A_v
    print 'CAUTION: Extinction law from Kim et al. is used.'
    wav_ext, k_lam = np.loadtxt(ROOT + '../database/extinction/extinction_law.txt', unpack=True)
    k_v = np.interp(0.550, wav_ext, k_lam)
    k = np.interp(1., wav_ext, k_lam)
    A_v = A_lam * (k_v / k)
    
    return A_v[:,:,0]


def get_stars_from_database(band, number, distance_range, ground, object_distance, seed=None):
    '''
    Selects magnitudes of fieldstars in IRAC and 2MASS from database 
    and assigns a random distance.
        
    Parameters
    ----------
        
    band : str
        Name of the detector band, here only ``'J_2MASS'``, ``'H_2MASS'``, ``'K_2MASS'``, 
        ``'IRAC1'``, ``'IRAC2'``, ``'IRAC3'``, ``'IRAC4'`` are allowed.
        
    number : int
        Number of selected stars, 289 is the maximum.
    
    distance_range : list
        List of distance of minimum and maximum distance in cm.
    
    ground : str
        ``'foreground'`` and ``'background'`` fieldstars to the modeled objects are possible.
        When ``'foreground'`` is enabled, e.g. for galaxies only stars distances in the foreground are
        produced. When ``'background'`` is enabled, 90% background and 10% foreground distances are 
        selected.
        
    object_distance : float
        Distance of modeled object in cm.
    
    seed : int, None
        Seed for random number generator. Set to a certain integer to reproduce the result.
        
        
    Returns
    -------
        
    mag : numpy.ndarray
        Magnitude of ``number`` stars in observed with ``band``.
        
    distance_stars : numpy.ndarray
        Distance of stars within ``distance_range``.
        
    '''
    
    # read stellar data file
    t = np.loadtxt(ROOT + '../database/fieldstars/field_stars.txt')
    
    # ensure that random numbers are the same every time
    if seed is not None:
        np.random.seed(seed)
    
    # maximum number of stars    
    if number > 289: 
        raise Exception('Only 289 are in the database file.')
        
    # avalable detectors
    bands = ['J_2MASS', 'H_2MASS', 'K_2MASS', 'IRAC1', 'IRAC2', 'IRAC3', 'IRAC4']
    if band not in bands:
        raise Exception('Only 2MASS and IRAC field stars are in the database.')
    
    # asign magnitudes
    mag_bands = {}
    for b in range(len(bands)):
        mag_bands[bands[b]] = t[:,b+2]
    
    # select magnitudes
    mag = mag_bands[band]
    np.random.shuffle(mag)
    mag = mag[:number]
    
    # forground or background stars  
    if ground not in ['foreground', 'background'] and isinstance(ground,float) is False:
        raise Exception('Stars added need to be foreground or background stars or fraction of distribution must be given.')
    if ground == 'foreground': frac = 1.
    if ground == 'background': frac = 0.
    else: frac = ground
    
    # select distances
    foreground_stars = np.random.uniform(distance_range[0], object_distance, int(number * frac))
    background_stars = np.random.uniform(distance_range[1], object_distance, number - int(number * frac))
    distance_stars = np.append(foreground_stars, background_stars)
        
    return mag, distance_stars
