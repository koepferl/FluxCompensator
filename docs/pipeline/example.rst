.. _label_example:

=============================
Example of a Simple Pipeline
=============================

Here you can see a small pipeline, which is reading in, deredding, change the resolution, convolving with PSF and filter and adding noise::

    import numpy as np
    
    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc
    
    from fluxcompensator.cube import *
    from fluxcompensator.psf import *
    from fluxcompensator.utils.resolution import *
    
    # read in from `Hyperion <http://www.hyperion-rt.org>`_
    m = ModelOutput('hyperion_output.rtout')
    array = m.get_image(group=0, inclination=0, distance=300*pc,
                        units='ergs/cm^2/s')
    
    # initial FluxCompensator array
    c = SyntheticCube(input_array=array, unit_out='ergs/cm^2/s',
                      name='test_cube')
    
    # dered with provided extinction law
    ext = c.extinction(A_v=20.)
    
    # change resolution to 10-times of the initial
    zoom = ext.change_resolution(new_resolution=10*ext.resolution['arcsec'],
                                 grid_plot=True)
    
    import fluxcompensator.database.missions as PSFs
        
    # call object from the psf database
    psf_object = getattr(PSFs, 'PACS1_PSF')
        
    # convolve with PSF
    psf = zoom.convolve_psf(psf_object)
    
    import fluxcompensator.database.missions as filters
        
    # call object from the filter database
    filter_input = getattr(filters, 'PACS1_FILTER')
    
    # convolve with filter
    filtered = psf.convolve_filter(filter_input, plot_rebin=None,
                                   plot_rebin_dpi=None)
    
    # add noise
    noise = filtered.add_noise(mu_noise=0, sigma_noise=5e-15, seed=2, diagnostics=None)

See Sections :ref:`label_post` and :ref:`label_output` to create SEDs, photometric fluxes and further outputs.
	
