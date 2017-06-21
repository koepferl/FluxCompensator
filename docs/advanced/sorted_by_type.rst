.. _label_sorted_by_type:

Application list of API
------------------------

Set-up
^^^^^^
* :class:`fluxcompensator.cube.SyntheticCube`
* :class:`fluxcompensator.image.SyntheticImage`
* :class:`fluxcompensator.sed.SyntheticSED`
* :class:`fluxcompensator.flux.SyntheticFlux`

Extinction
^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.extinction`
* :meth:`fluxcompensator.image.SyntheticImage.extinction`
* :meth:`fluxcompensator.sed.SyntheticSED.extinction`
* :meth:`fluxcompensator.flux.SyntheticFlux.extinction`

Resolution
^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.change_resolution`
* :meth:`fluxcompensator.image.SyntheticImage.change_resolution`

PSF convolution
^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.convolve_psf`
* :meth:`fluxcompensator.image.SyntheticImage.convolve_psf`

Filter convolution
^^^^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.convolve_filter`
* :meth:`fluxcompensator.sed.SyntheticSED.convolve_filter`

Noise
^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.add_noise`
* :meth:`fluxcompensator.image.SyntheticImage.add_noise`

Closer to real observations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
* :class:`fluxcompensator.interface.Interface2FITS`
* :meth:`fluxcompensator.interface.Interface2FITS.save2fits`
* :meth:`fluxcompensator.interface.Interface2FITS.add2observation`
* :class:`fluxcompensator.database.compact_pipeline.CompactPipeline`
* :meth:`fluxcompensator.image.SyntheticImage.add_to_observation`
* :meth:`fluxcompensator.image.SyntheticImage.add_field_stars`

Collapse to SED
^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.get_rough_sed`

Collapse to photometric flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.get_total_val`
* :meth:`fluxcompensator.image.SyntheticImage.get_total_val`
* :meth:`fluxcompensator.sed.SyntheticSED.get_total_val`

Convert to magnitudes
^^^^^^^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.flux.SyntheticFlux.magnitudes`

Plot image
^^^^^^^^^^
* :meth:`fluxcompensator.cube.SyntheticCube.plot_image`
* :meth:`fluxcompensator.image.SyntheticImage.plot_image`

Plot SED with total val and filter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.sed.SyntheticSED.plot_sed_multi_filter`
* :meth:`fluxcompensator.flux.SyntheticFlux.plot_sed_filter`

Additional classes
^^^^^^^^^^^^^^^^^^
* :class:`fluxcompensator.filter.Filter`
* :class:`fluxcompensator.psf.FilePSF`
* :class:`fluxcompensator.psf.FunctionPSF`
* :class:`fluxcompensator.psf.GaussianPSF`
* :class:`fluxcompensator.utils.pseudo.Pseudo`
* :class:`fluxcompensator.utils.one_filter.OneFilter`
* :class:`fluxcompensator.utils.plot.MakePlots`
* :class:`fluxcompensator.utils.plot_filters.PlotFilters`
* :class:`fluxcompensator.utils.resolution.ConservingZoom`
* :class:`fluxcompensator.utils.units.ConvertUnits`

Additional methods
^^^^^^^^^^^^^^^^^^
* :meth:`fluxcompensator.utils.tools.properties`
* :meth:`fluxcompensator.utils.tools.grid_units`
* :meth:`fluxcompensator.utils.tools.get_slices`
* :meth:`fluxcompensator.utils.tools.average_collapse`
* :meth:`fluxcompensator.utils.tools.central_wav`
* :meth:`fluxcompensator.utils.tools.get_projection`
* :meth:`fluxcompensator.utils.tools.where_is`
* :meth:`fluxcompensator.utils.tools.where_is_1D`