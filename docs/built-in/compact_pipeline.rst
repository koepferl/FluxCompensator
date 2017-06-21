.. _label_compact_pipeline:

==================
Compact Pipelines
==================

The FluxCompensator has several pre-constructed pipelines built-in, constructed with the class 

* :class:`fluxcompensator.database.compact_pipeline.CompactPipeline`. 

:class:`~fluxcompensator.database.compact_pipeline.CompactPipeline` objects harbor the filter & PSF information of several observational surveys, e.g. GLIMPSE, MIPSGal, HiGAL.

Built-in CompactPipeline objects
--------------------------------

* 2MASS survey :class:`~fluxcompensator.database.compact_pipeline.TWOMASSSURVEY_J`, :class:`~fluxcompensator.database.compact_pipeline.TWOMASSSURVEY_H`, :class:`~fluxcompensator.database.compact_pipeline.TWOMASSSURVEY_K`
* WISE survey :class:`~fluxcompensator.database.compact_pipeline.WISESURVEY_WISE1`, :class:`~fluxcompensator.database.compact_pipeline.WISESURVEY_WISE2`, :class:`~fluxcompensator.database.compact_pipeline.WISESURVEY_WISE3`, :class:`~fluxcompensator.database.compact_pipeline.WISESURVEY_WISE4`      
* GLIMPSE :class:`~fluxcompensator.database.compact_pipeline.GLIMPSE_IRAC1`, :class:`~fluxcompensator.database.compact_pipeline.GLIMPSE_IRAC2`, :class:`~fluxcompensator.database.compact_pipeline.GLIMPSE_IRAC3`, :class:`~fluxcompensator.database.compact_pipeline.GLIMPSE_IRAC4`      
* MIPSGal :class:`~fluxcompensator.database.compact_pipeline.MIPSGA_MIPS1`, :class:`~fluxcompensator.database.compact_pipeline.MIPSGAL_MIPS2`
* HiGAL :class:`~fluxcompensator.database.compact_pipeline.HIGAL_PACS1`, :class:`~fluxcompensator.database.compact_pipeline.HIGAL_PACS2`, :class:`~fluxcompensator.database.compact_pipeline.HIGAL_SPIRE1`, :class:`~fluxcompensator.database.compact_pipeline.HIGAL_SPIRE2`, :class:`~fluxcompensator.database.compact_pipeline.HIGAL_SPIRE3`

Building CompactPipelines
-------------------------

You can build a :class:`~fluxcompensator.database.compact_pipeline.CompactPipeline` in the following way::

    from fluxcompensator.database.compact_pipeline import CompactPipeline
    
    # setting up some filter)object and some psf_object
    filter_object = Filter(...)     # or Filter object from database
    psf_object = GaussianPSF(...)   # or PSF object from database, FilePSF(...), FunctionPSF(...)
    
    # initiate CompactPipeline
    compact_pipeline = CompactPipeline(PSF_object=psf_object, 
                       PSF_resolution=3., 
                       filter_object=filter_object)

.. note:: Here :class:`~fluxcompensator.psf.GaussianPSF` and :class:`~fluxcompensator.filter.Filter` are place holders some comparable ``filter_object`` or ``psf_object``. They can be replaced by the respected classes in the database or by self-constructed PSF and filter classes. To set up there objects see documentation of the :ref:`Filter Convolution <label_filter>` and the :ref:`PSF convolution <label_psf>`.

Pixel resolution can be given manually or extracted from the :ref:`label_database`. See the documentation of :ref:`label_database` for further instructions.


