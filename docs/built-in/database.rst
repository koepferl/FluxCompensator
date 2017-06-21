.. _label_database:

========
Database
========

The database is structured by ``missions``. Currently in contains

* 2MASS (called by ``twomass``)
* IRAS (called by ``iras``)
* WISE (called by ``wise``)
* SPITZER (called by ``spitzer``)
	* IRAC (called by ``irac``)
	* MIPS (called by ``mips``)
* HERSCHEL (called by ``herschel``)
	* PACS (called by ``pacs``)
	* SPIRE (called by ``spire``)

For the telescopes one can pass the ``diameter`` in cm with::

    from fluxcompensator.database.missions import *
    wise.diameter
    40.0
    
    from fluxcompensator.database.missions import wise
    wise.diameter
    40.0


Furthermore, the database in FluxCompensator contains the following objects always accessible for the user:

* :class:`~fluxcompensator.filter.Filter` (e.g. ``IRAC4``)
* :class:`~fluxcompensator.utils.one_filter.OneFilter` (e.g. ``IRAC4_PLOT``)
* :class:`~fluxcompensator.psf.FilePSF` (e.g. ``IRAC4_PSF``)

Additionally, detector relevant information such as

* PSF file resolution (e.g. ``IRAC4_PSF_RESOLUTION``) in arcsec per pixel (without oversampeling)
* zero-magnitude flux (e.g. ``IRAC4_ZERO``) in Jy

are passed, which might be needed for `flux-magnitude conversion <file:///Users/koepferl/simulated-observations-package/docs/_build/html/postprocessing/mag.html>`_ or convolving with the :ref:`label_GaussianPSF`::

    from fluxcompensator.database.missions import *
    IRAC4_PSF_RESOLUTION
    1.22
    IRAC4_ZERO
    64.9
    
    from fluxcompensator.database.missions import irac
    irac.IRAC4_PSF_RESOLUTION
    1.22
    irac.IRAC4_ZERO
    64.9
    
    from fluxcompensator.database.missions import spitzer
    spitzer.IRAC4_PSF_RESOLUTION
    1.22
    spitzer.IRAC4_ZERO
    64.9


Access Objects from Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A known object (e.g. ``IRAC4_FILTER``) can be passed in a script::

    import fluxcompensator.database.missions as filters
    filter_object = getattr(filters, 'IRAC4_FILTER')
    
    import fluxcompensator.database.missions as PSFs
    psf_object = getattr(PSFs, 'IRAC4_PSF')
    
    import fluxcompensator.database.missions as Plots
    plot_filter = getattr(Plots, 'IRAC4_FILTER_PLOT')

.. note:: That currently only for the missions in ``herschel`` and ``spitzer`` PSF files are available. Only if the PSF file is available the PSF file resolution (e.g. ``IRAC4_PSF_RESOLUTION``) is stored in the database. 


Access Attributes from Objects in Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Especially for the :class:`~fluxcompensator.filter.Filter` object it might be interesting to extract attributes from therein::

    from fluxcompensator.database.missions import *
    
    # name of Filter object
    IRAC4_FILTER.name
    'IRAC4'
    
    # central wavelength in microns
    IRAC4_FILTER.waf_0
    7.872
    
    # power law exponent alpha
    IRAC4_FILTER.alpha
    1.0
    
    # power law exponent beta
    IRAC4_FILTER.beta
    0.0	

Contents of Database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 1
   
   ../api/fluxcompensator.database.missions.twomass.rst
   ../api/fluxcompensator.database.missions.iras.rst
   ../api/fluxcompensator.database.missions.wise.rst
   ../api/fluxcompensator.database.missions.spitzer.rst
   ../api/fluxcompensator.database.missions.irac.rst
   ../api/fluxcompensator.database.missions.mips.rst
   ../api/fluxcompensator.database.missions.herschel.rst
   ../api/fluxcompensator.database.missions.pacs.rst
   ../api/fluxcompensator.database.missions.spire.rst
