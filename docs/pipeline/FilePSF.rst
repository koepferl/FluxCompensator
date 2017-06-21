.. _label_FilePSF:

=============
PSF File
=============

The most realistic way of creating a realistic synthetic observation is to convolve with a telescope specific PSF file provided by an observatory or self-constructed a PSF from real observations.

.. note:: If you use PSF files, you should first adjust the resolution of the image to the resolution of the mimiced detector. You can do this with :ref:`change_resolution <label_resolution>`.


Read Own File
^^^^^^^^^^^^^

If you have an actual PSF file ``my_psf.fits``, with the PSF from your telescope you want to "observe" your synthetic image with, you can read this file in the following way and construct a ``psf_object``::

	from fluxcompensator.psf import FilePSF
	
	# create PSF from my own file
	psf_object = FilePSF(psf_file='my_psf.fits', oversampled=4)
	
If there is no oversampeling in you file, just put ``oversampled=1``. 

Use the FluxCompensator's Built-in :ref:`label_database`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It might be, that you do not have an actual PSF file at hand. In this case you could convolve with a PSF provided by the FluxCompensator :ref:`label_database`. Here we have the publicly available PSF files of the telescopes, which are commonly used. For a reference list see the Appendix of Koepferl & Robitaille (subm. to ApJ). This :ref:`label_database` is far from being complete. If you have additional available files, please let me know. 

Currently PSFs from IRAC, MIPS, PACS and SPIRE are available here (e.g. MIPS for different temperatures). To construct the ``psf_object`` (e.g. of PACS1) add to your script::

    import fluxcompensator.database.missions as PSFs
    	
    # call object from the psf database
    psf_object = getattr(PSFs, 'PACS1_PSF')

The PSF objects in the :ref:`label_database` (e.g. ``PACS1_PSF``) can be called by using ``getattr`` and the ``str`` of the ``psf_object``.

Possible names of the attributes are::

	PSF	   oversampled
             
    IRAC1_PSF		 4
    IRAC2_PSF		 4
    IRAC3_PSF		 4
    IRAC4_PSF		 4

    MIPS1_PSF_10K	 5
    MIPS1_PSF_25K	 5
    MIPS1_PSF_50K	 5
    MIPS1_PSF_75K	 5
    MIPS1_PSF_100K	 5
    MIPS1_PSF_3000K	 5

    MIPS2_PSF_10K	 5
    MIPS2_PSF_25K	 5
    MIPS2_PSF_50K	 5
    MIPS2_PSF_75K	 5
    MIPS2_PSF_100K	 5
    MIPS2_PSF_3000K	 5

    MIPS3_PSF_10K	 5
    MIPS3_PSF_25K	 5
    MIPS3_PSF_50K	 5
    MIPS3_PSF_75K	 5
    MIPS3_PSF_100K	 5
    MIPS3_PSF_3000K	 5

    PACS1_PSF		10
    PACS2_PSF		10
    PACS3_PSF		10

    SPIRE1_PSF		10
    SPIRE2_PSF		10
    SPIRE3_PSF		10

For further information see:

* :class:`fluxcompensator.psf.FilePSF`