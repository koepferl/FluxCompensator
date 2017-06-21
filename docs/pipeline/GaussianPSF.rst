.. _label_GaussianPSF:

==============
Gaussian PSF
==============

It might be, that you do not have an actual PSF file of the PSF you want to "observe" your synthetic image with. In this case you could convolve with a Gaussian-shaped PSF. To first order an Airy disc can be approximated with a  Gaussian profile.

The ``psf_object`` is constructed in the following way::

    from fluxcompensator.psf import GaussianPSF
    
    # create Gaussian PSF
    psf_object = GaussianPSF(diameter=350.)

The ``diameter`` of the telescope is passed in cm. The standard deviation of the Gaussian ``sigma`` in units of pixel is calculate internally::

    sigma_psf = 0.44 * wav * 10.**(-4) / diameter / resolution['rad']

where ``wav`` and ``resolution['rad']`` are passed from the object :class:`~fluxcompensator.cube.SyntheticCube` or :class:`~fluxcompensator.image.SyntheticImage`.

For further information see:

* :class:`fluxcompensator.psf.GaussianPSF`
