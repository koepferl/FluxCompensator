.. _label_FunctionPSF:

================
PSF Function
================

It might be, that you do not have an actual PSF file of the PSF you want to "observe" your synthetic image with. In this case you could convolve with a 2D function::

	function(X, Y, wavelength)

To construct your own function in Python you can simply write something like this::

    def my_psf(X,Y,wavelength):
        # resolution in rad/pixel
        resolution = FC_object.FOV[0] / FC_object.distance / FC_object.pixel[0]
        
        # telescope diameter in cm
        D_telescope = 350.
        
        # standard deviation in pixel
        sigma = 0.44 * wavelength / D_telescope / resolution
        
        Z = np.exp(-(X**2 + Y**2)/(2 * sigma**2))
        
        return Z


.. note:: You can call attributes of the :ref:`FC_objects <label_objects>` like ``distance``, ``FOV``, etc. For more information see `FluxCompensator objects <file:///Users/koepferl/simulated-observations-package/docs/_build/html/FCobjects/objects.html>`_.

Now to construct the ``psf_object`` add to you script::

	from fluxcompensator.psf import FunctionPSF
	
	# create PSF from a function
	psf_object = FunctionPSF(psf_function=my_psf, width=32)

``width`` is the size in pixel of the created 2D PSF function grid . ``X`` and ``Y`` in the definition of the ``my_psf`` represent the x and y coordinates in an array centered around 0 to ``width/2``. ``X`` and ``Y`` are passed by :class:`~fluxcompensator.psf.FunctionPSF`. ``wavelength`` is passed by :meth:`fluxcompensator.cube.SyntheticCube.convolve_psf` or :meth:`fluxcompensator.image.SyntheticImage.convolve_psf`.

For further information see: 

* :class:`fluxcompensator.psf.FunctionPSF`
