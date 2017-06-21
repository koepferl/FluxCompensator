.. _label_mag:

=====================
Convert to Magnitudes
=====================

It is possible for any photometric scalar :class:`fluxcompensator.flux.SyntheticFlux` :ref:`FC_object <label_objects>` created with 

* :meth:`fluxcompensator.cube.SyntheticCube.get_total_val`
* :meth:`fluxcompensator.image.SyntheticImage.get_total_val`
* :meth:`fluxcompensator.sed.SyntheticSED.convolve_filter`
* :meth:`fluxcompensator.sed.SyntheticSED.get_total_val`

to convert to magnitudes with :meth:`fluxcompensator.flux.SyntheticFlux.magnitudes`, if the zero-point flux ``zero_point`` in Jy is available.

To do so (e.g. for IRAC4), add to your script::

    # convert to magnitudes
    FC_object_mag = FC_object.magnitudes(zero_point=64.9)

.. note:: The ``zero_point`` for common detectors can be accessed from the :ref:`label_database`.
     