.. _label_add_to_real_obs:

===============================================================
Blending Realistic Synthetic Observations in a Real Background
===============================================================

Once your 2D :ref:`FC_object <label_objects>` is a member of :class:`~fluxcompensator.image.SyntheticImage` (due to the :meth:`fluxcompensator.cube.SyntheticCube.convolve_filter` in the pipeline), you can combine it with a real observation (e.g. ``file.fits``). To add your realistic syntetic observation to a real observation and store it in a fits file, add to your script:: 

    # combine FC_object.val (2D) with a real observation
    FC_object.add_to_observation(fits_file='file.fits', name='add_model',   
                                 position_pix=[16., 32.], position_world=None,
                                 zero_edges=None)

The output file ``add_model.fits`` will be stored in the same directory where the script runs. The parametes

* ``position_pix`` : ``list`` with [x, y] postion of the model ``FC_object.val`` center in pixel coordinates in the observations ``file.fits``
* ``position_world`` : ``list`` with [x, y ] postion of the model ``FC_object.val`` center in world coordinates in the observations ``file.fits``. Only needed if ``position_pix = None``.
* ``zero_edges`` : If ``True`` edges of model ``FC_object.val`` are normalized to zero at the rim before adding to the observations ``file.fits``.