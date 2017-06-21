.. _label_fits:

====================
Creating Fits Files
====================

At any point in the :ref:`label_virtual_pipeline` you can store the ``val`` with

* 3D (x, y, wav) :class:`~fluxcompensator.cube.SyntheticCube`
* 2D (x, y) :class:`~fluxcompensator.image.SyntheticImage` 

dimension of the different :ref:`FC_objects <label_objects>`, respectivly, in a fits file (e.g. ``output.fits``). 

For a 3D :class:`~fluxcompensator.cube.SyntheticCube`, add to you script::

    from astropy.io import fits
    
    # store FC_object.val in fits file
    fits.writeto('output.fits', FC_object.val.swapaxes(0, 2).swapaxes(1, 2), 
                 clobber=True)
	
In this case you will find the file ``output.fits`` in the same directory as ``example.py``. 

 If you extend the example described in :ref:`label_cube`, the file will be exactly the same like the :download:`here <../media/output.fits>`.


For a 2D :class:`~fluxcompensator.image.SyntheticImage`, add to you script::

    from astropy.io import fits
    
    # store FC_object.val in fits file
    fits.writeto('output.fits', FC_object.val, clobber=True)
    