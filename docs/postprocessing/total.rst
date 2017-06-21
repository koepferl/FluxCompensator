.. _label_total:

===========================
Total Collapse Single Value
===========================

At some steps within the pipeline it might be interesting to collapse the :ref:`FC_objects <label_objects>`

* :class:`~fluxcompensator.cube.SyntheticCube`
* :class:`~fluxcompensator.sed.SyntheticSED`

with dimension

* 3D (wav, x, y)
* 1D (wav)

respectively, to the total photometric flux within some boundaries. This can be done with 

* :meth:`fluxcompensator.cube.SyntheticCube.get_total_val`
* :meth:`fluxcompensator.sed.SyntheticSED.get_total_val`. 

.. note:: This will not replace the :ref:`Filter Convolution <label_filter>` but it might help to see if the pipeline is working as you expected.

Especially after the :ref:`Filter Convolution <label_filter>` you might want to get the total photometric flux of a :ref:`FC_object <label_objects>` 

* :class:`~fluxcompensator.image.SyntheticImage` 

of dimension

* 2D (x, y)

with 

* :meth:`fluxcompensator.image.SyntheticImage.get_total_val`

without using boundaries.


Within Boundaries (3D & 1D to 0D)
----------------------------------

To do so, at any time you can collapse the ``val`` of :ref:`FC_objects <label_objects>` (3D or 1D) to a total photometric flux (0D) within the boundaries ``wav_1`` and ``wav_2`` in microns, with::

    # collapse FC_object.val (3D or 1D) within 30 and 60 microns      
    FC_object_tot = FC_object.get_total_val(wav_1=30., wav_2=60.)

.. note:: In order for this to work the :ref:`FC_objects <label_objects>` must have more than one entry in ``FC_object.wav``. This is why for :class:`~fluxcompensator.image.SyntheticImage` the method is different. For :class:`~fluxcompensator.flux.SyntheticFlux` this method makes no sense since the object is already an photometric flux of 0D.


Without Boundaries (2D to 0D)
-----------------------------

For :ref:`FC_objects <label_objects>` which are members of :class:`~fluxcompensator.image.SyntheticImage` (2D) only one ``FC_object.wav`` entry remains. Therefore, no boundaries need to be applied in the method :meth:`fluxcompensator.image.SyntheticImage.get_total_val`. Simply type::

    # collapse FC_object.val (2D)
    FC_object_tot = FC_object.get_total_val()

See Section :ref:`label_output` to create plots.
    
