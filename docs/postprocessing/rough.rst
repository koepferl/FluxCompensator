.. _label_rough:

========================
Collapse to a Rough SED
========================

Starting from a `3D ModelOutput <http://docs.hyperion-rt.org/en/stable/postprocessing/extracting_observables.html>`_ it might be interesting at some steps in the pipeline to collapse the 3D ``val`` of a :class:`~fluxcompensator.cube.SyntheticCube` (see :ref:`label_cube`) and extract a rough SED with :meth:`fluxcompensator.cube.SyntheticCube.get_rough_sed`. 

.. note:: This will not replace the `1D ModelOutput <http://docs.hyperion-rt.org/en/stable/postprocessing/extracting_observables.html>`_ but it might help to see if the pipeline is working.

To do so, at any time (as long the :ref:`FC_object <label_objects>` is still a :class:`~fluxcompensator.cube.SyntheticCube`) you can collapse it to a SED, with::

    # collapse 3D cube to rough SED
    rough = FC_object.get_rough_sed()

See Section :ref:`label_output` to create plots.