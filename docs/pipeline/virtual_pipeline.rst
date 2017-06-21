.. _label_virtual_pipeline:

================
Virtual Pipeline
================

Now that you have successfully initialized a FluxCompensator object it is time to transfrom it into a "realistic" synthetic observation. The following tools can be used to make the synthetic observation as "realistic" as possible. 

.. toctree::
   :maxdepth: 1
   
   extinction.rst
   resolution.rst
   psf.rst
   filter.rst
   noise.rst
   example.rst
   

.. note:: Some of the actions (e.g. :ref:`convolve_filter <label_filter>`, :ref:`get_total_val <label_total>`, :ref:`get_rough_sed <label_rough>`) change the dimension of the initial input. Therefore, not all operations can be interchanged. It is nesissary to physically understand the actions of every tool to avaid lots of error messages. The ``log`` of the FluxCompensator objects might help.

Examples Which Will Not Work
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. warning:: :ref:`Convolution with a filter <label_filter>` can only be performed once! 

For 3D :ref:`FC_objects <label_objects>` :class:`~fluxcompensator.cube.SyntheticCube` with the ``shape`` (x, y, wav)  (e.g. extracted from `ModelOutput.get_image <http://docs.hyperion-rt.org/en/stable/setup/setup_images.html>`_) after the :ref:`convolution with a filter <label_filter>` a 2D :class:`~fluxcompensator.image.SyntheticImage` of the ``shape`` (x, y) with ``wav`` at the filter wavelength remains. Convolving again with a filter or :ref:`collapsing to a SED <label_rough>` is not possible anymore.  For a :class:`~fluxcompensator.sed.SyntheticSED` object (e.g. extracted from `ModelOutput.get_sed <http://docs.hyperion-rt.org/en/stable/api/hyperion.model.ModelOutput.html?highlight=get_sed#hyperion.model.ModelOutput.get_sed>`_ or after :ref:`SyntheticCube.get_rough_sed <label_rough>`) when :ref:`convolving with a filter <label_filter>` , a scalar photometric flux remains which is a member of :class:`~fluxcompensator.flux.SyntheticFlux`. Now filter convoltion is not possible anymore. 

.. warning:: :ref:`Convolving with a PSF <label_psf>` can only be performed on 3D or 2D :ref:`FC_objects <label_objects>`. 

For :ref:`FC_objects <label_objects>` of 1D :class:`~fluxcompensator.sed.SyntheticSED` or 0D :class:`~fluxcompensator.flux.SyntheticFlux` it is not possible to :ref:`convolved with a PSF <label_psf>` which has a 2D shape.