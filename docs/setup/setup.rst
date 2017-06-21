.. _label_setup:

==========
Setting Up
==========

The FluxCompensator is written in a way, which makes it easy to write pipelines and automate them later. Therefore, it is better to write the code in scripts (e.g. ``example.py``) and run it by::

    Python example.py

In the `Hyperion image set-up  <http://docs.hyperion-rt.org/en/stable/setup/setup_images.html>`_ you can set up images and/ or SEDs as image groups. When post-processing them, depending on the image group, one has to choose between two initial setup of FluxCompensator.

.. toctree::
   :maxdepth: 1
   
   cube.rst
   sed.rst

When you want to make a realistic observation out of an image output from `Hyperion <http://www.hyperion-rt.org>`_ use :class:`~fluxcompensator.cube.SyntheticCube`. If you have the SED `Hyperion <http://www.hyperion-rt.org>`_ output, than use :class:`~fluxcompensator.sed.SyntheticSED` for setting up.

.. note:: In the process of producing a "realistic" observation, the 3D (x, y, wav) image group or the 1D (wav) sed group will change its dimensions. For instance convolving with a filter transforms 3D -> 2D and 1D -> scalar (photometric flux). Then automatically the class properties are changed to :class:`~fluxcompensator.image.SyntheticImage` or :class:`~fluxcompensator.flux.SyntheticFlux`, respectively.
