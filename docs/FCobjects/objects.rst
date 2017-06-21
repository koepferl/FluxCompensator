.. _label_objects:

===============================
FluxCompensator objects
===============================

The four objects passed in the FluxCompensator package. 

Dimensions
----------

All follow the same structure in attributes. They differ only in the attribute dimension of the physical property ``val`` and its corresponding wavelength ``wav``.

.. _making-a-table:

===============================================   =======================================   =======================================   =======================================
class                                             type                                      val dimension	                          wav dimension 
===============================================   =======================================   =======================================   =======================================
:class:`~fluxcompensator.cube.SyntheticCube`      3D cube                                   3D (x, y, wav)                            1D array
:class:`~fluxcompensator.image.SyntheticImage`    2D image                                  2D (x, y)                                 0D array
:class:`~fluxcompensator.sed.SyntheticSED`        1D sed                                    1D (wav)                                  1D array
:class:`~fluxcompensator.flux.SyntheticFlux`      0D val                                    0D (wav)                                  0D array
===============================================   =======================================   =======================================   =======================================

The `Hyperion <http://www.hyperion-rt.org>`_ output can be read by :class:`~fluxcompensator.cube.SyntheticCube` when we have a `3D Hyperion ModelOutput <http://docs.hyperion-rt.org/en/stable/setup/setup_images.html>`_ and by :class:`~fluxcompensator.sed.SyntheticSED` when we have a `1D Hyperion ModelOutput <http://docs.hyperion-rt.org/en/stable/setup/setup_images.html>`_. For further information see :ref:`label_setup`.
   
.. warning:: Within the :ref:`label_virtual_pipeline` and :ref:`label_post` and methods might change the dimension of the physical property ``val`` passed in the object and a new object type is returned. Also some methods do not work on certain dimension. 

Here are all possible combinations listed:

==============================================    ===========================================    =============================================
input FC_object                                   method documentation                           output FC_object                                 
==============================================    ===========================================    =============================================
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`extinction <label_extinction>`           :class:`~fluxcompensator.cube.SyntheticCube`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`change_resolution <label_resolution>`    :class:`~fluxcompensator.cube.SyntheticCube`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`convolve_psf <label_psf>`                :class:`~fluxcompensator.cube.SyntheticCube`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`convolve_filter <label_filter>`          :class:`~fluxcompensator.image.SyntheticImage`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`add_noise <label_noise>`                 :class:`~fluxcompensator.cube.SyntheticCube`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`get_rough_sed <label_rough>`             :class:`~fluxcompensator.sed.SyntheticSED`
:class:`~fluxcompensator.cube.SyntheticCube`      :ref:`get_total_val <label_total>`             :class:`~fluxcompensator.flux.SyntheticFlux`
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`extinction <label_extinction>`           :class:`~fluxcompensator.image.SyntheticImage`
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`change_resolution <label_resolution>`    :class:`~fluxcompensator.image.SyntheticImage`
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`convolve_psf <label_psf>`                :class:`~fluxcompensator.image.SyntheticImage`
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`convolve_filter <label_filter>`          ``ERROR``
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`add_noise <label_noise>`                 :class:`~fluxcompensator.image.SyntheticImage`
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`get_rough_sed <label_rough>`             ``ERROR``
:class:`~fluxcompensator.image.SyntheticImage`    :ref:`get_total_val <label_total>`             :class:`~fluxcompensator.flux.SyntheticFlux`
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`extinction <label_extinction>`           :class:`~fluxcompensator.sed.SyntheticSED`
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`change_resolution <label_resolution>`    ``ERROR``
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`convolve_psf <label_psf>`                ``ERROR``
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`convolve_filter <label_filter>`          :class:`~fluxcompensator.flux.SyntheticFlux`
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`add_noise <label_noise>`                 ``ERROR``
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`get_rough_sed <label_rough>`             ``ERROR``
:class:`~fluxcompensator.sed.SyntheticSED`        :ref:`get_total_val <label_total>`             :class:`~fluxcompensator.flux.SyntheticFlux`
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`extinction <label_extinction>`           :class:`~fluxcompensator.flux.SyntheticFlux`
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`change_resolution <label_resolution>`    ``ERROR``
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`convolve_psf <label_psf>`                ``ERROR``
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`convolve_filter <label_filter>`          ``ERROR``
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`add_noise <label_noise>`                 ``ERROR``
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`get_rough_sed <label_rough>`             ``ERROR``
:class:`~fluxcompensator.flux.SyntheticFlux`      :ref:`get_total_val <label_total>`             ``ERROR``
==============================================    ===========================================    =============================================

.. note:: Knowing the current object type and understanding the physical actions of the methods is essential. Attributes like ``log`` and ``stage`` come in handy. Here one can see the history and type of object.

Attributes & Properties
-----------------------

Attributes and properties of the objects can be called by adding to the script::

    # print wav attribute
    print c.wav
    
    # print resolution property in arcsec
    print c.resolution['arcsec']

Attributes like in ModelOutput of Hyperion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

wav : numpy.ndarray
    The wavelength of ``val`` entry in microns.

val : numpy.ndarray
    The physical property.

units : str
    Current units of ``val``.
    
distance : str
    Distance to the observed object in cm.
    
x_min : float 
    Physical offset from axis origin in FOV in cm.
    
x_max : float 
    Physical offset from axis origin in FOV in cm.

y_min : float
    Physical offset from axis origin in FOV in cm.

y_max : float 
    Physical offset from axis origin in FOV in cm.

lon_min : float
    Minimal longitudinal angle.
    
lon_max : float
    Maximal longitudinal angle.
    
lat_min : float
    Minimal latitudinal angle.
    
lat_max : float
    Maximal latitudinal angle.

pix_area_sr : float
    Pixel area per sr.
    
 
Attributes specific for the FluxCompensator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If ``input_array`` is already a :class:`~fluxcompensator.cube.SyntheticCube` object, the attributes are
passed. If input_array is not a :class:`~fluxcompensator.cube.SyntheticCube` object, :class:`~fluxcompensator.cube.SyntheticCube`
specific attributes are defined and then passed. 

unit_in : str
    Unit of ``val`` in ``input_array``. Valid options are:

        * ``'ergs/cm^2/s'``
        * ``'ergs/cm^2/s/Hz'``
        * ``'Jy'``
        * ``'mJy'``
        * ``'MJy/sr'``
      
grid_unit : float
    Physical unit of ``FOV`` axis in cm. Valid options are:
    
        * ``au`` in cm
        * ``pc`` in cm
        * ``kpc`` in cm
    
grid_unit_name
    Astronomical unit of ``FOV`` axis. Valid options are:

        * ``'au'``
        * ``'pc'``
        * ``'kpc'``
        
FOV : tuple
    Tuple ``FOV(x,y)`` of Field of View pixel entries.
    
        * pixel in x direction: ``FOV[0]``
        * pixel in y direction: ``FOV[1]``
        
name : str
    The name of the FluxCompensator object until another  
    input_array is called. The default is ``None``.

stage : str
    Gives current operation stage of SyntheticCube.
    E. g. ``'SyntheticCube: convolve_filter'``
    
log : list
    List of strings of the previous and current stages.
    
filter : dict
    Dictionary ``filter = {name, waf_0, waf_min, waf_max}`` 
    of the applied filter. 
    
        * name of filter:     ``filter['name']``
        * central wavelength: ``filter['waf_0']``
        * minimal wavelength:  ``filter['waf_min']``
        * maximal wavelength: ``filter['waf_max']``

Properties specific for the FluxCompensator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Properties are updated in the pipeline.

spacing_wav : float, ``None``        
	The property spacing_wav estimates the width of the logarithmic spaced wav entries.
        
pixel : tuple
	The property pixel is a tuple which resembles the current pixel in a
	value val. ``pixel(x,y)`` are calls as follows:
         
        * x: ``pixel[0]``
        * y: ``pixel[1]``
    
shape : tuple
	The property shape is a string, which resembles the current shape of
	the value val. 
        
        * scalar: ``'()'`` 
        * 1D:     ``'(wav)'`` 
        * 2D:     ``'(x, y)'`` 
        * 3D:     ``'(x, y , wav)'`` 

resolution : dict
	The property resolution tells you the current resolution. If we are already 
	in the :class:`~fluxcompensator.sed.SyntheticSED` or :class:`~fluxcompensator.flux.SyntheticFlux` dimension entries are considered as one large pixel.

            * resolution in arcsec per pixel : ``resolution['arcsec']``
            * resolution in rad per pixel    : ``resolution['rad']``

   