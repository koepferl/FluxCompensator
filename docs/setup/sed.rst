.. _label_sed:

-------------
SyntheticSED
-------------

To start with a `1D image group <http://docs.hyperion-rt.org/en/stable/setup/setup_images.html>`_ from the `Hyperion <http://www.hyperion-rt.org>`_ output ``hyperion.rtout`` your Python script ``example.py`` should contain::

    import numpy as np
    
    from hyperion.model import ModelOutput
    from hyperion.util.constants import pc
    from fluxcompensator.sed import *
    
    # read in from Hyperion
    m = ModelOutput('hyperion.rtout')
    array = m.get_sed(group=0, inclination=0, distance=300*pc,
                      units='ergs/cm^2/s')

Now your SED (from image ``group=0`` starting from 0, which contains a SED in this example) is scaled to ``300 pc`` in units of ``'ergs/cm^2/s'``. For further details see `ModelOutput <http://docs.hyperion-rt.org/en/stable/postprocessing/extracting_observables.html>`_.

To start with the FluxCompensator class :class:`~fluxcompensator.sed.SyntheticSED`, you simlpy write::

    # initial FluxCompensator array        
    FC_object = SyntheticSED(input_array=array, unit_out='ergs/cm^2/s', name='test_sed')

The output unit of FluxCompensator ``unit_out`` can be defined. Possible are all units like in `get_sed <http://docs.hyperion-rt.org/en/stable/api/hyperion.model.ModelOutput.html?highlight=get_sed#hyperion.model.ModelOutput.get_sed>`_. ``unit_out='ergs/cm^2/s'`` is the default. 
The name of the FluxCompensator run ``name`` can be defined with a ``str``. All outputs (e.g. plots) will start with this name.
