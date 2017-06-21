import pytest

import numpy as np
import os

from ..database.compact_pipeline import CompactPipeline, IRAC1, IRAC2, IRAC3, IRAC4, MIPS1, MIPS2, MIPS3, PACS1, PACS2, PACS3, SPIRE1, SPIRE2, SPIRE3, GLIMPSE_IRAC1, GLIMPSE_IRAC2, GLIMPSE_IRAC3, GLIMPSE_IRAC4, MIPSGAL_MIPS1, MIPSGAL_MIPS2, HIGAL_PACS1, HIGAL_PACS3, HIGAL_SPIRE1, HIGAL_SPIRE2, HIGAL_SPIRE3, WISESURVEY_WISE1, WISESURVEY_WISE2, WISESURVEY_WISE3, WISESURVEY_WISE4, TWOMASSSURVEY_J, TWOMASSSURVEY_H, TWOMASSSURVEY_K

def test_run():
    
    from fluxcompensator.psf import GaussianPSF
    psf_object = GaussianPSF(diameter=350.)
    
    import fluxcompensator.database.missions as filters
    filter_object = getattr(filters, 'PACS1_FILTER')
    
    pipe = CompactPipeline(PSF_object=psf_object, PSF_resolution=3., filter_object=filter_object)
    
def test_run_test_filter():
    
    assert IRAC1.filter_object.waf_0 == 3.550
    assert IRAC2.filter_object.waf_0 == 4.493
    assert IRAC3.filter_object.waf_0 == 5.731
    assert IRAC4.filter_object.waf_0 == 7.872
    
    assert MIPS1.filter_object.waf_0 == 23.68
    assert MIPS2.filter_object.waf_0 == 71.42
    assert MIPS3.filter_object.waf_0 == 155.9
    
    assert PACS1.filter_object.waf_0 == 70.
    assert PACS2.filter_object.waf_0 == 100.
    assert PACS3.filter_object.waf_0 == 160.

    assert SPIRE1.filter_object.waf_0 == 250.
    assert SPIRE2.filter_object.waf_0 == 350.
    assert SPIRE3.filter_object.waf_0 == 500.
