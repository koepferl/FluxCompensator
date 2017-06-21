import numpy as np

from .missions.irac import *
from .missions.mips import *
from .missions.pacs import *
from .missions.spire import *
from .missions.twomass import *
from .missions.wise import *
from .missions.twomass import diameter as twomass_dia
from .missions.wise import diameter as wise_dia


class CompactPipeline(object):
    def __init__(self, PSF_object, PSF_resolution, filter_object):

        '''
        Wraps the information of a PSF object, its resolution a corresponding Filter object of a virtual synthetic survey.
        
        Parameters
        ----------
        
        PSF_object : object
            
            You can use different objects:
            
            * GaussianPSF
            * FilePSF
            * FunctionPSF
            * pre-constructed PSF from database
            
        PSF_resolution : float, database
            
            Also PSF_resolution from database can be passed here.
            
        filter_object : object
            
            You can use several filters:
            
            * read-in filter
            * pre-constructed filter from database
            
        '''
        
        # store objects in dict
        self.PSF_object = PSF_object
        self.PSF_resolution = PSF_resolution
        self.filter_object = filter_object
        
IRAC1 = CompactPipeline(PSF_object=IRAC1_PSF, PSF_resolution=IRAC1_PSF_RESOLUTION, filter_object=IRAC1_FILTER)
IRAC2 = CompactPipeline(PSF_object=IRAC2_PSF, PSF_resolution=IRAC2_PSF_RESOLUTION, filter_object=IRAC2_FILTER)
IRAC3 = CompactPipeline(PSF_object=IRAC3_PSF, PSF_resolution=IRAC3_PSF_RESOLUTION, filter_object=IRAC3_FILTER)
IRAC4 = CompactPipeline(PSF_object=IRAC4_PSF, PSF_resolution=IRAC4_PSF_RESOLUTION, filter_object=IRAC4_FILTER)

MIPS1 = CompactPipeline(PSF_object=MIPS1_PSF_100K, PSF_resolution=MIPS1_PSF_RESOLUTION, filter_object=MIPS1_FILTER)
MIPS2 = CompactPipeline(PSF_object=MIPS2_PSF_100K, PSF_resolution=MIPS2_PSF_RESOLUTION, filter_object=MIPS2_FILTER)
MIPS3 = CompactPipeline(PSF_object=MIPS3_PSF_100K, PSF_resolution=MIPS3_PSF_RESOLUTION, filter_object=MIPS3_FILTER)

PACS1 = CompactPipeline(PSF_object=PACS1_PSF, PSF_resolution=PACS1_PSF_RESOLUTION, filter_object=PACS1_FILTER)
PACS2 = CompactPipeline(PSF_object=PACS2_PSF, PSF_resolution=PACS2_PSF_RESOLUTION, filter_object=PACS2_FILTER)
PACS3 = CompactPipeline(PSF_object=PACS3_PSF, PSF_resolution=PACS3_PSF_RESOLUTION, filter_object=PACS3_FILTER)

SPIRE1 = CompactPipeline(PSF_object=SPIRE1_PSF, PSF_resolution=SPIRE1_PSF_RESOLUTION, filter_object=SPIRE1_FILTER)
SPIRE2 = CompactPipeline(PSF_object=SPIRE2_PSF, PSF_resolution=SPIRE2_PSF_RESOLUTION, filter_object=SPIRE2_FILTER)
SPIRE3 = CompactPipeline(PSF_object=SPIRE3_PSF, PSF_resolution=SPIRE3_PSF_RESOLUTION, filter_object=SPIRE3_FILTER)

###########
# surveys #
###########

TWOMASSSURVEY_J = CompactPipeline(PSF_object=GaussianPSF(diameter=twomass_dia), PSF_resolution=1., filter_object=J_2MASS_FILTER)
TWOMASSSURVEY_H = CompactPipeline(PSF_object=GaussianPSF(diameter=twomass_dia), PSF_resolution=1., filter_object=H_2MASS_FILTER)
TWOMASSSURVEY_K = CompactPipeline(PSF_object=GaussianPSF(diameter=twomass_dia), PSF_resolution=1., filter_object=K_2MASS_FILTER)

WISESURVEY_WISE1 = CompactPipeline(PSF_object=GaussianPSF(diameter=wise_dia), PSF_resolution=1., filter_object=WISE1_FILTER)
WISESURVEY_WISE2 = CompactPipeline(PSF_object=GaussianPSF(diameter=wise_dia), PSF_resolution=1., filter_object=WISE2_FILTER)
WISESURVEY_WISE3 = CompactPipeline(PSF_object=GaussianPSF(diameter=wise_dia), PSF_resolution=1., filter_object=WISE3_FILTER)
WISESURVEY_WISE4 = CompactPipeline(PSF_object=GaussianPSF(diameter=wise_dia), PSF_resolution=1., filter_object=WISE4_FILTER)

# below like IRAC, MIPS, PACS, SPIRE

GLIMPSE_IRAC1 = CompactPipeline(PSF_object=IRAC1_PSF, PSF_resolution=IRAC1_PSF_RESOLUTION, filter_object=IRAC1_FILTER)
GLIMPSE_IRAC2 = CompactPipeline(PSF_object=IRAC2_PSF, PSF_resolution=IRAC2_PSF_RESOLUTION, filter_object=IRAC2_FILTER)
GLIMPSE_IRAC3 = CompactPipeline(PSF_object=IRAC3_PSF, PSF_resolution=IRAC3_PSF_RESOLUTION, filter_object=IRAC3_FILTER)
GLIMPSE_IRAC4 = CompactPipeline(PSF_object=IRAC4_PSF, PSF_resolution=IRAC4_PSF_RESOLUTION, filter_object=IRAC4_FILTER)

MIPSGAL_MIPS1 = CompactPipeline(PSF_object=MIPS1_PSF_100K, PSF_resolution=MIPS1_PSF_RESOLUTION, filter_object=MIPS1_FILTER)
MIPSGAL_MIPS2 = CompactPipeline(PSF_object=MIPS2_PSF_100K, PSF_resolution=MIPS2_PSF_RESOLUTION, filter_object=MIPS2_FILTER)

HIGAL_PACS1 = CompactPipeline(PSF_object=PACS1_PSF, PSF_resolution=PACS1_PSF_RESOLUTION, filter_object=PACS1_FILTER)
HIGAL_PACS3 = CompactPipeline(PSF_object=PACS3_PSF, PSF_resolution=PACS3_PSF_RESOLUTION, filter_object=PACS3_FILTER)

HIGAL_SPIRE1 = CompactPipeline(PSF_object=SPIRE1_PSF, PSF_resolution=SPIRE1_PSF_RESOLUTION, filter_object=SPIRE1_FILTER)
HIGAL_SPIRE2 = CompactPipeline(PSF_object=SPIRE2_PSF, PSF_resolution=SPIRE2_PSF_RESOLUTION, filter_object=SPIRE2_FILTER)
HIGAL_SPIRE3 = CompactPipeline(PSF_object=SPIRE3_PSF, PSF_resolution=SPIRE3_PSF_RESOLUTION, filter_object=SPIRE3_FILTER)
