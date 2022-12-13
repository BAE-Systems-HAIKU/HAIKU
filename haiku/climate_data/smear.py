import logging
from datetime import datetime
from typing import List, Tuple, Union

import numpy as np
import numpy.typing as npt
import copy
from haiku.climate_data.climate_data import ClimateData


class Smearer:
    """
       object to generate smeared ClimateData objects from ClimateData Object
       if smearing_factor is Float, apply guassian smearing to each data element where smearing_factor is the 1 sigma magnitude (relative smearing to data value)
       if smearing_factor is NDArray, must match dimensions of data; each element smeared by smearing factor array
    """
    @staticmethod
    def smear_climate_data_gaussian(data: ClimateData, smearing_factor:float=0.1) -> ClimateData:

        #apply guassian smearing to each climate_data data value
        #use smearing_factor as 1 sigma relative smear to apply
        smeared_data = copy.deepcopy(data)
        for i, value in enumerate(smeared_data.data):
            smeared_data.data[i]= value * (1 + np.random.normal(scale=smearing_factor))
        #make sure smearing doesn't go outside of training bounds
        smeared_data.data = np.clip(smeared_data.data,data.data.min(),data.data.max())
        return(smeared_data)
    

    @staticmethod
    def smear_climate_data_fully_mapped(data: ClimateData, smearing_factor:npt.NDArray) -> ClimateData:

        #apply absolute smearing from array for each data element
        #use smearing_factor as relative smear to apply (can vary with time)
        smeared_data = np.copy(data)
        for i, value in enumerate(smeared_data.data):
            smeared_data.data[i]= value * (1 + smearing_factor[i])
        return(smeared_data)

    @staticmethod
    def smear_climate_data_2d_mapped(data: ClimateData, smearing_factor:npt.NDArray) -> ClimateData:

        Exception("2d_mapped smearing is not yet implemented")
        
        #apply absolute smearing from array for each spatial grid element
        #use smearing_factor as relative smear to apply (independent of time)
        smeared_data = np.copy(data)
        for i, value in enumerate(smeared_data.data):
            smeared_data.data[i]= value * (1 + smearing_factor[i])
        return(smeared_data)
