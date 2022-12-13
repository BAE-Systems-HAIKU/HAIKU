"""
@classification
UNCLASSIFIED
@itar
ITAR CONTROLLED

@copyright
Copyright BAE Systems
Copyright (c) 2022 AIMdyn Inc.

Please reference the LICENSE.txt file included with this software
package for all terms and conditions related to this software
including license, rights and distribution.
"""

import numpy as np

def create_mask_dict(masked_latitudes, masked_longitudes, n_latitudes, n_longitudes):
    """Construct Mask object for use throughout HAIKU """
    # load mask
    #TODO: store this as an object/method of an object
    mask = {}
    mask['maskedLat'] = masked_latitudes
    mask['maskedLon'] = masked_longitudes
    # apply mask for plotting

    alphaMask = np.zeros(
        (n_longitudes,n_latitudes))
    alphaMask[mask['maskedLon'], mask['maskedLat']] = 1.0
    mask['alphaMask'] = alphaMask
    return(mask)
