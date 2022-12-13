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
import numpy.typing as npt
import datetime

from dateutil.relativedelta import relativedelta

#TODO: better to use datetime utils throughout rather than hacking central code to handle these functions that already exist there.

def date_array_for_range(date_min:datetime.datetime,
                         date_max:datetime.datetime,
                         delta:str = "month")->npt.NDArray[datetime.datetime]:
    assert (delta=="month"), "delta for date_int of type:" + str(delta) + "not yet supported; use month instead"
    #The +1 is to be inclusive of both start and stop dates
    n_dates = (date_max.year - date_min.year) * 12 + (date_max.month - date_min.month) + 1
    dates = np.zeros(n_dates,dtype=datetime.datetime)
    dates[0]=date_min

    for i in range(n_dates-1):
        dates[i+1]=dates[i]+relativedelta(months=+1)

    return(dates)

def date_int_array_for_range(date_min:int,
                             date_max:int,
                             delta:str = "month")->npt.NDArray[int]:
    assert (delta=="month"), "delta for date_int of type:" + str(delta) + "not yet supported; use month instead"
    n_dates = n_dates_in_range(date_min, date_max, delta)
    
    dates = np.zeros(n_dates,dtype=int)
    dates[0]=int(date_min)

    for i in range(n_dates-1):
        last_month = int(str(dates[i])[4:6])
        this_year  = int(str(dates[i])[0:4])
        this_day   = int(str(dates[i])[6:8])
        if last_month==12:
            this_year = this_year+1
            last_month=0
        this_month = last_month+1
        dates[i+1] = int(str(this_year).zfill(4)+str(this_month).zfill(2)+str(this_day).zfill(2))
    return(dates)
        
def n_dates_in_range(date_min:int, date_max:int, delta:str = "month")->npt.NDArray:
    """computes inclusive number of dates between two dates in integer format by delta step size"""
    
    assert (delta=="month"), "delta for date_int of type:" + str(delta) + "not yet supported; use month instead"
    month_diff = int(str(date_max)[4:6]) - \
                 int(str(date_min)[4:6])
    year_diff  = int(str(date_max)[0:4]) - \
                 int(str(date_min)[0:4])
    n_dates = year_diff * 12 + month_diff + 1
    return(n_dates)
