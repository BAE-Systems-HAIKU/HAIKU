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

from __future__ import annotations

import os
import pickle
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import copy

@dataclass
class ClimateData:

    #FIXME: currently climate data object from data/simulation has form:
    #[t,lon,lat], but from koopman prediciton: [lon,lat,t]
    #transpose() call works fine, but we should be consistent throughout
    #TODO: always store data as [x,y,t]
    data: npt.NDArray = np.empty((0, 0, 0))
    dates: npt.NDArray = np.empty(0)
    latitudes: npt.NDArray = np.empty(0)
    longitudes: npt.NDArray = np.empty(0)
    time_days: npt.NDArray = np.empty(0)
    description: str = ""

    # TODO: Refactor out the need for this
    date_int: npt.NDArray = np.empty(0)

    def __post_init__(self):
        self.date_int =\
            np.asarray([int(d.strftime("%Y%m%d")) for d in self.dates])

    def extend(self, o: ClimateData):
        """Concatenate with other climate data.

        If this object is not initialized, the other data
        will be used to initialize it.
        """
        # we either concatenate or initialize
        # based on if this object has any content
        if self.data.size == 0:
            # unintialized, so initialize
            self.data = o.data
            self.latitudes = o.latitudes
            self.longitudes = o.longitudes
            self.dates = o.dates
            self.description = o.description
            self.date_int = o.date_int
        else:
            # initialized, so concatenate
            # combine data along time axis
            self.data = np.concatenate((self.data, o.data), axis=0)

            # combine 1D arrays
            self.dates = np.concatenate((self.dates, o.dates), axis=0)
            self.date_int = np.concatenate((self.date_int, o.date_int), axis=0)

    def average(self, climate_data_list:list[ClimateData]):
        """averages several climate data together.  Assumes they have the same dimensions and metadata"""
        #reinitialize data object
        self.data.fill(0)
        for o in climate_data_list:
            self.data = self.data + o.data
        self.data= self.data/len(climate_data_list)

    def compute_sigmas(self, climate_data_list:list[ClimateData], width:float=.65):
        """computes pointwise statistical uncertainty over distribution of climate data objects.  
        Assumes they have the same dimensions and metadata
        width is the relative number of points to be within error bars, defaults to 2 sigma"""
        self.positive_band = np.zeros_like(self.data)
        self.negative_band = np.zeros_like(self.data)
        for x in range(len(self.data)):
            for y in range(len(self.data[x])):
                for z in range(len(self.data[x][y])):
                    positive_values=[]
                    negative_values=[]
                    for o in climate_data_list:
                        if o.data[x][y][z]>self.data[x][y][z]:
                            positive_values.append(o.data[x][y][z])
                        else:
                            negative_values.append(o.data[x][y][z])
                    #if no elements are present, set error to 0
                    #typically only the case where sea ice is not present (0)
                    if len(positive_values)==0:
                        positive_values.append(self.data[x][y][z])
                    if len(negative_values)==0:
                        negative_values.append(self.data[x][y][z])
                    positive_values.sort()
                    negative_values.sort(reverse=True)
                    self.positive_band[x][y][z]=positive_values[int(len(positive_values)*width)]
                    self.negative_band[x][y][z]=negative_values[int(len(negative_values)*width)]

    def return_rmse(self, data_target:npt.NDArray,use_errors=False)->npt.NDArray:
        #FIXME: assumes [lon, lat, t]

        if use_errors:
            error_rmse = [np.zeros(np.min((self.data.shape[2],data_target.shape[2]))),
                          np.zeros(np.min((self.data.shape[2],data_target.shape[2])))]
            for i in range(len(error_rmse[0])):
                error_rmse[0][i] = calcRMSE(self.negative_band[:,:,i],data_target[:,:,i])
                error_rmse[1][i] = calcRMSE(self.positive_band[:,:,i],data_target[:,:,i])
            return(error_rmse)
    
        rmse = np.zeros(np.min((self.data.shape[2],data_target.shape[2])))
        for i in range(len(rmse)):
            rmse[i] = calcRMSE(self.data[:,:,i],data_target[:,:,i])
        return(rmse)
            
    def return_spatial_correlation(self, data_target:npt.NDArray,use_errors=False)->npt.NDArray:
        #FIXME: assumes [lon, lat, t]
        if use_errors:
            error_spatial_correlation = [np.zeros(np.min((self.data.shape[2],data_target.shape[2]))),
                                         np.zeros(np.min((self.data.shape[2],data_target.shape[2])))]
            for i in range(len(error_spatial_correlation[0])):
                error_spatial_correlation[0][i] = calcSpatialCorr(
                    self.negative_band[:,:,i],data_target[:,:,i])
                error_spatial_correlation[1][i] = calcSpatialCorr(
                    self.positive_band[:,:,i],data_target[:,:,i])
            return(error_spatial_correlation)
    
        spatial_correlation = np.zeros(np.min((self.data.shape[2],data_target.shape[2])))
        for i in range(len(spatial_correlation)):
            spatial_correlation[i] = calcSpatialCorr(self.data[:,:,i],data_target[:,:,i])
        return(spatial_correlation)

    def return_coverage(self,use_errors=False)->npt.NDArray:
        #FIXME: assumes [lon, lat, t]
        if use_errors:
            error_coverage = [np.zeros(self.data.shape[2]),np.zeros(self.data.shape[2])]
            for i in range(len(error_coverage[0])):
                error_coverage[0][i] = calcCumulative(self.negative_band[:,:,i],
                                                      (self.longitudes,self.latitudes),"POLAR")
                error_coverage[1][i] = calcCumulative(self.positive_band[:,:,i],
                                                      (self.longitudes,self.latitudes),"POLAR")
            return(error_coverage)

        coverage = np.zeros(self.data.shape[2])
        for i in range(len(coverage)):
            coverage[i] = calcCumulative(self.data[:,:,i],(self.longitudes,self.latitudes),"POLAR")
        return(coverage)

    def return_climatological_mean(self,frequency="Monthly")->ClimateData:
        #TODO: implement other frequencies than monthly
        #FIXME: assumes [t, lon, lat]
        assert frequency == "Monthly", "only monthly frequency currently supported in climate_data.return_climatological_mean"

        climatological_mean = np.zeros((self.data.shape[0],self.data.shape[1], self.data.shape[2]))
        #extract the integer month from date integer
        #should be faster than using self.date.month
        monthInt = self.date_int // 100 % 100
        for imonth in range(12):
            #take the average of one month over the entire dataset
            #imonth +1 so that our months are indexed with Jan=0
            climatological_mean[monthInt == (imonth+1),:,:] = np.mean(
                self.data[monthInt == (imonth+1),:,:],axis=0)

        climatological_object = copy.deepcopy(self)
        climatological_object.data = climatological_mean
        return(climatological_object)
        
    def select_hemisphere(self, hemisphere:str,lat_bound:float=40):
        # trim data based on spatial region
        #FIXME: requires data shape [_,lat,_]
        if hemisphere == 'N':
            inds = np.argwhere(
                self.latitudes >= lat_bound)
            # flatten list of lists
            inds_list = [ind[0] for ind in inds]
            self.data = self.data[:, inds_list, :]
            self.latitudes = self.latitudes[inds_list]
        elif hemisphere == 'S':
            inds = np.argwhere(
                self.latitudes <= lat_bound)
            # flatten list of lists
            inds_list = [ind[0] for ind in inds]
            self.data = self.data[:, inds_list, :]
            self.latitudes = self.latitudes[inds_list]
        else:
            #POLAR
            pass
            
    def fill_arctic_top_hole(self,top_index:int=-7):
        # set nothern latitudes to 100% sea ice concentration
        # to deal with the pole hole; only relevant for CESM1 lat/lon.
        # this can be removed after the data
        # is interpolated over the pole hole
        #assumes data: [_,lat,_] shape
        self.data[:, top_index:,:]=100.0
            
    @property
    def filename(self) -> str:
        start_date = self.dates[0].strftime("%Y%m%d")
        end_date = self.dates[-1].strftime("%Y%m%d")
        return f"{self.description}_{start_date}_to_{end_date}"

    def serialize(self, filepath: str) -> None:
        if not filepath.endswith(".pkl"):
            filepath += ".pkl"
        try:
            os.makedirs(os.path.dirname(filepath))
        except FileExistsError:
            pass
            
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)

    @staticmethod
    def deserialize(filepath: str) -> ClimateData:
        with open(filepath, 'rb') as f:
            return pickle.load(f)

def calcRMSE(pred, true):
    """
    Calculate RMSE given predicted and true values. Both inputs are 2D gridded points.
    """
    rmse = np.sqrt( np.mean( np.square( pred - true ), axis=(0,1) ) )
    return rmse


def calcSpatialCorr(pred, true):
    """
    Calculate spatial correlation given predicted and true values. Both inputs are 2D gridded points.
    This needs some more refining before results are shown.     
    """
    scorr = np.corrcoef( pred.flatten(), true.flatten() )[0,1]
    return scorr


def calcCumulative(pred, grid:tuple, grid_type="N"):
    #maintain legacy input for N lon/lat grid
    #pred is 2-d spatial sea ice coverage np array
    #grid is either the latitude coordinates (w/ grid_type N)
    #or grid is a tuple of the (x,y) np arrays from the climate data object
    
    if grid_type=="N":
        return(calcCumulativeLat(pred,grid))
    #otherwise, assume polar with a grid that contains (x-centers, y-centers) in kilometers
    
    assert pred.shape[0] == len(grid[0]), f"prediction grid shape {pred.shape[0]} does not match grid {grid[0].shape}"
    assert pred.shape[1] == len(grid[1]), f"prediction grid shape {pred.shape[1]} does not match grid {grid[1].shape}"
       
    nx = len(grid[0])
    ny = len(grid[1])

    ice_area = 0
    
    miss_val = -9.99e8  # missing value over land grids
    conc_cutoff = 15    # only consider grids with more that 15% coverage of sea ice
    for ix in range(nx):
        for iy in range(ny):
            if( (pred[ix,iy]>=conc_cutoff) and (pred[ix,iy]!=miss_val) ):
                #grid is based on point centers, so this is slightly off
                #but the area doesn't vary much from point to point
                try:
                    dx = abs(grid[0][ix]-grid[0][ix+1])
                except:
                    dx = abs(grid[0][ix-1]-grid[0][ix])
                try:
                    dy = abs(grid[1][iy]-grid[1][iy+1])
                except:
                    dy = abs(grid[1][iy-1]-grid[1][iy])

                ice_area = ice_area + pred[ix,iy] * dx * dy

    # 100 is to convert from percent
    # 1000000 is to convert from sq m to million sq km.
    ice_area = ice_area / 100.0 / 1000000000000.0

    return ice_area # sea ice area


        
    
def calcCumulativeLat(pred, lat):
    """
    Calculate cumulative given 2D gridded points.
    #TODO: update this to automatically pull the correct grid
    """
    num_lat = pred.shape[1] #54  # number of latitudes
    num_lon = pred.shape[0] #288 # number of longitudes
    # lat_inc = 0.9375 # 180/192
    lat_inc = 0.94240837696 # 180/(192-1) 
    lon_inc = 1.25   # 360/288
    earth_rad = 6371.0  # radius of the Earth in km
    earth_deg = (2*np.pi*earth_rad)/360.  # 1 degree of the Earth in km
    grid_area = (earth_deg*lat_inc)*(earth_deg*lon_inc)  # area of each lat_inc by lon_inc grid in sq km
    x = 0
    y = 0
    miss_val = -9.99e8  # missing value over land grids
    conc_cutoff = 15    # only consider grids with more that 15% coverage of sea ice
    for ilon in range(num_lon):
        for ilat in range(num_lat):
            if( (pred[ilon,ilat]>=conc_cutoff) and (pred[ilon,ilat]!=miss_val) ):
                x = x + pred[ilon,ilat] * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area
                # x = x + 100 * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area
                y = y + 1.0 * np.cos( (lat[ilat]-lat_inc) /180*np.pi) * grid_area

    # 100 is to convert from percent
    # 1000000 is to convert to million sq km.
    x = x / 100.0 / 1000000.0
    y = y / 1000000.0

    return x # sea ice area
    #return y # sea ice extent 

