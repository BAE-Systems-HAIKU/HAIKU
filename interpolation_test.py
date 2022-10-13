from datetime import datetime
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from typing import List
from haiku.climate_data.climate_data import ClimateData
from haiku.climate_data.readNSIDC import readNSIDCSIC
from haiku.plotting.plotFunctions import plotSnapshot

from haiku.training.training_utilities import interpolate, interpolate_missing_months, months_difference

VMIN = 0
NUM_COLS = 2
TIMESTEPS = 5
VMAX = TIMESTEPS


def plot_matrix(matrix: npt.NDArray, col: int, num_cols: int):
    plt.subplot(1, num_cols, col)
    plt.title(f"T={col}")
    plt.imshow(
        matrix,
        vmin=VMIN,
        vmax=VMAX,
        cmap='Blues',
        aspect='auto',
        interpolation="none")

def plot_climate_data(self, climate_data: ClimateData, time_index: int,
                          display: bool, output_filepath: str = None):
        if climate_data.mask is not None:
            alpha_mask = climate_data.get_alpha_mask()
        else:
            alpha_mask = None

        formatted_date = datetime.strftime(
            climate_data.dates[time_index],
            "%Y-%m-%d"
        )
        title = f"{climate_data.name.value} at {formatted_date}"

        self.plot_geographic_data(
            climate_data.data[time_index],
            climate_data.lat,
            climate_data.lon,
            title,
            self.plotting_tools.get_label_from_data_type(
                climate_data.name
            ),
            display,
            alpha_mask,
            output_filepath
        )


def plot_geographic_data(data: npt.NDArray,
                         lat: npt.NDArray,
                         lon: npt.NDArray,
                         title: str,
                         colorbar_label: str,
                         display_to_screen: bool,
                         alpha_mask: npt.NDArray = None,
                         output_filepath: str = None):
    # clear anything existing on the axes
    plt.cla()

    # assume this is giving us (lat, lon)
    # array needs to be flipped vertically so display is correct
    data = np.flipud(data)

    # set bounds [left, right, bottom, top]
    extents = [
        lon.min(),
        lon.max(),
        lat.min(),
        lat.max()
    ]

    if alpha_mask is not None:
        img = plt.imshow(
            data,
            extent=extents,
            alpha=np.rot90(alpha_mask, axes=(0, 1)),
            aspect='auto',
            interpolation='nearest'
        )
    else:
        img = plt.imshow(data, extent=extents,
                         aspect='auto',
                         interpolation='nearest')

    cbar = plt.colorbar(img, fraction=0.018, pad=0.04)
    cbar.ax.set_ylabel(colorbar_label, rotation=90)

    plt.title(title)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")

    if output_filepath is not None:
        plt.savefig(output_filepath, bbox_inches='tight', pad_inches=.1)

    if not display_to_screen:
        plt.close()


def plot_climate_data(title: str, climate_data: ClimateData,
                      time_index: int, display: bool,
                      output_filepath: str = None):
    alpha_mask = None

    plot_geographic_data(
        climate_data.data[time_index],
        climate_data.latitudes,
        climate_data.longitudes,
        title,
        "Concentration",
        display,
        alpha_mask,
        output_filepath
    )


def test_dummy_data():
    y1 = np.ones((2, 2))
    y2 = np.ones((2, 2)) * 5
    x1 = 1
    x2 = 5
    x = 3
    y = interpolate(y1, y2, x1, x2, x)
    print(y)


def test_climate_data():
    # vars
    filepath_a = ("/home/haiku/repos/good/data/simdata/NSIDC/monthly/"
                  "north_regrid/"
                  "latlon_seaice_conc_monthly_nh_f08_198711_v03r01.nc")
    filepath_b = ("/home/haiku/repos/good/data/simdata/NSIDC/monthly/"
                  "north_regrid/"
                  "latlon_seaice_conc_monthly_nh_f08_198802_v03r01.nc")
    nsidc_var = "merged"

    # load data before
    SIC_a, dateInt_a, lat_a, lon_a, timeDays_a = readNSIDCSIC(
        filepath_a,
        nsidc_var,
        False
    )

    # load data after
    SIC_b, dateInt_b, lat_b, lon_b, timeDays_b = readNSIDCSIC(
        filepath_b,
        nsidc_var,
        False
    )

    # interpret time indices
    start_time = datetime.strptime(str(dateInt_a[0]), "%Y%m%d")
    end_time = datetime.strptime(str(dateInt_b[0]), "%Y%m%d")

    # we can interpret 0 as the starting time index
    # and elapsed_months as the ending time index
    # NOTE: This is only true since we look at data month by month
    elapsed_months = months_difference(start_time, end_time)

    missing_data = interpolate_missing_months(
        SIC_a,
        SIC_b,
        0,
        elapsed_months
    )

    for i, data in enumerate(missing_data):
        cd = ClimateData(
            data,
            dateInt_a,
            lat_a,
            lon_b,
            timeDays_a,
            "Interpolated Data"
            )

        plot_climate_data(
            "Interpolated Data",
            cd,
            0,
            False,
            f"{i}.png")



test_climate_data()
