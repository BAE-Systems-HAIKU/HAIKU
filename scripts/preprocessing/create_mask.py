import argparse
import os
import sys

from haiku.climate_data.mask import Mask
from netCDF4 import Dataset
import numpy as np

import matplotlib.pyplot as plt


def compute_mask(input_filepath: str, variable: str, minimum_value: float,
                 maximum_value: float, output_filepath: str) -> None:
    # extract the data from the dataset
    with Dataset(input_filepath) as dataset:
        dataset_vars = list(dataset.variables.keys())
        if variable is None:
            print("No variable specified.")
            print("Options: %s" % dataset_vars)
            print("Rerun with the `--variable` "
                  "flag and one of the above options.")
            sys.exit(1)

        if variable not in dataset.variables:
            raise ValueError(
                "Specified variable %s not found in extracted data. "
                "Valid variables: %s" % (variable, dataset_vars))

        print("Loading variable %s from dataset" % variable)
        mapping = np.squeeze(np.array(dataset.variables[variable]))

    # create map mask from data
    # contains all DESIRED points
    M = np.where((mapping >= minimum_value) & (mapping <= maximum_value))

    # cast to numpy array
    M = np.asarray(M)

    # plot
    print("Plotting...")
    fig, axes = plt.subplots(1, 2)

    # unmasked map
    axes[0].set_title("Unmasked Map")
    axes[0].set_xlabel("Longitude")
    axes[0].set_ylabel("Latitude")
    axes[0].imshow(np.flipud(mapping))

    # enforce masking pixel value
    mapping[M[0], M[1]] = -100

    # data before interpolation
    temp_plot = np.ma.masked_where((mapping == -100), mapping)
    cmap = plt.cm.get_cmap("viridis").copy()
    cmap.set_bad(color='red')
    plt.imshow(np.flipud(temp_plot), cmap=cmap)
    axes[1].set_title("Masked Map")
    axes[1].set_xlabel("Longitude")
    axes[1].imshow(np.flipud(temp_plot), cmap=cmap)

    plot_filename = output_filepath + ".png"
    plot_filepath =\
        os.path.join(os.path.dirname(output_filepath), plot_filename)
    plt.savefig(plot_filepath)
    plt.close()

    print("Serializing to %s" % output_filepath)
    mask = Mask(M)
    mask.serialize(output_filepath)


def arg_parse():
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("input_filepath", type=str,
                        help="Path to dataset to use for mask computation")
    parser.add_argument("output_filepath", type=str,
                        help="Where to save the mask")
    parser.add_argument("minimum_value", type=float,
                        help="Values lower than this will "
                             "not be included in the mask")
    parser.add_argument("maximum_value", type=float,
                        help="Values greater than this will "
                             "not be included in the mask")
    parser.add_argument("--variable", type=str, default=None,
                        help=("Variable to use for map. "
                              "Only applicable if input dataset contains "
                              "more than one variable of interest. "
                              "If not specified, program will output "
                              "the list of options to choose from."))
    args = parser.parse_args()

    # verify args
    assert os.path.exists(args.input_filepath), \
        f"{args.input_filepath} not found"
    assert args.minimum_value <= args.maximum_value, \
        ("Minimum not less than maximum! "
         f"{args.minimum_value} !< {args.maximum_value}")
    return(args)

if __name__ == "__main__":
    args = arg_parse()
    # setup output directory
    # explicitly set output directory if use has only specified filename
    if os.path.dirname(args.output_filepath) == "":
        args.output_filepath = os.path.join(os.curdir, args.output_filepath)
    os.makedirs(os.path.dirname(args.output_filepath), exist_ok=True)

    # process
    compute_mask(
        args.input_filepath,
        args.variable,
        args.minimum_value,
        args.maximum_value,
        args.output_filepath)
