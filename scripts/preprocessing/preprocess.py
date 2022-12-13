import argparse
import logging
import os
import re

from haiku.climate_data.preprocessor import Preprocessor


def configure_logging(logging_level: str, filename: str):
    """Create and format the global logging logger with console/file logging.

    logging_level can be: CRITICAL, ERROR, WARNING, INFO, or DEBUG
    Note: Can also use logging enums for this value (e.g. logging.INFO)
    """
    # grab the global logging logger and set the logging level
    logger = logging.getLogger()

    # ensure any other handlers are removed
    for handler in logger.handlers:
        logger.removeHandler(handler)

    logger.setLevel(logging_level)

    # create formatter to be used by the global logger
    formatter = logging.Formatter("[%(asctime)s.%(msecs)03d] - "
                                  "%(filename)s:%(funcName)s:%(levelname)s - "
                                  "%(message)s", "%Y-%m-%d %H:%M:%S")

    # we create a stream handler for console logging
    # and a file handler for file logging
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    fh = logging.FileHandler(filename)
    fh.setFormatter(formatter)

    # finally, we add the console logging handler to the logger
    logger.addHandler(ch)
    logger.addHandler(fh)

def main(args):

    # verify args
    assert os.path.exists(args.input_directory), \
        f"Input directory not found: {args.input_directory}"
    os.makedirs(args.output_directory, exist_ok=True)

    # configure logging (put log in same directory as script)
    log_filepath =\
        os.path.join(os.path.dirname(__file__), "haiku_preprocessing.log")
    configure_logging(logging.INFO, log_filepath)

    # accumulate list of directories to process
    # if recursive is set, then all subdirectories will be added as well
    directories = []
    for dirpath, _, _ in os.walk(args.input_directory):
        directories.append(dirpath)
        if not args.recursive:
            break
    logging.info("Directories to process: %s", directories)

    # preprocess each directory
    preprocessor = Preprocessor(args.reference_gridfile)
    for i, directory in enumerate(directories):
        # see if there is a pattern match
        match = re.search(args.pattern, directory)
        if match is None:
            logging.warning("Skipping directory as pattern does not match: %s",
                            directory)
            continue
        logging.info("Processing %s", directory)
        preprocessor.preprocess(
            directory,
            args.output_directory,
            args.only_plot_first
        )

def parse_args():
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("input_directory", type=str,
                        help="Directory containing data files")
    parser.add_argument("reference_gridfile", type=str,
                        help="Filepath to grid description textfile")
    parser.add_argument("output_directory", type=str,
                        help="Directory to store serialized output")
    parser.add_argument("--recursive", "-r", action="store_true",
                        help="If specified, recursively process directories "
                             "within the input directory")
    parser.add_argument("--pattern", type=str, default="",
                        help="Optionally only include directories "
                             "with a pattern match")
    parser.add_argument("--only_plot_first", action="store_true",
                        help="If specified, will only plot the first timestep "
                             "of preprocessed data")
    return(parser.parse_args())

if __name__ == "__main__":
    args = parse_args()
    main(args)
