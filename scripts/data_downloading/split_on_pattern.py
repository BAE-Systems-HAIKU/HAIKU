"""Separate files into unique subdirectories based on a pattern.

Extracts matching element of filepath to use as subdirectory name,
then moves the file into that subdirectory.

Useful for separating datasets lumped together but with distinguishing
features in their filename (e.g., ensemble runs, different hemispheres)
"""

import argparse
import os
import re
import shutil


if __name__ == "__main__":
    # cli args
    parser = argparse.ArgumentParser()
    parser.add_argument("input_directory", type=str,
                        help="Directory containing datasets")
    parser.add_argument("--recursive", "-r", action="store_true",
                        help="Recursively process subdirectories")
    parser.add_argument("--pattern", type=str, default="\.\d{3}\.",
                        help="Optional pattern if different than default")
    args = parser.parse_args()

    # verify
    assert os.path.exists(args.input_directory)

    # accumulate list of directories to process
    # if recursive is set, then all subdirectories will be added as well
    directories = []
    for dirpath, _, _ in os.walk(args.input_directory):
        directories.append(dirpath)
        if not args.recursive:
            break

    # accumulate filepaths per directory
    for directory in directories:
        pattern_match_map = {}
        for filename in os.listdir(directory):
            # identify ensemble
            match = re.search(args.pattern, filename)
            if match is None:
                print("SKIPPING: No pattern match for %s" % filename)
                continue
            if len(match.groups()) > 1:
                print("SKIPPING: More than one matching groups %s"
                      % match.groups())
                continue
            pattern_match = match.group(0).replace(".", "")
            print("Pattern Match: %s" % pattern_match)

            # add filepath to ensemble map
            filepath = os.path.join(directory, filename)
            if pattern_match not in pattern_match_map:
                pattern_match_map[pattern_match] = []
            pattern_match_map[pattern_match].append(filepath)

        # move ensembles to subdirectories
        for pattern_match, filepaths in pattern_match_map.items():
            # create subdirectory
            directory_filepath =\
                os.path.join(directory, pattern_match)
            os.makedirs(directory_filepath, exist_ok=True)

            # move files
            for filepath in filepaths:
                filename = os.path.basename(filepath)
                new_filepath = os.path.join(directory_filepath, filename)
                shutil.move(filepath, new_filepath)
