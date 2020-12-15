# -*- coding: utf-8 -*-

"""
Remove nuclei segmented from a 2-d projection.
"""

import os
import argparse
import time
import datetime
import sys
import re

import bigfish.stack as stack
import bigfish.segmentation as segmentation
from utils import Logger


if __name__ == "__main__":
    print()

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_directory",
                        help="Path of the output directory.",
                        type=str,
                        default="/Users/arthur/output/2019_racha")
    parser.add_argument("--log_directory",
                        help="Path of the log directory.",
                        type=str,
                        default="/Users/arthur/output/2019_racha/log")

    # initialize parameters
    args = parser.parse_args()
    output_directory = args.output_directory
    projection_directory = os.path.join(output_directory,
                                        "nuc_projection_focus")
    mask_directory = os.path.join(output_directory, "nuc_mask_1")
    removed_directory = os.path.join(output_directory, "nuc_without_mask_1")
    log_directory = args.log_directory
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(projection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(projection_directory))
    if not os.path.isdir(mask_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_directory))
    if not os.path.isdir(removed_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(removed_directory))
    if not os.path.isdir(log_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(log_directory))

    # initialize logging
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M:%S")
    log_file = os.path.join(log_directory, "log_removing_nuclei_" + date)
    sys.stdout = Logger(log_file)

    print("Running {0} file...".format(os.path.basename(__file__)), "\n")
    start_time = time.time()

    print("Output directory: {0}".format(output_directory))
    print("Projection directory: {0}".format(projection_directory))
    print("Mask directory: {0}".format(mask_directory))
    print("Removed nuclei directory: {0}".format(removed_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Removing nuclei...")

    # start analysis
    projections = os.listdir(projection_directory)
    pattern = "[0-9A-Za-z]*_[A-Za-z]*_[a-z]*_[A-Za-z]*_[0-9A-Za-z]*_[0-9]*_[0-9]"

    nb_images = 0
    for filename_png in projections:
        if not re.match(pattern, filename_png):
            print("\t '{0}' is not a valid file. Skipped."
                  .format(filename_png))
            continue

        # filename
        filename = str(filename_png.split(".")[0])
        filename_base = str(filename.split("_")[:-1])
        i = int(filename.split("_")[-1])
        print("\t", filename)

        # projection
        path = os.path.join(projection_directory, filename_png)
        nuc_projected = stack.read_image(path)

        # mask
        path = os.path.join(mask_directory, filename + ".tiff")
        mask = stack.read_image(path)

        # remove segmented nuclei
        unsegmented_nuclei = segmentation.remove_segmented_nuc(nuc_projected,
                                                               mask)
        path = os.path.join(removed_directory, filename_png)
        stack.save_image(unsegmented_nuclei, path)

        nb_images += 1

    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
