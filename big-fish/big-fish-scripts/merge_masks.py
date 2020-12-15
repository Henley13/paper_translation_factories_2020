# -*- coding: utf-8 -*-

"""
Merge nuclei masks.
"""

import os
import argparse
import time
import datetime
import sys
import re

import bigfish.stack as stack
import bigfish.segmentation as segmentation
import bigfish.plot as plot

import numpy as np

from skimage.morphology import remove_small_objects, remove_small_holes

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
    mask_directory_1 = os.path.join(output_directory, "nuc_mask_1")
    mask_directory_2 = os.path.join(output_directory, "nuc_mask_2")
    mask_directory = os.path.join(output_directory, "nuc_mask")
    mask_directory_original = os.path.join(output_directory,
                                           "nuc_mask_original")
    plot_directory_original = os.path.join(output_directory, "plot",
                                           "segmentation_nuc_original")
    plot_directory = os.path.join(output_directory, "plot", "segmentation_nuc")
    log_directory = args.log_directory
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(projection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(projection_directory))
    if not os.path.isdir(mask_directory_1):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_directory_1))
    if not os.path.isdir(mask_directory_2):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_directory_2))
    if not os.path.isdir(mask_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_directory))
    if not os.path.isdir(mask_directory_original):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_directory_original))
    if not os.path.isdir(plot_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory))
    if not os.path.isdir(plot_directory_original):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory_original))
    if not os.path.isdir(log_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(log_directory))

    # initialize logging
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M:%S")
    log_file = os.path.join(log_directory, "log_merge_masks_" + date)
    sys.stdout = Logger(log_file)

    print("Running {0} file...".format(os.path.basename(__file__)), "\n")
    start_time = time.time()

    print("Output directory: {0}".format(output_directory))
    print("Projection directory 1: {0}".format(projection_directory))
    print("Mask directory 1: {0}".format(mask_directory_1))
    print("Mask directory 2: {0}".format(mask_directory_2))
    print("Mask directory: {0}".format(mask_directory))
    print("Mask directory original: {0}".format(mask_directory_original))
    print("Plot directory original: {0}".format(plot_directory_original))
    print("Plot directory: {0}".format(plot_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Merging masks...")

    # start analysis
    masks_1 = os.listdir(mask_directory_1)
    pattern = "[0-9A-Za-z]*_[A-Za-z]*_[a-z]*_[A-Za-z]*_[0-9A-Za-z]*_[0-9]*_[0-9]"

    nb_images = 0
    for filename_tiff in masks_1:
        if not re.match(pattern, filename_tiff):
            print("\t '{0}' is not a valid file. Skipped."
                  .format(filename_tiff))
            continue

        # filename
        filename = str(filename_tiff.split(".")[0])
        filename_base = str(filename.split("_")[:-1])
        i = int(filename.split("_")[-1])
        print("\t", filename)

        # nuclei projection
        path = os.path.join(projection_directory, filename + ".png")
        nuc_focus = stack.read_image(path)

        # mask 1
        path = os.path.join(mask_directory_1, filename_tiff)
        mask_1 = stack.read_image(path)

        # mask 2
        path = os.path.join(mask_directory_2, filename_tiff)
        mask_2 = stack.read_image(path)

        # mask
        mask = segmentation.merge_labels(mask_1, mask_2)

        # postprocess mask and mask_1
        mask_1 = segmentation.dilate_erode_labels(mask_1)
        mask = segmentation.dilate_erode_labels(mask)

        # save mask and mask_1
        path = os.path.join(mask_directory_original, filename + ".png")
        stack.save_image(mask_1, path)
        path = os.path.join(mask_directory, filename + ".png")
        stack.save_image(mask, path)

        # plot
        path = os.path.join(plot_directory_original, filename + ".png")
        plot.plot_segmentation(nuc_focus, mask_1, rescale=False,
                               title=filename, framesize=(15, 5),
                               remove_frame=True, path_output=path, ext="png",
                               show=False)
        path = os.path.join(plot_directory, filename + ".png")
        plot.plot_segmentation(nuc_focus, mask, rescale=False, title=filename,
                               framesize=(15, 5), remove_frame=True,
                               path_output=path, ext="png", show=False)

        nb_images += 1

    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
