# -*- coding: utf-8 -*-

"""
Filter FISH images with LoG and background removal.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
import bigfish.detection as detection
from utils import Logger
from loader import (get_metadata_directory, generate_filename_base,
                    images_generator)


if __name__ == "__main__":
    print()

    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("experience_directory",
                        help="Name of the experience directory.",
                        type=str)
    parser.add_argument("--base_directory",
                        help="Path of the data directory.",
                        type=str,
                        default="/Users/arthur/data/2019_racha")
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
    experience_directory = args.experience_directory
    base_directory = args.base_directory
    output_directory = args.output_directory
    output_directory_log = os.path.join(output_directory, "cyt_filter_log")
    output_directory_background = os.path.join(output_directory,
                                               "cyt_filter_background")
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(output_directory_log):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory_log))
    if not os.path.isdir(output_directory_background):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory_background))
    if not os.path.isdir(log_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(log_directory))

    # initialize logging
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M:%S")
    log_file = os.path.join(
        log_directory, "log" + "_" + experience_directory)
    sys.stdout = Logger(log_file)

    print("Running {0} file...".format(os.path.basename(__file__)), "\n")
    start_time = time.time()

    print("Data directory: {0}".format(base_directory))
    print("Experience directory name: {0}".format(experience_directory))
    print("Output directory: {0}".format(output_directory))
    print("LoG filter directory: {0}".format(output_directory_log))
    print("Background filter directory: {0}"
          .format(output_directory_background))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Filter images...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory)

    nb_images = 0
    for i, image in enumerate(generator):
        filename = filename_base + "_" + str(i)
        print("\t", image.shape, image.dtype, filename)

        # cyt
        cyt = image[0, 1, :, :, :]

        # get sigma
        sigma_z, sigma_yx = detection.get_sigma(resolution_z=300,
                                                resolution_yx=103,
                                                psf_z=350, psf_yx=150)
        sigma_log = (sigma_z, sigma_yx, sigma_yx)
        sigma_background = (sigma_z*5, sigma_yx*5, sigma_yx*5)

        # LoG filter
        cyt_filtered_log = stack.log_filter(cyt, sigma_log, keep_dtype=True)
        path = os.path.join(output_directory_log, filename + ".tiff")
        stack.save_image(cyt_filtered_log, path)

        # background filter
        cyt_filtered_background = stack.remove_background_gaussian(
            image=cyt, sigma=sigma_background)
        path = os.path.join(output_directory_background, filename + ".tiff")
        stack.save_image(cyt_filtered_background, path)

        nb_images += 1

    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
