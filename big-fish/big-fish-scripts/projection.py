# -*- coding: utf-8 -*-

"""
Project cytoplasm and nuclei in 2-d.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
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
    output_directory_cyt_mip = os.path.join(output_directory,
                                            "cyt_projection_mip")
    output_directory_cyt_focus = os.path.join(output_directory,
                                              "cyt_projection_focus")
    output_directory_nuc_focus = os.path.join(output_directory,
                                              "nuc_projection_focus")
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(output_directory_cyt_mip):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory_cyt_mip))
    if not os.path.isdir(output_directory_cyt_focus):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory_cyt_focus))
    if not os.path.isdir(output_directory_nuc_focus):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory_nuc_focus))
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
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Projecting images...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory)

    nb_images = 0
    for i, image in enumerate(generator):
        filename = filename_base + "_" + str(i)
        print("\t", image.shape, image.dtype, filename)

        # cyt and nuc
        nuc = image[0, 0, :, :, :]
        cyt = image[0, 1, :, :, :]

        # projections
        nuc_focus = stack.focus_projection_fast(nuc,
                                                proportion=0.7,
                                                neighborhood_size=7)
        nuc_focus = stack.rescale(nuc_focus, channel_to_stretch=0)
        cyt_focus = stack.focus_projection_fast(cyt,
                                                proportion=0.75,
                                                neighborhood_size=7)
        cyt_in_focus = stack.in_focus_selection(cyt,
                                                proportion=0.80,
                                                neighborhood_size=30)
        cyt_mip = stack.maximum_projection(cyt_in_focus)

        # save projections
        path = os.path.join(output_directory_nuc_focus, filename + ".png")
        stack.save_image(nuc_focus, path)
        path = os.path.join(output_directory_cyt_focus, filename + ".png")
        stack.save_image(cyt_focus, path)
        path = os.path.join(output_directory_cyt_mip, filename + ".png")
        stack.save_image(cyt_mip, path)

        nb_images += 1

    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
