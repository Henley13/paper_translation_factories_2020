# -*- coding: utf-8 -*-

"""
Detect spots in the FISH image.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
import numpy as np

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
    parser.add_argument("threshold",
                        help="Detection threshold.",
                        type=int)
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
    threshold = args.threshold
    base_directory = args.base_directory
    output_directory = args.output_directory
    projection_cyt_directory = os.path.join(output_directory,
                                            "cyt_projection_mip")
    log_filter_directory = os.path.join(output_directory, "cyt_filter_log")
    spot_directory = os.path.join(output_directory, "spot_detection")
    plot_directory = os.path.join(output_directory, "plot", "spot_detection")
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(log_filter_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(log_filter_directory))
    if not os.path.isdir(spot_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(spot_directory))
    if not os.path.isdir(plot_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory))
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
    print("Detection threshold: {0}".format(threshold))
    print("Output directory: {0}".format(output_directory))
    print("MIP directory: {0}".format(projection_cyt_directory))
    print("LoG filter directory: {0}".format(log_filter_directory))
    print("Spot directory: {0}".format(spot_directory))
    print("Plot directory: {0}".format(plot_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Detecting spots...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory,
                                 return_image=False)

    # get sigma
    sigma_z, sigma_yx = detection.get_sigma(resolution_z=300,
                                            resolution_yx=103,
                                            psf_z=350,
                                            psf_yx=150)
    sigma = (sigma_z, sigma_yx, sigma_yx)

    nb_images = 0
    for i, _ in enumerate(generator):
        filename = filename_base + "_" + str(i)
        print("\t", filename)

        # LoG filter
        path = os.path.join(log_filter_directory, filename + ".tiff")
        cyt_filtered_log = stack.read_image(path)

        # cyt maximum projection
        path = os.path.join(projection_cyt_directory, filename + ".png")
        cyt_mip = stack.read_image(path)
        cyt_mip_contrast = stack.rescale(cyt_mip, channel_to_stretch=0)

        # detect spot
        mask_lm = detection.local_maximum_detection(cyt_filtered_log,
                                                    minimum_distance=2)
        spots, radius, _ = detection.spots_thresholding(image=cyt_filtered_log,
                                                        sigma=sigma,
                                                        mask_lm=mask_lm,
                                                        threshold=threshold)

        # save detected spots
        path = os.path.join(spot_directory, filename)
        np.savez(path, spots=spots, radius=radius)

        # save plot detection
        path = os.path.join(plot_directory, filename + ".png")
        plot.plot_spot_detection(tensor=cyt_mip_contrast,
                                 spots=spots,
                                 radius_yx=radius[-1],
                                 rescale=False,
                                 framesize=(15, 5),
                                 remove_frame=True,
                                 path_output=path,
                                 ext="png",
                                 show=False)

        nb_images += 1

    print()
    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
