# -*- coding: utf-8 -*-

"""
Detect foci in the FISH image.
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

    # input-output directories
    decomposition_directory = os.path.join(output_directory,
                                           "cluster_decomposition")
    foci_directory = os.path.join(output_directory, "foci_detection")
    cyt_projection_directory = os.path.join(output_directory,
                                            "cyt_projection_mip")

    # plot directory
    plot_directory_foci = os.path.join(output_directory, "plot",
                                       "foci_detection")

    # check directories exist
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(decomposition_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(decomposition_directory))
    if not os.path.isdir(foci_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(foci_directory))
    if not os.path.isdir(cyt_projection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cyt_projection_directory))
    if not os.path.isdir(plot_directory_foci):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory_foci))
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
    print("Decomposition directory: {0}".format(decomposition_directory))
    print("Foci directory: {0}".format(foci_directory))
    print("Cytoplasm projection directory: {0}"
          .format(cyt_projection_directory))
    print("Plot decomposition directory: {0}".format(plot_directory_foci))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Detecting foci...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory,
                                 return_image=False)

    nb_images = 0
    for i, _ in enumerate(generator):
        filename = filename_base + "_" + str(i)
        print("\t", filename)

        # spots
        path = os.path.join(decomposition_directory, filename + ".npz")
        data = np.load(path)
        spots_out_cluster = data["spots_out_cluster"]
        spots_in_cluster = data["spots_in_cluster"]
        clusters = data["clusters"]
        radius_spots = data["radius_spots"]

        # cytoplasm maximum projection
        path = os.path.join(cyt_projection_directory, filename + ".png")
        cyt_mip = stack.read_image(path)
        cyt_mip_contrast = stack.rescale(cyt_mip, channel_to_stretch=0)

        # detect foci
        spots = np.concatenate((spots_out_cluster, spots_in_cluster[:, :3]),
                               axis=0)
        clustered_spots = detection.cluster_spots(spots=spots,
                                                  resolution_z=300,
                                                  resolution_yx=103,
                                                  radius=350,
                                                  nb_min_spots=5)
        foci = detection.extract_foci(clustered_spots=clustered_spots)

        # save foci
        path = os.path.join(foci_directory, filename)
        np.savez(path,
                 clustered_spots=clustered_spots,
                 foci=foci)

        # plot foci
        path = os.path.join(plot_directory_foci, filename + ".png")
        plot.plot_foci_detection(cyt_mip_contrast,
                                 spots=None,
                                 foci=foci,
                                 radius_spots_yx=radius_spots[-1],
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
