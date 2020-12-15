# -*- coding: utf-8 -*-

"""
Decompose cluster of spots in the FISH image.
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
    spot_directory = os.path.join(output_directory, "spot_detection")
    background_filter_directory = os.path.join(output_directory,
                                               "cyt_filter_background")
    cyt_projection_directory = os.path.join(output_directory,
                                            "cyt_projection_mip")
    plot_directory_decomposition = os.path.join(output_directory, "plot",
                                                "cluster_decomposition")
    plot_directory_reference = os.path.join(output_directory, "plot",
                                            "reference_spot")
    decomposition_directory = os.path.join(output_directory,
                                           "cluster_decomposition")
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(spot_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(spot_directory))
    if not os.path.isdir(background_filter_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(background_filter_directory))
    if not os.path.isdir(cyt_projection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cyt_projection_directory))
    if not os.path.isdir(decomposition_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(decomposition_directory))
    if not os.path.isdir(plot_directory_decomposition):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory_decomposition))
    if not os.path.isdir(plot_directory_reference):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_directory_reference))
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
    print("Spot directory: {0}".format(spot_directory))
    print("Background filter directory: {0}"
          .format(background_filter_directory))
    print("Cytoplasm projection directory: {0}"
          .format(cyt_projection_directory))
    print("Plot decomposition directory: {0}"
          .format(plot_directory_decomposition))
    print("Plot reference directory: {0}".format(plot_directory_reference))
    print("Decomposition directory: {0}".format(decomposition_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Decomposing clusters...")

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
        path = os.path.join(spot_directory, filename + ".npz")
        data = np.load(path)
        spots = data["spots"]
        radius = tuple(data["radius"])
        radius_yx = radius[-1]

        # cytoplasm maximum projection
        path = os.path.join(cyt_projection_directory, filename + ".png")
        cyt_mip = stack.read_image(path)
        cyt_mip_contrast = stack.rescale(cyt_mip, channel_to_stretch=0)

        # background filter
        path = os.path.join(background_filter_directory, filename + ".tiff")
        cyt_no_background = stack.read_image(path)

        # cluster decomposition
        (spots_out_cluster,
         spots_in_cluster,
         clusters,
         reference_spot) = detection.run_decomposition(
            cyt_no_background,
            spots,
            radius,
            min_area=2,
            resolution_z=300,
            resolution_yx=103,
            psf_z=350,
            psf_yx=150)

        # save decomposition
        path = os.path.join(decomposition_directory, filename)
        np.savez(path,
                 spots_out_cluster=spots_out_cluster,
                 spots_in_cluster=spots_in_cluster,
                 clusters=clusters,
                 reference_spot=reference_spot,
                 radius_spots=radius)

        # save plot reference spot
        path = os.path.join(plot_directory_reference, filename + ".png")
        image_reference = []
        reference_spot_ = reference_spot / reference_spot.max()
        for j in range(reference_spot_.shape[0]):
            image_reference.append(reference_spot_[j])
        plot.plot_images(image_reference,
                         framesize=(15, int(5 * reference_spot_.shape[0] / 3)),
                         remove_frame=True,
                         path_output=path,
                         ext="png",
                         show=False)

        # save plot decomposition
        all_spots = np.concatenate(
            (spots_out_cluster, spots_in_cluster[:, :3]),
            axis=0)
        path = os.path.join(plot_directory_decomposition, filename + ".png")
        plot.plot_spot_detection(tensor=cyt_mip_contrast,
                                 spots=all_spots,
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
