# -*- coding: utf-8 -*-

"""
Segment cells from projected FISH images.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
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
                        help="Segmentation threshold.",
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
    mask_nuc_directory = os.path.join(output_directory, "nuc_mask")
    mask_cyt_directory = os.path.join(output_directory, "cyt_mask")
    projection_cyt_directory = os.path.join(output_directory,
                                            "cyt_projection_focus")
    plot_directory = os.path.join(output_directory, "plot", "segmentation")
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(mask_nuc_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_nuc_directory))
    if not os.path.isdir(mask_cyt_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(mask_cyt_directory))
    if not os.path.isdir(projection_cyt_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(projection_cyt_directory))
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
    print("Segmentation threshold: {0}".format(threshold))
    print("Output directory: {0}".format(output_directory))
    print("Mask nuclei directory: {0}".format(mask_nuc_directory))
    print("Mask cytoplasm directory: {0}".format(mask_cyt_directory))
    print("Projection cytoplasm directory: {0}"
          .format(projection_cyt_directory))
    print("Plot directory: {0}".format(plot_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Segmenting images...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory,
                                 return_image=False)

    nb_images = 0
    for i, _ in enumerate(generator):
        filename = filename_base + "_" + str(i)
        print("\t", filename)

        # cyt focus projection
        path = os.path.join(projection_cyt_directory, filename + ".png")
        cyt_projected = stack.read_image(path)
        cyt_projected_contrast = stack.rescale(cyt_projected,
                                               channel_to_stretch=0)

        # nuclei labelled
        path = os.path.join(mask_nuc_directory, filename + ".png")
        nuc_labelled = stack.read_image(path)

        # compute binary mask
        mask = segmentation.build_cyt_binary_mask(cyt_projected,
                                                  threshold=threshold)
        mask[nuc_labelled > 0] = True

        # compute relief
        relief = segmentation.build_cyt_relief(cyt_projected,
                                               nuc_labelled=nuc_labelled,
                                               mask_cyt=mask,
                                               alpha=0.99)

        # segment cytoplasm
        cyt_mask = segmentation.cyt_watershed(relief=relief,
                                              nuc_labelled=nuc_labelled,
                                              mask=mask,
                                              smooth=7)

        # save cytoplasm mask
        path = os.path.join(mask_cyt_directory, filename + ".png")
        stack.save_image(cyt_mask, path)

        # save plot segmentation
        path = os.path.join(plot_directory, filename + ".png")
        plot.plot_segmentation_boundary(tensor=cyt_projected_contrast,
                                        mask_nuc=nuc_labelled,
                                        mask_cyt=cyt_mask,
                                        framesize=(10, 10),
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
