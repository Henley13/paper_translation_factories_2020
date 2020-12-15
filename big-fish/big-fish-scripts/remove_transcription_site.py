# -*- coding: utf-8 -*-

"""
Remove transcription sites in the FISH image.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
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
    nuc_mask_directory = os.path.join(output_directory, "nuc_mask")
    foci_directory = os.path.join(output_directory, "foci_detection")
    transcription_site_directory = os.path.join(output_directory,
                                                "transcription_site_removal")

    # check directories exist
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))
    if not os.path.isdir(nuc_mask_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(nuc_mask_directory))
    if not os.path.isdir(foci_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(foci_directory))
    if not os.path.isdir(transcription_site_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(transcription_site_directory))
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
    print("Nuclei mask directory: {0}".format(nuc_mask_directory))
    print("Foci directory: {0}".format(foci_directory))
    print("Transcription site removal directory: {0}"
          .format(transcription_site_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Removing transcription sites...")

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
        path = os.path.join(foci_directory, filename + ".npz")
        data = np.load(path)
        clustered_spots = data["clustered_spots"]
        foci = data["foci"]

        # nuclei masks
        path = os.path.join(nuc_mask_directory, filename + ".png")
        mask_nuc = stack.read_image(path)
        nuc = mask_nuc > 0

        # spots out of foci and inside foci
        spots_out_foci = clustered_spots.copy()
        spots_out_foci = spots_out_foci[spots_out_foci[:, 3] == -1, :]
        spots_in_foci = clustered_spots.copy()
        spots_in_foci = spots_in_foci[spots_in_foci[:, 3] != -1, :]

        # remove foci inside nuclei
        spots_in_foci_cleaned, foci_cleaned = stack.remove_transcription_site(
            mask_nuc=nuc,
            spots_in_foci=spots_in_foci,
            foci=foci)

        # save transcription site-free coordinates
        path = os.path.join(transcription_site_directory, filename)
        np.savez(path,
                 spots_out_foci=spots_out_foci,
                 spots_in_foci=spots_in_foci_cleaned,
                 foci=foci_cleaned)

        nb_images += 1

    print()
    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
