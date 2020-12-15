# -*- coding: utf-8 -*-

"""
Compute spatial features for every detected cell.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.classification as classification

import numpy as np
import pandas as pd
from utils import Logger
from loader import get_metadata_directory, generate_filename_base


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
    log_directory = args.log_directory

    # input-output directories
    cells_directory = os.path.join(output_directory, "individual_cell")
    dataframe_extraction_directory = os.path.join(output_directory,
                                                  "dataframes",
                                                  "extract_cells")
    dataframe_features_directory = os.path.join(output_directory,
                                                "dataframes", "features")

    # check directories exist
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))

    if not os.path.isdir(cells_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cells_directory))
    if not os.path.isdir(dataframe_extraction_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(dataframe_extraction_directory))
    if not os.path.isdir(dataframe_features_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(dataframe_features_directory))

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

    print("Cells directory: {0}".format(cells_directory))
    print("Dataframes extraction directory: {0}"
          .format(dataframe_extraction_directory))
    print("Dataframes features directory: {0}"
          .format(dataframe_features_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_fov' \n")

    print("Computing features...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_experience = generate_filename_base(experience)

    # load cell extraction dataframe for the related experience
    path = os.path.join(dataframe_extraction_directory,
                        filename_experience + ".csv")
    df_extraction = pd.read_csv(path, sep=';', encoding="utf-8")
    print("\t", "{0}: input dataframe {1}"
          .format(filename_experience, df_extraction.shape))

    # initialize dataframe features
    features_name = classification.get_features_name(
        names_features_distance=True,
        names_features_intranuclear=True,
        names_features_protrusion=True,
        names_features_dispersion=True,
        names_features_topography=True,
        names_features_foci=True,
        names_features_area=True)
    columns = ["cell"] + features_name

    # compute features
    filename_cells = []
    features = []
    for i, row in df_extraction.iterrows():
        # get cell data
        filename_cell = row["cell"]
        filename_cells.append(filename_cell)
        path = os.path.join(cells_directory, filename_cell + ".npz")
        data = np.load(path)
        cyt_coord = data["cyt"]
        nuc_coord = data["nuc"]
        rna_coord = data["rna"]
        rna_coord = np.unique(rna_coord, axis=0)

        # compute features
        features_cell = classification.get_features(
            cyt_coord,
            nuc_coord,
            rna_coord,
            compute_distance=True,
            compute_intranuclear=True,
            compute_protrusion=True,
            compute_dispersion=True,
            compute_topography=True,
            compute_foci=True,
            compute_area=True)
        features.append(features_cell)

    # gather results
    features = np.stack(features, axis=0)
    df_features = pd.DataFrame(features, columns=features_name)
    df_features["cell"] = filename_cells
    df_features = df_features[columns]

    # save
    path = os.path.join(dataframe_features_directory,
                        filename_experience + ".csv")
    df_features.reset_index(drop=True, inplace=True)
    df_features.to_csv(path, sep=';', header=True, index=False,
                       encoding="utf-8")
    print("Shape of the experience dataframe: {0}".format(df_features.shape))
    print()

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
