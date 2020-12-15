# -*- coding: utf-8 -*-

"""
Merge dataframe returned from each experience during the step of features
computation.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.classification as classification
import pandas as pd

from utils import Logger
from loader import (get_experiences, generate_filename_base)


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
    log_directory = args.log_directory

    # input-output directories
    dataframe_directory = os.path.join(output_directory, "dataframes")
    dataframe_features_directory = os.path.join(output_directory,
                                                "dataframes",
                                                "features")

    # check directories exist
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))

    if not os.path.isdir(dataframe_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(dataframe_directory))
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
        log_directory, "log" + "_merge_features_dataframe")
    sys.stdout = Logger(log_file)

    print("Running {0} file...".format(os.path.basename(__file__)), "\n")
    start_time = time.time()

    print("Output directory: {0}".format(output_directory))
    print("Dataframes directory: {0}".format(dataframe_directory))
    print("Dataframes features directory: {0}"
          .format(dataframe_features_directory))
    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Merging dataframes saved during features computation...")

    # initialize dataframe
    features_name = classification.get_features_name(
        names_features_distance=True,
        names_features_intranuclear=True,
        names_features_protrusion=True,
        names_features_dispersion=True,
        names_features_topography=True,
        names_features_foci=True,
        names_features_area=True)
    columns = ["cell"] + features_name
    df_features = pd.DataFrame(columns=columns)

    nb_dataframes = 0
    experiences = get_experiences()
    for experience in experiences:

        # load dataframe
        filename = generate_filename_base(experience)
        path = os.path.join(dataframe_features_directory, filename + ".csv")
        if not os.path.exists(path):
            print("\t", filename, "does not exist")
            continue
        df_experience = pd.read_csv(path, sep=';', encoding="utf-8")
        print("\t", "{0}: {1}".format(filename, df_experience.shape))

        # concatenate dataframes
        df_features = pd.concat([df_features, df_experience])

        nb_dataframes += 1

    # save features dataframe
    path = os.path.join(dataframe_directory, "features.csv")
    df_features.reset_index(drop=True, inplace=True)
    df_features.to_csv(path,
                       sep=';',
                       header=True,
                       index=False,
                       encoding="utf-8")
    print()
    print("Shape of the final features dataframe: {0}"
          .format(df_features.shape))
    print("Columns:")
    columns_name = df_features.columns
    for col in columns_name:
        print("\t", col)

    print()
    print("Done ({0} dataframes)!".format(nb_dataframes), "\n")

    # merge and save final dataframe
    path = os.path.join(dataframe_directory, "cells.csv")
    df_cell = pd.read_csv(path, sep=';', encoding="utf-8")
    df = pd.merge(left=df_cell, right=df_features, how='inner', on=["cell"])
    path = os.path.join(dataframe_directory, "merged.csv")
    df.reset_index(drop=True, inplace=True)
    df.to_csv(path, sep=';', header=True, index=False, encoding="utf-8")
    print()
    print("Shape of the merged dataframe: {0}".format(df.shape))
    print("Columns:")
    columns_name = df.columns
    for col in columns_name:
        print("\t", col)

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
