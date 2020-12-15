# -*- coding: utf-8 -*-

"""
Extract individual cells results (segmentation and detection) from field of
view analysis.
"""

import os
import argparse
import time
import datetime
import sys

import bigfish.stack as stack
import bigfish.plot as plot

import numpy as np
import pandas as pd
from utils import Logger
from loader import (get_metadata_directory, generate_filename_base,
                    images_generator, image_to_dismiss)


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
    nuc_mask_directory = os.path.join(output_directory, "nuc_mask")
    cyt_mask_directory = os.path.join(output_directory, "cyt_mask")
    cyt_projection_directory = os.path.join(output_directory,
                                            "cyt_projection_mip")
    detection_directory = os.path.join(output_directory,
                                       "transcription_site_removal")
    cells_directory = os.path.join(output_directory, "individual_cell")
    dataframe_directory = os.path.join(output_directory, "dataframes",
                                       "extract_cells")

    # plot directories
    plot_individual_cell_directory = os.path.join(output_directory, "plot",
                                                  "individual_cell")

    # check directories exist
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(output_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(output_directory))

    if not os.path.isdir(nuc_mask_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(nuc_mask_directory))
    if not os.path.isdir(cyt_mask_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cyt_mask_directory))
    if not os.path.isdir(cyt_projection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cyt_projection_directory))
    if not os.path.isdir(detection_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(detection_directory))
    if not os.path.isdir(cells_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(cells_directory))
    if not os.path.isdir(dataframe_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(dataframe_directory))
    if not os.path.isdir(plot_individual_cell_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(plot_individual_cell_directory))

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
    print("Cytoplasm mask directory: {0}".format(cyt_mask_directory))
    print("Cytoplasm projection directory: {0}"
          .format(cyt_projection_directory))
    print("Detection directory: {0}".format(detection_directory))
    print("Cells directory: {0}".format(cells_directory))
    print("Dataframes directory: {0}".format(dataframe_directory))
    print("Plot individual cell directory: {0}"
          .format(plot_individual_cell_directory))

    print("Log directory: {0}".format(log_directory))
    print("Log file: {0}".format(log_file.split("/")[-1]))
    print("Date: {0}".format(date), "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Extracting cells...")

    # start analysis
    experience = get_metadata_directory(experience_directory)
    filename_base = generate_filename_base(experience)
    generator = images_generator(base_directory, experience_directory,
                                 return_image=False)

    # initialize dataframe experience
    df_experience = pd.DataFrame({"id_cell": [],
                                  "cell": [],
                                  "image": [],
                                  "experience": [],
                                  "gene": [],
                                  "author": [],
                                  "puromycin": [],
                                  "drug": [],
                                  "paper": [],
                                  "nb_rna": [],
                                  "nb_rna_in_foci": [],
                                  "nb_rna_out_foci": [],
                                  "nb_foci": []})

    nb_images = 0
    for i, _ in enumerate(generator):
        filename = filename_base + "_" + str(i)
        gene = experience["gene"]
        author = experience["author"]
        puro = experience["puro"]
        drug = experience["other"]
        paper = experience["paper"]
        if image_to_dismiss(experience_directory, i):
            print("\t", filename, "DISMISSED")
            continue
        print("\t", filename)

        # masks
        path = os.path.join(nuc_mask_directory, filename + ".png")
        mask_nuc = stack.read_image(path)
        path = os.path.join(cyt_mask_directory, filename + ".png")
        mask_cyt = stack.read_image(path)

        # cytoplasm maximum projection
        path = os.path.join(cyt_projection_directory, filename + ".png")
        cyt_mip = stack.read_image(path)

        # spots and foci
        path = os.path.join(detection_directory, filename + ".npz")
        data = np.load(path)
        spots_out_foci = data["spots_out_foci"]
        spots_in_foci = data["spots_in_foci"]
        foci = data["foci"]

        # extract coordinates
        results = stack.extract_coordinates_image(cyt_labelled=mask_cyt,
                                                  nuc_labelled=mask_nuc,
                                                  spots_out=spots_out_foci,
                                                  spots_in=spots_in_foci,
                                                  foci=foci)
        if len(results) == 0:
            print("\t \t no cell extracted !")
            continue

        # save individual cells results
        id_cells = []
        nb_rna_cells = []
        nb_rna_in_foci_cells = []
        nb_rna_out_foci_cells = []
        nb_foci_cells = []
        filename_cells = []
        for i_cell, results_cell in enumerate(results):

            # get results per cell
            filename_cell = filename + "_cell_{0}".format(i_cell)
            cyt_coord, nuc_coord, rna_coord, cell_foci, frame = results_cell
            nb_rna = len(rna_coord)
            nb_rna_in_foci = len(rna_coord[rna_coord[:, 3] != -1, :])
            nb_rna_out_foci = len(rna_coord[rna_coord[:, 3] == -1, :])
            nb_foci = len(cell_foci)

            # save result for the dataframes
            filename_cells.append(filename_cell)
            id_cells.append(i_cell)
            nb_rna_cells.append(nb_rna)
            nb_rna_in_foci_cells.append(nb_rna_in_foci)
            nb_rna_out_foci_cells.append(nb_rna_out_foci)
            nb_foci_cells.append(nb_foci)

            # save result per cell
            (min_y, min_x, max_y, max_x) = frame
            image_cyt = cyt_mip[min_y: max_y, min_x: max_x]
            mask_cyt_cell = mask_cyt[min_y: max_y, min_x: max_x]
            mask_nuc_cell = mask_nuc[min_y: max_y, min_x: max_x]
            path = os.path.join(cells_directory, filename_cell)
            np.savez(path,
                     cyt=cyt_coord,
                     nuc=nuc_coord,
                     rna=rna_coord,
                     foci=cell_foci,
                     cell=frame,
                     image_cyt=image_cyt,
                     mask_cyt=mask_cyt_cell,
                     mask_nuc=mask_nuc_cell)

            # plot individual cell
            path = os.path.join(plot_individual_cell_directory,
                                filename_cell + ".png")
            plot.plot_cell(cyt_coord=cyt_coord,
                           nuc_coord=nuc_coord,
                           rna_coord=rna_coord,
                           foci_coord=cell_foci,
                           image_cyt=image_cyt,
                           mask_cyt=mask_cyt_cell,
                           mask_nuc=mask_nuc_cell,
                           count_rna=True,
                           title=filename_cell,
                           remove_frame=False,
                           rescale=True,
                           framesize=(15, 10),
                           path_output=path,
                           ext="png",
                           show=False)

        # save information in a dataframe
        gene_cells = [gene] * len(id_cells)
        author_cells = [author] * len(id_cells)
        puro_cells = [puro] * len(id_cells)
        drug_cells = [drug] * len(id_cells)
        paper_cells = [paper] * len(id_cells)
        image_cells = [filename] * len(id_cells)
        experience_cells = [filename_base] * len(id_cells)
        df_image = pd.DataFrame({"id_cell": id_cells,
                                 "cell": filename_cells,
                                 "image": image_cells,
                                 "experience": experience_cells,
                                 "gene": gene_cells,
                                 "author": author_cells,
                                 "puromycin": puro_cells,
                                 "drug": drug_cells,
                                 "paper": paper_cells,
                                 "nb_rna": nb_rna_cells,
                                 "nb_rna_in_foci": nb_rna_in_foci_cells,
                                 "nb_rna_out_foci": nb_rna_out_foci_cells,
                                 "nb_foci": nb_foci_cells})
        df_experience = pd.concat([df_experience, df_image])

        nb_images += 1

    # save dataframe experience
    path = os.path.join(dataframe_directory, filename_base + ".csv")
    df_experience.reset_index(drop=True, inplace=True)
    df_experience.to_csv(path,
                         sep=';',
                         header=True,
                         index=False,
                         encoding="utf-8")
    print()
    print("Shape of the experience dataframe: {0}".format(df_experience.shape))

    print()
    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
