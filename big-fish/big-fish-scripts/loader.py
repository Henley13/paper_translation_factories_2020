# -*- coding: utf-8 -*-

"""
Parse files and read images.
"""

import os
import argparse
import time
import sys
import datetime

import numpy as np

import bigfish.stack as stack
from utils import Logger


def _all_experiences():
    """List all the experiences made.

    :return:
    experiments : List[dict]
        A list of dictionaries with metadata about experiences (gene, author,
        paper, usage of puromycin, other manipulation, batch, directory).

    """
    experiments = [

        # Aubin's data

        {"gene": "DYNC1H1",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w1_bac_dync1h1_RachaAuto"},

        {"gene": "BUB1",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w4_bac_bub1_RachaPlus"},

        {"gene": "ATP6A2",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w3_bac_atp6ap2_Racha"},

        {"gene": "SPEN",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w1_bac_spen_Racha"},

        {"gene": "KIF1C",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w3_Kif1c_ctrl"},

        {"gene": "RAB13",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w1_Rab13_ctrl"},

        {"gene": "KIF20B",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w7_bac_kif20b_Racha"},

        {"gene": "MYO18A",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w8_bac_myo18a_Racha"},

        {"gene": "PAK2",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w5_bac_pak2_Racha"},

        {"gene": "CEP192",
         "author": None,
         "puro": False,
         "paper": "aubin",
         "other": None,
         "batch": 1,
         "directory": "w11_bac_cep192_Abed"},

        # Racha's data

        {"gene": "AKAP1",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AKAP1_Already_Analyzed_by_Florian_First_paper_version"},

        {"gene": "AKAP1",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AKAP1_4092_060319_Adham"},

        {"gene": "AP1S2",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AP1S2_3311_180219_Adham"},

        {"gene": "FLNA",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "FLNA_Edouard"},

        {"gene": "CRKL",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CRKL_5551_170419_Adham"},

        {"gene": "CRKL",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CRKL_Edouard_1"},

        {"gene": "CRKL",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "CRKL_Edouard_2"},

        {"gene": "CEP170P1",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CEP170P1_Edouard_1"},

        {"gene": "CEP170P1",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "CEP170P1_Edouard_2"},

        {"gene": "CEP170P1",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 3,
         "directory": "CEP170P1_Edouard_3"},

        {"gene": "AURKA",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AURKA_2476_240119_Adham"},

        {"gene": "AURKA",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AURKA_Edouard_1"},

        {"gene": "AURKA",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "AURKA_Edouard_2"},

        {"gene": "AURKA",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 3,
         "directory": "AURKA_Edouard_3"},

        {"gene": "AKAP9",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AKAP9_Already_Analyzed_by_Florian_First_paper_version"},

        {"gene": "MYSNE2",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "mSYNE2_4430_240119_Adham"},

        {"gene": "PLEC",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "PLEC_8150_240119_Adham"},

        {"gene": "HSP90B1",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "HSP90B1_3509_240119_Adham"},

        {"gene": "BUB1",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "BUB1_2334_240119_Adham"},

        {"gene": "HMMR",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "HMMR_6087_220119_Adham"},

        {"gene": "ASPM",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "ASPM_3346_220119_Adham"},

        {"gene": "CTNNB1",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CTNNB1_2019_EB"},

        {"gene": "KIF4A",
         "author": "edouard",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "KIF4A_2019_EB"},

        {"gene": "KIF5B",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "KIF5B"},

        {"gene": "MYH3H",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "MYH3_2019"},

        {"gene": "DYNLL2",
         "author": "adham",
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "DYNLL2_3137_230419_Adham"},

        {"gene": "KIF1C",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "KIF1C"},

        # Puromycin experiments

        {"gene": "KIF4A",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "KIF4A_puro"},

        {"gene": "AURKA",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AURKA_puro_1"},

        {"gene": "AURKA",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "AURKA_puro_2"},

        {"gene": "AURKA",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 3,
         "directory": "AURKA_puro_3"},

        {"gene": "DYNC1H1",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "DYNC1H1_puro"},

        {"gene": "BUB1",
         "author": "adham",
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "BUB1_2334_230419_puro_Adham"},

        {"gene": "CTNNB1",
         "author": "adham",
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CTNNB1_3843_210419_puro_Adham"},

        {"gene": "CTNNB1",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "CTNNB1_Puro_1h"},

        {"gene": "KIF1C",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "KIF1C_puro_1h_"},

        {"gene": "ASPM",
         "author": "adham",
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "ASPM_3346_220419_puro_Adham"},

        {"gene": "HMMR",
         "author": "adham",
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "HMMR_6087_170419_puro_Adham"},

        {"gene": "MYH3H",
         "author": "adham",
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "MYH3_2279_230419_puro_Adham"},

        {"gene": "MYH3H",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "MYH3_Puro"},

        {"gene": "KIF1C",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 2,
         "directory": "KIF1C_puro"},

        {"gene": "AKAP1",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AKAP1_puro_1h"},

        {"gene": "AKAP9",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "Akap9_puro_1h"},

        {"gene": "AP1S2",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "AP1S2_puro_1h"},

        {"gene": "HSP90B1",
         "author": None,
         "puro": True,
         "paper": "racha",
         "other": None,
         "batch": 1,
         "directory": "HSP90_puro_1h"},

        # Others

        {"gene": "CTNNB1",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": "LG007",
         "batch": 1,
         "directory": "293T_CTNNB1_007_2h"},

        {"gene": "CTNNB1",
         "author": None,
         "puro": False,
         "paper": "racha",
         "other": "DMSO",
         "batch": 1,
         "directory": "293T_CTNNB1_contDMSO"}

    ]

    return experiments


def _get_gene(gene, experiences="all"):
    """Filter experiences involving a specific gene.

    :param gene: str
        Name of the gene for which we want to filter the experiences.
    :param experiences: List[dict]
        A list of dictionaries with metadata about experiences we need to
        filter.

    :return:

    experiences_to_keep : List[dict]
        A list of dictionaries with metadata about experiences involving the
        requested gene (gene, author, paper, usage of puromycin, other
        manipulation, batch, directory).

    """
    if experiences == "all":
        experiences_to_filter = _all_experiences()
    else:
        experiences_to_filter = experiences

    experiences_to_keep = []
    for experience in experiences_to_filter:
        if experience["gene"] == gene:
            experiences_to_keep.append(experience)

    return experiences_to_keep


def _get_paper(paper, experiences="all"):
    """Filter experience used for a specific paper

    :param paper: str
        Paper targeted.
    :param experiences: List[dict]
        A list of dictionaries with metadata about experiences we need to
        filter.

    :return:

    experiences_to_keep : List[dict]
        A list of dictionaries with metadata about experiences involving the
        requested paper (gene, author, paper, usage of puromycin, other
        manipulation, batch, directory).

    """
    if experiences == "all":
        experiences_to_filter = _all_experiences()
    else:
        experiences_to_filter = experiences

    experiences_to_keep = []
    for experience in experiences_to_filter:
        if experience["paper"] == paper:
            experiences_to_keep.append(experience)

    return experiences_to_keep


def _get_puro(puro, experiences="all"):
    """Filter experience with the use of puromycin.

    :param puro: bool
        Use of puromycin or not.
    :param experiences: List[dict]
        A list of dictionaries with metadata about experiences we need to
        filter.

    :return:

    experiences_to_keep : List[dict]
        A list of dictionaries with metadata about experiences involving
        puromycin or not (gene, author, paper, usage of puromycin, other
        manipulation, batch, directory).

    """
    if experiences == "all":
        experiences_to_filter = _all_experiences()
    else:
        experiences_to_filter = experiences

    experiences_to_keep = []
    for experience in experiences_to_filter:
        if experience["puro"] is puro:
            experiences_to_keep.append(experience)

    return experiences_to_keep


def get_experiences(gene=None, paper=None, puro=None, experiences="all"):
    """Filter experiences.

    :param gene: str
        Name of the gene for which we want to filter the experiences.
    :param paper: str
        Paper targeted.
    :param puro: bool
        Use of puromycin or not.
    :param experiences: List[dict]
        A list of dictionaries with metadata about experiences we need to
        filter.

    :return:

    experiences_to_keep : List[dict]
        A list of dictionaries with metadata about experiences involving
        puromycin or not (gene, author, paper, usage of puromycin, other
        manipulation, batch, directory).

    """
    if experiences == "all":
        experiences_to_keep = _all_experiences()
    else:
        experiences_to_keep = experiences

    if gene is not None:
        experiences_to_keep = _get_gene(gene, experiences_to_keep)

    if paper is not None:
        experiences_to_keep = _get_paper(paper, experiences_to_keep)

    if puro is not None:
        experiences_to_keep = _get_puro(puro, experiences_to_keep)

    return experiences_to_keep


def get_directories(gene=None, paper=None, puro=None, experiences="all"):
    """Get list of directories name.

    :param gene: str
        Name of the gene for which we want to filter the experiences.
    :param paper: str
        Paper targeted.
    :param puro: bool
        Use of puromycin or not.
    :param experiences: List[dict]
        A list of dictionaries with metadata about experiences we need to
        filter.

    :return:

    directories : List[str]
        A list of directories name.

    """
    experiences = get_experiences(gene, paper, puro, experiences)
    directories = []
    for experience in experiences:
        experience_directory = experience["directory"]
        directories.append(experience_directory)

    return directories


def get_metadata_directory(data_directory):
    """Get metadata about a specific directory.

    :param data_directory: str
        Targeted directory.

    :return:

    experience : dict
        Dictionary with metadata bout the directory (gene, author, use of
        puromycin, paper where the data is analyzed, other use of drugs, batch,
        directory).

    """
    experiences = get_experiences()
    for experience in experiences:
        experience_directory = experience["directory"]
        if experience_directory == data_directory:
            return experience

    return None


def _map_aubin(base_directory, directory):
    """Get a map to read the data from a specific directory.

    :param base_directory: str
        Path of the main data directory.
    :param directory: str
        Name of the directory to inspect.

    :return:

    data_map : List[tuple]
        A list of tuple with a recipe and the path of the directory.

    """
    if directory == "w1_bac_dync1h1_RachaAuto":
        opt = "w1_bac_dync1h1_296"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10', 'p11', 'p12', 'p13', 'p14', 'p15', 'p16', 'p17', 'p18',
                'p19', 'p20', 'p21', 'p22', 'p24', 'p25', 'p26', 'p27',
                'p28', 'p29', 'p30', 'p31', 'p32']
        c = ["DAPI", "Cy3"]

    elif directory == "w4_bac_bub1_RachaPlus":
        opt = "w4_bac_bub1_2334"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10', 'p11']
        c = ["dapi", "Cy3"]

    elif directory == "w3_bac_atp6ap2_Racha":
        opt = "w3_bac_atp6ap2_3304"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10']
        c = ["dapi", "Cy3"]

    elif directory == "w1_bac_spen_Racha":
        opt = "w1_bac_spen_3451"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10', 'p11']
        c = ["dapi", "Cy3"]

    elif directory == "w3_Kif1c_ctrl":
        opt = "w3_HeLa_Kif1c_2"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10', 'p11']
        c = ["dapi", "Cy3"]

    elif directory == "w1_Rab13_ctrl":
        opt = "w1_HeLa_Rab13_1"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10', 'p11', 'p12']
        c = ["dapi", "Cy3"]

    elif directory == "w7_bac_kif20b_Racha":
        opt = "w7_bac_kif20b_5120"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10']
        c = ["dapi", "cy3"]

    elif directory == "w8_bac_myo18a_Racha":
        opt = "w8_bac_myo18a_3107"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10']
        c = ["dapi", "cy3"]

    elif directory == "w5_bac_pak2_Racha":
        opt = "w5_bac_pak2_4061"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10']
        c = ["dapi", "Cy3"]

    elif directory == "w11_bac_cep192_Abed":
        opt = "w11_bac_cep192_2229"
        fovs = ['p01', 'p02', 'p03', 'p04', 'p05', 'p06', 'p07', 'p08', 'p09',
                'p10']
        c = ["DAPI", "Cy3"]

    else:
        return None

    recipe = {"fov": fovs,
              "c": c,
              "ext": "TIF",
              "opt": opt,
              "pattern": "opt_fov_c.ext"}

    data_map = [(recipe, os.path.join(base_directory, directory))]

    return data_map


def _map_adham(base_directory, directory):
    """Get a map to read the data from a specific directory.

    :param base_directory: str
        Path of the main data directory.
    :param directory: str
        Name of the directory to inspect.

    :return:

    data_map : List[tuple]
        A list of tuple with a recipe and the path of the directory.

    """
    fovs = ["P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09",
            "P10"]

    if directory == "AKAP1_4092_060319_Adham":
        opt = "AKAP1_4092_060319_Adham"

    elif directory == "AP1S2_3311_180219_Adham":
        opt = "AP1S2_3311_180219_Adham"

    elif directory == "CRKL_5551_170419_Adham":
        opt = "CRKL_5551_170419_Adham"

    elif directory == "AURKA_2476_240119_Adham":
        opt = "AURKA_2476_240119_Adham"

    elif directory == "mSYNE2_4430_240119_Adham":
        fovs = ["P01", "P02", "P03", "P05", "P06", "P08", "P09", "P10"]
        opt = "mSYNE2_4430_240119_Adham"

    elif directory == "PLEC_8150_240119_Adham":
        fovs = ["P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09"]
        opt = "PLEC_8150_240119_Adham"

    elif directory == "HSP90B1_3509_240119_Adham":
        opt = "HSP90B1_3509_240119_Adham"

    elif directory == "BUB1_2334_240119_Adham":
        opt = "BUB1_2334_240119_Adham"

    elif directory == "HMMR_6087_220119_Adham":
        opt = "HMMR_6087_220119_Adham"

    elif directory == "ASPM_3346_220119_Adham":
        opt = "ASPM_3346_220119_Adham"

    elif directory == "DYNLL2_3137_230419_Adham":
        opt = "DYNLL2_3137_230419_Adham"

    elif directory == "BUB1_2334_230419_puro_Adham":
        opt = "BUB1_2334_230419_puro_Adham"

    elif directory == "CTNNB1_3843_210419_puro_Adham":
        opt = "CTNNB1_3843_210419_puro_Adham"

    elif directory == "ASPM_3346_220419_puro_Adham":
        opt = "ASPM_3346_220419_puro_Adham"

    elif directory == "HMMR_6087_170419_puro_Adham":
        opt = "HMMR_6087_170419_puro_Adham"

    elif directory == "MYH3_2279_230419_puro_Adham":
        opt = "MYH3_2279_230419_puro_Adham"

    else:
        return None

    recipe = {"fov": fovs,
              "c": ["w1DAPI", "w2DsRed"],
              "ext": "TIF",
              "opt": opt,
              "pattern": "opt_fov_c.ext"}
    data_map = [(recipe, os.path.join(base_directory, directory))]

    return data_map


def _map_other(base_directory, directory):
    """Get a map to read the data from a specific directory.

    :param base_directory: str
        Path of the main data directory.
    :param directory: str
        Name of the directory to inspect.

    :return:

    data_map : List[tuple]
        A list of tuple with a recipe and the path of the directory.

    """
    if directory == "AKAP1_Already_Analyzed_by_Florian_First_paper_version":
        recipe = {"fov": ["1", "2", "3", "4", "5"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "FLNA_Edouard":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CRKL_Edouard_1":
        recipe = {"fov": ["1", "2", "3"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CRKL_Edouard_2":
        recipe = {"fov": ["1", "2", "3", "4"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CEP170P1_Edouard_1":
        recipe = {"fov": ["1", "2"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CEP170P1_Edouard_2":
        recipe = {"fov": ["1", "2"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CEP170P1_Edouard_3":
        recipe = {"fov": ["1", "2", "3", "4", "5"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_Edouard_1":
        recipe = {"fov": ["1", "2"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_Edouard_2":
        recipe = {"fov": ["1"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_Edouard_3":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7", "8"],
                  "c": ["Hoechst", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AKAP9_Already_Analyzed_by_Florian_First_paper_version":
        recipe = {"fov": ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9",
                          "s10"],
                  "c": ["w3DAPI", "w1DsRed"],
                  "ext": "TIF",
                  "opt": "BAC4231_AKAP9_Endogene_NT_1",
                  "pattern": "opt_c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CTNNB1_2019_EB":
        recipe = {"fov": ["1", "2", "3", "4", "10", "11", "12", "13", "14",
                          "15"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "KIF4A_2019_EB":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                          "11", "12", "13", "14", "15", "16", "17", "18", "19",
                          "20"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory, "KIF4A"))]

    elif directory == "KIF5B":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "MYH3_2019":
        recipe = {"fov": ["1", "3", "4", "5", "6", "7", "8", "9", "10", "11",
                          "12"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "293T_CTNNB1_007_2h":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                          "11", "12", "13", "14", "15"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "293T_CTNNB1_contDMSO":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7", "8", "9"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "KIF4A_puro":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_puro_1":
        recipe = {"fov": ["1", "2", "3"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_puro_2":
        recipe = {"fov": ["2", "3", "4", "5", "6"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AURKA_puro_3":
        recipe = {"fov": ["7", "8", "9", "10"],
                  "c": ["Dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "DYNC1H1_puro":
        recipe = {"fov": ["1", "2", "3", "4"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "CTNNB1_Puro_1h":
        recipe = {"fov": ["1", "2", "3"],
                  "c": ["Hoechst", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "KIF1C_puro_1h_":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "MYH3_Puro":
        recipe = {"fov": ["1", "2"],
                  "c": ["Hoechst", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "KIF1C_puro":
        recipe = {"fov": ["", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                          "10", "11", "12", "13", "14", "15", "16", "17", "18",
                          "19", "20", "21", "22", "23", "24", "25", "26", "27",
                          "28", "29", "30"],
                  "c": ["w3Hoechst", "w1Cy3"],
                  "opt": "HeLaKYOTO-FISHendo-LNA-KIF1C-PURO",
                  "ext": "TIF",
                  "pattern": "opt_fov_c.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "KIF1C":
        recipe = {"fov": ["", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                          "10", "11", "12", "13", "14", "15", "16", "17", "18",
                          "19", "20", "21", "22", "23", "24", "25", "26", "27",
                          "28", "29", "30"],
                  "c": ["w3Hoechst", "w1Cy3"],
                  "opt": "HeLaKYOTO-FISHendo-LNA-KIF1C",
                  "ext": "TIF",
                  "pattern": "opt_fov_c.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AKAP1_puro_1h":
        recipe = {"fov": ["1", "2", "4", "6", "7"],
                  "c": ["dapi", "Cy3"],
                  "ext": "tif",
                  "pattern": "c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "Akap9_puro_1h":
        recipe = {"fov": ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9",
                          "s10"],
                  "c": ["w3DAPI", "w1DsRed"],
                  "opt": "BAC4231_AKAP9_Endogne_puro30min_1",
                  "ext": "tif",
                  "pattern": "opt_c_fov.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "AP1S2_puro_1h":
        recipe = {"fov": ["1", "2", "3", "4", "5", "6"],
                  "c": ["w1Hoechst", "w3Cy3"],
                  "opt": "AP1S2_puro60",
                  "ext": "TIF",
                  "pattern": "opt_fov_c.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    elif directory == "HSP90_puro_1h":
        recipe = {"fov": ["", "1", "2", "3"],
                  "c": ["w3DAPI", "w1DsRed"],
                  "opt": "HSP90B1_puro60",
                  "ext": "TIF",
                  "pattern": "opt_fov_c.ext"}
        data_map = [(recipe, os.path.join(base_directory, directory))]

    else:
        return None

    return data_map


def map_folder_recipe(base_directory, directory):
    """Get a map to read the data from a specific directory.

    :param base_directory: str
        Path of the main data directory.
    :param directory: str
        Name of the directory to inspect.

    :return:

    data_map : List[tuple]
        A list of tuple with a recipe and the path of the directory.

    """
    if directory in ['w1_bac_dync1h1_RachaAuto', 'w4_bac_bub1_RachaPlus',
                     'w3_bac_atp6ap2_Racha', 'w1_bac_spen_Racha',
                     'w3_Kif1c_ctrl', 'w1_Rab13_ctrl', 'w7_bac_kif20b_Racha',
                     'w8_bac_myo18a_Racha', 'w5_bac_pak2_Racha',
                     'w11_bac_cep192_Abed']:
        data_map = _map_aubin(base_directory, directory)

    elif directory in ['AKAP1_4092_060319_Adham', 'AP1S2_3311_180219_Adham',
                       'CRKL_5551_170419_Adham', 'AURKA_2476_240119_Adham',
                       'mSYNE2_4430_240119_Adham', 'PLEC_8150_240119_Adham',
                       'HSP90B1_3509_240119_Adham', 'BUB1_2334_240119_Adham',
                       'HMMR_6087_220119_Adham', 'ASPM_3346_220119_Adham',
                       'DYNLL2_3137_230419_Adham',
                       'BUB1_2334_230419_puro_Adham',
                       'CTNNB1_3843_210419_puro_Adham',
                       'ASPM_3346_220419_puro_Adham',
                       'HMMR_6087_170419_puro_Adham',
                       'MYH3_2279_230419_puro_Adham']:
        data_map = _map_adham(base_directory, directory)

    elif directory in ["AKAP1_Already_Analyzed_by_Florian_First_paper_version",
                       "FLNA_Edouard", "CRKL_Edouard_1", "CRKL_Edouard_2",
                       "CEP170P1_Edouard_1", "CEP170P1_Edouard_2",
                       "CEP170P1_Edouard_3", "AURKA_Edouard_1",
                       "AURKA_Edouard_2", "AURKA_Edouard_3", "CTNNB1_2019_EB",
                       "KIF4A_2019_EB",
                       "AKAP9_Already_Analyzed_by_Florian_First_paper_version",
                       "KIF5B", "MYH3_2019", "293T_CTNNB1_007_2h",
                       "293T_CTNNB1_contDMSO", "KIF4A_puro", "AURKA_puro_1",
                       "AURKA_puro_2", "AURKA_puro_3", "DYNC1H1_puro",
                       "CTNNB1_Puro_1h", "KIF1C_puro_1h_", "MYH3_Puro",
                       "KIF1C_puro", "KIF1C", "AKAP1_puro_1h", "Akap9_puro_1h",
                       "AP1S2_puro_1h", "HSP90_puro_1h"]:
        data_map = _map_other(base_directory, directory)

    else:
        data_map = None

    return data_map


def check_files_exist(base_directory):
    """Check every files described in recipes exist.

    :param base_directory: str
        Path of the main data directory.

    :return:

    """
    experiences = get_experiences()
    for experience in experiences:
        experience_directory = experience["directory"]
        data_map = map_folder_recipe(base_directory, experience_directory)
        if data_map is None:
            continue
        for (recipe, directory) in data_map:
            stack.check_recipe(recipe, directory)

    return


def generate_filename_base(experience):
    """Generate the base of the filename (without the id of the fov).

    :param experience: dict
        Dictionary with metadata bout the directory (gene, author, use of
        puromycin, paper where the data is analyzed, other use of drugs, batch,
        directory).

    :return:

    filename_base : str
        Base of the filename used to save the output.

    """
    gene = experience["gene"]
    author = experience["author"]
    if author is None:
        author_str = "unknown"
    else:
        author_str = author
    puro = experience["puro"]
    if puro:
        puro_str = "puro"
    else:
        puro_str = "nopuro"
    paper = experience["paper"]
    other = experience["other"]
    if other is None:
        other_str = "nodrug"
    else:
        other_str = other
    batch = experience["batch"]
    filename_base = "{0}_{1}_{2}_{3}_{4}_{5}".format(
        gene, author_str, puro_str, paper, other_str, batch)

    return filename_base


def images_generator(base_directory, data_directory, return_image=True):
    """Generate images from every input directories.

    :param base_directory: str
        Path of the main data directory.
    :param data_directory: str
        Input directory to parse.
    :param return_image: bool
        Load the image.

    :return:

    image : np.ndarray, np.uint
        Image with shape (r, c, z, y, x).

    """
    data_map = map_folder_recipe(base_directory, data_directory)
    for recipe, input_folder in data_map:
        nb_fov = stack.utils.count_nb_fov(recipe)
        for i_fov in range(nb_fov):
            if return_image:
                if _bad_input(input_folder, i_fov):
                    image = _load_bad_input(recipe, input_folder, i_fov)
                else:
                    pass
                    image = stack.build_stack(recipe, input_folder,
                                              input_dimension=3,
                                              i_fov=i_fov,
                                              check=True)
            else:
                image = None

            yield image


def _bad_input(input_directory, i_fov):
    """Check if the input files are good or not.

    :param input_directory: str
        Full path of the input directory.
    :param i_fov: int
        Index of the fov to generate.

    :return: bool
        Is the image directly readable or not ?

    """
    if "CTNNB1_2019_EB" in input_directory and i_fov in [3]:
        return True

    elif "AURKA_Edouard_2" in input_directory and i_fov in [0]:
        return True

    elif ("AURKA_Edouard_3" in input_directory
          and i_fov in [0, 1, 2, 3, 4, 5, 6, 7]):
        return True

    elif "CTNNB1_Puro_1h" in input_directory and i_fov in [0, 1, 2]:
        return True

    elif "AURKA_puro_2" in input_directory and i_fov in [1, 2, 3, 4]:
        return True

    elif "AURKA_puro_3" in input_directory and i_fov in [0, 1, 2, 3]:
        return True

    elif "MYH3_Puro" in input_directory and i_fov in [0, 1]:
        return True

    elif "w7_bac_kif20b_Racha" in input_directory and i_fov in [6]:
        return True

    else:
        return False


def _load_bad_input(recipe, input_directory, i_fov):
    """Load an image badly saved (MIP most of the times).

    :param recipe: dict
        Map the images according to their field of view, their round,
        their channel and their spatial dimensions. Only contain the keys
        'pattern', 'fov', 'r', 'c', 'z', 'ext' or 'opt'.
    :param input_directory: str
        Full path of the input directory.
    :param i_fov: int
        Index of the fov to generate.

    :return:

    image : np.ndarray, np.uint
        Image with shape (r, c, z, y, x).

    """
    recipe = stack.utils.fit_recipe(recipe)
    if "CTNNB1_2019_EB" in input_directory and i_fov in [3]:
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=0)
        nuc = stack.read_image(path)
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=1)
        cyt = stack.read_image(path)
        nb_z = cyt.shape[0]
        nuc = nuc[12, :, :]
        nuc = np.stack([nuc] * nb_z, axis=0)
        image = np.stack([nuc, cyt], axis=0)
        image = image[np.newaxis, ...]

    elif "w7_bac_kif20b_Racha" in input_directory and i_fov in [6]:
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=0)
        nuc = stack.read_image(path)
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=1)
        cyt = stack.read_image(path)
        nb_z = cyt.shape[0]
        nuc = nuc[17, :, :]
        nuc = np.stack([nuc] * nb_z, axis=0)
        image = np.stack([nuc, cyt], axis=0)
        image = image[np.newaxis, ...]

    else:
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=0)
        nuc = stack.read_image(path)
        path = stack.utils.get_path_from_recipe(recipe, input_directory,
                                                fov=i_fov, c=1)
        cyt = stack.read_image(path)
        nb_z = cyt.shape[0]
        nuc = np.stack([nuc] * nb_z, axis=0)
        image = np.stack([nuc, cyt], axis=0)
        image = image[np.newaxis, ...]
    stack.check_array(image,
                      ndim=5,
                      dtype=[np.uint8, np.uint16])

    return image


def image_to_dismiss(input_directory, i_fov):
    """Check if the image should be analyzed or not.

    :param input_directory: str
        Full path of the input directory.
    :param i_fov: int
        Index of the fov to generate.

    :return: bool
        Is the image dismissed ?

    """
    if ("w1_bac_dync1h1_RachaAuto" in input_directory
            and i_fov in [15, 16, 17, 24, 25, 26, 27]):
        return True

    elif "w4_bac_bub1_RachaPlus" in input_directory and i_fov in [0]:
        return True

    elif "w1_bac_spen_Racha" in input_directory and i_fov in [9]:
        return True

    if "w3_Kif1c_ctrl" in input_directory and i_fov in [1, 2, 8, 9, 10]:
        return True

    if "w7_bac_kif20b_Racha" in input_directory and i_fov in [0, 1]:
        return True

    if "w8_bac_myo18a_Racha" in input_directory and i_fov in [5, 9]:
        return True

    if "w5_bac_pak2_Racha" in input_directory and i_fov in [3, 4]:
        return True

    if "w11_bac_cep192_Abed" in input_directory and i_fov in [6]:
        return True

    if ("CRKL_5551_170419_Adham" in input_directory
            and i_fov in [ 0, 1, 2, 7, 8]):
        return True

    if "CRKL_Edouard_1" in input_directory and i_fov in [1, 2]:
        return True

    if "CRKL_Edouard_2" in input_directory and i_fov in [0]:
        return True

    if "CEP170P1_Edouard_1" in input_directory and i_fov in [1]:
        return True

    if "CEP170P1_Edouard_3" in input_directory and i_fov in [0]:
        return True

    if "AURKA_2476_240119_Adham" in input_directory and i_fov in [4]:
        return True

    if "AURKA_Edouard_3" in input_directory and i_fov in [2]:
        return True

    if "PLEC_8150_240119_Adham" in input_directory and i_fov in [6]:
        return True

    if "ASPM_3346_220119_Adham" in input_directory and i_fov in [1, 2]:
        return True

    if ("KIF4A_2019_EB" in input_directory
            and i_fov in [0, 1, 2, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]):
        return True

    if "KIF5B" in input_directory and i_fov in [0, 1, 3]:
        return True

    if "MYH3_2019" in input_directory and i_fov in [1, 2, 3, 4, 9, 10]:
        return True

    if "DYNC1H1_puro" in input_directory and i_fov in [0]:
        return True

    if "CTNNB1_3843_210419_puro_Adham" in input_directory and i_fov in [4]:
        return True

    if "KIF1C_puro_1h_" in input_directory and i_fov in [0]:
        return True

    if "ASPM_3346_220419_puro_Adham" in input_directory and i_fov in [8]:
        return True

    if "MYH3_2279_230419_puro_Adham" in input_directory and i_fov in [0]:
        return True

    else:
        return False


if __name__ == "__main__":
    print()

    parser = argparse.ArgumentParser()
    parser.add_argument("--base_directory",
                        help="Path of the data directory.",
                        type=str,
                        default="/Users/arthur/data/2019_racha")
    parser.add_argument("--log_directory",
                        help="Path of the log directory.",
                        type=str,
                        default="/Users/arthur/output/2019_racha/log")
    args = parser.parse_args()
    base_directory = args.base_directory
    log_directory = args.log_directory
    if not os.path.isdir(base_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(base_directory))
    if not os.path.isdir(log_directory):
        raise ValueError("Directory does not exist: {0}"
                         .format(log_directory))

    now = datetime.datetime.now()
    date = now.strftime("%Y%m%d%H%M%S")
    log_file = os.path.join(log_directory, "log_loader" + "_" + date)
    sys.stdout = Logger(log_file)

    print("Running {0} file...".format(os.path.basename(__file__)), "\n")
    start_time = time.time()

    print("Data directory: {0}".format(base_directory))
    print("Log directory: {0}".format(log_directory), "\n")

    print("Checking if every files exist and recipes...")
    check_files_exist(base_directory)
    print("Done !", "\n")

    print("Files are saved with the pattern "
          "'gene_author_puromycin_paper_drug_batch_fov' \n")

    print("Loading images...")
    nb_images = 0
    directories = get_directories()
    for data_directory in directories:
        print("\t directory {0}".format(data_directory))
        if data_directory is None:
            continue
        experience = get_metadata_directory(data_directory)
        filename_base = generate_filename_base(experience)
        generator = images_generator(base_directory, data_directory)
        for i, image in enumerate(generator):
            filename = filename_base + "_" + str(i)
            print("\t\t", image.shape, image.dtype, filename)
            nb_images += 1
    print("Done ({0} images)!".format(nb_images), "\n")

    end_time = time.time()
    duration = int(round((end_time - start_time) / 60))
    print("Duration: {0} minutes.".format(duration))
