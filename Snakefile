"""
Get stats for raw sequencing data in each sample.
"""

import collections
import os
import pysam
import re

import pandas as pd
import numpy as np

import ccsstatlib.plot
import ccsstatlib.stats

import matplotlib.pyplot as plt

localrules: subread_stats_get_tables


###############
# Definitions #
###############

# Get config
FOFN_FILE_NAME = config.get('fofn', None)

SAMPLE_NAME = config.get('sample', None)

if FOFN_FILE_NAME is None:
    raise RuntimeError('Missing FOFN file in config (list of input PacBio BAM or BAX files). Example: snakemake ... --config fofn=/path/to/sample.fofn')

if SAMPLE_NAME is None:
    if FOFN_FILE_NAME.lower().endswith('.fofn'):
        SAMPLE_NAME = os.path.basename(FOFN_FILE_NAME).rsplit('.', 1)[0]
    else:
        raise RuntimeError(
            'Sample name must be set (e.g. snakemake ... --config fofn=/path/to/sample.fofn sample=NAME). If the '
            'sample name is not set, the FOFN file name must end in .fofn, and the sample name will be the FOFN file '
            'name without the .fofn extension.'
        )


#######################
### Read FOFN Table ###
#######################


def get_cell_dict(as_list=False):
    """
    Get a list of input files per cell. Returns a dictionary where keys are cell names and values are a list of
    input BAM files for that cell.

    :param as_list: If `True`, return a list of cells in the order they were found in the FOFN.
    """

    with open(FOFN_FILE_NAME, 'r') as in_file:

        cell_dict = collections.defaultdict(list)
        cell_list = list()

        line_count = 0

        for line in in_file:

            line_count += 1

            line = line.strip()

            if not line or line.startswith('#'):
                continue

            ext_match = re.match(r'.*(\.bam|\.fastq|\.fastq\.gz)$', line, re.IGNORECASE)

            if ext_match is None:
                raise RuntimeError(
                    'Unrecognized file type in FOFN on line {}: Must end with ".bam", ".fastq", or ".fastq.gz": {}'.format(
                        line_count,
                        line
                    )
                )

            cell_name = os.path.basename(line)[:-len(ext_match[1])]

            if not as_list:
                if cell_name in cell_dict:
                    raise RuntimeError('Duplicate cell on line {}: {} ("{}" conflicts with previous entry "{}")'.format(
                        line_count,
                        cell_name,
                        os.path.basename(line),
                        os.path.basename(cell_dict[cell_name])
                    ))

                cell_dict[cell_name] = line

            else:
                if cell_name in cell_list:
                    raise RuntimeError('Duplicate cell name on line {}: {}'.format(line_count, cell_name))

                cell_list.append(cell_name)


    return cell_dict if not as_list else cell_list


#########
# Rules #
#########

#
# Default target
#

# ccs_stats_get_tables
#
# Get all tables.
rule ccs_stats_get_tables:
    input:
        tab_sample=f'{SAMPLE_NAME}/sample_summary.tsv.gz',
        tab_cell=f'{SAMPLE_NAME}/cell_summary.tsv.gz',
        png_cdf=f'{SAMPLE_NAME}/plot/cdf_cell.png',
        png_qv=f'{SAMPLE_NAME}/plot/hist-line_cell_qv.png',
        png_passes=f'{SAMPLE_NAME}/plot/hist-line_cell_passes.png',
        png_den=f'{SAMPLE_NAME}/plot/density_cell_len.png'


include: 'rules/ccs_stats.snakefile'
include: 'rules/plot.snakefile'
