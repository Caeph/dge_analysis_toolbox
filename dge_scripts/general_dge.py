import pandas as pd
import numpy as np
from scipy import stats
import rpy2
import rpy2.robjects as robjects

robjects.r('library(BiocManager)')
# TODO make sure every pkg is installed


def __subset_matrix(parameters, treatment_group, control_group):
    samples_in_treatment = parameters.sample_file[parameters.sample_file["groupID"] == treatment_group]["sampleID"].values
    samples_in_control = parameters.sample_file[parameters.sample_file["groupID"] == control_group]["sampleID"].values
    relevant_count_matrix_subset = parameters.count_matrix[["gene_ID", *samples_in_treatment, *samples_in_control]]
    return relevant_count_matrix_subset

def perform_dge(parameters):
    to_compare = parameters.contrasts
    for i, row in to_compare.iterrows():
        treatment_group = row["treatment"]
        control_group = row['control']

        matrix_subset = __subset_matrix(parameters, treatment_group, control_group)

        print(f"Running analysis with treatment group {treatment_group} and control group {control_group}")


    #deseq_results = __perform_analysis_with_deseq(parameters)
    #edger_results = __perform_analysis_with_deseq(parameters)


def __perform_analysis_with_deseq(parameters):
    robjects.r('library(DESeq2)')
    results = ...
    return results


def __perform_analysis_with_edgeR(parameters):
    robjects.r('library(edgeR)')
    results = ...
    return results
