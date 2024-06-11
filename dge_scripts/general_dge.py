import pandas as pd
import numpy as np
from scipy import stats
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

pandas2ri.activate()

ro.r('library(BiocManager)')
# Load required R packages
deseq2 = importr('DESeq2')
edgeR = importr('edgeR')


# TODO make sure every pkg is installed


def __subset_matrix(parameters, treatment_group, control_group):
    samples_in_treatment = parameters.sample_file[parameters.sample_file["groupID"] == treatment_group][
        "sampleID"].values
    samples_in_control = parameters.sample_file[parameters.sample_file["groupID"] == control_group]["sampleID"].values
    relevant_count_matrix_subset = parameters.count_matrix[[*samples_in_treatment, *samples_in_control]]
    relevant_count_matrix_subset.index = parameters.count_matrix['gene_ID']
    return relevant_count_matrix_subset, samples_in_treatment, samples_in_control


def _filter_matrix_on_cpm_threshold(rnaseq_matrix, cpm_genecount_min_threshold=2):
    r_rnaseqMatrix = pandas2ri.py2rpy(rnaseq_matrix)
    normed_r_rnaseqMatrix = edgeR.cpm(r_rnaseqMatrix)
    normed_r_rnaseqMatrix = pd.DataFrame(normed_r_rnaseqMatrix,
                                         columns=rnaseq_matrix.columns,
                                         index=rnaseq_matrix.index
                                         )

    # filter the ORIGINAL matrix -- select only gene with sufficient normed read count
    normed_readcount_for_genes = (normed_r_rnaseqMatrix.values > 1).sum(axis=1)
    filtered_r_rnaseqMatrix = rnaseq_matrix.loc[normed_readcount_for_genes >= cpm_genecount_min_threshold]

    return filtered_r_rnaseqMatrix


def perform_dge(parameters):
    to_compare = parameters.contrasts
    for i, row in to_compare.iterrows():
        treatment_group = row["treatment"]
        control_group = row['control']
        print(f"Running analysis with treatment group {treatment_group} and control group {control_group}")
        matrix_subset, samples_in_treatment, samples_in_control = __subset_matrix(parameters,
                                                                                  treatment_group,
                                                                                  control_group)
        matrix_subset = matrix_subset.round(decimals=0)
        prepped_matrix_subset = _filter_matrix_on_cpm_threshold(
            matrix_subset)  # the bottom read count threshold could be adjusted here

        deseq_results = __perform_analysis_with_deseq(prepped_matrix_subset,
                                                      samples_in_treatment,
                                                      samples_in_control,
                                                      treatment_group,
                                                      control_group,
                                                      parameters)
        gene_annotations = __annotate_gene_names(parameters, deseq_results.index)


        break


def __annotate_gene_names(parameters, gene_col):
    ...
    # TODO
    return gene_col


def __perform_analysis_with_deseq(prepped_matrix_subset, samples_in_treatment, samples_in_control,
                                  treatment_name, control_name,
                                  parameters
                                  ):
    # tutorial says counts in matrix should NOT be normalized
    conditions = pd.DataFrame({"conditions": [*[treatment_name for _ in range(len(samples_in_treatment))],
                                              *[control_name for _ in range(len(samples_in_control))]],
                               },
                              index=[*samples_in_treatment, *samples_in_control]
                              )
    # just to be sure, reorder the columns in prepped matrix subset
    prepped_matrix_subset = prepped_matrix_subset[[*samples_in_treatment, *samples_in_control]]

    # do the analysis in R
    processing_func = ro.r('''
    function(rnaseqMatrix, conditions) {
      library(DESeq2)
      ddsData <- DESeqDataSetFromMatrix(countData = rnaseqMatrix,
                              colData = conditions,
                              design = ~ conditions)
      dds = DESeq(ddsData)
      return(dds)
    }
    ''')
    # TODO uklidit to tak, aby vysledek byl dataframe

    colData = pandas2ri.py2rpy(conditions)
    rnaseqMatrix = pandas2ri.py2rpy(prepped_matrix_subset)
    dds = processing_func(rnaseqMatrix, colData)
    results_dds = ro.r('''
    function(dds, treatment, control) {
        output = results(dds, c("conditions", treatment, control)) 
        return(as.data.frame(output))
    }
    ''')(dds, treatment_name, control_name)
    results_dds = pandas2ri.rpy2py(results_dds)

    counts = pd.DataFrame(ro.r('counts')(dds, normalized=True),
                          columns=[*samples_in_treatment, *samples_in_control],
                          index=prepped_matrix_subset.index)

    results_dds["baseMean_treatment"] = counts[samples_in_treatment].T.mean()
    results_dds['baseMean_control'] = counts[samples_in_control].T.mean()
    results_dds['padj'] = results_dds['padj'].fillna(1)

    return results_dds.sort_values(by="pvalue", ascending=True)


def __perform_analysis_with_edgeR(parameters, prepped_matrix_subset):
    results = ...
    return results
