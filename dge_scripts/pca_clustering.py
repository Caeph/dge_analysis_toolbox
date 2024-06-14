import pandas as pd
from sklearn.decomposition import PCA
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


def pca_cluster_on_deseq(all_deseq_results, parameters):
    deseq_dds, deseq_results, padj_deseq_results, filtered_deseq_results = all_deseq_results
    deseq2 = importr('DESeq2')
    rbase = importr('base')
    rld = deseq2.rlog(deseq_dds, blind=True)
    rld_mat = ro.r('assay')(rld)
    columns = rbase.colnames(rld)
    rld_mat_df = pd.DataFrame(rld_mat, index=deseq_results.index, columns=columns)

    # pca
    pca = PCA(n_components=2)
    pca_fitted = pd.DataFrame(pca.fit_transform(rld_mat_df.T), index=columns)  # because transposed

    # clustering
    correlation_matrix = rld_mat_df.corr()
    return pca_fitted, pca.explained_variance_ratio_, correlation_matrix


def pca_cluster_on_edger(all_edger_results, parameters):
    edgeR_results, padj_edger_results, filtered_edger_results = all_edger_results
