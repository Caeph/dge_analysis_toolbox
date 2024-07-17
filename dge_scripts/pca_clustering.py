import pandas as pd
from sklearn.decomposition import PCA
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


def pca_cluster_on_deseq(all_deseq_results, parameters, pca_ntop=500):
    deseq_dds, deseq_results = all_deseq_results
    deseq2 = importr('DESeq2')
    rbase = importr('base')
    rld = deseq2.rlog(deseq_dds, blind=True)
    rld_mat = ro.r('assay')(rld)
    columns = rbase.colnames(rld)
    rld_mat_df = pd.DataFrame(rld_mat, index=deseq_results.index, columns=columns)

    # clustering
    correlation_matrix = rld_mat_df.corr()

    # pca
    pca_results = ro.r("""
        function(object, intgroup="conditions", ntop=500, pcsToUse=1:2) {
  rv <- rowVars(assay(object), useNames=TRUE)

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  
  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(V1=pca$x[,pcsToUse[1]],
                  V2=pca$x[,pcsToUse[2]],
                  group=group, intgroup.df, name=colnames(object))
  colnames(d)[1:2] <- pcs
  return(c(percentVar[pcsToUse], d))
}
        """
    )(rld)
    pc1_explained_ratio, pc2_explained_ratio = pca_results[0][0], pca_results[1][0]
    levels = {i:level for i, level in enumerate(list(ro.r['levels'](pca_results[5])))}
    groups = pandas2ri.rpy2py(pca_results[5]) - 1
    groups = [levels[val] for val in groups]

    pca_fitted = pd.DataFrame(
        {
            "PCA component 1": pca_results[2],
            "PCA component 2": pca_results[3],
            # "groupID": groups,
            "sampleID": pca_results[6]
        }
    )

    return pca_fitted, [pc1_explained_ratio, pc2_explained_ratio], correlation_matrix


def pca_cluster_on_edger(all_edger_results, parameters):
    edgeR_results, padj_edger_results, filtered_edger_results = all_edger_results
    raise NotImplementedError("PCA on edger is not implemented")
