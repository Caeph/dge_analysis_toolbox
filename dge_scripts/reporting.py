import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns

"""
plot_MA = function(featureNames, logCounts, logFoldChange, FDR, xlab="logCounts", ylab="logFC", title="MA plot", pch=20, top_gene_labels_show=20) {
  plot(logCounts, logFoldChange, col=ifelse(FDR<=0.05, "red", "black"), xlab=xlab, ylab=ylab, main=title, pch=pch);
  text(logCounts[1:top_gene_labels_show], logFoldChange[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
}

plot_Volcano = function(featureNames, logFoldChange, FDR, xlab="logFC", ylab="-1*log10(FDR)", title="Basic Volcano plot", pch=20, top_gene_labels_show=20) {
  plot(logFoldChange, -1*log10(FDR), col=ifelse(FDR<=0.05, "red", ifelse(abs(logFoldChange)>=1 & FDR<=0.05, "red", "black")), xlab=xlab, ylab=ylab, main=title, pch=pch);
  text(logFoldChange[1:top_gene_labels_show], (-1*log10(FDR))[1:top_gene_labels_show], labels=featureNames[1:top_gene_labels_show], cex= 0.7, pos=3)
}
"""


def plot_clustering(correlation_matrix, title, linkage='average'):
    sns.clustermap(correlation_matrix, method=linkage, cmap='vlag', linewidths=.5)
    plt.title(f'{title}\n{linkage} linkage hier. clustering')

    # TODO how to return figure
