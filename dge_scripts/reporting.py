import os

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


def plot_clustering(correlation_matrix, tag, out_dir, parameters, linkage='average'):
    print(f"Plotting hierarchical clustering with {linkage}...")
    sns.clustermap(correlation_matrix, method=linkage, cmap='viridis', linewidths=.5)
    # plt.title(f'{tag}\n{linkage} linkage hier. clustering')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_hierarchical_clustering.png")
    plt.savefig(filepath)
    plt.close()
    return filepath


def plot_pca(fitted_pca, pca_explained_ratio, tag, out_dir, parameters):
    print(f"Plotting PCA...")
    pc_labels = [f"PC{i+1}: explains {np.round(val*100, decimals=1)}% variance" for i, val in enumerate(
        pca_explained_ratio)]

    pc1_label, pc2_label = pc_labels[0], pc_labels[1]

    fitted_pca = pd.merge(fitted_pca, parameters.sample_file, left_index=True, right_on='sampleID', how='inner')
    p1 = sns.scatterplot(x=0,  # Horizontal axis
                         y=1,  # Vertical axis
                         data=fitted_pca,  # Data source
                         size=10,
                         hue="groupID",
                         legend=False
                         )

    for line in range(0, fitted_pca.shape[0]):
        p1.text(fitted_pca[0][line] + 0.1, fitted_pca[1][line],
                fitted_pca["sampleID"][line], horizontalalignment='left',
                size='medium', color='black', weight='semibold')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_pca.png")

    plt.title(f'{tag}\nPCA (explaines {np.round(np.sum(pca_explained_ratio)*100, decimals=1)}% variance)')
    plt.xlabel(pc1_label)
    plt.ylabel(pc2_label)

    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()
    return filepath