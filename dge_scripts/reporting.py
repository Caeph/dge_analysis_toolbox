import os

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def plot_clustering(correlation_matrix, tag, out_dir, parameters, linkage='average'):
    print(f"Plotting hierarchical clustering with {linkage}: {tag}...")
    sns.clustermap(correlation_matrix, method=linkage, cmap='viridis', linewidths=.5)
    # plt.title(f'{tag}\n{linkage} linkage hier. clustering')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_hierarchical_clustering.png")
    plt.savefig(filepath)
    plt.close()
    return filepath


def plot_pca(fitted_pca, pca_explained_ratio, tag, out_dir, parameters):
    print(f"Plotting PCA: {tag}...")
    pc_labels = [f"PC{i + 1}: explains {np.round(val * 100, decimals=1)}% variance" for i, val in enumerate(
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

    # annotation
    for line in range(0, fitted_pca.shape[0]):
        p1.text(fitted_pca[0].iloc[line] + 0.1, fitted_pca[1].iloc[line],
                fitted_pca["sampleID"][line], horizontalalignment='left',
                size='medium', color='black', weight='semibold')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_pca.png")

    plt.title(f'{tag}\nPCA (explaines {np.round(np.sum(pca_explained_ratio) * 100, decimals=1)}% variance)')
    plt.xlabel(pc1_label)
    plt.ylabel(pc2_label)

    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()
    return filepath


def plot_volcano(results_df, tag, out_dir, parameters, advanced=False):
    # log2(FC) on x-axis
    # -1*log10(padj) on y-axis
    # color in red those with -1*log10(padj) > thr
    if advanced:
        desc = "advanced"
    else:
        desc = "basic"
    print(f"Plotting {desc} volcano plot: {tag}")

    results_df['log10_padj*(-1)'] = - np.log10(results_df['padj'])
    # define hues
    if advanced:
        results_df["gene_color"] = "low_padj_low_fcc"
        results_df.loc[(results_df['padj'] < parameters.padj_alpha), "gene_color"] = "low_fc"
        # results_df['log2FoldChange']
        negative_mask = (results_df['log2FoldChange'] < -parameters.fold_change_threshold) & (
                    results_df['padj'] < parameters.padj_alpha)
        results_df.loc[negative_mask, "gene_color"] = "negative_fc"
        positive_mask = (results_df['log2FoldChange'] > parameters.fold_change_threshold) & (
                    results_df['padj'] < parameters.padj_alpha)
        results_df.loc[positive_mask, "gene_color"] = "positive_fc"
        palette = {
            "low_padj_low_fcc": 'grey',
            "low_fc": 'black',
            "positive_fc": "green",
            "negative_fc": "red"
        }
        plottitle = "advanced"
    else:
        results_df["gene_color"] = (results_df['padj'] < parameters.padj_alpha)
        palette = {True: "red", False: "black"}
        plottitle = "basic"

    p1 = sns.scatterplot(data=results_df,
                         x='log2FoldChange',
                         y='log10_padj*(-1)',
                         hue='gene_color',
                         palette=palette,
                         size=1,
                         legend=False,
                         linewidth=0
                         )

    if advanced:
        plt.axvline(-parameters.fold_change_threshold, linestyle='--')
        plt.axvline(parameters.fold_change_threshold, linestyle='--')
        plt.axhline(parameters.padj_alpha, linestyle='--')

    # annotation of outliers
    to_annotate = results_df.head(parameters.annotate_extremes_no)  # this will work -- the dataframe is already ordered
    x_col, y_col = 'log2FoldChange', 'log10_padj*(-1)'
    for line in range(0, to_annotate.shape[0]):
        p1.text(to_annotate[x_col].iloc[line] + 0.1, to_annotate[y_col].iloc[line],
                to_annotate.index[line], horizontalalignment='left',
                size='small', color='black')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_{plottitle}_volcano_plot.png")
    plt.title(f'{tag}\n{plottitle} volcano plot')
    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()
    return filepath


def plot_MA(results_df, tag, out_dir, counts_col_name, parameters):
    # log counts on x-axis
    # log FC on y-axis
    # color: FDR < thr
    print(f"Plotting MA plot: {tag}...")
    results_df["gene_color"] = (results_df['padj'] < parameters.padj_alpha)
    results_df["log2Counts"] = np.log2(results_df[counts_col_name])
    palette = {True: "red", False: "black"}
    p1 = sns.scatterplot(data=results_df,
                    x='log2Counts',
                    y='log2FoldChange',
                    hue='gene_color',
                    palette=palette,
                    size=1,
                    legend=False,
                    linewidth=0
                    )
    plt.axvline(0, linestyle='--')
    plt.axhline(0, linestyle='--')

    # annotation of outliers
    to_annotate = results_df.head(parameters.annotate_extremes_no)  # this will work -- the dataframe is already ordered
    x_col, y_col = 'log2FoldChange', 'log10_padj*(-1)'
    for line in range(0, to_annotate.shape[0]):
        p1.text(to_annotate[x_col].iloc[line] + 0.1, to_annotate[y_col].iloc[line],
                to_annotate.index[line], horizontalalignment='left',
                size='small', color='black')

    filetag = tag.replace(' ', '_')  # just to be safe
    filepath = os.path.join(out_dir, f"{filetag}_MA_plot.png")
    plt.title(f'{tag}\nMA plot')
    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()
    return filepath
