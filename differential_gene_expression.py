import argparse
import datetime
import os
from pathlib import Path

import pandas as pd

from dge_scripts.general_dge import perform_dge
from dge_scripts.pca_clustering import pca_cluster_on_deseq, pca_cluster_on_edger
import dge_scripts.reporting as reporting

parser = argparse.ArgumentParser()

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("--count_matrix", type=str,
                           help="""Path to tab delimited file with the following structure:  first columns contains GeneIDs, 
                    and each other column represents read counts (not normalized) for each sample.
                    """,
                           required=True
                           )
requiredNamed.add_argument("--sample_file", type=str,
                           help="""
                    A path to a tab delimited file with groupID in first column and sampleID in the second.
                    """,
                           required=True
                           )
requiredNamed.add_argument("--contrasts", type=str,
                           help="""
                    A path to a tab delimited file indicating which groups should be compared in a pairwise manner. 
                    Keep always the same structure: 1st column - treatment, 2nd column - control
                    """,
                           required=True
                           )

# optional parameters
parser.add_argument("--organism", default="mus", type=str,
                    help="""
                    Organism to work on. Supported: mus (Mus musculus) | homo (Homo sapiens). 
                    Default: homo.
                    """
                    )
parser.add_argument("--padj_alpha", default=0.05, type=float,
                    help="""
                    FDR maximum threshold for results selection.
                    """)
parser.add_argument("--fold_change_threshold", default=1, type=float,
                    help="""
                    Minimal log fold change threshold
                    """)
parser.add_argument("--gene_annotation_resource", default=None, type=str,
                    help="""
                    Resource to use for gene annotation. If undefined, no annotation is performed.
                    Supported options: ENSEMBL
                    """)
parser.add_argument("--output_directory_path", default=None, type=str,
                    help="""
                    Path to store the directory. If the directory exists, the run will fail.
                    """
                    )

available_info_on_organisms = {
    "mus": dict(organism_name="Mus musculus",
                database="org.Mm.eg.db",
                organism_shortcut='mmu'),
    "homo": dict(organism_name="Homo sapiens",
                 database="org.Hs.eg.db",
                 organism_shortcut='hsa')
}


class DGE_parameters:
    def __init__(self, args):
        self.count_matrix_filename = args.count_matrix
        self.count_matrix = pd.read_csv(args.count_matrix, sep='\t')

        self.sample_filename = args.sample_file
        self.sample_file = pd.read_csv(args.sample_file, sep='\t', header=None)
        self.sample_file.columns = ["groupID", "sampleID"]

        self.contrasts_filename = args.contrasts
        self.contrasts = pd.read_csv(args.contrasts, sep='\t', header=None)
        self.contrasts.columns = ["treatment", "control"]

        # optional params
        self.organism_info = available_info_on_organisms[args.organism]
        self.padj_alpha = args.padj_alpha
        self.fold_change_threshold = args.fold_change_threshold
        self.gene_annotation_resource = args.gene_annotation_resource

        # no of extreme genes to show annotation for
        self.annotate_extremes_no = 20

        # if new parameter is added, it should be added to these lists too
        # they serve for reporting, no practical measure though
        self.__optional_parameters = [self.organism_info,
                                      self.padj_alpha,
                                      self.fold_change_threshold,
                                      self.gene_annotation_resource,
                                      self.annotate_extremes_no]
        self.__optional_parameter_labels = ["organism info",
                                            "FDR maximum threshold",
                                            "fold change minimum threshold",
                                            "gene annotation resource",
                                            "no of extreme genes to annotate"]

        # output prep
        if args.output_directory_path is None:
            now = datetime.datetime.now().strftime("%d-%m-%Y-%H:%M:%S")
            build = ["output-dge-analysis", now]
            for param, val in zip(["padj", "fc"],
                                  [self.padj_alpha, self.fold_change_threshold]):
                build.append(f"{param}={val}")
            args.output_directory_path = "_".join(build)

        self.output_dir = args.output_directory_path
        if Path(self.output_dir).exists():
            raise FileExistsError(f"The path {self.output_dir} already exists.")
        os.makedirs(self.output_dir, exist_ok=False)

        # recap params to the output
        with open(os.path.join(self.output_dir, "parameters.tsv"), mode='w') as writer:
            print("Parameter\tvalue\tinfo", file=writer)
            for filename, label in zip(
                    [self.count_matrix_filename, self.sample_filename, self.contrasts_filename],
                    ["count matrix", "sample file", "contrasts"]):
                print(f"{label}\t{filename}\trequired", file=writer)
            for value, label in zip(self.__optional_parameters, self.__optional_parameter_labels):
                stringified_values = str(value).replace("\n", " ")
                print(f"{label}\t{stringified_values}\t", file=writer)

    def report(self):
        # required -- filenames and other stuff
        print("Input files overview:")
        for filename, dataframe, label in zip(
                [self.count_matrix_filename, self.sample_filename, self.contrasts_filename],
                [self.count_matrix, self.sample_file, self.contrasts],
                ["count matrix", "sample file", "contrasts"]
        ):
            dfsize = len(dataframe)
            columns = ", ".join(dataframe.columns)
            print(f"{label}: taken from {filename}, it has {dfsize} data rows")
            print(f"detected columns: {columns}")
            print("--")

        print()
        print("Other parameters:")
        for value, label in zip(self.__optional_parameters, self.__optional_parameter_labels):
            stringified_values = str(value).replace("\n", " ")
            print(f"{label}: {stringified_values}")
        print("--")
        print()


def venn_diagrams(all_edger_results, all_deseq_results, out_dir, parameters):
    venn_tags = []  # tag, path

    edger_df, deseq_df = all_edger_results[-1], all_deseq_results[-1]

    # fdr filter
    edger_df_sub = edger_df[edger_df["padj"] < parameters.padj_alpha]
    tag_edger, tag_deseq = f"edgeR padj<{parameters.padj_alpha}", f"DESeq2 padj<{parameters.padj_alpha}"
    deseq_df_sub = deseq_df[deseq_df['padj'] < parameters.padj_alpha]
    filepath = reporting.plot_venn_diagram(edger_df_sub,
                                           deseq_df_sub,
                                           tag_edger,
                                           tag_deseq,
                                           out_dir,
                                           parameters)
    venn_tags.append(
        [f"Venn {tag_edger} vs. {tag_deseq}", filepath]
    )

    # fdr and fold change filter
    edger_df_sub = edger_df_sub[edger_df_sub["log2FoldChange"] > parameters.fold_change_threshold]
    deseq_df_sub = deseq_df_sub[deseq_df_sub["log2FoldChange"] > parameters.fold_change_threshold]
    tag_edger = f"edgeR padj<{parameters.padj_alpha} FC>{parameters.fold_change_threshold}"
    tag_deseq = f"DESeq2 padj<{parameters.padj_alpha} FC>{parameters.fold_change_threshold}"
    filepath = reporting.plot_venn_diagram(edger_df_sub,
                                           deseq_df_sub,
                                           tag_edger,
                                           tag_deseq,
                                           out_dir,
                                           parameters)
    venn_tags.append(
        [f"Venn {tag_edger} vs. {tag_deseq}", filepath]
    )

    return venn_tags


def main(args):
    parameters = DGE_parameters(args)
    parameters.report()

    print("Running general DGE analysis...")
    full_results = perform_dge(parameters)
    figures_dir = os.path.join(parameters.output_dir, "figures")
    os.makedirs(figures_dir)

    for treatment_id, control_id, all_deseq_results, all_edger_results in full_results:
        deseq_pca_fitted, deseq_pca_info, deseq_correlation_matrix = pca_cluster_on_deseq(all_deseq_results,
                                                                                          parameters)

        # plot the plots
        plot_files = []
        current_figures_dir = os.path.join(figures_dir, f"{treatment_id}__vs__{control_id}".replace(" ", "_"))
        os.makedirs(current_figures_dir)

        deseq_tag = f"DESeq2_DGE_{treatment_id}_vs_{control_id}"
        edger_tag = f"edgeR_DGE_{treatment_id}_vs_{control_id}"

        # update counts:
        deseq_df = all_deseq_results[-1]
        deseq_df["baseMean+1"] = deseq_df["baseMean"] + 1

        # volcano and MA
        for df, tag, counts_col in zip([all_deseq_results[-1], all_edger_results[-1]],
                                       [deseq_tag, edger_tag],
                                       ["baseMean+1", "logCPM"]
                                       ):
            # basic volcano
            path = reporting.plot_volcano(df,
                                          tag,
                                          current_figures_dir,
                                          parameters)
            plot_files.append((f"{tag} -- basic volcano", path))

            # advanced volcano
            path = reporting.plot_volcano(df,
                                          tag,
                                          current_figures_dir,
                                          parameters,
                                          advanced=True)
            plot_files.append((f"{tag} -- advanced volcano", path))

            # MA plot
            path = reporting.plot_MA(df,
                                     tag,
                                     current_figures_dir,
                                     counts_col,
                                     parameters)
            plot_files.append((f"{tag} -- MA plot", path))

        # clustering
        path = reporting.plot_clustering(deseq_correlation_matrix,
                                         deseq_tag,
                                         current_figures_dir,
                                         parameters)
        plot_files.append((f"{deseq_tag} -- hierarchical clustering", path))

        # PCA
        path = reporting.plot_pca(deseq_pca_fitted,
                                  deseq_pca_info,
                                  deseq_tag,
                                  current_figures_dir,
                                  parameters)
        plot_files.append((f"{deseq_tag} -- pca", path))

        # venn diagrams with comparison
        venn_files = venn_diagrams(all_edger_results,
                                   all_deseq_results,
                                   current_figures_dir,
                                   parameters)
        plot_files.extend(venn_files)

    # \plot the plots


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
