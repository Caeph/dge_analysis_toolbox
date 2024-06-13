import argparse
import datetime
import os
from pathlib import Path

import pandas as pd

from dge_scripts.general_dge import perform_dge
from dge_scripts.pca_clustering import pca_cluster_on_deseq, pca_cluster_on_edger

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
                    A path to a tab delimited file indicating which groups should be comparised in a pairwise manner. 
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

        # if new parameter is added, it should be added to these lists too
        # they serve for reporting, no practical measure though
        self.__optional_parameters = [self.organism_info,
                                      self.padj_alpha,
                                      self.fold_change_threshold,
                                      self.gene_annotation_resource]
        self.__optional_parameter_labels = ["organism info",
                                            "FDR maximum threshold",
                                            "fold change minimum threshold",
                                            "gene annotation resource"]

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


def main(args):
    parameters = DGE_parameters(args)
    parameters.report()

    print("Running general DGE analysis...")
    all_deseq_results, all_edger_results = perform_dge(parameters)
    deseq_pca_fitted, deseq_correlation_matrix = pca_cluster_on_deseq(all_deseq_results, parameters)

    # prep pca and hierarchical clustering


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
