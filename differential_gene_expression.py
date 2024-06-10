import argparse
import pandas as pd

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
parser.add_argument("--organism", default="homo", type=str,
                    help="""
                    Organism to work on. Supported: mus (Mus musculus) | homo (Homo sapiens). 
                    Default: homo.
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
        self.count_matrix = pd.read_csv(args.count_matrix)

        self.sample_filename = args.sample_file
        self.sample_file = pd.read_csv(args.sample_file, sep='\t')

        self.contrasts_filename = args.contrasts
        self.contrasts = pd.read_csv(args.contrasts, sep='\t')

        self.organism_info = available_info_on_organisms[args.organism]

        # if new parameter is added, it should be added to these lists too
        self.optional_parameters = [self.organism_info]
        self.parameter_labels = ["organism_info"]

    def report(self):
        # required -- filenames and other stuff
        print("Required files:")
        for filename, dataframe, label in zip(
            [self.count_matrix, self.sample_filename, self.contrasts_filename]
        ):
            dfsize = len(dataframe)
            columns = ", ".join(dataframe.columns.values)
            print(f"{label}: taken from {filename}, it has {dfsize} data rows")
            print(f"detected columns: {columns}")
            print("--")

        print()
        print("Other parameters")
        for value, label in zip():
            stringified_values = str(value).replace("\n", " ")
            print(f"{label}: {stringified_values}")



def report_setup(args):
    parameters = DGE_parameters(args)
    parameters.report()


def main(args):
    report_setup(args)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    args = parser.parse_args()
    main(args)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
