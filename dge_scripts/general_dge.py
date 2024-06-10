import pandas as pd
import numpy as np
from scipy import stats
import rpy2
import rpy2.robjects as robjects

robjects.r('library(BiocManager)')
# TODO make sure everything is installed
robjects.r('library(BiocManager)')


def perform_dge(parameters):
    ...


def __perform_analysis_with_deseq(parameters):
    robjects.r('library(DESeq2)')


def __perform_analysis_with_edgeR(parameters):
    robjects.r('library(edgeR)')
