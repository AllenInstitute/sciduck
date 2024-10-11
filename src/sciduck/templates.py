import numpy as np
import pandas as pd
from anndata import AnnData
import warnings

##
def allen_institute_basic(adata: AnnData):
    """
    Adds Allen Institute basic QC standards to the AnnData object. Does not perform the filtering.

    Args:
        adata: The AnnData object where constraints will be added.
    """
    ## Filter cells/nuclei on UMI and gene count thresholds, showing the default values.
    add_range_constraint(adata, column="counts", gt=2000, lt=100000)
    add_range_constraint(adata, column="genes", gt=1000, lt=13000)

    ## Filter cells/nuclei on mitochondrial gene expression, showing the default values.
    add_range_constraint(adata, "doublet_score", lt = 0.3)
    add_range_constraint(adata, "pct_counts_mt", lt = 3.0)
    add_range_constraint(adata, "GEX_Reads_mapped_confidently_to_genome", gt = 0.0)
    add_range_constraint(adata, "GEX_Reads_mapped_to_genome", gt = 0.0)
    add_range_constraint(adata, "GEX_Reads_with_TSO", lt = 1.0)

    ## Neuron / Non-Neuron QC constraints
    add_range_constraint(adata, "genes", gt = 2000, subset = "Class", subset_values = ['Excitatory', 'Inhibitory'])
    add_range_constraint(adata, "genes", gt = 1000, subset = "Class", subset_values = ['Astrocytes', 'Oligodendrocytes', 'Microglia', 'Endothelial', 'Pericytes'])

    return adata
