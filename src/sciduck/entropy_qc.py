from anndata import AnnData
import scanpy as sc
from scipy.spatial import KDTree
from scipy.stats import entropy
import numpy as np
import pandas as pd

def cell_entropy(adata: AnnData,
                    annotation_columns: list,
                    nearest_neighbors: int = 15, 
                    dim: str = "X_scVI"
                ) -> AnnData:
    """
    Compute entropy (mixing) of annotations within a cells local neighborhood. 

    Args:
        adata: AnnData object with `annotations`
        annotation_columns: Cell level annotations.
        nearest_neighbors: Number of nearest neighbors.
        dim: Dimensionality reduction
        
    Returns:
        Returns the updated AnnData object.
    """
    ## Build nearest neighbor tree for fast lookup
    print("Building nearest neighbor tree.")
    nnTree = KDTree(adata.obsm[dim])
    nearest_dist, nearest_ind = nnTree.query(adata.obsm[dim], k=nearest_neighbors)
    ##
    for anno in annotation_columns:
        print("Computing entropy on: " + anno)
        adata.obs[anno + "_entropy"] = -1 ## Initialize with a value outside range of entropy.
        for cell in range(0, adata.shape[0]):
            nearest_neighbors = nearest_ind[cell,:]
            adata.obs.loc[adata.obs.index[cell], anno + "_entropy"] = entropy(adata.obs.loc[adata.obs.index[nearest_neighbors],anno].value_counts()/len(nearest_neighbors))
    return adata

def cluster_entropy(adata: AnnData,
                    group_by: str,
                    annotation_columns: list
                ) -> AnnData:
    """
    Compute entropy (mixing) of annotations within a cell set (group_by). 

    Args:
        adata: AnnData object with `annotations`
        group_by: AnnData.obs column that groups cells.
        annotation_columns: Cell level annotations.
        
    Returns:
        Returns the updated AnnData object.
    """
    ##
    for anno in annotation_columns:
        print("Computing entropy on: " + anno)
        adata.obs[anno + "_entropy"] = adata.obs.groupby(group_by)[anno].value_counts(normalize=True).groupby(group_by).apply(lambda x: entropy(x))
    return adata

def cluster_annotation_entropy(adata: AnnData,
                               cluster_col: str,
                               label_cols: list
                              ) -> AnnData:
    """
    Compute per-cluster entropy for each label column individually,
    and the joint entropy across all label columns for each cluster.
    
    The computed values are added to `adata.obs` for all cells,
    using the following keys:
        - For each label: <label>_cluster_entropy
        - For joint entropy: joint_entropy

    Args:
        adata (AnnData): Annotated data matrix with `.obs` containing both
                         clustering labels and categorical annotation columns.
        cluster_col (str): The name of the column in `adata.obs` that defines clusters (e.g., 'leiden').
        label_cols (list of str): List of one or more columns in `adata.obs` for which to compute entropy.
                                  These should be categorical annotations (e.g., 'batch', 'cell_type').

    Returns:
        AnnData: The input `adata` object with new `.obs` columns added:
                 - `<label>_cluster_entropy`: entropy of the label distribution within each cluster
                 - `joint_entropy`: entropy of the joint label distribution within each cluster
    """
    print(f"Computing cluster-wise annotation entropy for labels: {label_cols}")
    cluster_groups = adata.obs.groupby(cluster_col)

    for cluster, group in cluster_groups:
        cluster_index = group.index  # Index of all cells in the current cluster
        stats = {}

        # Compute entropy for each label separately
        for label in label_cols:
            label_counts = group[label].value_counts(normalize=True)
            stats[f"{label}_cluster_entropy"] = entropy(label_counts)

        # Compute joint entropy across all label combinations
        if len(label_cols) >= 2:
            # Create a single string per row that represents the combination of labels
            joint_labels = group[label_cols].astype(str).agg("_".join, axis=1)
            joint_counts = joint_labels.value_counts(normalize=True)
            stats["joint_entropy"] = entropy(joint_counts)

        # Assign entropy values to each cell in the cluster
        for key, val in stats.items():
            adata.obs.loc[cluster_index, key] = val

    return adata

    from anndata import AnnData

def extract_entropy_df(adata: AnnData,
                       cluster_col: list,
                       label_cols: list,
    ) -> pd.DataFrame:
    """
    Extract cluster-level entropy values from .obs of an AnnData object.

    Args:
        adata: AnnData with entropy columns added by `cluster_annotation_entropy`.
        cluster_col: Column in .obs identifying clusters (e.g., 'leiden').
        label_cols: List of label columns used in the entropy calculation.

    Returns:
        pd.DataFrame with one row per cluster and columns for each entropy metric.
    """
    entropy_cols = [f"{label}_cluster_entropy" for label in label_cols]

    if len(label_cols) >= 2:
        entropy_cols.append("joint_entropy")

    return (
        adata.obs[[cluster_col] + entropy_cols]
        .drop_duplicates(subset=cluster_col)
        .set_index(cluster_col)
        .sort_index()
    
