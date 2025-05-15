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
    Compute entropy (mixing) of annotations within a cell set (group_by),
    and assign that entropy to all cells in the group.
    """
    for anno in annotation_columns:
        print("Computing entropy on:", anno)
        
        # First compute per-group entropy
        group_entropies = (
            adata.obs.groupby(group_by)[anno]
            .value_counts(normalize=True)
            .groupby(level=0)
            .apply(lambda x: entropy(x.values))
        )
        
        # Assign back to cells
        adata.obs[f"{anno}_entropy"] = adata.obs[group_by].map(group_entropies)

    return adata

def cluster_annotation_entropy(adata: AnnData,
                               cluster_col: str,
                               label_cols: list,
                               joint_entropy_key: str = "joint_entropy"
                              ) -> AnnData:
    """
    Compute per-cluster entropy for each label column individually (via `cluster_entropy`),
    and the joint entropy across all label columns for each cluster.
    Adds per-cell values to adata.obs.
    """
    print(f"Computing cluster-wise annotation entropy for labels: {label_cols}")
    
    # Step 1: Compute individual entropies using helper
    adata = cluster_entropy(adata, group_by=cluster_col, annotation_columns=label_cols)

    # Step 2: Compute joint entropy
    if len(label_cols) >= 2:
        joint_entropies = (
            adata.obs.groupby(cluster_col)[label_cols]
            .apply(lambda df: entropy(df.astype(str).agg('_'.join, axis=1).value_counts(normalize=True).values))
        )
        adata.obs[joint_entropy_key] = adata.obs[cluster_col].map(joint_entropies)

    return adata

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
    
