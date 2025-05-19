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

def cluster_annotation_entropy(
    adata: AnnData,
    cluster_col: str,
    label_cols: list,
    joint_entropy_key: str = "joint_entropy",
    exclude_label: str = "unknown",
    threshold: float = 0.1,
    debug: bool = False
) -> AnnData:
    """
    Compute per-cluster entropy for each label column individually using `cluster_entropy`,
    and compute entropy based on disagreement between the two labels,
    excluding cells where either label is `exclude_label`.
    Also adds a boolean column `inconsistent_label` which is True if all entropies >= threshold.
    """
    assert len(label_cols) == 2, "Exactly two label columns required"

    # Step 1: compute per-label entropy
    adata = cluster_entropy(adata, group_by=cluster_col, annotation_columns=label_cols)

    # Step 2: compute joint (disagreement) entropy
    joint_entropies = {}
    for cluster, group in tqdm(adata.obs.groupby(cluster_col), desc="Computing mismatch entropy"):
        labels = group[label_cols].astype(str).apply(lambda x: x.str.strip().str.lower())

        valid_mask = (labels[label_cols[0]] != exclude_label) & (labels[label_cols[1]] != exclude_label)
        valid_labels = labels[valid_mask]

        if valid_labels.empty:
            joint_entropies[cluster] = 0.0
            continue

        mismatches = valid_labels[valid_labels[label_cols[0]] != valid_labels[label_cols[1]]]
        disagreement_fraction = len(mismatches) / len(valid_labels)
        joint_entropies[cluster] = disagreement_fraction

        if debug:
            print(f"\n--- Cluster {cluster} mismatch debug ---")
            print(f"Valid cells: {len(valid_labels)}")
            print(f"Mismatches: {len(mismatches)}")
            print(f"Disagreement fraction: {disagreement_fraction:.4f}")
            print("Top mismatch pairs:")
            print(
                mismatches.value_counts()
                .rename("count")
                .reset_index()
                .head()
            )

    adata.obs[joint_entropy_key] = adata.obs[cluster_col].map(joint_entropies)

    # Step 3: define inconsistent_label flag
    col1_entropy = adata.obs[f"{label_cols[0]}_entropy"]
    col2_entropy = adata.obs[f"{label_cols[1]}_entropy"]
    joint_entropy = adata.obs[joint_entropy_key]

    inconsistent_flag = (
        (col1_entropy >= threshold) &
        (col2_entropy >= threshold) &
        (joint_entropy >= threshold)
    )

    adata.obs["inconsistent_label"] = inconsistent_flag

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
    entropy_cols = [f"{label}_entropy" for label in label_cols]

    if len(label_cols) >= 2:
        entropy_cols.append("joint_entropy")

    return (
        adata.obs[[cluster_col] + entropy_cols]
        .drop_duplicates(subset=cluster_col)
        .set_index(cluster_col)
        .sort_index()
    
