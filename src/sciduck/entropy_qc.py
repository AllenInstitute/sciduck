from anndata import AnnData
import scanpy as sc
from scipy.stats import entropy
import numpy as np

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
    for anno in annotations:
        print("Computing entropy on: " + anno)
        adata.obs[anno + "_entropy"] = -1 ## Initialize with a value outside range of entropy.
        for cell in range(0, adata.shape[0]):
            nearest_neighbors = nearest_ind[cell,:]
            adata.obs.loc[adata.obs.index[cell], anno + "_entropy"] = scipy.stats.entropy(adata.obs.loc[adata.obs.index[nearest_neighbors],anno].value_counts()/len(nearest_neighbors))
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
    for anno in annotations:
        print("Computing entropy on: " + anno)
        adata.obs[anno + "_entropy"] = adata.obs.groupby(group_by)[anno].value_counts(normalize=True).groupby(group_by).apply(lambda x: entropy(x))
