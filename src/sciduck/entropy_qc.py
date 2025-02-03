from anndata import AnnData
import scanpy as sc
import numpy as np

def cluster_entropy_qc_metric(adata: AnnData,
                                cluster_column: str,
                                annotation_column: str
                            ) -> list:
    """
    Compute entropy (mixing) of an annotation within a pre-defined cluster.

    Args:
        adata: Anndata object.
        cluster_column: Column name in `obs` containing cluster labels.
        annotation_columns: Column name in `obs` containing annotations.

    Returns:
        A list of cluster entropy quality control metrics.
    """
    ##
    print("Computing cluster entropy on: " + annotation_column)
    ## Initialize with a value outside range of entropy.
    qc_entropy = [-1] * adata.shape[0]
    ## Compute entropy on annotation column for each cluster
    for cluster in np.unique(adata.obs["cluster_column"]):
        qc_entropy[adata.obs["cluster_column"] == cluster] = \
            scipy.stats.entropy(adata.obs.loc[adata.obs["cluster_column"] == cluster, annotation_column].value_counts() / sum(adata.obs["cluster_column"] == cluster))
    ##
    return qc_entropy

def cell_entropy_qc_metric(adata: AnnData,
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
