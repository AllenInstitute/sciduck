from anndata import AnnData
import scanpy as sc
import numpy as np
import warnings

def filter_on_counts_genes(adata: AnnData, 
                            min_counts: int = 2000, 
                            max_counts: int = 100000, 
                            min_genes: int = 1000, 
                            max_genes: int = 13000, 
                            inplace: bool = False
                        ) -> AnnData | None:
    """
    Filter samples based on number of counts (UMIs) and genes detected.

    Args:
        adata: Anndata object.
        min_counts: Minimum counts detected.
        max_counts: Maximum counts detected.
        min_genes: Minimum genes detected.
        max_genes: Maximum genes detected.
        inplace: Update the adata object in place or return unmodified object with keeper_cells flagged.

    Returns:
        Returns either AnnData | None
    """
    if isinstance(adata, AnnData):
        ## Filter cells based on counts and genes detected
        if 'keeper_cells' not in adata.obs.columns:
            adata.obs["keeper_cells"] = [True] * adata.shape[0]
        
        ## Filter cells based on counts detected
        adata.obs["keeper_cells"] &= sc.pp.filter_cells(adata, min_counts = min_counts, inplace=False)[0]
        adata.obs["keeper_cells"] &= sc.pp.filter_cells(adata, max_counts = max_counts, inplace=False)[0]
        
        ## Filter cells based on genes detected
        adata.obs["keeper_cells"] &= sc.pp.filter_cells(adata, min_genes = min_genes, inplace=False)[0]
        adata.obs["keeper_cells"] &= sc.pp.filter_cells(adata, max_genes = max_genes, inplace=False)[0]
        
        ## 
        if inplace:
            adata._inplace_subset_obs(adata.obs["keeper_cells"])
            return None
        else:
            return adata

def filter_on_precomputed_metrics(adata: AnnData, 
                                    doublet_score: float = 0.3, 
                                    pct_counts_mt: float = 3.0, 
                                    GEX_Reads_mapped_confidently_to_genome: float = 0.0, 
                                    GEX_Reads_mapped_to_genome: float = 0.0, 
                                    GEX_Reads_with_TSO: float = 1.0, 
                                    inplace: bool = False,
                                ) -> AnnData | None:
    """
    Filter samples based on precomputed quality control metrics.

    Args:
        adata: Anndata object.
        doublet_score: Maximum doublet score.
        pct_counts_mt: Maximum percentage of counts in mitochondrial genes.
        GEX_Reads_mapped_confidently_to_genome: Minimum percentage of confidently mapped reads. There is no pre-defined good practice for threshold, user must specify.
        GEX_Reads_mapped_to_genome: Minimum percentage of reads mapped to genome. There is no pre-defined good practice for threshold, user must specify.
        GEX_Reads_with_TSO: Maximum percentage of reads with TSO, per cell. There is no pre-defined good practice for threshold, user must specify.
        inplace: Update the adata object in place or return unmodified object with keeper_cells flagged.

    Returns:
        Returns either AnnData | None
    """
    if isinstance(adata, AnnData):
        ## Filter cells based on counts and genes detected
        if 'keeper_cells' not in adata.obs.columns:
            adata.obs["keeper_cells"] = [True] * adata.shape[0]

        ## Gather function parameters
        filter_params = locals(); del filter_params['adata']; del filter_params['inplace']
        
        ## Check that columns are present in the adata object and skip if not present.
        for qc_metric in list(filter_params.keys()):
            if qc_metric not in adata.obs.columns:
                warnings.warn(f"A {qc_metric} column is missing in the adata object. Skipping this metric.", UserWarning)
                del filter_params[qc_metric]

        ## Filter cells based on qc metrics
        if "doublet_score" in filter_params.keys():
            adata.obs["keeper_cells"] &= adata.obs["doublet_score"] < doublet_score
        if "pct_counts_mt" in filter_params.keys():
            adata.obs["keeper_cells"] &= adata.obs["pct_counts_mt"] < pct_counts_mt
        if "GEX_Reads_mapped_confidently_to_genome" in filter_params.keys():
            adata.obs["keeper_cells"] &= adata.obs["GEX_Reads_mapped_confidently_to_genome"] > GEX_Reads_mapped_confidently_to_genome
        if "GEX_Reads_mapped_to_genome" in filter_params.keys():
            adata.obs["keeper_cells"] &= adata.obs["GEX_Reads_mapped_to_genome"] > GEX_Reads_mapped_to_genome
        if "GEX_Reads_with_TSO" in filter_params.keys():
            adata.obs["keeper_cells"] &= adata.obs["GEX_Reads_with_TSO"] < GEX_Reads_with_TSO

        ##
        if inplace:
            adata._inplace_subset_obs(adata.obs["keeper_cells"])
            return None
        else:
            return adata

def filter_utilizing_coarse_labels(adata: AnnData, 
                                    coarse_label_column: str,
                                    coarse_label_map: dict = {"Neurons": ["Excitatory", "Inhibitory"], 
                                                                "Non-Neurons": ["Astrocytes", "Oligodendrocytes", "Microglia", "Endothelial", "Pericytes"]},
                                    coarse_label_gene_threshold: dict = {"Neurons": 2000, "Non-Neurons": 1000},
                                    inplace: bool = False
                                ) -> AnnData | None:
    """
    Filter samples based on coarse label specific thresholds for genes detected.

    Args:
        adata: Anndata object.
        coarse_label_column: Column name in adata.obs containing coarse labeling identifying neuron and non-neuronal cell types.
        coarse_label_map: Coarse cell type labels to define specific filtering on.
        coarse_label_gene_threshold: Minimum genes detected for each coarse label.
        inplace: Update the adata object in place or return unmodified object with keeper_cells flagged.

    Returns:
        Returns either AnnData | None
    """
    if isinstance(adata, AnnData):
        ## Filter cells based on counts and genes detected
        if 'keeper_cells' not in adata.obs.columns:
            adata.obs["keeper_cells"] = [True] * adata.shape[0]
        ## Filter specific quality thresholds for each coarse label
        for coarse_label, coarse_labels in coarse_label_map.items():
            adata.obs["keeper_cells"] &= ~adata.obs[coarse_label_column].isin(coarse_labels) | sc.pp.filter_cells(adata, min_genes = coarse_label_threshold[coarse_label], inplace=False)[0]
        ##
        if inplace:
            adata._inplace_subset_obs(adata.obs["keeper_cells"])
            return None
        else:
            return adata

# def filter_on_cluster_metrics(adata: AnnData,
#                                 cluster_column: str,
#                                 annotation_columns: list,
#                                 annotation_thresholds: dict, ## HOW TO ENCODE DIRECTION, E,G. > 0.3 doublet score or < 0.9 TSO.
#                             ) -> AnnData | None:
#     """
#     Filter samples based on entropy of quality control metrics or metadata for each cluster.
#
#     Args:
#         adata: type AnnData: Anndata object.
#         cluster_column: Column name in adata.obs containing cluster labels.
#         annotation_columns: Column name in adata.obs containing entropy values.
#         annotation_thresholds: Minimum entropy values for each annotation being considered.
#        
#     Returns:
#         Returns either AnnData | None
#     """
#     if isinstance(adata, AnnData):
#         ## Filter cells based on counts and genes detected
#         if 'keeper_cells' not in adata.obs.columns:
#             adata.obs["keeper_cells"] = [True] * adata.shape[0]
#         ## Compute cluster medians for each annotation
#         for anno in annotation_columns:
#             try:
#                 adata.obs[anno + "_cluster_median"] = adata.obs.groupby("cluster")[anno].median().round(4)
#             except:
#                 warnings.warn(f"Unable to compute cluster median for {anno}. Skipping this metric.", UserWarning)
#                 continue
#         ## Apply filtering based on user provided thresholds
#         for anno, threshold in annotation_thresholds.items():
#             adata.obs["keeper_cells"] &= adata.obs[anno + "_entropy"] > threshold
#         ##
#         if inplace:
#             adata._inplace_subset_obs(adata.obs["keeper_cells"])
#             return None
#         else:
#             return adata
