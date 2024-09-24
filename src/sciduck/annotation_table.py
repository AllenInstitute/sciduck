from anndata import AnnData
import scanpy as sc
import numpy as np

def build_annotation_table(adata: AnnData, 
                            categorical_annotations: list = ["donor_name", "load_name", "roi"],
                            numeric_annotations: list = ["doublet_score", "pct_counts_mt"],
                            min_percent: float = 0.05,
                            annotation_alerts: dict = {"donor_name": 0.90, "load_name": 0.95, "roi": 0.95},
                            add_dominant_library_info: bool = True,
                            add_entropy_info: bool = True,
                            mapping_summary: dict = {}):
    """
    Build a standardized table of annotations describing each cluster to be used for taxonomy development or communication.
    
    Args:
        adata (.h5ad): Anndata object with `annotations` in adata.obs.
        annotations (list of str): Cluster summary annotations specific to data.
        mapping_summary (dict): A dictionary to store the mapping summary. User can pass in an existing dictionary to append to if desired.
    
    Returns:
        dict: A mapping_summary containing character and numeric summaries along with donor/lib/roi composition alerts.
    """
    ## Annotations summaries
    for anno in annotations:
        freq = adata.obs.groupby(["cluster", anno])[anno].size()
        mapping_summary[anno] = {}
        for cluster in np.unique(adata.obs.cluster):
            cluster_anno = freq.loc[freq.index.get_level_values("cluster") == cluster].sort_values(ascending=False)
            cluster_anno = np.round(cluster_anno / cluster_anno.sum(), 2)
            record = ""
            for level in cluster_anno.index.get_level_values(anno):
                proportion = cluster_anno.loc[cluster_anno.index.get_level_values(anno) == level].values[0]
                if proportion > min_percent: ## Only show the annotation term when its percentage > 0.05
                    record += level + "(" + str(proportion) + ")" if level == cluster_anno.index.get_level_values(anno)[0] else "| " + level + "(" + str(proportion) + ")"
            mapping_summary[anno][cluster] = record
            if anno in annotation_alerts.keys():
                alert_status = "Balanced" if cluster_anno.max() < annotation_alerts[anno] else "Cautious"
                mapping_summary[anno + "_composition_alert"][cluster] = alert_status

    ## Calculate the median values for basic qc stats
    mapping_summary["cluster_size"] = adata.obs.groupby("cluster").size()
    mapping_summary["gene.counts"] = adata.obs.groupby("cluster")["n_genes_by_counts"].median().round(4)
    mapping_summary["umi.counts"] = adata.obs.groupby("cluster")["n_counts"].median().round(4)
    mapping_summary["umi_gene_ratio"] = ratio.round(4)

    ## Calculate the median values for numeric annotations each cluster
    for col in numeric_annotations:
        mapping_summary[f'{col}_median'] = adata.obs.groupby('cluster')[col].median().round(4)
    
    return mapping_summary

def add_dominant_library_info(adata: AnnData,
                                mapping_summary: dict
                            ) -> dict:
    """Add dominant library information to the mapping summary.
    Args:
        adata (.h5ad): Anndata object with `annotations` in adata.obs.
        mapping_summary (dict): A pre-existing mapping summary dictionary to add information on dominant library.
    
    Returns:
        dict: A mapping_summary containing dominant library information.

    """

    ## Calculate the total number of cells for each library
    library_freq_total = adata.obs.groupby("load_name")["load_name"].size()
    mapping_summary["dominant_library_percentage"] = {}
    mapping_summary["dominant_library_info"] = {}

    ## Add dominant library information to the mapping summary for each cluster
    for cluster in np.unique(adata.obs.cluster):

        ## Calculate the library composition for each cluster
        cluster_library_freq = adata.obs[adata.obs["cluster"] == cluster].groupby("load_name").size()
        cluster_library_normalized = np.round(cluster_library_freq / cluster_library_freq.sum(), 4)      

        ## Information for the dominant library in the cluster
        dominant_library = cluster_library_freq.idxmax()
        dominant_percentage = cluster_library_normalized[dominant_library]
        total_cells_in_dominant_library = library_freq_total[dominant_library]

        ## Update the mapping summary
        mapping_summary["dominant_library_percentage"][cluster] = dominant_percentage
        mapping_summary["dominant_library_info"][cluster] = f"{dominant_library} ({dominant_percentage})% | {total_cells_in_dominant_library}# cells."

    return mapping_summary