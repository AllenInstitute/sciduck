from anndata import AnnData
import scanpy as sc
import numpy as np
import pandas as pd

def build_annotation_table(adata: AnnData, 
                            group_by = "cluster",
                            categorical_annotations: list = ["donor_name", "load_name", "roi"],
                            numeric_annotations: list = ["doublet_score", "pct_counts_mt"],
                            min_percent: float = 0.05,
                            annotation_alerts: dict = {"donor_name": 0.90, "load_name": 0.95, "roi": 0.95},
                            mapping_summary: dict = {}) -> dict:
    """
    Build a standardized table of annotations describing each cluster to be used for taxonomy development or communication.
    
    Args:
        adata: Anndata object with `annotations` in `obs`.
        group_by: Anndata.obs to use for grouping cells.
        categorical_annotations: List of categorical metadata to include.
        numeric_annotations: List of numeric metadata to include.
        min_percent: Minimal percentrage to print annotation.
        mapping_summary: A dictionary to store the mapping summary. User can pass in an existing dictionary to append to if desired.
    
    Returns:
        A mapping_summary containing character and numeric summaries along with donor/lib/roi composition alerts.
    """
    ## Annotations summaries
    for anno in categorical_annotations:
        print("Processing annotation: ", anno)
        freq = adata.obs.groupby([group_by, anno])[anno].size()
        mapping_summary[anno] = {}
        if anno in annotation_alerts.keys() and f"{anno}_composition_alert" not in mapping_summary:
            mapping_summary[f"{anno}_composition_alert"] = {}
        for cluster in np.unique(adata.obs[group_by]):
            cluster_anno = freq.loc[freq.index.get_level_values(group_by) == cluster].sort_values(ascending=False)
            cluster_anno = np.round(cluster_anno / cluster_anno.sum(), 2)
            record = ""
            for level in cluster_anno.index.get_level_values(anno):
                proportion = cluster_anno.loc[cluster_anno.index.get_level_values(anno) == level].values[0]
                if proportion > min_percent: ## Only show the annotation term when its percentage > 0.05
                    record += level + "(" + str(proportion) + ")" if level == cluster_anno.index.get_level_values(anno)[0] else "| " + level + "(" + str(proportion) + ")" 
            mapping_summary[anno][cluster] = record
            if anno in annotation_alerts.keys():
                alert_status = "Balanced" if cluster_anno.max() < annotation_alerts[anno] else "Cautious"
                mapping_summary[f"{anno}_composition_alert"][cluster] = alert_status

    ## Calculate the median values for basic qc stats
    mapping_summary[str(group_by) + "_size"] = adata.obs.groupby(group_by).size()

    ## Calculate the median values for numeric annotations each cluster
    for col in numeric_annotations:
        print("Processing annotation: ", col)
        mapping_summary[f'{col}_median'] = adata.obs.groupby(group_by)[col].median().round(4)

    return pd.DataFrame.from_dict(mapping_summary)

def add_dominant_library_info(adata: AnnData, 
                    mapping_summary: dict,
                    library_metadata_column: str = "load_name",
                    ) -> dict:
    """
    Add dominant library information to the mapping summary.

    Args:
        adata: Anndata object with `annotations` in `obs`.
        mapping_summary: A pre-existing mapping summary dictionary to add information on dominant library.
        library_metadata_column: Defaults to `load_name` and describes the sequencing batch.
    
    Returns:
        A mapping_summary containing dominant library information.
    """

    ## Calculate the total number of cells for each library
    library_freq_total = adata.obs.groupby(library_metadata_column)[library_metadata_column].size()
    mapping_summary["dominant_library_percentage"] = {}
    mapping_summary["dominant_library_info"] = {}

    ## Add dominant library information to the mapping summary for each cluster
    for cluster in np.unique(adata.obs.cluster):

        ## Calculate the library composition for each cluster
        cluster_library_freq = adata.obs[adata.obs["cluster"] == cluster].groupby(library_metadata_column).size()
        cluster_library_normalized = np.round(cluster_library_freq / cluster_library_freq.sum(), 4)      

        ## Information for the dominant library in the cluster
        dominant_library = cluster_library_freq.idxmax()
        dominant_percentage = cluster_library_normalized[dominant_library]
        total_cells_in_dominant_library = library_freq_total[dominant_library]

        ## Update the mapping summary
        mapping_summary["dominant_library_percentage"][cluster] = dominant_percentage
        mapping_summary["dominant_library_info"][cluster] = f"{dominant_library} ({dominant_percentage})% | {total_cells_in_dominant_library}# cells."

    return mapping_summary
