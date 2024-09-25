import numpy as np
import pandas as pd
from anndata import AnnData
import warnings

def _check_adata(adata, column):
    if not isinstance(adata, AnnData):
        warnings.warn('adata does not appear to be an anndata object', UserWarning)
    if 'qc_constraints' not in adata.uns.keys():
        adata.uns['qc_constraints'] = {}
    if column not in adata.uns['qc_constraints'].keys():
        adata.uns['qc_constraints'][column] = {}

def add_range_constraint(adata, column, gt=None, lt=None, subset=None, subset_values=None):
    """
    Add a range constraint for a specific column in an AnnData object. This function allows filtering
    cells where values in a column fall within the specified range [gt, lt]. Open intervals are supported 
    by leaving 'gt' or 'lt' as None. Subsetting based on an additional column will apply the condition to
    the subset and retain the rest of the cells.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object where constraints will be added under adata.uns['qc_constraints'].
    column : str
        The column in 'adata.obs' for which the range constraint will be applied.
    gt : float or None, optional
        The lower bound of the range (inclusive). If None, the lower bound is open.
    lt : float or None, optional
        The upper bound of the range (inclusive). If None, the upper bound is open.
    subset : str or None, optional
        The column name in 'adata.obs' to use for subsetting the data.
    subset_values : list or None, optional
        A list of values to subset on, applied to 'subset' column.
    """
    constraint = {
        'gt': gt,
        'lt': lt,
        'subset': {},
        'exclude': [],
        'groupby': None
    }

    if subset and subset_values:
        constraint['subset'][subset] = subset_values

    _check_adata(adata, column)
    
    # Assign a unique key for each constraint
    constraint_key = f'constraint_{len(adata.uns["qc_constraints"][column]) + 1}'
    adata.uns['qc_constraints'][column][constraint_key] = constraint

def add_exclude_constraint(adata, column, exclude_values, subset=None, subset_values=None):
    """
    Add an exclude constraint for a specific column in an AnnData object. This function filters
    out cells where the column value is in the 'exclude_values' list. Subsetting based on 
    additional conditions is also supported.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object where constraints will be added.
    column : str
        The column in 'adata.obs' for which the exclude constraint will be applied.
    exclude_values : list
        A list of values that should be excluded from 'column'.
    subset : str or None, optional
        The column name in 'adata.obs' to use for subsetting the data.
    subset_values : list or None, optional
        A list of values to subset on, applied to 'subset' column.
    """
    
    constraint = {
        'gt': None,
        'lt': None,
        'subset': {},
        'exclude': exclude_values,
        'groupby': None
    }

    if subset and subset_values:
        constraint['subset'][subset] = subset_values

    _check_adata(adata, column)
    
    # Assign a unique key for each constraint
    constraint_key = f'constraint_{len(adata.uns["qc_constraints"][column]) + 1}'
    adata.uns['qc_constraints'][column][constraint_key] = constraint

def add_group_level_constraint(adata, column, groupby, gt=None, lt=None, agg_func='mean'):
    """
    Add a group-level constraint for a specific column in an AnnData object. Cells are grouped 
    by the 'groupby' column, and an aggregation function (e.g., 'mean') is applied to each group 
    for the specified column. Cells in groups that meet the 'gt' and 'lt' conditions are kept.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object where constraints will be added.
    column : str
        The column in 'adata.obs' for which the group-level constraint will be applied.
    groupby : str
        The column in 'adata.obs' used to group the cells (e.g., cluster ID).
    gt : float or None, optional
        The lower bound for the aggregated group value. If None, the lower bound is open.
    lt : float or None, optional
        The upper bound for the aggregated group value. If None, the upper bound is open.
    agg_func : str, optional
        The aggregation function to apply to the groups. Can be 'mean', 'sum', 'std', or 'median'.
        Default is 'mean'.
    """
    constraint = {
        'gt': gt,
        'lt': lt,
        'subset': {},
        'exclude': [],
        'groupby': groupby,
        'agg_func': agg_func
    }

    _check_adata(adata, column)
    
    # Assign a unique key for each constraint
    constraint_key = f'constraint_{len(adata.uns["qc_constraints"][column]) + 1}'
    adata.uns['qc_constraints'][column][constraint_key] = constraint

def apply_constraints(adata, inplace=False):
    """
    Apply all quality control (QC) constraints stored in 'adata.uns["qc_constraints"]' to the 'adata.obs' 
    DataFrame. This function handles both individual cell-level constraints (range, exclude) and 
    group-level constraints (using a 'groupby' column). The resulting 'keeper_cells' column in 'adata.obs'
    will indicate which cells passed the constraints.

    Parameters:
    -----------
    adata : AnnData
        The AnnData object where constraints are applied.
    inplace : bool, optional
        If True, only the cells passing the constraints will be retained in 'adata'.
        If False (default), the 'keeper_cells' column will be added to 'adata.obs' to indicate 
        which cells passed the filtering.
    
    Returns:
    --------
    AnnData
        The filtered AnnData object with the 'keeper_cells' column added to 'adata.obs'.
        If 'inplace=True', the returned object will only contain the filtered cells.
    """
    if isinstance(adata, AnnData):
        keeper_cells = pd.Series([True] * adata.n_obs, index=adata.obs.index)
        agg_funcs = {'sum': np.sum, 'mean': np.mean, 'std': np.std, 'median': np.median}

        for col, constraints_dict in adata.uns['qc_constraints'].items():
            if col not in adata.obs.columns:
                warnings.warn(f'Column {col} not in obs, skipping this metric', UserWarning)
                continue

            column_data = adata.obs[col]
            col_keeper = pd.Series([True] * adata.n_obs, index=adata.obs.index)

            # Iterate over the dictionary of constraints
            for constraint_key, constraints in constraints_dict.items():
                gt = pd.Series([True] * adata.n_obs, index=adata.obs.index)
                lt = pd.Series([True] * adata.n_obs, index=adata.obs.index)
                exclude = pd.Series([True] * adata.n_obs, index=adata.obs.index)
                subset = pd.Series([True] * adata.n_obs, index=adata.obs.index)

                # Apply group-level constraint
                if 'groupby' in constraints and constraints['groupby'] is not None:
                    groupby = constraints['groupby']
                    agg_func = agg_funcs[constraints['agg_func']]

                    if groupby in adata.obs:
                        grouped = adata.obs.groupby(groupby)[col].agg(agg_func)

                        group_gt = pd.Series([True] * len(grouped), index=grouped.index)
                        group_lt = pd.Series([True] * len(grouped), index=grouped.index)

                        if constraints['gt'] is not None:
                            group_gt = grouped >= constraints['gt']

                        if constraints['lt'] is not None:
                            group_lt = grouped <= constraints['lt']

                        group_keeper = group_gt & group_lt
                        groups_to_keep = group_keeper[group_keeper].index
                        cur_sub = adata.obs[groupby].isin(groups_to_keep)
                        col_keeper &= cur_sub
                else:
                    if constraints['gt'] is not None:
                        gt = column_data >= constraints['gt']

                    if constraints['lt'] is not None:
                        lt = column_data <= constraints['lt']

                    if constraints['exclude']:
                        exclude = ~column_data.isin(constraints['exclude'])

                    if constraints['subset'].items():
                        for subset_col, subset_values in constraints['subset'].items():
                            subset_in = adata.obs[subset_col].isin(subset_values) if subset_col in adata.obs else pd.Series([True] * adata.n_obs, index=adata.obs.index)
                            cur_sub = (gt & lt & exclude)
                            cur_sub |= ~subset_in
                            col_keeper &= cur_sub
                    else:
                        col_keeper &= (gt & lt & exclude)

                keeper_cells &= col_keeper

        adata.obs['keeper_cells'] = keeper_cells

        if inplace:
            adata._inplace_subset_obs(adata.obs['keeper_cells'])

        return adata

