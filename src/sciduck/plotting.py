import seaborn as sns
import matplotlib.pyplot as plt
from typing import Sequence, Literal
from anndata import AnnData

def violin(
        adata: AnnData,
        y: str,
        groupby: str|None = None,
        log: bool = True,
        scale: Literal["area", "count", "width"] = "width",
        order: Sequence[str] | None = None,
        xlabel: str = "",
        size: int = 4,
        x_label_rotation: float = 0,
        ylines: Sequence[float] = []
        ):
    """
    Violin plot.

    Parameters
    ----------
    adata
        Annotated data matrix.
    y
        columns in `.obs` to plot for y axis.
    groupby
        The key of the observation grouping to use for x axis.
    log
        Use log scale for y.
    scale
        The method used to scale the width of each violin.
        If 'width' (the default), each violin will have the same width.
        If 'area', each violin will have the same area.
        If 'count', a violin’s width corresponds to the number of observations.
    order
        Order of x.
    xlabel
        Label of the x axis. Defaults to `groupby`.
    size
        Size of the x labels.
    x_label_rotation
        Rotation of xtick labels.
    ylines
        Add horizontal lines to the plot.

    Returns
    -------
    Tuple[Axes, Axes]
        A tuple containing two matplotlib Axes objects:
        - The first Axes object is the bar plot showing the count of observations for each group.
        - The second Axes object is the violin plot showing the distribution of values in `y` for each group.
    Examples
    --------
    import sciduck as sd
    # Call the function
    ax1, ax2=sd.pl.violin(adata, groupby='AIT21_Class', y='n_genes_by_counts', x_label_rotation=20, size=12, ylines=[1000,2000,3000])
    # Customize further if needed
    ax2.set_title("Violin Plot")

    """

    count_values = adata.obs[groupby].value_counts().reset_index()
    count_values.columns = [groupby, 'cell_counts']
    if order is None:
        order = count_values[groupby].tolist()

    palette = sns.color_palette("husl", len(count_values))

    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.5, 2])  # Adjust ratios as needed


    # Create the bar plot on top of the violin plot
    ax1 = fig.add_subplot(gs[0, 0])
    sns.barplot(x=groupby, y='cell_counts', data=count_values, palette=palette, ax=ax1, order=order)
    ax1.set_xlabel('')  # Remove x-axis label
    ax1.set_xticks([])  # Remove x-ticks
    ax1.set_xticklabels([])  # Remove x-tick labels
    if log:
        ax1.set_yscale('log')

    # Create the violin plot
    ax2 = fig.add_subplot(gs[1, 0])
    sns.violinplot(x=groupby, y=y, data=adata.obs, inner=None, palette=palette, ax=ax2, order=order, density_norm=scale)
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=x_label_rotation, ha='right', size=size)  # Rotate x labels
    if log:
        ax2.set_yscale('log')
    if xlabel !="":
        ax2.set_xlabel(xlabel)
    for y in ylines:  # add horizontal lines
        ax2.axhline(y=y, color='red', linestyle='--', linewidth=1)

    plt.tight_layout()

    return ax1, ax2