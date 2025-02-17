# sciduck: single cell data quality control toolkit

`sciduck` provides a transparent framework for single cell/nuclei quality control and annotation. The package was built following **anndata** and **scanpy** to easily fit into related workflows. It provides tranparency to single cell analysis pipelines for reproducible data analysis and annotation.  The Python-based implementation efficiently deals with
datasets of more than one million cells.

<div align="center">
<img
    src="/logo/sciduck_DALLEv0.png"
    width="200"
>
</div>

## Tutorials:

[Quality Control](https://github.com/AllenInstitute/sciduck/blob/main/tutorials/standard_workflow.ipynb)

In the `Quality Control` tutorial we showcase how to perform various range and group contraints.

[Plotting](https://github.com/AllenInstitute/sciduck/blob/main/tutorials/sciduck_plotting_example.ipynb)

In the `Plotting` tutorial we show how to use the sciduck plotting functions to identify quality control cuttoffs.

[Annotation Table (Taxonomy)](https://github.com/AllenInstitute/sciduck/blob/main/tutorials/annotation_table.ipynb)

In the `Annotation Table` tutorial we show how to build an annotation table for use in taxonomy development. Typically what is hosted on google sheets at the Allen Institute.

## Documentation

[Read the docs](https://sciduck.readthedocs.io/en/latest/sciduck.html)

## Installation

`pip install sciduck`
