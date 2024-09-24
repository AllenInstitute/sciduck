# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../'))

# The full version, including alpha/beta/rc tags
release = '00.00.01'

# -- General configuration ------------------------------------------------

extensions = [
    "nbsphinx",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.linkcode",
    "sphinx_autodoc_typehints",
    "sphinx_plotly_directive",
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

myst_enable_extensions = [
    "amsmath",
    #"colon_fence",
    #"deflist",
    "dollarmath",
    #"fieldlist",
    #"html_admonition",
    #"html_image",
    #"linkify",
    #"replacements",
    #"smartquotes",
    #"strikethrough",
    #"substitution",
    #"tasklist",
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'
# autodoc_default_flags = ['members']
napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False

intersphinx_mapping = {
    "anndata": ("https://anndata.readthedocs.io/en/stable/", None),
    "ipython": ("https://ipython.readthedocs.io/en/stable/", None),
    "matplotlib": ("https://matplotlib.org/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "python": ("https://docs.python.org/3", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/reference/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
    "scanpy": ("https://scanpy.readthedocs.io/en/stable/", None),
}

smv_branch_whitelist = r'main'  # Include all branches

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'sciduck'
copyright = '2024, Nelson Johansen'
author = 'Nelson Johansen'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# pip install sphinx-rtd-theme; https://sphinx-themes.org/sample-sites/sphinx-rtd-theme/

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
