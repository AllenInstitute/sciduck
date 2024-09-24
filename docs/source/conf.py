# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import pathlib
sys.path.insert(0, os.path.abspath('../../src/'))

# -- General configuration ------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
]

# Generate the API documentation when building
# autosummary_generate = True
# autodoc_member_order = 'bysource'

napoleon_google_docstring = True
napoleon_numpy_docstring = False

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

# Autodoc settings
autodoc_default_options = {
    'members': True,
    'undoc-members': False,  # Hide undocumented members
    'private-members': False,  # Hide private members (_member)
    'special-members': False,  # Hide special members (__member__)
    'inherited-members': False,  # Hide inherited members
    'show-inheritance': True,
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
# pip install sphinx-rtd-theme; https://sphinx-themes.org/sample-sites/sphinx-rtd-theme/

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
