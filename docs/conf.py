# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'data-driven-representative-models-of-biobased-complex-fluids'
copyright = '2025, Daniel J. York, Isaac Vidal-Daza, Francisco J. Martin-Martinez'
author = 'Daniel J. York, Isaac Vidal-Daza, Francisco J. Martin-Martinez'
release = '2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',  # Optional for Google/NumPy docstring support
]

# Paths and templates
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# HTML output
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    "collapse_navigation": False,  # Ensures the sidebar is expandable
    "sticky_navigation": True,     # Keeps the sidebar visible while scrolling
    "navigation_depth": 5,         # Allows deeper levels of nesting
}

html_static_path = ['_static']

master_doc = 'index'

# Autodoc default options (you can modify this as needed)
autodoc_default_options = {
    'members': True,            # Include members (functions, classes, etc.)
    'undoc-members': True,      # Include members without docstrings
    'show-inheritance': True,   # Show inheritance in class documentation
}
