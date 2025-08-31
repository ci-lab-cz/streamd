import os
import sys

# Ensure project root is on sys.path
sys.path.insert(0, os.path.abspath('..'))

project = 'StreaMD'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'myst_parser',
]
autosummary_generate = True
html_theme = 'sphinx_rtd_theme'
# Include both docstring sections
napoleon_google_docstring = False
napoleon_numpy_docstring = True
