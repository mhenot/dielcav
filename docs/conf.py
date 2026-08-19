# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "dielcav"
copyright = "2026, Marceau Hénot"
author = "Marceau Hénot"
release = "0.1"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "nbsphinx",
    "nbsphinx_link",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

nbsphinx_execute = "never"

autodoc_member_order = "bysource"

html_theme = "alabaster"
