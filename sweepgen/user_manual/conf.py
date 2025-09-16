# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "waLBerla SweepGen"
copyright = "2025, Frederik Hennig"
author = "Frederik Hennig"

from sweepgen import __version__

version = __version__

language = "en"
default_role = "any"
pygments_style = "sphinx"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.extlinks",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "attrs_inline",
    "attrs_block",
]

autodoc_member_order = "bysource"
autodoc_typehints = "description"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3.8", None),
    "numpy": ("https://docs.scipy.org/doc/numpy/", None),
    "sympy": ("https://docs.sympy.org/latest/", None),
    "pystencils": ("https://pycodegen.pages.i10git.cs.fau.de/docs/pystencils/2.0dev/", None),
}

BRANCH = "codegen-update"  # FIXME

extlinks = {
    "walberla-example-app": (
        f"https://i10git.cs.fau.de/walberla/walberla/-/tree/{BRANCH}/apps/examples/%s",
        "%s"
    )
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_theme_options = {
    "logo": {
        "image_light": "_static/sweepgen_light.svg",
        "image_dark": "_static/sweepgen_dark.svg",
    }
}
html_favicon = "_static/favicon.svg"
