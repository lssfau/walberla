#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import sphinx_rtd_theme

sys.path.append(os.path.abspath('.'))


html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_theme = 'sphinx_rtd_theme'


extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'doxylink',
]

doxylink_baseurl = "http://www.walberla.net/doxygen"

source_suffix = '.rst'
master_doc = 'index'

project = 'waLBerla'
copyright = '2020, LSS waLBerla Team'

version = ''
release = ''

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

htmlhelp_basename = 'waLBerladoc'


intersphinx_mapping = {'python': ('http://docs.python.org/3.8', None),
                       'numpy': ('http://docs.scipy.org/doc/numpy/', None),
                       'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
                       'matplotlib': ('http://matplotlib.sourceforge.net/', None)}
