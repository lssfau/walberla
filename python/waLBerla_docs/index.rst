Documentation of waLBerla's Python interface
============================================

Quick Start / Tutorials:
------------------------

This is the documentation of waLBerla`s Python coupling.
It enables to use the core functionality of waLBerla from Python. A primary advantage of
the Python coupling is that it allows having a view on waLBerla`s data as
NumPy arrays. With this, it is easily possible to provide pre and postprocessing
routines using powerful packages from Python. Furthermore, it is possible to set up entire
simulations just from Python. This feature is primarily used in the code generation framework `pystencils <https://pycodegen.pages.i10git.cs.fau.de/pystencils/>`_
within its `parallel data handling <https://pycodegen.pages.i10git.cs.fau.de/pystencils/notebooks/03_tutorial_datahandling.html>`_ which entirely builds upon waLBerla`s Python coupling.
waLBerla`s Python bindings are built with `pybind11 <https://pybind11.readthedocs.io/en/stable/#>`_, a lightweight header-only package which is shipped as a submodule.
Thus, there is no need to install any additional software.


Installation:
------------------------

To install waLberla as Python package in your path run the following command in your build folder ::

   make pythonModuleInstall


API Documentation:
------------------


.. toctree::
   :maxdepth: 2

   modules/blockforest
   modules/core
   modules/field
   modules/plot
   modules/tools
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

