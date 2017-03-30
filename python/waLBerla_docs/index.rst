Documentation of waLBerla's Python interface
============================================

Quick Start / Tutorials:
------------------------

You can quickly try out waLBerla's Python interface in a hosted IPython notebook without
having to install anything:

http://demo.walberla.net

This site contains interactive tutorials, illustrating how to set up a 2D lattice Boltzmann
simulation with waLBerla.


Installation with conda:
------------------------

To run waLBerla on your own machine the simplest way to get going is the installation via
the `conda package manager <http://conda.pydata.org>`_::

   conda install --channel lssfau walberla


Run in docker:
--------------

Docker is a lightweight virtualization solution. We provide a docker image that
contains the same environment as hosted on http://demo.walberla.net.
With this image you can run and develop waLBerla simulations on your own machine without having to manually
install the dependencies. All you need is a running installation of  `Docker <www.docker.com>`_.
Run the waLBerla image with the following command and navigate in your browser to http://localhost:8888 ::

   docker run -it -p 8888:8888 walberla/runenv-ubuntu-python


API Documentation:
------------------


.. toctree::
   :maxdepth: 2

   modules/blockforest
   modules/core
   modules/field
   modules/geometry
   modules/lbm
   modules/plot
   modules/tools
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

