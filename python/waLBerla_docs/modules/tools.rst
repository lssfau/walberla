************
tools module
************

Jobscript Generation
====================

.. automodule:: waLBerla.tools.jobscripts
   :members:
   
   
Sqlite3 Helpers
===============

.. autofunction:: waLBerla.tools.sqlitedb.storeSingle

.. autofunction:: waLBerla.tools.sqlitedb.storeMultiple

.. autofunction:: waLBerla.tools.sqlitedb.checkAndUpdateSchema

.. autofunction:: waLBerla.tools.sqlitedb.mergeSqliteFiles


PRM Files
=========

.. automodule:: waLBerla.tools.config
   :members:


LBM Unit Conversion
===================

.. autoclass:: waLBerla.tools.lbm_unitconversion.PintUnitConverter
   :members:
   :undoc-members:
   :special-members: __init__


.. autofunction:: waLBerla.tools.lbm_unitconversion.extractLatticeFactors

.. autofunction:: waLBerla.tools.lbm_unitconversion.computeLatticeFactors


HTML Reports
============

.. automodule:: waLBerla.tools.report
   :members:

.. autofunction:: waLBerla.tools.report.runWebserver

.. autofunction:: waLBerla.tools.report.cliFrontend

.. autofunction:: waLBerla.tools.report.generate
   
.. autofunction:: waLBerla.tools.report.setupFlaskApp
   