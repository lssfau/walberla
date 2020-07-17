"""This module creates HTML reports based on the Flask and Jinja2 packages.

.. attention::
    This module is not actively developed any more - consider using IPython Notebook instead

It provides custom template markup for generating matplotlib graphs with data from a sqlite3 database.
To use this module one first has to write the simulation results to a sqlite3 database e.g. using waLBerla.tools.sqlite

Then create a Jinja2 HTML template. Here is an example using base templates provided by this module::

    {% extends "waLBerla/bootstrap_report.html" %}
    {% block content %}
    <div class="container">
      <div class="text-center">
        <h1> MySetup </h1>

        {% matplotlib %}
            # Some SQL query that returns a two column result
            q = "SELECT capillaryNr,shapeFactor FROM runs WHERE yCells=160 AND zCells=160 AND bubbleDiameter=60 ORDER BY capillaryNr"  # noqa: E501
            # plt is matplotlib.pyplot extended with the custom 'dbplot' function
            plt.dbplot( q, label="sim shapefactor", marker="o" )
            plt.legend( loc='center')
        {% endmatplotlib %}
      </div>
    </div>

Then use one of the functions provided by this module to create an HTML page from this template.
Either once or run a small webserver that reloads the template if it was changed.
"""

from __future__ import print_function, absolute_import, division, unicode_literals

from waLBerla.tools.report.report import setupFlaskApp, runWebserver, generate, cliFrontend

__all__ = ['setupFlaskApp', 'runWebserver', 'generate', 'cliFrontend']
