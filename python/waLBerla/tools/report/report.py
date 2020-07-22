from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import io
import sqlite3
from glob import glob
from functools import partial
import matplotlib.pyplot as plt
from matplotlib import rc

from jinja2 import Environment, PackageLoader, PrefixLoader, FileSystemLoader, ChoiceLoader
from jinja2 import nodes
from jinja2.ext import Extension

try:
    import mpld3

    mpld3_available = True
except ImportError:
    mpld3_available = False

rc('font', size='8.0')
# matplotlib.style.use('ggplot')


# mpl3d is a library for displaying interactive plots in html pages
# however some features like legends in bar plots are not supported yet
# thus this features is by default disabled
use_mpld3 = False

write_pdfs = True
pdf_output_dir = "pdfs"


def plot_to_html(name=None):
    """Converts a matplotlib figure, such that it can be displayed in an html page
        either using an embedded svg tag (default)
        or by using mpld3
    """
    result = ""
    if mpld3_available and use_mpld3:
        html = mpld3.fig_to_html(plt.gcf())
        result = html
    else:
        figStr = io.StringIO()
        plt.savefig(figStr, format='svg', bbox_inches='tight')
        result = figStr.getvalue()

    if write_pdfs and name:
        if not os.path.exists(pdf_output_dir):
            os.makedirs(pdf_output_dir)

        filename = os.path.join(pdf_output_dir, name + '.pdf')
        filepath = os.path.join(filename)
        plt.savefig(filepath, bbox_inches='tight')
        result = '<a href="%s" > %s </a>' % (filename, result)

    plt.clf()

    return result


def dbQuery(database_getter, query):
    """Executes an SQL query on the given sqlite database handle
       Returns tuple of two elements, first element is a list of column names,
       the second element contains the data (row by row)"""
    database_handle = database_getter()
    c = database_handle.cursor()
    c.execute(query)
    column_names = []
    for column_description in c.description:
        column_names += [column_description[0]]
    return column_names, c.fetchall()


def numpy_arr_from_db(database_getter, query):
    """Executes sqlite query and returns result as numpy arrays
        Returns 2-tuple: xVals, [yVals0, yVals1, ... ] (depending on number of result columns
                         yVals* are numpy arrays, xVals can also be of other type (f.e. strings)
    """
    column_names, data = dbQuery(database_getter, query)
    xVals = [data[i][0] for i in range(len(data))]
    yArrays = []
    if data:
        for j in range(1, len(data[0])):
            yVals = np.array([data[i][j] for i in range(len(data))])
            yArrays.append(yVals)

    return xVals, yArrays


def dbPlot(database_getter, query, **kwargs):
    column_names, data = dbQuery(database_getter, query)

    xVals, yArrays = numpy_arr_from_db(database_getter, query)
    for j, yVals in enumerate(yArrays):
        for value in yVals:
            assert (isinstance(value, (int, float)))
        if 'label' in kwargs:
            plt.plot(xVals, yVals, **kwargs)
        else:
            plt.plot(xVals, yVals, label=column_names[j], **kwargs)


def remove_indentation(input):
    input = input.replace('\t', '    ')
    indentation = min([len(s) - len(s.lstrip()) for s in input.splitlines(False) if len(s.lstrip()) > 0])
    return "\n".join([s[indentation:] for s in input.splitlines(False)])


class MatplotlibExtension(Extension):
    # a set of names that trigger the extension.
    tags = set(['matplotlib', 'plot'])

    def __init__(self, environment):
        super(MatplotlibExtension, self).__init__(environment)

        # add the defaults to the environment
        environment.extend(
            plt_object=plt,
            output_folder=".",
            plot_funcs={},
            vars={},
            database_file="database.sqlite"
        )

    def parse(self, parser):
        lineno = next(parser.stream).lineno
        args = []
        if parser.stream.current.type != 'block_end':
            lineno = parser.stream.current.lineno
            args.append(parser.parse_expression())
        else:
            args.append(nodes.Const(None))

        body = parser.parse_statements(['name:endmatplotlib', 'name:endplot', 'name:end_matplotlib', 'name:end_plot'],
                                       drop_needle=True)
        args.append(nodes.ContextReference())
        args.append(nodes.Name('i', 'load'))
        return nodes.CallBlock(self.call_method('_execute_matplotlib', args),
                               [], [], body).set_lineno(lineno)

    def _execute_matplotlib(self, name, context, i, caller):
        tag_body = caller()
        plt.dbplot = partial(dbPlot, self.environment.database_getter)
        context = context.get_all()

        context.update({'plt': plt,
                        'vars': self.environment.vars,
                        'numpy_arr_from_db': partial(numpy_arr_from_db, self.environment.database_getter),
                        'db_query': partial(dbQuery, self.environment.database_getter),
                        'np': np})

        context.update(self.environment.plot_funcs)
        context.update({'i': i})
        try:
            exec(remove_indentation(tag_body), context)
        except Exception as e:
            print(remove_indentation(tag_body))
            raise e

        return plot_to_html(name)


# User Functions


def setupFlaskApp(app, context={}, plot_funcs={}, database_file="database.sqlite"):
    """Call this function to set up your custom flask app. This is an advanced function,
       probably you want to use runWebserver instead"""
    from flask import g

    def get_db():
        db = getattr(g, '_database', None)
        if db is None:
            db = g._database = sqlite3.connect(database_file)
        return db

    app.jinja_loader = ChoiceLoader([FileSystemLoader('.'),
                                     PrefixLoader(
                                         {'waLBerla': PackageLoader('waLBerla', 'tools', 'report', 'templates')})])

    app.jinja_env.add_extension(MatplotlibExtension)

    app.jinja_env.database_file = database_file
    app.jinja_env.database_getter = get_db
    app.jinja_env.plot_funcs = plot_funcs


def runWebserver(context={}, plot_funcs={}, database_file="database.sqlite",
                 open_browser=False, debug=True):
    """Runs a small local webserver using the flask module ( pip3 install flask ) serving the report.
       When refreshing in the browser the report is updated."""
    from flask import Flask
    from flask import render_template, send_from_directory

    context['database_file'] = database_file

    app = Flask(__name__)
    app.debug = True
    setupFlaskApp(app, context=context, plot_funcs=plot_funcs, database_file=database_file)

    @app.route("/")
    def main_route():
        reports = [n[2:-5] for n in glob("t_*.html")]
        return render_template("waLBerla/report_overview.html", reports=reports, **context)

    @app.route("/<template_name>")
    def report_route(template_name):
        template_name = "t_" + template_name + ".html"
        return render_template(template_name, **context)

    @app.route("/pdfs/<path:filename>")
    def pdf_route(filename):
        return send_from_directory(os.path.join(os.getcwd(), pdf_output_dir), filename)

    if open_browser:
        import webbrowser
        webbrowser.open('http://127.0.0.1:5000/')

    app.run(debug=debug)


def generate(context={}, plot_funcs={}, database_file="database.sqlite", input_output_list=None, open_browser=False):
    """Generates a html report. Uses jinja2 templating engine with custom
       matplotlib tag, which inserts matplotlib figures as svg graphics"""

    def get_db():
        if get_db.databaseHandle is None:
            get_db.databaseHandle = sqlite3.connect(database_file)
        return get_db.databaseHandle

    get_db.databaseHandle = None

    context['database_file'] = database_file

    loaders = ChoiceLoader([FileSystemLoader('.'),
                            PrefixLoader({'waLBerla': PackageLoader('waLBerla', 'tools', 'report', 'templates')})])

    env = Environment(loader=loaders, extensions=[MatplotlibExtension])
    env.database_file = database_file
    env.database_getter = get_db
    env.plot_funcs = plot_funcs

    if input_output_list is None:
        input_output_list = [(input, input[2:]) for input in glob("t_*.html")]
        context.update({'reports': [r[1] for r in input_output_list]})
        input_output_list.append(('waLBerla/report_overview.html', 'index.html'))

    for e in input_output_list:
        template = env.get_template(e[0])
        with open(e[1], "w") as f:
            print("-- Generating " + e[1])
            f.write(template.render(**context))

    if open_browser:
        import webbrowser
        webbrowser.open(os.path.join(os.getcwd(), 'index.html'))


def cliFrontend(context={}):
    """Command line interface for above functions"""
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--interactive", action="store_true",
                        help="Run a webserver to interactively develop your reports. (flask required)")
    parser.add_argument("-d", "--database", default="database.sqlite", help="Database file (default: database.sqlite)")
    parser.add_argument("-b", "--open_browser", action="store_true", help="Open Browser (on")

    args = parser.parse_args()
    if args.interactive:
        runWebserver(context=context, database_file=args.database, open_browser=args.open_browser)
    else:
        generate(context=context, database_file=args.database, open_browser=args.open_browser)
