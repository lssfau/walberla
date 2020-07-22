# from docutils.parsers.rst import Directive
# from sphinx.addnodes import download_reference
# from sphinx.writers.html import HTMLTranslator
from sphinx.domains import Domain
from docutils import nodes


class DoxyDomain(Domain):
    name = "doxylink"


def generate_doxygen_link(type, value):
    assert (type == "class" or type == "struct" or type == "namespace" or type == "file")
    value = value.replace("_", "__")
    value = value.replace(":", "_1")
    value = value.replace("/", "_2")
    value = value.replace(".", "_8")

    if type == "file":
        type = ""

    result = type + value + ".html"
    return result


def doxylink_class(*args, **kwargs):
    return doxygenlink_role("class", *args, **kwargs)


def doxylink_struct(*args, **kwargs):
    return doxygenlink_role("struct", *args, **kwargs)


def doxylink_namespace(*args, **kwargs):
    return doxygenlink_role("namespace", *args, **kwargs)


def doxylink_file(*args, **kwargs):
    return doxygenlink_role("file", *args, **kwargs)


def doxygenlink_role(type, name, rawtext, text, lineno, inliner, options={}, content=[]):
    try:
        app = inliner.document.settings.env.app
        doxylink = generate_doxygen_link(type, text)
        node = make_link_node(rawtext, app, text, doxylink, options)
        return [node], []
    except ValueError:
        msg = inliner.reporter.error(
            "Error parsing doxylink. It has to be of form [class,struct,namespace]:identifierName")
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]


def make_link_node(rawtext, app, text, doxylink, options):
    try:
        base = app.config.doxylink_baseurl
        if not base:
            raise AttributeError
    except AttributeError as err:
        raise ValueError('doxygenlink_base_url configuration value is not set (%s)' % str(err))
    #
    slash = '/' if base[-1] != '/' else ''
    ref = base + slash + '/' + doxylink
    node = nodes.reference(rawtext, text, refuri=ref, reftitle='(C++ Documentation)')
    return node


def setup(app):
    app.add_domain(DoxyDomain)

    app.add_role_to_domain('doxylink', 'class', doxylink_class)
    app.add_role_to_domain('doxylink', 'struct', doxylink_struct)
    app.add_role_to_domain('doxylink', 'namespace', doxylink_namespace)
    app.add_role_to_domain('doxylink', 'file', doxylink_file)

    app.add_config_value('doxylink_baseurl', None, 'env')
    return
