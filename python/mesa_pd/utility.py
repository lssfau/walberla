# -*- coding: utf-8 -*-

from jinja2 import Environment, FileSystemLoader
import os


class TerminalColor:
    DEFAULT = '\033[0m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'


def find(f, seq):
    """Return first item in sequence where f(item) == True."""
    for item in seq:
        if f(item):
            return item
    return None


def cap_first(s):
    return s[0].capitalize() + s[1:]


def get_jinja_environment():
    dirname = os.path.dirname(__file__)
    env = Environment(loader=FileSystemLoader(dirname + '/templates'))
    env.filters['capFirst'] = cap_first
    return env


def generate_file(path, template, context={}, filename=None):
    if filename is None:
        filename = template.replace(".templ", "")
    env = get_jinja_environment()
    print(f"generating: {(path / filename)}")
    with open(path / filename, "wb") as fout:
        content = env.get_template(template).render(context)
        fout.write(content.encode('utf8'))
