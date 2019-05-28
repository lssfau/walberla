# -*- coding: utf-8 -*-

from jinja2 import Environment, FileSystemLoader
import os

class TerminalColor:
   DEFAULT = '\033[0m'
   GREEN   = '\033[92m'
   YELLOW  = '\033[93m'
   RED     = '\033[91m'

def find(f, seq):
   """Return first item in sequence where f(item) == True."""
   for item in seq:
      if f(item):
         return item
   return None

def capFirst(s):
   return s[0].capitalize() + s[1:]

def getJinjaEnvironment():
   dirname = os.path.dirname(__file__)
   env = Environment(loader=FileSystemLoader(dirname + '/templates'))
   env.filters['capFirst'] = capFirst
   return env

def checkInterface(path, template, accessor):
   with open(path + template, 'r') as myfile:
      data = myfile.read()
   for prop in accessor.properties:
      for func in prop.getFunctionNames():
         if not (func in data):
            raise RuntimeError("{} required but not used in kernel ({})".format(func, template))

def generateFile(path, template, context = {}, filename=None):
   if filename==None:
      filename = template.replace(".templ", "")
   env = getJinjaEnvironment()
   print("generating: " + path + filename)
   fout = open(path + filename, "wb")
   content = env.get_template(template).render(context)
   fout.write(content.encode('utf8'))
   fout.close()
