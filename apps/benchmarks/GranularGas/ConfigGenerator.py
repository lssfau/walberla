# -*- coding: utf-8 -*-

from jinja2 import Environment, FileSystemLoader
import os

class Parameter:
   def __init__(self, name, type, defValue=""):
      """Propery of a data strcuture

      Parameters
      ----------
      name : str
         name of the property
      type : str
         type of the property
      defValue : str
         default value the property should be initialized with
      """

      self.name     = name
      self.type     = type
      self.defValue = defValue

   def __str__(self):
      return "name: {}, type: {}, defValue: {}".format(self.name, self.type, self.defValue)

class Config:
   def __init__(self):
      self.parameters = []

   def parameterExists(self, name):
      for v in self.parameters:
         if v.name==name:
            return True
      return False

   def addParameter(self, name, type, defValue):
      if self.parameterExists( name ):
         print("parameters already added: " + name)
      else:
         self.parameters.append( Parameter(name, type, defValue) )

   def generateFile(self, template):
      context = dict()
      context["parameters"]    = self.parameters

      path = ""
      filename = template.replace(".templ", "")
      dirname = os.path.dirname(__file__)
      env = Environment(loader=FileSystemLoader(dirname))
      print("generating: " + path + filename)
      fout = open(path + filename, "wb")
      content = env.get_template(template).render(context)
      fout.write(content.encode('utf8'))
      fout.close()

   def generate(self):
      print("="*90)
      print("Config File:")
      print("")
      print("{0: <30}{1: <30}{2: <30}".format("Name", "Type", "Def. Value"))
      print("="*90)
      for param in self.parameters:
         print("{0: <30.29}{1: <30.29}{2: <30.29}".format(param.name, param.type, param.defValue))
      print("="*90)

      self.generateFile("Parameters.templ.h")
      self.generateFile("Parameters.templ.cpp")
