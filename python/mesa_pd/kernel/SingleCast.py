# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class SingleCast:
   def __init__(self, shapes):
      self.shapes = shapes
      self.accessor = Accessor()
      self.accessor.require("shape", "BaseShape*", access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["interface"] = self.accessor.properties
      context["shapes"]    = self.shapes
      generateFile(path, 'kernel/SingleCast.templ.h', context)
