# -*- coding: utf-8 -*-

from ..utility import generateFile

class ShapePackUnpack:
   def __init__(self, shapes):
      self.shapes = shapes
   def generate(self, path):
      context = dict()
      context["shapes"] = self.shapes
      generateFile(path, 'mpi/ShapePackUnpack.templ.h', context)
