# -*- coding: utf-8 -*-

import numpy as np
from ..utility import generateFile

class ShapeStorage:
   def __init__(self, p, shapes):
      p.addProperty("shapeID",         "size_t",                  defValue="",          syncMode="COPY")
      p.addProperty("rotation",        "walberla::mesa_pd::Rot3", defValue="",          syncMode="ALWAYS")
      p.addProperty("angularVelocity", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
      p.addProperty("torque",          "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

      self.shapes          = shapes

   def generate(self, path):
      context = dict()
      context["shapes"]          = self.shapes

      generateFile(path, 'data/ShapeStorage.templ.h', context)
