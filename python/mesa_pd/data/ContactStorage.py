# -*- coding: utf-8 -*-

import numpy as np
from ..Container import Container
from ..utility import generateFile

class ContactStorage(Container):
   def __init__(self):
      super().__init__()
      self.addProperty("uid",                "walberla::id_t",           defValue = "walberla::id_t(-1)",   syncMode="NEVER")

   def generate(self, path):
      self.unrollDimension()

      print("="*90)
      print("Creating ContactStorage Datastructure:")
      print("")
      print("{0: <20}{1: <30}{2: <20}{3: <10}".format("Type", "Name", "Def. Value", "SyncMode"))
      print("="*90)
      for prop in self.properties:
         print("{0: <20.19}{1: <30.29}{2: <20.19}{3: <10.9}".format(prop.type, prop.name, prop.defValue, prop.syncMode))
      print("="*90)

      context = dict()
      context["includes"]    = self.includes
      context["properties"]  = self.properties

      generateFile(path, 'data/ContactStorage.templ.h', context, filename='data/ContactStorage.h')
      generateFile(path, 'data/ContactAccessor.templ.h', context, filename='data/ContactAccessor.h')
