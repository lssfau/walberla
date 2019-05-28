# -*- coding: utf-8 -*-

import numpy as np
from ..Container import Container
from ..utility import generateFile

class ContactHistory(Container):
   def __init__(self):
      super().__init__()

   def generate(self, path):
      self.unrollDimension()

      print("="*90)
      print("Creating ContactHistory Datastructure:")
      print("")
      print("{0: <20}{1: <30}{2: <20}{3: <10}".format("Type", "Name", "Def. Value", "SyncMode"))
      print("="*90)
      for prop in self.properties:
         print("{0: <20.19}{1: <30.29}{2: <20.19}{3: <10.9}".format(prop.type, prop.name, prop.defValue, prop.syncMode))
      print("="*90)

      context = dict()
      context["includes"]    = self.includes
      context["properties"]  = self.properties

      generateFile(path, 'data/ContactHistory.templ.h', context, filename='data/ContactHistory.h')
      generateFile(path, 'mpi/notifications/ContactHistoryNotification.templ.h', context)
