# -*- coding: utf-8 -*-

from .Property import Property, unrollDimension
from .utility import find, TerminalColor

class Container:
   def __init__(self):
      """Base class for a container which manages includes and properties
      """

      self.properties = []
      self.includes   = []

   def addProperty(self, name, type, access="grs", defValue="", syncMode = "ALWAYS", dim=1):
      prop = find( lambda x : x.name == name, self.properties )
      if (prop == None):
         #print(TerminalColor.GREEN + "creating particle property: {}".format(name) + TerminalColor.DEFAULT)
         self.properties.append( Property(name, type, access=access, defValue=defValue, syncMode=syncMode, dim=dim) )
      else:
         if not (prop.type == type and prop.name == name and prop.defValue == defValue and prop.dim == dim):
            raise RuntimeError(TerminalColor.RED + "property definition differs from previous one:\nPREVIOUS {}\nNEW {}".format(prop, Property(name, type, defValue=defValue, syncMode=syncMode, dim=dim)) + TerminalColor.DEFAULT)
         print(TerminalColor.YELLOW + "reusing particle property: {}".format(name) + TerminalColor.DEFAULT)

   def addInclude(self, include):
      if (include in self.includes):
         print(TerminalColor.YELLOW + "reusing particle include: {}".format(include) + TerminalColor.DEFAULT)
      else:
         #print(TerminalColor.GREEN + "creating particle include: {}".format(include) + TerminalColor.DEFAULT)
         self.includes.append(include)

   def unrollDimension(self):
      self.properties = unrollDimension(self.properties)
