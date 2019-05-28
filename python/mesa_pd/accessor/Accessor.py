# -*- coding: utf-8 -*-

from ..Property import Property
from ..utility import find, TerminalColor, generateFile

class Accessor:
   def __init__(self):
      self.properties = []

   def require(self, name, type, access):
      """requires that a certain property is accessible

      Parameters
      ----------
      name : str
         name of the property requested
      type : str
         type of the requested property
      access : str
         'g' for getter (getName)
         'r' for reference (getNameRef)
         's' for setter (setName)
         any combination is possible

      Example
      -------
      require("position", "walberla::mesa_pd::Vec3", "sg")
      """

      prop = find( lambda x : x.name == name, self.properties )
      if (prop == None):
         #print(TerminalColor.GREEN + "[{}] creating particle property: {}".format(info, name) + TerminalColor.DEFAULT)
         self.properties.append( Property(name, type, access = access) )
      else:
         if not (prop.type == type):
            raise RuntimeError(TerminalColor.RED + "requirement definition differs from previous one:\n{} {}\n{} {}".format(name, type, prop.name, prop.type) + TerminalColor.DEFAULT)
         foo = "".join(sorted(access + prop.access))
         prop.access = ''.join([foo[i] for i in range(len(foo)-1) if foo[i+1]!= foo[i]]+[foo[-1]])


   def mergeRequirements(self, accessor):
      for req in accessor.properties:
         self.require( req.name, req.type, req.access )

   def printSummary(self):
      print("="*90)
      print("Requirements for Accessor:")
      print("")
      print("{0: <30}{1: <30}{2: <30}".format("Name", "Type", "Access"))
      print("="*90)
      for gs in self.properties:
         print("{0: <30.29}{1: <30.29}{2: <30.29}".format(gs.name, gs.type, gs.access))
      print("="*90)
