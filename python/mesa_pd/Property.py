# -*- coding: utf-8 -*-

from .utility import TerminalColor

class Property:
   def __init__(self, name, type, access="grs", defValue="", syncMode = "ALWAYS", dim=1):
      """Propery of a data strcuture

      Parameters
      ----------
      name : str
         name of the property
      type : str
         type of the property
      access : str
         'g' for getter (getName)
         'r' for reference (getNameRef)
         's' for setter (setName)
         any combination is possible
      defValue : str
         default value the property should be initialized with
      syncMode : str
         'NEVER', this property does not have to be synced
         'COPY', this property must be synced on creation
         'MIGRATION', this property must be synced when the ownership changes
         'ALWAYS', this property has to be synced in every iteration
      dim : int
         dimensions of the property
      """

      #sort access specifier and remove duplicates
      foo = "".join(sorted(access))
      access = ''.join([foo[i] for i in range(len(foo)-1) if foo[i+1]!= foo[i]]+[foo[-1]])

      for acc in access:
         if not (acc in ["g","s","r"]):
            raise RuntimeError("{} is not a valid access specifier in {}".format(acc, access))

      if (not syncMode in ["NEVER", "COPY", "MIGRATION", "ALWAYS"]):
         raise RuntimeError(TerminalColor.RED + "{} is no valid sync for property: {}".format(syncMode, name) + TerminalColor.DEFAULT)

      if (dim < 1):
         raise RuntimeError(TerminalColor.RED + "dimension has to be >=1: {}".format(dim) + TerminalColor.DEFAULT)

      self.name     = name
      self.type     = type
      self.access   = access
      self.defValue = defValue
      self.syncMode = syncMode
      self.dim      = dim

   def __str__(self):
      return "name: {}, type: {}, access: {}, defValue: {}, syncMode: {}, dim: {}".format(self.name, self.type, self.access, self.defValue, self.syncMode, self.dim)

   def getAccessFunction(self):
      """Returns a list of accessor function names
      """

      funcs = []
      if 'g' in self.access:
         funcs.append("get" + self.name.capitalize())
      if 's' in self.access:
         funcs.append("set" + self.name.capitalize())
      if 'r' in self.access:
         funcs.append("get" + self.name.capitalize() + "Ref")
      return funcs

def unrollDimension(props):
   """Unrolls all more dimensional properties into one dimensional properties

   Iterates over all elements. Copies all one dimensional properties.
   More dimensional properties get split into one dimensional properties with added suffix.

   Parameters
   ----------
   props : list
      list of properties to be unrolled

   Returns
   -------
      list of unrolled properties
   """

   unrolled = []
   for prop in props:
      if (prop.dim == 1):
         unrolled.append(Property(prop.name, prop.type, prop.access, prop.defValue, prop.syncMode, prop.dim))
      else:
         for d in range(prop.dim):
            unrolled.append(Property("{}{}".format(prop.name,d), prop.type, prop.access, prop.defValue, prop.syncMode, 1))
   return unrolled
