# -*- coding: utf-8 -*-

from ..utility import generateFile

class ReduceProperty:
   def generate(self, path):
      generateFile(path, 'mpi/ReduceProperty.templ.h')
