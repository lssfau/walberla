# -*- coding: utf-8 -*-

from ..utility import generateFile

class ReduceContactHistory:
   def generate(self, path):
      generateFile(path, 'mpi/ReduceContactHistory.templ.h')
