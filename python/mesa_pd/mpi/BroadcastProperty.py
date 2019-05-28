# -*- coding: utf-8 -*-

from ..utility import generateFile

class BroadcastProperty:
   def generate(self, path):
      generateFile(path, 'mpi/BroadcastProperty.templ.h')
