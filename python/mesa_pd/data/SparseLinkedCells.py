# -*- coding: utf-8 -*-

import numpy as np
from ..utility import generateFile

class SparseLinkedCells:
   def generate(self, path):
      generateFile(path, 'data/SparseLinkedCells.templ.h')
