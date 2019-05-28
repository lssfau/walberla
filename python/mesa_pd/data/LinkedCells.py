# -*- coding: utf-8 -*-

import numpy as np
from ..utility import generateFile

class LinkedCells:
   def generate(self, path):
      generateFile(path, 'data/LinkedCells.templ.h')
