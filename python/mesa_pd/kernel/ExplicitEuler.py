# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class ExplicitEuler:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("position",        "walberla::mesa_pd::Vec3", access="gs")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3", access="gs")
      self.accessor.require("invMass",         "walberla::real_t",        access="g" )
      self.accessor.require("force",           "walberla::mesa_pd::Vec3", access="gs" )
      self.accessor.require("flags",           "walberla::mesa_pd::data::particle_flags::FlagT", access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["interface"] = self.accessor.properties
      generateFile(path, 'kernel/ExplicitEuler.templ.h', context)

      context["InterfaceTestName"] = "ExplicitEulerInterfaceCheck"
      context["KernelInclude"] = "kernel/ExplicitEuler.h"
      context["ExplicitInstantiation"] = "template void kernel::ExplicitEuler::operator()(const size_t p_idx1, Accessor& ac) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/ExplicitEulerInterfaceCheck.cpp')
