# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
from mesa_pd.utility import generateFile

class ExplicitEulerWithShape:
   def __init__(self):
      self.accessor = Accessor()
      self.accessor.require("position",        "walberla::mesa_pd::Vec3", access="gs")
      self.accessor.require("linearVelocity",  "walberla::mesa_pd::Vec3", access="gs")
      self.accessor.require("invMass",         "walberla::real_t",        access="g" )
      self.accessor.require("force",           "walberla::mesa_pd::Vec3", access="gs" )
      self.accessor.require("rotation",        "walberla::mesa_pd::Rot3", access="gs")
      self.accessor.require("angularVelocity", "walberla::mesa_pd::Vec3", access="gs")
      self.accessor.require("invInertiaBF",    "walberla::mesa_pd::Mat3", access="g" )
      self.accessor.require("torque",          "walberla::mesa_pd::Vec3", access="gs" )
      self.accessor.require("flags",           "walberla::mesa_pd::data::particle_flags::FlagT", access="g")

   def getRequirements(self):
      return self.accessor

   def generate(self, path):
      context = dict()
      context["interface"] = self.accessor.properties
      generateFile(path, 'kernel/ExplicitEulerWithShape.templ.h', context)

      context["InterfaceTestName"] = "ExplicitEulerWithShapeInterfaceCheck"
      context["KernelInclude"] = "kernel/ExplicitEulerWithShape.h"
      context["ExplicitInstantiation"] = "template void kernel::ExplicitEulerWithShape::operator()(const size_t p_idx1, Accessor& ac) const;"
      generateFile(path, 'tests/CheckInterface.templ.cpp', context, '../../tests/mesa_pd/kernel/interfaces/ExplicitEulerWithShapeInterfaceCheck.cpp')
