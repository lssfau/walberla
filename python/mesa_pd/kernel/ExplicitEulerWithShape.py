# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class ExplicitEulerWithShape:
    def __init__(self):
        self.context = dict()
        self.context['interface'] = [create_access("position", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("invMass", "walberla::real_t", access="g"),
                                     create_access("force", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("rotation", "walberla::mesa_pd::Rot3", access="gs"),
                                     create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("invInertiaBF", "walberla::mesa_pd::Mat3", access="g"),
                                     create_access("torque", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g")]

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/ExplicitEulerWithShape.templ.h', ctx)

        ctx["InterfaceTestName"] = "ExplicitEulerWithShapeInterfaceCheck"
        ctx["KernelInclude"] = "kernel/ExplicitEulerWithShape.h"
        ctx["ExplicitInstantiation"] = "template void kernel::ExplicitEulerWithShape::operator()(const size_t p_idx1, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/ExplicitEulerWithShapeInterfaceCheck.cpp')
