# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class ExplicitEuler:
    def __init__(self):
        self.context = dict()
        self.context['interface'] = [create_access("position", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("invMass", "walberla::real_t", access="g"),
                                     create_access("force", "walberla::mesa_pd::Vec3", access="gs"),
                                     create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g")]

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/ExplicitEuler.templ.h', ctx)

        ctx["InterfaceTestName"] = "ExplicitEulerInterfaceCheck"
        ctx["KernelInclude"] = "kernel/ExplicitEuler.h"
        ctx["ExplicitInstantiation"] = "template void kernel::ExplicitEuler::operator()(const size_t p_idx1, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/ExplicitEulerInterfaceCheck.cpp')
