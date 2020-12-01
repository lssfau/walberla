# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class SemiImplicitEuler:
    def __init__(self, integrate_rotation=True):
        self.context = {'bIntegrateRotation': integrate_rotation,
                        'interface': [
                            create_access("position", "walberla::mesa_pd::Vec3", access="gs"),
                            create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gs"),
                            create_access("invMass", "walberla::real_t", access="g"),
                            create_access("force", "walberla::mesa_pd::Vec3", access="gs"),
                            create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g")
                        ]
                        }

        if integrate_rotation:
            self.context['interface'].append(create_access("rotation", "walberla::mesa_pd::Rot3", access="gs"))
            self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gs"))
            self.context['interface'].append(create_access("invInertiaBF", "walberla::mesa_pd::Mat3", access="g"))
            self.context['interface'].append(create_access("torque", "walberla::mesa_pd::Vec3", access="gs"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/SemiImplicitEuler.templ.h', ctx)

        ctx["InterfaceTestName"] = "SemiImplicitEulerInterfaceCheck"
        ctx["KernelInclude"] = "kernel/SemiImplicitEuler.h"
        ctx["SemiImplicitEulerInstantiation"] = \
            "template void kernel::SemiImplicitEuler::operator()(const size_t p_idx1, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/SemiImplicitEulerInterfaceCheck.cpp')
