# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class PFCDamping:
    def __init__(self, integrate_rotation=True):
        self.context = {'bIntegrateRotation': integrate_rotation,
                        'interface': [
                            create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gs"),
                            create_access("force", "walberla::mesa_pd::Vec3", access="gs"),
                            create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g")
                        ]
                        }

        if integrate_rotation:
            self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gs"))
            self.context['interface'].append(create_access("torque", "walberla::mesa_pd::Vec3", access="gs"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/PFCDamping.templ.h', ctx)

        ctx["InterfaceTestName"] = "PFCDampingInterfaceCheck"
        ctx["KernelInclude"] = "kernel/PFCDamping.h"
        ctx["ExplicitInstantiation"] = \
            "template void kernel::PFCDamping::operator()(const size_t p_idx1, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/PFCDampingInterfaceCheck.cpp')
