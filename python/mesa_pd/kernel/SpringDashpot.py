# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class SpringDashpot:
    def __init__(self):
        self.context = {'interface': []}
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("force", "walberla::mesa_pd::Vec3", access="r"))
        self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("torque", "walberla::mesa_pd::Vec3", access="r"))
        self.context['interface'].append(create_access("type", "uint_t", access="g"))
        self.context['interface'].append(
            create_access("contactHistory", "std::map<walberla::id_t, walberla::mesa_pd::Vec3>", access="gs"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        ctx["parameters"] = ["stiffness", "dampingN", "dampingT", "friction"]

        generate_file(module['module_path'], 'kernel/SpringDashpot.templ.h', ctx)

        ctx["InterfaceTestName"] = "SpringDashpotInterfaceCheck"
        ctx["KernelInclude"] = "kernel/SpringDashpot.h"
        ctx["ExplicitInstantiation"] = \
            "template void kernel::SpringDashpot::operator()(" \
            "const size_t p_idx1, " \
            "const size_t p_idx2, " \
            "Accessor& ac, " \
            "const Vec3& contactPoint, " \
            "const Vec3& contactNormal, " \
            "const real_t& penetrationDepth) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/SpringDashpotInterfaceCheck.cpp')
