# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class ForceLJ:
    def __init__(self):
        self.context = {
            'interface': [
                create_access("position", "walberla::mesa_pd::Vec3", access="g"),
                create_access("force", "walberla::mesa_pd::Vec3", access="r"),
                create_access("type", "uint_t", access="g")
            ]
        }

    def generate(self, module):
        ctx = {
            'module': module,
            **self.context,
            "parameters": ["epsilon", "sigma"]
        }

        generate_file(module['module_path'], 'kernel/ForceLJ.templ.h', ctx)

        ctx["InterfaceTestName"] = "ForceLJInterfaceCheck"
        ctx["KernelInclude"] = "kernel/ForceLJ.h"
        ctx["ExplicitInstantiation"] = \
            "template void kernel::ForceLJ::operator()(const size_t p_idx1, const size_t p_idx2, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/ForceLJInterfaceCheck.cpp')
