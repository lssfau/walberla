# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class HeatConduction:
    def __init__(self):
        self.context = {
            'interface': [
                create_access("temperature", "walberla::real_t", access="g"),
                create_access("heatFlux", "walberla::real_t", access="gsr"),
                create_access("type", "uint_t", access="g")
            ]
        }

    def generate(self, module):
        ctx = {'module': module,
               **self.context,
               "parameters": ["conductance"]}

        generate_file(module['module_path'], 'kernel/HeatConduction.templ.h', ctx)

        ctx["InterfaceTestName"] = "HeatConductionInterfaceCheck"
        ctx["KernelInclude"] = "kernel/HeatConduction.h"
        ctx["ExplicitInstantiation"] = \
            "template void kernel::HeatConduction::operator()(" \
            "const size_t p_idx1, " \
            "const size_t p_idx2, " \
            "Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/HeatConductionInterfaceCheck.cpp')
