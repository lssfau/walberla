# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class TemperatureIntegration:
    def __init__(self):
        self.context = {
            'interface': [
                create_access("temperature", "walberla::real_t", access="gs"),
                create_access("heatFlux", "walberla::real_t", access="gs"),
                create_access("invMass", "walberla::real_t", access="g"),
                create_access("type", "uint_t", access="g")
            ]
        }

    def generate(self, module):
        ctx = {
            'module': module,
            **self.context,
            "parameters": ["invSpecificHeat"]
        }
        generate_file(module['module_path'], 'kernel/TemperatureIntegration.templ.h', ctx)

        ctx["InterfaceTestName"] = "TemperatureIntegrationInterfaceCheck"
        ctx["KernelInclude"] = "kernel/TemperatureIntegration.h"
        ctx[
            "ExplicitInstantiation"] = \
            "template void kernel::TemperatureIntegration::operator()(const size_t p_idx1, Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/TemperatureIntegrationInterfaceCheck.cpp')
