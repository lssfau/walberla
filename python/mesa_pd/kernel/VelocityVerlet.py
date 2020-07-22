# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class VelocityVerlet:
    def __init__(self, integrate_rotation=True):
        self.context = {'bIntegrateRotation': integrate_rotation, 'interface': []}
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="gs"))
        self.context['interface'].append(create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gs"))
        self.context['interface'].append(create_access("invMass", "walberla::real_t", access="g"))
        self.context['interface'].append(create_access("force", "walberla::mesa_pd::Vec3", access="gs"))
        self.context['interface'].append(create_access("oldForce", "walberla::mesa_pd::Vec3", access="gs"))

        if integrate_rotation:
            self.context['interface'].append(create_access("rotation", "walberla::mesa_pd::Rot3", access="gs"))
            self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gs"))
            self.context['interface'].append(create_access("invInertiaBF", "walberla::mesa_pd::Mat3", access="g"))
            self.context['interface'].append(create_access("torque", "walberla::mesa_pd::Vec3", access="gs"))
            self.context['interface'].append(create_access("oldTorque", "walberla::mesa_pd::Vec3", access="gs"))

        self.context['interface'].append(
            create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g"))

    def generate(self, module):
        ctx = {'module': module, **self.context}

        generate_file(module['module_path'], 'kernel/VelocityVerlet.templ.h', ctx)

        ctx["InterfaceTestName"] = "VelocityVerletInterfaceCheck"
        ctx["KernelInclude"] = "kernel/VelocityVerlet.h"
        ctx["ExplicitInstantiation"] = \
            "template void kernel::VelocityVerletPreForceUpdate::operator()(" \
            "const size_t p_idx1, " \
            "Accessor& ac) const;\n" + \
            "template void kernel::VelocityVerletPostForceUpdate::operator()(" \
            "const size_t p_idx1, " \
            "Accessor& ac) const;"
        generate_file(module['test_path'], 'tests/CheckInterface.templ.cpp', ctx,
                      'kernel/interfaces/VelocityVerletInterfaceCheck.cpp')
