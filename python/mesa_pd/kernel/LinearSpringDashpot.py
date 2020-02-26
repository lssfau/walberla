# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


class LinearSpringDashpot:
    def __init__(self):
        self.context = {'interface': []}
        self.context['interface'].append(create_access("uid", "walberla::id_t", access="g"))
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
        ctx["parameters"] = ["stiffnessN", "stiffnessT", "dampingN", "dampingT", "frictionCoefficientStatic",
                             "frictionCoefficientDynamic"]

        generate_file(module['module_path'], 'kernel/LinearSpringDashpot.templ.h', ctx)
