# -*- coding: utf-8 -*-

from mesa_pd.utility import generate_file


class SpringDashpotSpring:
    def generate(self, module):
        ctx = {'module': module}
        ctx["parameters"] = ["stiffnessN", "dampingN", "stiffnessT", "coefficientOfFriction"]

        generate_file(module['module_path'], 'kernel/SpringDashpotSpring.templ.h', ctx)
