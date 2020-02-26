# -*- coding: utf-8 -*-

from mesa_pd.accessor import create_access
from mesa_pd.utility import generate_file


def create_property(name, type, defValue=""):
    """
    Parameters
    ----------
    name : str
       name of the property
    type : str
       type of the property
    defValue : str
       default value the property should be initialized with
    """

    return {'name': name, 'type': type, 'defValue': defValue}


class IntegrateParticlesHCSITS():
    def __init__(self):
        self.context = {'properties': [], 'interface': []}

        self.context['properties'].append(create_property("speedLimiterActive", "bool", defValue="false"))
        self.context['properties'].append(create_property("speedLimitFactor", "real_t", defValue="real_t(1.0)"))

        self.context['interface'].append(create_access("uid", "walberla::id_t", access="g"))
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("rotation", "walberla::mesa_pd::Rot3", access="gr"))
        self.context['interface'].append(create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("dv", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("dw", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(
            create_access("flags", "walberla::mesa_pd::data::particle_flags::FlagT", access="g"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/IntegrateParticlesHCSITS.templ.h', ctx)
