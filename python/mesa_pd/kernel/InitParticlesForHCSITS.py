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


class InitParticlesForHCSITS:
    def __init__(self):
        self.context = {'properties': [], 'interface': []}

        self.context['properties'].append(
            create_property("globalAcceleration", "walberla::mesa_pd::Vec3", defValue="0"))

        self.context['interface'].append(create_access("uid", "walberla::id_t", access="g"))
        self.context['interface'].append(create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("invMass", "walberla::real_t", access="g"))
        self.context['interface'].append(create_access("inertia", "walberla::mesa_pd::Mat3", access="g"))
        self.context['interface'].append(create_access("invInertia", "walberla::mesa_pd::Mat3", access="g"))
        self.context['interface'].append(create_access("dv", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("dw", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("torque", "walberla::mesa_pd::Vec3", access="r"))
        self.context['interface'].append(create_access("force", "walberla::mesa_pd::Vec3", access="r"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/InitParticlesForHCSITS.templ.h', ctx)
