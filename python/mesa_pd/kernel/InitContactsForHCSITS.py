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


class InitContactsForHCSITS:
    def __init__(self):
        self.context = {'properties': [], 'interface': []}

        self.context['properties'].append(create_property("erp", "real_t", defValue="real_t(0.8)"))
        self.context['properties'].append(create_property("maximumPenetration", "real_t", defValue="0"))

        self.context['interface'].append(create_access("uid", "walberla::id_t", access="g"))
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("invInertia", "walberla::mesa_pd::Mat3", access="g"))
        self.context['interface'].append(create_access("invMass", "walberla::real_t", access="g"))

    def getRequirements(self):
        return self.paccessor

    def generate(self, module):
        ctx = {'module': module, **self.context}
        ctx["material_parameters"] = ["friction"]
        generate_file(module['module_path'], 'kernel/InitContactsForHCSITS.templ.h', ctx)
