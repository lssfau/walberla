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


class HCSITSRelaxationStep():
    def __init__(self):
        self.context = {'properties': [], 'interface': []}

        self.context['properties'].append(create_property("maxSubIterations", "size_t", defValue="20"))
        self.context['properties'].append(
            create_property("relaxationModel", "RelaxationModel", defValue="InelasticFrictionlessContact"))
        self.context['properties'].append(create_property("deltaMax", "real_t", defValue="0"))
        self.context['properties'].append(create_property("cor", "real_t", defValue="real_t(0.2)"))

        self.context['interface'].append(create_access("uid", "walberla::id_t", access="g"))
        self.context['interface'].append(create_access("position", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("linearVelocity", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("angularVelocity", "walberla::mesa_pd::Vec3", access="g"))
        self.context['interface'].append(create_access("invMass", "walberla::real_t", access="g"))
        self.context['interface'].append(create_access("invInertia", "walberla::mesa_pd::Mat3", access="g"))
        self.context['interface'].append(create_access("dv", "walberla::mesa_pd::Vec3", access="gr"))
        self.context['interface'].append(create_access("dw", "walberla::mesa_pd::Vec3", access="gr"))

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'], 'kernel/HCSITSRelaxationStep.templ.h', ctx)
