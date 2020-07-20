# -*- coding: utf-8 -*-

from ..utility import TerminalColor, find, generate_file


def create_contact_storage_property(name, type, defValue=""):
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


class ContactStorage():
    def __init__(self):
        self.context = {'includes': [], 'properties': []}
        self.add_property("uid", "walberla::id_t", defValue="walberla::id_t(-1)")

    def add_property(self, name, type, defValue=""):
        prop = find(lambda x: x['name'] == name, self.context['properties'])
        if (prop is None):
            # print(f"{TerminalColor.GREEN} creating property: {name} {TerminalColor.DEFAULT}")
            self.context['properties'].append(create_contact_storage_property(name, type, defValue=defValue))
        else:
            if not (prop['type'] == type and prop['name'] == name and prop['defValue'] == defValue):
                new_prop = create_contact_storage_property(name, type, defValue=defValue)
                raise RuntimeError(
                    f"{TerminalColor.RED} property definition differs from previous one:\n"
                    f"PREVIOUS {prop}\n"
                    f"NEW {new_prop} {TerminalColor.DEFAULT}")
            print(f"{TerminalColor.YELLOW} reusing property: {name} {TerminalColor.DEFAULT}")

    def add_include(self, include):
        if (include in self.context['includes']):
            print(f"{TerminalColor.YELLOW} reusing include: {include} {TerminalColor.DEFAULT}")
        else:
            # print(f"{TerminalColor.GREEN} creating include: {include} {TerminalColor.DEFAULT}")
            self.context['includes'].append(include)

    def print(self):
        print("=" * 90)
        print("Creating ContactStorage Datastructure:")
        print("")
        print("{0: <20}{1: <30}{2: <20}{3: <10}".format("Type", "Name", "Def. Value", "SyncMode"))
        print("=" * 90)
        for prop in self.properties:
            print("{0: <20.19}{1: <30.29}{2: <20.19}{3: <10.9}".format(prop.type, prop.name, prop.defValue,
                                                                       prop.syncMode))
        print("=" * 90)

    def generate(self, module):
        ctx = {'module': module}
        ctx.update(self.context)

        generate_file(module['module_path'], 'data/ContactStorage.templ.h', ctx, filename='data/ContactStorage.h')
        generate_file(module['module_path'], 'data/ContactAccessor.templ.h', ctx, filename='data/ContactAccessor.h')
