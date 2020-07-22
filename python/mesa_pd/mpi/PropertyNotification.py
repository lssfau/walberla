# -*- coding: utf-8 -*-

from ..utility import find, generate_file, TerminalColor


class PropertyNotification:
    def __init__(self, name):
        self.context = {
            'name': name,
            'properties': []
        }

    def add_property(self, name, type, reset_value):
        prop = find(lambda x: x['name'] == name, self.context['properties'])
        if (prop is None):
            self.context['properties'].append({'name': name, 'type': type, 'resetValue': reset_value})
        else:
            if not (prop['type'] == type and prop['name'] == name and prop['resetValue'] == reset_value):
                new_prop = {'name': name, 'type': type, 'resetValue': reset_value}
                raise RuntimeError(
                    f"{TerminalColor.RED} property definition differs from previous one:\n"
                    f"PREVIOUS {prop}\n"
                    f"NEW {new_prop} {TerminalColor.DEFAULT}")
            print(f"{TerminalColor.YELLOW} reusing property: {name} {TerminalColor.DEFAULT}")

    def generate(self, module):
        ctx = {'module': module, **self.context}
        generate_file(module['module_path'],
                      'mpi/notifications/PropertyNotification.templ.h',
                      ctx,
                      f'mpi/notifications/{self.context["name"]}.h')
