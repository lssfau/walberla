# -*- coding: utf-8 -*-

from pathlib import Path
import shutil
import os


class Module:
    def __init__(self, path, module_name='mesa_pd'):
        """Propery of a data strcuture

        Parameters
        ----------
        path : str
           path to the root waLBerla folder
        module_name : str
           name of the generated module
        """

        self.context = {}
        self.context['path'] = Path(path).resolve()
        self.context['name'] = module_name
        self.context['module_path'] = self.context['path'] / 'src' / self.context['name']
        self.context['test_path']   = self.context['path'] / 'tests' / self.context['name']

        self.components = []

    def add(self, component):
        self.components.append(component)
        return component

    def rename(self):
        for root, dirnames, filenames in os.walk(self.context['module_path']):
            for filename in filenames:
                filedata = None
                print(f'{root}/{filename}')
                with open(f'{root}/{filename}', 'r') as file:
                    filedata = file.read()

                filedata = filedata.replace('mesa_pd', self.context['name'])

                with open(f'{root}/{filename}', 'w') as file:
                    file.write(filedata)

    def generate(self, folder_check=True):
        if (folder_check):
            print(f"This operation will overwrite the content of: {self.context['module_path']}")
            answer = input("Continue? (y to confirm)")
            if (answer != "y"):
                return

        mesa_pd_folder = (Path(__file__).parents[2] / 'src' / 'mesa_pd').resolve()
        if (mesa_pd_folder != self.context['module_path']):
            shutil.rmtree(self.context['module_path'])
            shutil.copytree(mesa_pd_folder, self.context['module_path'])

        for d in self.components:
            d.generate(self.context)

        self.rename()