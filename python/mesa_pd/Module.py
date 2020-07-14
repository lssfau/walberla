# -*- coding: utf-8 -*-

from pathlib import Path
import shutil


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

        path = Path(path).resolve()
        self.context = {
            'path': path,
            'name': module_name,
            'module_path': path / 'src' / module_name,
            'test_path': path / 'tests' / module_name,
            'enableOpenMP': False
        }

        self.components = []

    def add(self, component):
        self.components.append(component)
        return component

    def enable_openmp(self, enabled):
        self.context['enableOpenMP'] = enabled

    def rename(self):
        for filename in (f for f in self.context['module_path'].glob('**/*') if f.is_file()):
            filedata = None
            # print(f'renaming module name: {filename}')
            with open(filename, encoding="utf-8") as fin:
                filedata = fin.read()

            filedata = filedata.replace('mesa_pd', self.context['name'])

            with open(filename, 'w', encoding="utf-8") as fout:
                fout.write(filedata)

    def generate(self, folder_check=True):
        print(f"This operation will overwrite the content of: {self.context['module_path']}")
        if (folder_check):
            answer = input("Continue? (y to confirm)")
            if (answer != "y"):
                return

        mesa_pd_folder = (Path(__file__).parents[2] / 'src' / 'mesa_pd').resolve()
        if (mesa_pd_folder != self.context['module_path']):
            if not self.context['module_path'].exists():
                self.context['module_path'].mkdir(parents=True)
            shutil.rmtree(self.context['module_path'])
            shutil.copytree(mesa_pd_folder, self.context['module_path'])

        for d in self.components:
            d.generate(self.context)

        self.rename()
