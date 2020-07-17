#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def is_walberla_root(walberla_root):
    return (walberla_root / 'src' / 'walberla.h').exists()


def trace_dependencies(modules, base_modules):
    deps = base_modules.copy()
    for mod in base_modules:
        deps |= modules.get(mod, set())
    if deps == base_modules:
        return deps
    else:
        return trace_dependencies(modules, deps)


def get_module_dependencies(module_dir):
    with open(module_dir / 'CMakeLists.txt', 'r') as fin:
        m = re.search(r'DEPENDS([\w\s]*)', fin.read())
        if m is not None:
            stripped = {x.strip() for x in m.group(1).split(' ') if x != ''}
            return (module_dir.name, {x for x in stripped if (x not in ['', 'BUILD_ONLY_IF_FOUND'])})
    return (module_dir.name, set())


def get_dependency_graph(walberla_root):
    modules_dir = walberla_root / 'src'
    modules = (get_module_dependencies(x) for x in modules_dir.iterdir() if x.is_dir())
    return {i: v for i, v in modules}


def color_dependencies(fout, dependencies, base_module):
    for dep in trace_dependencies(dependencies, {base_module}):
        if (dep == base_module):
            fout.write(f'  {dep}[fillcolor=limegreen, style="filled"];\n')
            continue
        if (dep in dependencies[base_module]):
            fout.write(f'  {dep}[fillcolor=gold, style="filled"];\n')
            continue
        fout.write(f'  {dep}[fillcolor=brown1, style="filled"];\n')


def write_dependency_graph(fout, dependencies):
    for module, deps in dependencies.items():
        for dep in deps:
            fout.write(f'  {dep} -> {module};\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates a dot file containing the waLBerla dependency graph.')
    parser.add_argument('-d', default='..', help='waLBerla root directory', type=Path)
    parser.add_argument('-o', default='module_dependencies.dot', help='output dot file')
    parser.add_argument('--trace', default=None, help='Colors all direct and indirect dependencies of TRACE.')
    args = parser.parse_args()

    if not is_walberla_root(args.d):
        raise ValueError(f'{args.d} is not a waLBerla root directory!')

    dependencies = get_dependency_graph(args.d)

    with open(args.o, 'w') as fout:
        fout.write('digraph module_dependencies {\n')
        write_dependency_graph(fout, dependencies)
        if args.trace is not None:
            color_dependencies(fout, dependencies, args.trace)
        fout.write('}\n')
