import argparse

#!/usr/bin/env python3

import os
import re


def create_cmake_lists (folder, module_name, top_level=True, deps = []):
    header_and_source = []
    subfolder = []
    newline = '\n    '
    add_sub = '\nadd_subdirectory( '
    for entry in os.scandir(folder):
        if entry.name == 'all.h':
            continue
        if entry.is_file():
            if entry.name.endswith('.h') or entry.name.endswith('.cpp'):
                header_and_source += [entry.name]
        else:
            if entry.name != 'doc' and entry.name != 'CMakeFiles':
                subfolder += [entry.name]
            create_cmake_lists(folder + '/' + entry.name, module_name, False)

    print(subfolder)
    if not header_and_source:
        return

    if top_level:
        with open(folder + '/CMakeListsRefactor.txt', 'w') as f:
            content = f"""add_library( {module_name} )
target_link_libraries( {module_name} PUBLIC {' '.join(x for x in deps)} )
target_sources( {module_name}
    PRIVATE
    {newline.join(x for x in header_and_source)}     
    )

add_subdirectory( {add_sub.join(x + ' )' for x in subfolder)}
"""
            f.write(content)
    else:
        with open(folder + '/CMakeLists.txt', 'w') as f:
            content = f"""target_sources( {module_name}
    PRIVATE
    {newline.join(x for x in header_and_source)}     
    )
"""
            f.write(content)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refactor CMakeLists.txt file')
    parser.add_argument('folder', type=str, help='Folder to be refactored; a CMakeListsRefactor.txt file will be created in each subfolder')
    args = parser.parse_args()

    print(args.folder)
    with open(args.folder + '/CMakeLists.txt', 'r') as f:
        depsRaw = re.findall(r'DEPENDS (.*)\)', f.read(), re.DOTALL)
    deps = []
    if depsRaw:
        deps = depsRaw[0].split()

    create_cmake_lists( args.folder, args.folder.split('/')[-1], True, deps )

