#!/usr/bin/env python3
from __future__ import print_function
# from sys import argv
from os import path
import argparse

usage = """
    Example usage:
    Assume you have a waLBerla app in ~/code/app and the waLBerla sources in ~/code/wlb
    Add the following lines to your .bashrc:

        createShortcuts.py ~/code/wlb wl-     ~/code/wlb  > ~/.waLBerla_shortcuts
        createShortcuts.py ~/code/app myapp-  ~/code/app --app --build ~/build/app  --walberla_source ~/code/wlb >> ~/.waLBerla_shortcuts
        source file_to_source

    Then the following aliases are available:
        wl-s     cd ~/code/wlb
        wl-ss    cd ~/code/wlb/src
        wl-st    cd ~/code/wlb/tests
        wl-sa    cd ~/code/wlb/apps

        app-s    cd ~/code/app
        app-sa   cd ~/code/app/apps
        app-b    cd ~/build/app
        app-ba   cd ~/build/app/apps

        app-make
        app-ccmake
        app-git
"""

parser = argparse.ArgumentParser()
parser.add_argument('prefix', help='Shortcut Prefix')
parser.add_argument('src', help='Source directory')
parser.add_argument('--build', help='Build directory')
parser.add_argument('--cmake_def', help='CMake Definitions')
parser.add_argument('--walberla_source', help='waLBerla source directory')
parser.add_argument('--app', help='Generate shortcut for apps', action="store_true")
parser.add_argument('--no_ccache', help='CMake macros do not use ccache', action="store_true")
parser.add_argument('--shell', help='Supported shells bash and fish', default="bash")
args = parser.parse_args()

if args.app:
    if not args.prefix:
        print('When creating shortcuts for a waLBerla app you have to define a prefix (--prefix)')
        exit(1)
    if not args.walberla_source:
        print('When creating shortcuts for a waLBerla app you have to define'
              'the path to waLBerla sources (--walberla_source)')
        exit(1)

# Assume that a "folderComplete.py" script is located in the same directory with this script
dir_of_current_script = path.dirname(path.realpath(__file__))
complete_script = path.join(dir_of_current_script, 'folderComplete.py')

cd_shortcut_rule = dict()

cd_shortcut_rule['bash'] = """
{short_command_name}() {{
  cd {dir}/$1
}}
_{short_command_name}() {{
  cur="${{COMP_WORDS[COMP_CWORD]}}"
  COMPREPLY=(`{complete_script} {dir} $cur`)
}}
complete -o filenames  -o nospace -F _{short_command_name} {short_command_name}

"""

cd_shortcut_rule['fish'] = """
function {short_command_name}
    cd {dir}
end
complete -f --command {short_command_name} --arguments '(__fish_complete_directories {dir} )'
"""


def create_cd_shortcut_rule(short_command_name, dir, shell):
    return cd_shortcut_rule[shell].format(complete_script=complete_script,
                                          short_command_name=short_command_name,
                                          dir=dir)


########################################################################################################################
# CMake Aliases
########################################################################################################################

build_aliases = dict()
build_aliases['bash'] = """
alias {prefix}cmake-gcc='CC="{ccache}gcc" CXX="{ccache}g++" cmake {definitions} {source_dir}'
alias {prefix}cmake-intel='CC="{ccache}icc" CXX="{ccache}icpc" cmake {definitions} {source_dir}'
alias {prefix}cmake-clang='CC="{ccache}clang -Qunused-arguments -fcolor-diagnostics" CXX="{ccache}clang++ -Qunused-arguments -fcolor-diagnostics" cmake {definitions} {source_dir}'
alias {prefix}cmake-clang-nocolor='CC="{ccache}clang -Qunused-arguments" CXX="{ccache}clang++ -Qunused-arguments" cmake {definitions} {source_dir}'

alias {prefix}cmake='cmake {build_dir}'
alias {prefix}ccmake='ccmake {build_dir}'
alias {prefix}cmake-gui='cmake-gui {build_dir}'
alias {prefix}make='make -f {build_dir}/apps/Makefile'
"""

build_aliases['fish'] = """
function {prefix}cmake-gcc
    set CC  {ccache}gcc
    set CXX {ccache}g++
    cmake {definitions} {source_dir}
end
function {prefix}cmake-intel
    set CC  {ccache}icc
    set CXX {ccache}icpc
    cmake {definitions} {source_dir}
end
function {prefix}cmake-clang
    set CC  {ccache}clang -Qunused-arguments -fcolor-diagnostics
    set CXX {ccache}clang++ -Qunused-arguments -fcolor-diagnostics
    cmake {definitions} {source_dir}
end
function {prefix}cmake-clang-nocolor
    set CC  {ccache}clang -Qunused-arguments
    set CXX {ccache}clang++ -Qunused-arguments
    cmake {definitions} {source_dir}
end

alias {prefix}cmake='cmake {build_dir}'
alias {prefix}ccmake='ccmake {build_dir}'
alias {prefix}cmake-gui='cmake-gui {build_dir}'
alias {prefix}make='make -f {build_dir}/apps/Makefile'
"""


def create_build_aliases(prefix, source_dir, build_dir, definitions, shell):
    if args.no_ccache:
        ccache = ""
    else:
        ccache = "ccache "
    return build_aliases[shell].format(prefix=prefix,
                                       source_dir=source_dir,
                                       build_dir=build_dir,
                                       definitions=definitions,
                                       ccache=ccache)


source_aliases = """
alias {prefix}git='git --git-dir {source_dir}/.git --work-tree {source_dir}'
alias {prefix}git-buildscript='git --git-dir {source_dir}/.git --work-tree {source_dir} push origin master:personal/$USER/buildscript -f'
"""


def create_source_aliases(prefix, source_dir):
    return source_aliases.format(prefix=prefix, source_dir=source_dir)


cmake_defs = args.cmake_def
if args.app:
    cmake_defs += " -DWALBERLA_DIR=" + args.walberla_source

print(create_cd_shortcut_rule(args.prefix + "s", args.src, args.shell))
print(create_cd_shortcut_rule(args.prefix + "ss", path.join(args.src, "src"), args.shell))
print(create_cd_shortcut_rule(args.prefix + "st", path.join(args.src, "tests"), args.shell))
print(create_cd_shortcut_rule(args.prefix + "sa", path.join(args.src, "apps"), args.shell))
print(create_source_aliases(args.prefix, args.src))

if args.build:
    print(create_build_aliases(args.prefix, args.src, args.build, cmake_defs, args.shell))
    print(create_cd_shortcut_rule(args.prefix + "b", args.build, args.shell))
    print(create_cd_shortcut_rule(args.prefix + "bs", path.join(args.build, "src"), args.shell))
    print(create_cd_shortcut_rule(args.prefix + "bt", path.join(args.build, "tests"), args.shell))
    print(create_cd_shortcut_rule(args.prefix + "ba", path.join(args.build, "apps"), args.shell))
