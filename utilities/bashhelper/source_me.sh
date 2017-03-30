###########################################################################################
#############  waLBerla Bash shortcuts                         ############################
###########################################################################################
#
# This script sets up your local environment for waLBerla developing
#
# Features:
#    - shortcuts for switching between source and build directory
#    - cmake commands that automatically set up ccache
#    - commands to open web tools (buildbot, trac)
#
# Setup:
#   Insert the following into your .bashrc file and adapt your paths
#   export WALBERLA_SOURCE_DIR=/path/to/walberla/sources
#   export WALBERLA_BUILD_DIR=/path/to/build/directory
#   export WALBERLA_CMAKE_DEFS="-DCMAKE_BUILD_TYPE=Debug -DWALBERLA_ENABLE_GUI=ON"
#   source $WALBERLA_SOURCE_DIR/utilities/bashhelper/source_me.sh
#  
# Shortcuts: ( all cd* commands take as argument a relative path )
#   cdb   - cd to build directory 
#   cdbt  - cd to build/tests
#   cds   - cd to source directory
#   cdss  - cd to sourceDir/src
#   cdst  - cd to sourceDir/tests
#   sbs   - switch to corresponding path in build and source directory
#
#
# Fast Access to web tools:
#   open-buildbot
#   open-trac 
#
# CMake and ccache:  go to your (empty) build directory and then type
#   cmake-gcc or cmake-clang to 
#
###########################################################################################



WSD=$WALBERLA_SOURCE_DIR
WBD=$WALBERLA_BUILD_DIR
CMAKE_DEF=$WALBERLA_CMAKE_DEFS

alias cmake-gcc="CC=\"ccache gcc\" CXX=\"ccache g++\" cmake $CMAKE_DEF $WSD"
alias cmake-intel="CC=\"ccache icc\" CXX=\"ccache icpc\" cmake $CMAKE_DEF $WSD"
alias cmake-clang='CC="ccache clang -Qunused-arguments -fcolor-diagnostics" CXX="ccache clang++ -Qunused-arguments -fcolor-diagnostics" cmake $CMAKE_DEF $WSD'
alias cmake-clang-nocolor='CC="ccache clang -Qunused-arguments" CXX="ccache clang++ -Qunused-arguments" cmake $CMAKE_DEF $WSD'
alias sbs='cd `$WSD/utilities/bashhelper/correspondingDirs.py $WSD $WBD`'

export COMPLETE_SCRIPT="$WSD/utilities/bashhelper/folderComplete.py"


# Python modules for gdb pretty printer

export PYTHONPATH=$PYTHONPATH:$WSD/python
export PATH=$PATH:$WSD/python/waLBerla/source_checker

# Web
alias open-trac='firefox https://www10.informatik.uni-erlangen.de/trac/walberla &'
alias open-buildbot='firefox https://www10.informatik.uni-erlangen.de/buildbot/waterfall &'



# -------------------------------- Shortcuts -------------------------------------------

# Build directory
cdb() {
  cd $WBD/$1
}
_cdb() {
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=(`$COMPLETE_SCRIPT $WBD $cur`)
}
complete -o filenames  -o nospace -F _cdb cdb


# Build test directory
cdbt() {
  cd $WBD/tests/$1
}
_cdbt() {
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=(`$COMPLETE_SCRIPT $WBD/tests $cur`)
}
complete -o filenames  -o nospace -F _cdbt cdbt


# Source dir
cds() {
  cd $WSD/$1
}
_cds() {
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=(`$COMPLETE_SCRIPT $WSD $cur`)
}
complete -o filenames  -o nospace -F _cds cds


# Source-test dir
cdst() {
  cd $WSD/tests/$1
}
_cdst() {
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=(`$COMPLETE_SCRIPT $WSD/tests $cur`)
}
complete -o filenames  -o nospace -F _cdst cdst


# Source-source dir
cdss() {
  cd $WSD/src/$1
}
_cdss() {
  cur="${COMP_WORDS[COMP_CWORD]}"
  COMPREPLY=(`$COMPLETE_SCRIPT $WSD/src $cur`)
}
complete -o filenames  -o nospace -F _cdss cdss

