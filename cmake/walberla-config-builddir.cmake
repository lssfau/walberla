# This file is copied (and configured) to the build directory root.
# Also a export(PACKAGE waLBerla ) is done to register the package at the
# cmake registry.
# 
# When a different project needs waLBerla, it calls find_package (waLBerla) 
# which first locates the build folder, using the cmake package registry, 
# and the automatically executes this file. 


# Set the source and binary dir
set ( walberla_SOURCE_DIR @walberla_SOURCE_DIR@ )
set ( walberla_BINARY_DIR @walberla_BINARY_DIR@ )
