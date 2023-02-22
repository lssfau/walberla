waLBerla_link_files_to_builddir( *.prm )

waLBerla_add_executable(NAME    DeformationField
                        FILES   DeformationField.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    SingleVortex
                        FILES   SingleVortex.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)

waLBerla_add_executable(NAME    ZalesakDisk
                        FILES   ZalesakDisk.cpp
                        DEPENDS blockforest boundary core domain_decomposition field lbm postprocessing timeloop vtk)