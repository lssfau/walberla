add_library( lbm_generated)

target_link_libraries( lbm_generated
        PUBLIC
        blockforest
        boundary
        communication
        core
        domain_decomposition
        field
        geometry
        gui
        stencil
        timeloop
        vtk
        )

add_subdirectory( boundary )
add_subdirectory( communication )
add_subdirectory( gpu )
add_subdirectory( evaluation )
add_subdirectory( field )
add_subdirectory( refinement )
add_subdirectory( storage_specification )
add_subdirectory( sweep_collection )