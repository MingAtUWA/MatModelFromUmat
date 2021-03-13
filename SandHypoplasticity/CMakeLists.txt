#===================================
#
# Sand Hypoplasticity Project
#
#===================================
add_library(
    SandHypoplasticity SHARED
    # 1
    sand_hypoplasticity_integration.f90
    #umat_hypoplasticity.for
    # 2
    # rewrite umat for efficency
    # sand_hypoplasticity_integration_dir.f90
    umat_hypoplasticity.f90
    )

target_link_libraries(
    SandHypoplasticity PUBLIC
    FortranUtils
    )
