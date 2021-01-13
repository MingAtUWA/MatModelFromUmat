#===================================
#
# Sand Hypoplasticity Project
#
#===================================
add_library(
    SandHypoplasticity SHARED
    sand_hypoplasticity_integration.f90
    umat_hypoplasticity.for
    )

target_link_libraries(
    SandHypoplasticity PUBLIC
    FortranUtils
    )
