! wrapper of umat_hypoplasticity() subroutines
subroutine sand_hypoplasticity_integration(stress, ddsdde, statev, dstran, props)
    !DEC$ ATTRIBUTES DLLEXPORT::sand_hypoplasticity_integration  
    implicit none

    ! variables not used
    ! ndi: number of direct stress components
    integer, parameter::ndi = 3
    ! nshr: number of shear stress components
    integer, parameter::nshr = 3
    ! ntens = ndi + nshr: size of stress or strain components
    integer, parameter::ntens = 6
    ! stress
    real(8), intent(inout)::stress(ntens)
    ! ddsdde: tangential stiffness matrix
    real(8), intent(out)::ddsdde(ntens, ntens)
    ! total strain increment
    real(8), intent(in)::dstran(ntens)
    
    ! nstatv: number of state variables
    integer, parameter::nstatv = 13
    ! statev: statev variables
    real(8), intent(inout)::statev(nstatv)
    
    ! nprop: number of material properties
    integer, parameter::nprops = 16
    ! props: material properties
    real(8), intent(in)::props(nprops)
    
    !cmname: user-defined material name
    character(len=80)::cmname = 'SandHypoplasticity'
    ! strain
    real(8)::stran(ntens)
    data stran /0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
    ! sse: specific elastic strain energy
    real(8)::sse = 0.0
    ! spd: plastic dissipation
    real(8)::spd = 0.0
    ! scd: ¡°creep¡± dissipation
    real(8)::scd = 0.0
    ! ddsddt(ntens): variation of the stress increments with respect to the temperature.
    real(8)::ddsddt(ntens)
    ! rpl: volumetric heat generation per unit time by mechanical working of the material
    real(8)::rpl = 0.0
    ! drpldt: variation of rpl with respect to temperature
    real(8)::drpldt = 0.0
    !drplde(ntens): variation of RPL with respect to strain increments
    real(8)::drplde(ntens)
    ! temp: temperature at start of increment
    real(8)::temp = 0.0
    ! dtemp: increment of temperature
    real(8)::dtemp = 0.0
    ! celent: characteristic element length
    real(8)::celent = 0.0
    ! predef: predefined field interpolated from nodes at the start of the increment
    real(8)::predef(1)
    ! dpred: increments of predefined field
    real(8)::dpred(1)
    ! coords: coordinates of this point
    real(8)::coords(3)
    ! drot(3,3): increment of rigid body rotation
    real(8)::drot(3, 3)
    ! dfgrd0(3,3): deformation gradient at the beginning of this increment
    real(8)::dfgrd0(3, 3)
    ! dfgrd1(3,3): deformation gradient at the end of the increment
    real(8)::dfgrd1(3, 3)
    ! noel: element number
    integer::noel = 0
    ! npt: integration point number
    integer::npt = 0
    ! layer: layer number (for composite shells and layered solids)
    integer::layer = 0
    ! kspt: section point number within current layer
    integer::kspt = 0
    ! time(2): 1. step time at the beginning of current increment
    !          2. total time at the beginning of current increment
    real(8)::time(2)
    ! dtime: time increment
    real(8)::dtime = 1.0
    ! pnewdt: ratio of suggested new time increment, new dt = pnewdt * dtime
    real(8)::pnewdt = 0.0
    ! kstep: step number
    integer::kstep = 0
    ! kinc: increment number
    integer::kinc = 0
    
    statev(13) = dtime
    
    call umat_hypoplasticity(stress, statev, ddsdde,    &
            sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
            stran, dstran, time, dtime,                 &
            temp, dtemp, predef, dpred, cmname,         &
            ndi, nshr, ntens, nstatv, props, nprops,    &
            coords, drot, pnewdt,                       &
            celent, dfgrd0, dfgrd1, noel, npt, layer,   &
            kspt, kstep, kinc)
    
end subroutine sand_hypoplasticity_integration