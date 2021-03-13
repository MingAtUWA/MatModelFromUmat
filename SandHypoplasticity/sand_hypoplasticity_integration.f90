! wrapper of umat_hypoplasticity() subroutines
recursive subroutine sand_hypoplasticity_integration(stress, &
    ddsdde, statev, dstran, props, error)
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
    
    integer, intent(out)::error
    
    integer, save::aaa
    
    !cmname: user-defined material name
    character*80, parameter::cmname = 'SandHypoplasticity'
    ! strain
    real(8), parameter::stran(ntens) = 0.0
    !data stran /0.0, 0.0, 0.0, 0.0, 0.0, 0.0/
    ! sse: specific elastic strain energy
    real(8), parameter::sse = 0.0
    ! spd: plastic dissipation
    real(8), parameter::spd = 0.0
    ! scd: ¡°creep¡± dissipation
    real(8), parameter::scd = 0.0
    ! ddsddt(ntens): variation of the stress increments with respect to the temperature.
    real(8), parameter::ddsddt(ntens) = 0.0
    ! rpl: volumetric heat generation per unit time by mechanical working of the material
    real(8), parameter::rpl = 0.0
    ! drpldt: variation of rpl with respect to temperature
    real(8), parameter::drpldt = 0.0
    !drplde(ntens): variation of RPL with respect to strain increments
    real(8), parameter::drplde(ntens) = 0.0
    ! temp: temperature at start of increment
    real(8), parameter::temp = 0.0
    ! dtemp: increment of temperature
    real(8), parameter::dtemp = 0.0
    ! celent: characteristic element length
    real(8), parameter::celent = 0.0
    ! predef: predefined field interpolated from nodes at the start of the increment
    real(8), parameter::predef(1) = 0.0
    ! dpred: increments of predefined field
    real(8), parameter::dpred(1) = 0.0
    ! coords: coordinates of this point
    real(8), parameter::coords(3) = 0.0
    ! drot(3,3): increment of rigid body rotation
    real(8), parameter::drot(3, 3) = 0.0
    ! dfgrd0(3,3): deformation gradient at the beginning of this increment
    real(8), parameter::dfgrd0(3, 3) = 0.0
    ! dfgrd1(3,3): deformation gradient at the end of the increment
    real(8), parameter::dfgrd1(3, 3) = 0.0
    ! noel: element number
    integer, parameter::noel = 0
    ! npt: integration point number
    integer, parameter::npt = 0
    ! layer: layer number (for composite shells and layered solids)
    integer, parameter::layer = 0
    ! kspt: section point number within current layer
    integer, parameter::kspt = 0
    ! time(2): 1. step time at the beginning of current increment
    !          2. total time at the beginning of current increment
    real(8), parameter::time(2) = 0.0
    ! dtime: time increment
    real(8), parameter::dtime = 1.0
    ! pnewdt: ratio of suggested new time increment, new dt = pnewdt * dtime
    real(8), parameter::pnewdt = 0.0
    ! kstep: step number
    integer, parameter::kstep = 0
    ! kinc: increment number
    integer, parameter::kinc = 0
    
    call umat_hypoplasticity(stress, statev, ddsdde,    &
            sse, spd, scd, rpl, ddsddt, drplde, drpldt, &
            stran, dstran, time, dtime,                 &
            temp, dtemp, predef, dpred, cmname,         &
            ndi, nshr, ntens, nstatv, props, nprops,    &
            coords, drot, pnewdt,                       &
            celent, dfgrd0, dfgrd1, noel, npt, layer,   &
            kspt, kstep, kinc)
    
    ! no error output
    error = 0
    
end subroutine sand_hypoplasticity_integration