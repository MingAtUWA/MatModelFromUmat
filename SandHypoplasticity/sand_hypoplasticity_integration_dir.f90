! simplified version of umat_hypoplasticity()
! statev(13) need to set as initial integration time step
subroutine sand_hypoplasticity_integration(stress, &
    ddsdde, statev, dstran, props, error)
    !DEC$ ATTRIBUTES DLLEXPORT::sand_hypoplasticity_integration  
    implicit none
    ! stress
    real(8), intent(inout)::stress(6)
    ! ddsdde: tangential stiffness matrix
    real(8), intent(out)::ddsdde(6, 6)
    ! statev: statev variables
    real(8), intent(inout)::statev(16)
    ! total strain increment
    real(8), intent(in)::dstran(6)
    ! props: material properties
    real(8), intent(in)::props(13)
    integer, intent(out)::error  
    
    !   1. nasvdim    = maximum number of additional state variables
    !   2. tolintT    = prescribed error tolerance for the adaptive 
    !                     substepping scheme
    !   3. maxnint    = maximum number of time substeps allowed.
    !                     If the limit is exceeded abaqus is forced to reduce 
    !                     the overall time step size (cut-back) 
    !   4. DTmin      = minimum substeps size allowed.
    !                     If the limit is exceeded abaqus is forced to reduce 
    !                     the overall time step size (cut-back)
    !   5. perturb    = perturbation parameter for numerical computation of Jacobian matrices
    !   6. nfasv      = number of first additional state variable in statev field 
    
    logical::elprsw
    integer::i, nfev, nparms, inittension
    real(8)::dot_vect_h, theta, dtsub
    real(8)::sig_n(6), sig_np1(6), DDtan(6,6), pore
    real(8)::deps_np1(6), depsv_np1, norm_deps
    real(8)::norm_deps2, qq, cos3t, I1, I2, I3, norm_D2, norm_D
    ! elastic parameter in tensile state
    real(8)::youngel, nuel
    
    real(8), parameter::dtime = 1.0
    integer, parameter::nprops = 16
    integer, parameter::ntens = 6 ! stress, strain components number
    integer, parameter::nasvdim = 15
    integer, parameter::nydim = 6 + nasvdim
    integer, parameter::nasv = 7 ! additional state variable number
    integer, parameter::nyact = 6 + nasv ! nyact <= nydim
    integer, parameter::nfasv = 1
    integer, parameter::maxnint = 10000
    integer, parameter::maxninttest = 1000
    real(8), parameter::tolintT = 1.0d-2
    real(8), parameter::tolintTtest = 1.0d-1
    real(8), parameter::DTmin = 1.0d-17
    real(8), parameter::perturb = 1.0d-5
    real(8)::parms(nprops);
    ! additional state variables
    real(8)::asv(nasvdim)
    ! stresses, additional state variables
    real(8)::y(nydim), y_n(nydim), dy(nydim)
    
    ! error:
    ! 0 no problem in time integration
    ! 3 in tensile state
    ! >= 10 severe error, terminate calculation
    ! 11 too many iteration required
    ! 12 substep too small
    ! 13 parameter incorrect, fb not defined
    error = 0

    ! check material parameters and move them to array parms(nparms)
    call check_parms_h(props, nprops, parms, nparms, error)

    ! print informations about time integration, useful when problems occur for debug purpose
    elprsw = .false.

    ! suggested time substep size and initial excess pore pressure
    dtsub = statev(13)
    if ((dtsub <= 0.0) .or. (dtsub > dtime)) then
        dtsub = dtime
    endif
    
    pore = -statev(8)
      
    ! init void ratio
    if (statev(7) < 0.001) then 
        statev(7) = props(16)
    endif

    ! additional state variables
    do i = 1, nasv
        asv(i) = statev(i - 1 + nfasv)
    enddo

    ! volume strain increment and effective stress
    do i = 1, 6        
        sig_n(i) = 0.0
        deps_np1(i) = 0.0
    enddo
    call move_sig_h(stress, ntens, pore, sig_n)
    call move_eps_h(dstran, ntens, deps_np1, depsv_np1)

    ! incremental strain 
    norm_D2 = dot_vect_h(2, deps_np1, deps_np1, 6)
    norm_D = sqrt(norm_D2)

    call iniy_h(y, nydim, nasv, ntens, sig_n, asv)
    call push_h(y, y_n, nydim)

    nfev = 0
    
    ! check whether the initial state is tensile
    ! inittension = 0 not tensile; 1 tensile
    inittension = 0
    call check_RKF_h(inittension, y, nyact, nasv, parms, nparms)
    if (inittension == 0) then
    
        if (norm_D == 0.0) then
            ! no stress update
            do i = 1, nyact
                y(i) = y_n(i)
            enddo
        else
            ! normal RKF23 integration, main integration
            call rkf23_update_h(y, nyact, nasv, dtsub, tolintT, maxnint, &
                DTmin, deps_np1, parms, nparms, nfev, elprsw, dtime, error)
            if (error == 3) then ! may be close to tensile state
                ! don't update stress and state variables
                do i = 1, nyact        
                    y(i) = y_n(i)
                end do
            elseif (error >= 10) then ! fatal error
                return
            endif
        endif
        
        ! compute ddsdde
        call perturbate_h(y_n, y, nyact, nasv, dtsub, tolintT, maxnint, DTmin, &
            deps_np1, parms, nparms, nfev, elprsw, theta, ntens, DDtan, dtime, error)
    
    else ! in tensile state, elasticity
        youngel = 100.0
        nuel = 0.48
	    call calc_elasti_h(y, nyact, nasv, dtsub, tolintT, maxnint, DTmin, deps_np1, &
            parms, nparms, nfev, elprsw, dtime, DDtan, youngel, nuel, error)
    endif
    
    ! update dtsub and nfev
    if (dtsub < 0.0) then 
        dtsub = 0.0
    else if (dtsub > dtime) then 
        dtsub = dtime
    endif
    statev(13) = dtsub
    statev(10) = dfloat(nfev)
      
    ! convert solution (stress + consistent tangent matrix) to abaqus format
    ! update pore pressure and total stresses 
    call solout_h(stress, ntens, asv, nasv, ddsdde, y, nydim, pore, &
            depsv_np1, parms, nparms, DDtan)
      
    ! updated additional state variables to abaqus statev vector
    do i = 1, nasv
        statev(i - 1 + nfasv) = asv(i) 
    enddo
    
    ! transfer additional information to statev vector
    do i = 1, 6
        sig_np1(i) = y(i) ! intergranular strain
    enddo
    statev(8) = -pore
    statev(9) = -(sig_np1(1) + sig_np1(2) + sig_np1(3)) / 3.0

    ! update statev variables
    if (inittension == 0) then
        call calc_statev_h(sig_np1, statev, parms, nparms, nasv, nasvdim, deps_np1)
    endif

    return
endsubroutine sand_hypoplasticity_integration