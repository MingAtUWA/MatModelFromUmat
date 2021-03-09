! Based on RKF23 implementation by D. Masin
subroutine umat_hypoplasticity(stress,statev,ddsdde,sse,spd,scd, &
    rpl,ddsddt,drplde,drpldt, &
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname, &
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc) 
    implicit none

    character*80::cmname
    integer::ntens, ndi, nshr, nstatv, nprops, noel, npt, layer, kspt, kstep, kinc, inittension
    
    real(8)::stress(ntens), statev(nstatv), &
    ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens), &
    stran(ntens), dstran(ntens), time(2), predef(1), dpred(1), &
    props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
    real(8)::sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent

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
    integer::i, error, nfev, nparms

    real(8)::dot_vect_h
    real(8)::parms(nprops), theta, dtsub
    real(8)::sig_n(6), sig_np1(6), DDtan(6,6), pore
    real(8)::deps_np1(6), depsv_np1, norm_deps
    real(8)::norm_deps2, qq, cos3t, I1, I2, I3, norm_D2, norm_D
    ! elastic parameter in tensile state
    real(8)::youngel, nuel
    
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

    ! additional state variables
    real(8)::asv(nasvdim)

    ! stresses, additional state variables
    real(8)::y(nydim), y_n(nydim), dy(nydim)
    
    ! error:
    ! 0 no problem in time integration
    ! 1 problems in evaluation of the time rate, (e.g. undefined stress state), reduce time integration substeps
    ! 3 problems in time integration, reduce abaqus load increment (cut-back)
    ! 10 severe error, terminate calculation
    error = 0

    ! check material parameters and move them to array parms(nparms)
    call check_parms_h(props, nprops, parms, nparms, error)

    ! print informations about time integration, useful when problems occur for debug purpose
    elprsw = .false.

    ! suggested time substep size and initial excess pore pressure
    dtsub = statev(13)
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
    do i=1,6        
        sig_n(i)=0
        deps_np1(i)=0
    enddo
    call move_sig_h(stress,ntens,pore,sig_n)
    call move_eps_h(dstran,ntens,deps_np1,depsv_np1)

    ! incremental strain 
    norm_D2 = dot_vect_h(2, deps_np1, deps_np1, 6)
    norm_D = sqrt(norm_D2)

    call iniy_h(y, nydim, nasv, ntens,sig_n,asv)
    call push_h(y, y_n, nydim)

    ! check whether the initial state is not tensile
    ! inittension = 0 not tensile; 1 tensile
    inittension = 0
    call check_RKF_h(inittension, y, nyact, nasv, parms, nparms)

    ! local integration using adaptive RKF23 method, consistent Jacobian and error estimation
    if ((dtsub <= 0.0) .or. (dtsub > dtime)) then
        dtsub = dtime
    endif

    nfev = 0
    if (inittension == 0) then
      
        if(norm_D == 0.0) then
            ! no stress update
            do i = 1, nyact
                y(i) = y_n(i)
            enddo
        else ! testing == 0
            ! normal RKF23 integration, main integration
            ! error here
            call rkf23_update_h(y, nyact, nasv, dtsub, tolintT, maxnint, &
                DTmin, deps_np1, parms, nparms, nfev, elprsw, dtime, error)
            ! check integration error
            if (error == 3) then ! may be close to tensile state
                !call wrista_h(1, y, nydim, deps_np1, dtime, coords, statev, &
                !    nstatv, parms, nparms, noel, npt, ndi, nshr, kstep, kinc)
                do i = 1, nyact        
                    y(i) = y_n(i)
                enddo
            elseif (error == 10) then ! fatal error
                return
            endif
        endif
        
        ! compute ddsdde
        call perturbate_h(y_n,y,nyact,nasv,dtsub,tolintT,maxnint,DTmin, &
            deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DDtan, dtime,error)
    
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
      
    ! convert solution (stress + cons. tangent) to abaqus format
    ! update pore pressure and total stresses 
    call solout_h(stress, ntens, asv, nasv, ddsdde, y, nydim, pore, &
            depsv_np1, parms, nparms, DDtan)
      
    ! updated additional state variables to abaqus statev vector
    do i = 1, nasv
        statev(i - 1 + nfasv) = asv(i) 
    enddo
      
    ! transfer additional information to statev vector
    do i = 1,6
        sig_np1(i) = y(i) ! intergranular strain
    enddo
    statev(8) = -pore
    statev(9) = -(sig_np1(1) + sig_np1(2) + sig_np1(3)) / 3.0

    ! update statev variables
    if (inittension == 0) then
        call calc_statev_h(sig_np1, statev, parms, nparms, nasv, nasvdim, deps_np1)
    endif

    return
endsubroutine umat_hypoplasticity 

! copy props(nprops) to parms(nparms)
subroutine check_parms_h(props, nprops, parms, nparms, error)
    implicit none
    integer::nprops,nparms,i,error
    real(8)::props(nprops),parms(nprops)
    real(8)::phi_deg,phi,hs,en,ed0,ec0,ei0,alpha,beta
    real(8)::m_R,m_T,r_uc,beta_r,chi,bulk_w,p_t
    
    nparms = nprops
    do i = 1, nprops
        parms(i) = props(i)
    enddo

    ! recover material parameters
    phi_deg = parms(1)
    phi = phi_deg * 3.14159265358979323846 / 180.0
    parms(1) = phi
    p_t = parms(2)
    hs = parms(3)
    en = parms(4)
    ed0 = parms(5)
    ec0 = parms(6)
    ei0 = parms(7)
    alpha = parms(8)
    beta = parms(9)
    m_R = parms(10) 
    m_T = parms(11)
    r_uc = parms(12)
    beta_r = parms(13)
    chi = parms(14)
    bulk_w = parms(15)        
    
    return
endsubroutine

! finds vectors F_sigma and F_q in F(y)
subroutine get_F_sig_q_h(sig, q, nasv, parms, nparms, deps, F_sig, F_q, error)
    implicit none
    real(8)::dot_vect_h  
    real(8)::sig(6),q(nasv),parms(nparms),deps(6)
    real(8)::MM(6,6),HH(nasv,6),F_sig(6),F_q(nasv)
    real(8)::LL(6,6),NN(6),norm_D,norm_D2
    integer::nparms, nasv, ii, istrain, error

    ! compute tangent operators
    if(parms(10) <= 0.5) then
        istrain = 0 
    else 
        istrain = 1
    endif

    call get_tan_h(deps, sig, q, nasv, parms, nparms, &
            MM, HH, LL, NN, istrain, error)

    ! compute F_sig = MM * deps
    if (istrain == 1) then
        call matmul_h(MM, deps, F_sig, 6, 6, 1)
    else 
        call matmul_h(LL, deps, F_sig, 6, 6, 1)
		norm_D2 = dot_vect_h(2, deps, deps, 6)
		norm_D = dsqrt(norm_D2)
        do ii = 1, 6
            F_sig(ii) = F_sig(ii) + NN(ii) * norm_D
        enddo
    endif

    ! compute F_q = HH * deps
    call matmul_h(HH, deps, F_q, nasv, 6, 1)
    
    return
endsubroutine

! computes matrices M and H for Masin hypoplastic model for clays with intergranular strains
subroutine get_tan_h(deps, sig, q, nasv, parms, nparms, &
    MM, HH, LL, NN, istrain, error)
!stress and strain convention: tension and extension positive
    implicit none
    integer::nparms,nasv,i,j,error
    real(8)::dot_vect_h
    real(8)::sig(6),q(nasv),parms(nparms),deps(6)
    real(8)::eta(6),eta_dev(6),del(6),void,sig_star(6)
    real(8)::eta_del(6),eta_delta(6),eta_eps(6)
    real(8)::norm_del,norm_del2,norm_deps,norm_deps2,eta_dn2
    real(8)::pp,qq,cos3t,I1,I2,I3,tanpsi
    real(8)::a,a2,FF,fd,fs
    real(8)::num,den,aF,Fa2,eta_n2,norm_m,norm_m2
    real(8)::II(6,6),IU(6,6)
    real(8)::MM(6,6),HH(nasv,6),LL(6,6),NN(6),AA(6,6),m(6)
    integer istrain
    real(8)::m_dir(6),m_dir1(6),Leta(6),H_del(6,6),H_e(6)
    real(8)::load,rho
    real(8)::temp1,temp2,temp3,temp4
    real(8)::phi,hs,en,ed0,ec0,ei0,alpha,beta,r_uc
    real(8)::m_R,m_T,beta_r,chi,bulk_w,p_t,sinphi,sinphi2
    real(8)::ec,ed,ei,bauer,fb,fe,az
    
    real(8), parameter::sqrt2 = 1.414213562
    real(8), parameter::twosqrt2 = 2.828427125
    real(8), parameter::sqrt3 = 1.732050808
    real(8), parameter::sqrt6 = 2.449489743
    real(8), parameter::ln2m1 = 1.442695041 !one/dlog(two)
    real(8), parameter::tiny = 1.0d-17
    
    do i = 1, 6
        do j = 1, 6
            MM(i,j) = 0.0
            LL(i,j) = 0.0
            II(i,j) = 0.0
            IU(i,j) = 0.0
            H_del(i,j) = 0.0
        enddo
        eta_del(i) = 0.0
        eta_delta(i) = 0.0
        eta_eps(i) = 0.0
    enddo

    do i = 1, nasv
        do j = 1, 6
            HH(i,j) = 0.0
        enddo
    enddo

    ! fourth order identity tensors in Voigt notation
    II(1,1) = 1.0
    II(2,2) = 1.0
    II(3,3) = 1.0
    II(4,4) = 0.5
    II(5,5) = 0.5
    II(6,6) = 0.5
    
    IU(1,1) = 1.0
    IU(2,2) = 1.0
    IU(3,3) = 1.0
    IU(4,4) = 1.0
    IU(5,5) = 1.0
    IU(6,6) = 1.0

    ! recover material parameters
    phi	  = parms(1)
    hs    = parms(3)
    en    = parms(4)
    ed0   = parms(5)
    ec0   = parms(6)
    ei0   = parms(7)
    alpha = parms(8)
    beta  = parms(9)
    m_R = parms(10) 
    m_T = parms(11)
    r_uc = parms(12)
    beta_r = parms(13)
    chi = parms(14)
    bulk_w = parms(15)
    p_t = parms(2)

    sinphi = dsin(phi)
    sinphi2 = sinphi * sinphi
    
    ! recover internal state variables
    del(1)=q(1)
    del(2)=q(2)
    del(3)=q(3)
    del(4)=q(4)
    del(5)=q(5)
    del(6)=q(6)
    void=q(7)
    
    ! axis translation due to cohesion (p_t>0)
    sig_star(1) = sig(1) - p_t
    sig_star(2) = sig(2) - p_t
    sig_star(3) = sig(3) - p_t
    sig_star(4) = sig(4)
    sig_star(5) = sig(5)
    sig_star(6) = sig(6)

    ! strain increment and intergranular strain directions
    norm_deps2 = dot_vect_h(2,deps,deps,6)
    norm_del2 = dot_vect_h(2,del,del,6)
    norm_deps = dsqrt(norm_deps2)
    norm_del = dsqrt(norm_del2)
    
    if (norm_del >= tiny) then
        do i=1,6
            eta_del(i) = del(i) / norm_del
        enddo
    endif

    eta_delta(1) = eta_del(1)
    eta_delta(2) = eta_del(2)
    eta_delta(3) = eta_del(3)
    eta_delta(4) = 0.5 * eta_del(4)
    eta_delta(5) = 0.5 * eta_del(5)
    eta_delta(6) = 0.5 * eta_del(6)

    if (norm_deps >= tiny) then
        do i = 1, 6
            eta_eps(i) = deps(i) / norm_deps
        enddo
    endif

    ! auxiliary stress tensors
    call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
    ! T_cap
    eta(1) = sig_star(1) / I1
    eta(2) = sig_star(2) / I1
    eta(3) = sig_star(3) / I1
    eta(4) = sig_star(4) / I1
    eta(5) = sig_star(5) / I1
    eta(6) = sig_star(6) / I1   
    ! T_cap_dot
    eta_dev(1) = eta(1) - 1.0/3.0
    eta_dev(2) = eta(2) - 1.0/3.0
    eta_dev(3) = eta(3) - 1.0/3.0
    eta_dev(4) = eta(4)
    eta_dev(5) = eta(5)
    eta_dev(6) = eta(6)

    ! functions a and F
    eta_dn2 = dot_vect_h(1,eta_dev,eta_dev,6)
    tanpsi = sqrt3 * dsqrt(eta_dn2)
    temp1 = 1.0/8.0 * tanpsi * tanpsi + (2.0 - tanpsi * tanpsi)/(2.0 + sqrt2*tanpsi*cos3t)
    temp2 = tanpsi / twosqrt2

    a = sqrt3 * (3.0 - sin(phi)) / (twosqrt2 * sin(phi))
    a2 = a * a
    FF = dsqrt(temp1) - temp2

    ! barotropy and pyknotropy functions
    bauer = dexp(-(-I1/hs)**en)
    ed = ed0 * bauer
    ec = ec0 * bauer
    ei = ei0 * bauer

    ! fb
    temp1 = 3.0 + a * a - a * sqrt3 * ((ei0-ed0)/(ec0-ed0))**alpha
    if (temp1 < 0.0) stop 'factor fb not defined'
    fb = hs/en/temp1*(1.0 + ei)/ei*(ei0/ec0)**beta*(-I1/hs)**(1.0-en)
    !fb = hs/en/temp1*(one+ei)/ei*(-I1/hs)**(one-en)
      
    ! fe
    fe = (ec/void)**beta
    !fe = (ei/void)**beta
      
    fs = fb * fe

    ! fd
    if (void >= ed) then
        fd = ((void - ed) / (ec - ed))**alpha
    else
        fd = 0.0
    endif

    ! tensor L
    eta_n2 = dot_vect_h(1,eta,eta,6)
    do i = 1,6
        do j = 1, 6
            LL(i,j) = (II(i,j)*FF*FF + a2*eta(i)*eta(j)) / eta_n2
        enddo
    enddo
    
    ! tensor NN
    do i = 1, 6
        NN(i) = FF * a * (eta(i)+eta_dev(i)) / eta_n2
    enddo
    
    ! BEGIN intergranular STRAIN
    if(istrain == 1) then

        ! loading function
        load = dot_vect_h(2,eta_del,eta_eps,6)
        ! intergranular strain--related tensors
        rho = norm_del / r_uc
        if (rho > 1.0) then
            rho = 1.0
        endif

        call matmul_h(LL,eta_del,Leta,6,6,1)
          
        ! tangent stiffness M(sig,q,eta_eps)
        temp1 = ((rho**chi) * m_T + (1.0 - rho**chi) * m_R) * fs

        if (load > 0.0) then ! loading
          
            temp2 = (rho**chi) * (1.0 - m_T) * fs
            temp3 = (rho**chi) * fs * fd
            do i = 1, 6
                do j = 1, 6
                    AA(i,j) = temp2 * Leta(i) * eta_delta(j) + temp3 * NN(i) * eta_delta(j)
                    MM(i,j) = temp1 * LL(i,j) + AA(i,j)
                enddo
            enddo
              
        else ! reverse loading
            temp4 = (rho**chi)*(m_R-m_T)*fs
            do i = 1, 6
                do j = 1, 6
                    AA(i,j) = temp4*Leta(i)*eta_delta(j)
                    MM(i,j) = temp1*LL(i,j)+AA(i,j)
                enddo
            enddo
        endif

        ! intergranular strain evolution function
        ! NOTE: H_del transforms a strain-like vector into a strain-like vector
        ! eta_del(i) instead of eta_delta(i), I = 6x6 unit matrix
        if (load > 0.0) then
            do i = 1, 6
                do j = 1, 6
                    H_del(i,j) = IU(i,j) - (rho**beta_r) * eta_del(i) * eta_delta(j)
                enddo
            enddo
        else
            do i = 1, 6
                H_del(i,i) = 1.0
            enddo
        endif

        ! void ratio evolution function (tension positive)
        do i = 1, 6 
            if (i <= 3) then
                H_e(i) = 1.0 + void
            else
                H_e(i) = 0.0
            endif
        enddo

        ! assemble hardening matrix
        ! for clay hypoplasticity only, can be ignored for sand
        do i = 1, nasv
            if (i <= 6) then
                do j = 1, 6
                    HH(i,j) = H_del(i,j)
                enddo
            else
                do j = 1, 6
                    HH(i,j)=H_e(j)
                    enddo
            endif
        enddo
        
    ! not accouting for intergranular strain
    else if (istrain == 0) then
      
        do i = 1, 6 
            if (i <= 3) then
                H_e(i) = 1.0 + void
            else
                H_e(i) = 0.0
            endif
        enddo
          
        do i = 1, nasv
            if (i <= 6) then
                do j = 1, 6
                    HH(i,j) = 0.0
                enddo
            else
                do j = 1, 6
                    HH(i,j) = H_e(j)
                enddo
            endif
        enddo
        
    endif ! end istrain/noistrain switch   

    do i = 1, 6
        do j = 1, 6
            LL(i,j) = LL(i,j) * fs
        enddo
        NN(i) = NN(i) * fs * fd
    enddo        
    
    return
endsubroutine

! initializes the vector of state variables
subroutine iniy_h(y, nydim, nasv, ntens, sig, qq)
    implicit none
    integer i, nydim, nasv, ntens
    real(8)::y(nydim), qq(nasv), sig(ntens)
    do i=1, nydim
        y(i) = 0.0
    enddo
    do i=1, ntens
        y(i) = sig(i)
    enddo
    ! additional state variables
    do i=1,nasv
        y(6+i) = qq(i)
    enddo
    return
endsubroutine

! calculate invariants of strain tensor
subroutine inv_eps_h(eps, eps_v, eps_s, sin3t)
    implicit none
    integer::i
    real(8)::eps(6), edev(6), edev2(6), ev3, tredev3
    real(8)::eps_v, eps_s, sin3t, norm2, numer, denom
    ! volumetric strain
    eps_v = eps(1)+eps(2)+eps(3)
    ev3 = eps_v / 3.0
    ! deviator strain
    edev(1) = eps(1) - ev3
    edev(2) = eps(2) - ev3
    edev(3) = eps(3) - ev3
    edev(4) = eps(4) * 0.5
    edev(5) = eps(5) * 0.5
    edev(6) = eps(6) * 0.5
    ! second invariant
    norm2 = edev(1) * edev(1) + edev(2) * edev(2) + edev(3) * edev(3) &
        + 2.0 * (edev(4) * edev(4) + edev(5) * edev(5) + edev(6) * edev(6))
    eps_s = dsqrt(norm2 * 2.0 / 3.0)
    ! components of (edev_ij)(edev_jk)
    edev2(1) = edev(1) * edev(1) + edev(4) * edev(4) + edev(5) * edev(5)
    edev2(2) = edev(4) * edev(4) + edev(2) * edev(2) + edev(6) * edev(6)
    edev2(3) = edev(6) * edev(6) + edev(5) * edev(5) + edev(3) * edev(3)
    edev2(4) = 2.0 * (edev(1) * edev(4) + edev(4) * edev(2) + edev(6) * edev(5))
    edev2(5) = 2.0 * (edev(5) * edev(1) + edev(6) * edev(4) + edev(3) * edev(5))
    edev2(6) = 2.0 * (edev(4) * edev(5) + edev(2) * edev(6) + edev(6) * edev(3))
    ! Lode angle
    if(eps_s == 0.0) then 
        sin3t = -1.0           
    else
        tredev3 = 0.0
        do i = 1, 6
            tredev3 = tredev3 + edev(i) * edev2(i)
        enddo
        numer = dsqrt(6.0d0) * tredev3
        denom = (dsqrt(norm2))**3
        sin3t = numer / denom
        if (dabs(sin3t) > 1.0) then
            sin3t = sin3t / dabs(sin3t)
        endif
    endif
    return
endsubroutine

! calculate invariants of stress tensor
subroutine inv_sig_h(sig, pp, qq, cos3t, I1, I2, I3)
! NOTE: Voigt notation is used with the following index conversion
! 11 -> 1, 22 -> 2, 33 -> 3, 12 -> 4, 13 -> 5, 23 -> 6
    implicit none
    real(8)::sig(6), sdev(6), eta(6), eta_d(6), eta_d2(6)
    real(8)::xmin1,xmin2,xmin3
    real(8)::tretadev3,pp,qq,cos3t,I1,I2,I3
    real(8)::norm2,norm2sig,norm2eta,numer,denom
    real(8)::dot_vect_h
    real(8), parameter::sqrt6 = 2.449489743
    real(8), parameter::tiny = 1.0d-18
    
    ! trace and mean stress
    I1 = sig(1) + sig(2) + sig(3)
    pp = I1 / 3.0
    ! deviator stress
    sdev(1) = sig(1) - pp
    sdev(2) = sig(2) - pp
    sdev(3) = sig(3) - pp
    sdev(4) = sig(4)
    sdev(5) = sig(5)
    sdev(6) = sig(6)
    
    ! normalized stress and dev. normalized stress
    if(I1 /= 0.0) then
        eta(1) = sig(1)/I1
        eta(2) = sig(2)/I1
        eta(3) = sig(3)/I1
        eta(4) = sig(4)/I1
        eta(5) = sig(5)/I1
        eta(6) = sig(6)/I1
    else
        eta(1) = sig(1) / tiny
        eta(2) = sig(2) / tiny
        eta(3) = sig(3) / tiny
        eta(4) = sig(4) / tiny
        eta(5) = sig(5) / tiny
        eta(6) = sig(6) / tiny        
    endif
    eta_d(1) = eta(1) - 1.0/3.0
    eta_d(2) = eta(2) - 1.0/3.0
    eta_d(3) = eta(3) - 1.0/3.0
    eta_d(4) = eta(4)
    eta_d(5) = eta(5)
    eta_d(6) = eta(6)
    ! second invariants
    norm2 = dot_vect_h(1,sdev,sdev,6)
    norm2sig = dot_vect_h(1,sig,sig,6)
    norm2eta = dot_vect_h(1,eta_d,eta_d,6)
    qq = dsqrt(1.5 * norm2)
    I2 = 0.5 * (norm2sig - I1*I1)
    ! components of (eta_d_ij)(eta_d_jk)
    eta_d2(1) = eta_d(1)*eta_d(1)+eta_d(4)*eta_d(4)+eta_d(5)*eta_d(5)
    eta_d2(2) = eta_d(4)*eta_d(4)+eta_d(2)*eta_d(2)+eta_d(6)*eta_d(6)
    eta_d2(3) = eta_d(6)*eta_d(6)+eta_d(5)*eta_d(5)+eta_d(3)*eta_d(3)
    eta_d2(4) = eta_d(1)*eta_d(4)+eta_d(4)*eta_d(2)+eta_d(6)*eta_d(5)
    eta_d2(5) = eta_d(5)*eta_d(1)+eta_d(6)*eta_d(4)+eta_d(3)*eta_d(5)
    eta_d2(6) = eta_d(4)*eta_d(5)+eta_d(2)*eta_d(6)+eta_d(6)*eta_d(3)        
    ! Lode angle
    if(norm2eta < tiny) then 
        cos3t = -1.0            
    else
        tretadev3 = dot_vect_h(1,eta_d,eta_d2,6)
        numer = -sqrt6*tretadev3
        denom=(dsqrt(norm2eta))**3
        cos3t=numer/denom
        if(dabs(cos3t) > 1.0) then
            cos3t = cos3t / dabs(cos3t)
        endif
    endif 
    ! determinant
    xmin1=sig(2)*sig(3)-sig(6)*sig(6)
    xmin2=sig(4)*sig(3)-sig(6)*sig(5)
    xmin3=sig(4)*sig(6)-sig(5)*sig(2)
    I3=sig(1)*xmin1-sig(4)*xmin2+sig(5)*xmin3
    return
endsubroutine

! calculate consistent tengential stiffness matrix
subroutine perturbate_h(y_n,y_np1,n,nasv,dtsub,err_tol,maxnint, DTmin, &
    deps_np1,parms,nparms,nfev,elprsw,theta,ntens,DD,dtime, error)
    implicit none
    logical::elprsw
    integer::ntens, jj, kk, i
    integer::n, nasv, nparms, nfev
    integer::maxnint, error, istrain
    real(8)::y_n(n), y_np1(n), y_star(n), parms(nparms)
    real(8)::dtsub, err_tol, DTmin, dtime
    real(8)::theta, sig(6), q(nasv)
    real(8)::deps_np1(6), deps_star(6)
    real(8)::dsig(6), DD(6,6), HHtmp(nasv,6)
    real(8)::LL(6,6), NN(6)
          
    ! initialize DD and y_star
    if (parms(10) <= 0.5) then
        istrain = 0 
    else 
        istrain = 1
    endif

    do kk = 1, 6
        do jj = 1, 6
            DD(kk, jj) = 0.0
        enddo
    enddo
      
    do i = 1, 6
        sig(i) = y_n(i)
    enddo
      
    do i = 1, nasv
        q(i) = y_n(6+i)
    enddo
        
    call push_h(y_n,y_star,n)

    if (error /= 10) then
        call get_tan_h(deps_np1, sig, q, nasv, parms, nparms, &
                    DD, HHtmp, LL, NN, istrain, error)                
    endif
    
    if (istrain == 0) then
        do kk = 1, 6
            do jj = 1, 6
                DD(kk, jj) = LL(kk, jj)
            enddo
        enddo
    else
        do kk=1,6
            do jj=1,6
                DD(kk, jj) = parms(10) * LL(kk, jj)
            enddo
        enddo
    endif
    
    return
endsubroutine  

! calculate coefficient kRK from current state y and strain increment deps
! after Masin's hypoplastic model for clays with intergranular strains
subroutine rhs_h(y, ny, nasv, parms, nparms, deps, kRK, nfev, error)
    implicit none
    integer error,ny,nparms,nasv,i,nfev
    real(8)::y(ny),kRK(ny),parms(nparms),deps(6)
    real(8)::sig(6),q(nasv)
    real(8)::F_sig(6),F_q(nasv)
    
    ! update counter for the number of f(y) evaluations
    nfev = nfev + 1
    
    ! initialize kRK
    do i = 1, ny
        kRK(i) = 0.0
    enddo

    ! recover current state variables (sig,q)                   
    do i = 1, 6
        sig(i) = y(i)
    enddo

    do i = 1, nasv
        q(i) = y(6 + i)
    enddo
      
    ! build F_sig(6) and F_q(nasv) vectors
    call get_F_sig_q_h(sig, q, nasv, parms, nparms, deps, F_sig, F_q, error)
      
    if(error == 10) return

    do i = 1, 6
        kRK(i) = F_sig(i)
    enddo 
    write(*, *) kRK(1), kRK(2), kRK(3), kRK(4), kRK(5), kRK(6)
    
    do i = 1, nasv
        kRK(6 + i) = F_q(i)
    enddo                   

    return
endsubroutine

! numerical solution of y'=f(y) using self-apdative explicit RKF23 scheme
subroutine rkf23_update_h(y, n, nasv, dtsub, err_tol, maxnint, DTmin, &
    deps_np1, parms, nparms, nfev, elprsw, dtime, error)
    implicit none
    logical elprsw
    integer n,nasv,nparms,i,ksubst,nfev
    integer maxnint,error,error_RKF
    real(8)::y(n), parms(nparms), dtsub, err_tol, DTmin
    real(8)::deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
    real(8)::y_hat(n)
    real(8)::T_k,DT_k,dtime
    real(8)::kRK_1(n),kRK_2(n),kRK_3(n)
    real(8)::norm_R,S_hull, temp
    
    error_RKF = 0
    T_k = 0.0      
    DT_k = dtsub / dtime
    ksubst = 0
    nfev = 0
    do i = 1, n
        y_k(i) = y(i)
    enddo
        
    ! start substepping 
    do while (T_k < (1.0 - 1.0d-3)) 
        ksubst = ksubst + 1
        if (ksubst > maxnint) then ! too many substep required
            error = 11
            return
        endif          
        
        write(*, *) 'y_k: ', y_k(1), y_k(2), y_k(3), y_k(4), y_k(5), y_k(6)
        call check_RKF_h(error_RKF,y_k,n,nasv,parms,nparms)
        if (error_RKF == 1) then 
            error = 3
            return
        endif
        call rhs_h(y_k,n,nasv,parms,nparms,deps_np1,kRK_1, nfev, error)
        if (error == 10) return
        
        ! cal y_2
        temp = 0.5 * DT_k
        do i = 1, n
            y_2(i) = y_k(i) + temp * kRK_1(i)
        enddo
        write(*, *) 'y_2', y_2(1), y_2(2), y_2(3), y_2(4), y_2(5), y_2(6)
        
        call check_RKF_h(error_RKF, y_2, n, nasv, parms, nparms)
        if (error_RKF == 1) then 
            error = 3
            return
        endif
        call rhs_h(y_2,n,nasv,parms,nparms,deps_np1,kRK_2, nfev,error)
        if(error == 10) return
                                     
        ! cal y_3
        do i = 1, n
            y_3(i) = y_k(i) - DT_k * kRK_1(i) + 2.0 * DT_k * kRK_2(i)
        enddo
        write(*, *) 'y_3', y_3(1), y_3(2), y_3(3), y_3(4), y_3(5), y_3(6)

        call check_RKF_h(error_RKF, y_3, n, nasv, parms, nparms)
        if (error_RKF == 1) then 
            error=3
            return
        endif
        call rhs_h(y_3, n, nasv, parms,nparms,deps_np1,kRK_3, nfev,error)
        if (error==10) return
                      
        ! approx. solutions of 2nd (y_til) and 3rd (y_hat) order
        do i = 1, n        
            y_til(i) = y_k(i) + DT_k * kRK_2(i)
            y_hat(i) = y_k(i) + DT_k * (1.0/6.0 * kRK_1(i) + 2.0/3.0 * kRK_2(i) + 1.0/6.0 * kRK_3(i))
        enddo

        ! local error estimate
        call norm_res_h(y_til,y_hat,n,nasv,norm_R)
        ! check if output y_hat can be used as an input into the next step
        call check_RKF_h(error_RKF,y_hat,n,nasv,parms,nparms)

        if (error_RKF /= 0) then
            error = 3
	        return
        endif

        ! time step size estimator according to Hull  	
        if (norm_R /= 0.0) then
            S_hull = 0.9 * DT_k * (err_tol/norm_R)**(1.0/3.0)
        else
            S_hull = 1
        endif

        if (norm_R < err_tol) then ! substep accepted                    
            ! update y_k and T_k and estimate new substep size DT_k
            do i = 1, n        
                y_k(i) = y_hat(i)
            enddo
            T_k = T_k + DT_k
            DT_k = min(4.0 * DT_k, S_hull) ! adjust DT_k
            dtsub = DT_k * dtime
            DT_k = min(1.0 - T_k, DT_k)       
        else ! substep not accepted, recompute with smaller substep
            DT_k = max(DT_k / 4.0, S_hull)
            ! check for minimum step size
            if (DT_k < DTmin) then ! substep too small
                error = 12
                return
            endif
        endif
    enddo
    
    ! complete updating y
    do i = 1,n
        y(i) = y_k(i)
    enddo
    return
endsubroutine

! checks if RKF23 solout vector y is OK for hypoplasticity
subroutine check_RKF_h(error_RKF, y, ny, nasv, parms, nparms)
    implicit none
    integer::error_RKF,ny,nasv,i,nparms,testnan,iopt
    real(8)::y(ny),parms(nparms)
    real(8)::sig(6),pmean,sig_star(6)
    real(8)::xN1(3),xN2(3),xN3(3),S(3),P,Q,tmin, p_t
    real(8), parameter::minstress = 0.9
    
    ! mean stress
    p_t = parms(2)
    do i = 1, 6
        sig(i) = y(i)
    enddo
    sig_star(1) = sig(1) - p_t
    sig_star(2) = sig(2) - p_t
    sig_star(3) = sig(3) - p_t
    sig_star(4) = sig(4)
    sig_star(5) = sig(5)
    sig_star(6) = sig(6)
    write (*, *) 's11-s33: ', y(1), y(2), y(3)
    write (*, *) 's12-s31: ', y(4), y(5), y(6)
    pmean = -(sig_star(1) + sig_star(2) + sig_star(3)) / 3.0    
    if (pmean <= minstress) then ! tension
        error_RKF = 1
    endif

    ! get minimum principal stress
    iopt = 0
    call PrnSig_h(iopt, sig_star, xN1, xN2, xN3, S(1),S(2),S(3), P, Q)
    tmin = -S(1)
    if (tmin >= -S(2)) tmin = -S(2)
    if (tmin >= -S(3)) tmin = -S(3)
    if (tmin <= minstress) then ! tension
        error_RKF = 1
    endif
    
    return
endsubroutine

! ========================== Result formatting =======================

! copy the vector of state variables to umat output
subroutine solout_h(stress, ntens, asv, nasv, ddsdde, y, &
    nydim, pore, depsv_np1, parms, nparms, DD)
! NOTE: solid mechanics convention for stress and strain components
!       pore is always positive in compression
    implicit none
    integer nydim,nasv,nparms,ntens,i,j
    real(8)::y(nydim),asv(nasv),stress(ntens)
    real(8)::ddsdde(ntens,ntens),DD(6,6)
    real(8)::parms(nparms),bulk_w,pore,depsv_np1 

    ! update excess pore pressure (if undrained conditions)
    ! compression as positive
    bulk_w = parms(15)
    pore = pore - bulk_w * depsv_np1

    ! updated total stresses (effective stresses stored in y(1:6))
    do i = 1, ntens
        if (i<=3) then
            stress(i) = y(i) - pore
        else
            stress(i) = y(i)
        endif
    enddo

    ! additional state variables (first 6 components are intergranular strains)
    do i = 1, nasv
        asv(i) = y(6+i)
    enddo

    ! consistent tangent stiffness
    do j = 1, ntens
        do i = 1, ntens
            ddsdde(i,j) = DD(i,j)      
        enddo
    enddo

    do j = 1, 3
        do i = 1, 3
            ddsdde(i,j) = ddsdde(i,j) + bulk_w        
        enddo
    enddo
      
    return
endsubroutine

! computes additional state variables for postprocessing
subroutine calc_statev_h(stress, statev, parms, nparms, nasv, nasvdim, deps)
    implicit none
    logical elprsw
    integer ntens,jj,kk,i
    integer n,nasv,nparms,nfev,nasvdim
    integer maxnint,error

    real(8)::parms(nparms),dot_vect_h
    real(8)::stress(6),statev(nasvdim)
    real(8)::deps(6),tmax,tmin
    real(8)::MM(6,6),HHtmp(nasv,6)
    real(8)::LL(6,6),NN(6)
    integer istrain
    real(8)::iopt
    real(8)::I1,I2,I3,cos3t,pp,qq
    real(8)::sin2phi,sinphi,sig_star(6),p_t
    real(8)::norm_del,norm_del2,del(6)

    ! calc phimob (statev 11) from Matsuoka-Nakai YS
    p_t = parms(2)
    do i = 1, 3
        sig_star(i) = stress(i) - p_t
    enddo
    do i = 4, 6
        sig_star(i) = stress(i)
    enddo
    call inv_sig_h(sig_star,pp,qq,cos3t,I1,I2,I3)
      
    if(I3 /= 0) then
        sin2phi = (9.0d0+I1*I2/I3)/(1.0d0+I1*I2/I3)
    else 
    sin2phi = 0.0
    endif
      
    if(sin2phi < 0) then
        sin2phi = 0.0
    endif 
    if(sin2phi > 1) then
        sin2phi = 1.0
    endif 
    sinphi = sqrt(sin2phi)
      
    statev(11) = asin(sinphi)*180.0d0/3.141592d0

    ! cal normalized length of intergranular strain rho (statev 12)
    if (parms(10) <= 0.5) then
        istrain = 0 
    else 
        istrain = 1
    endif

    if(istrain == 1) then
        
        do i = 1, 6
            del(i) = statev(i)
        enddo       
        
        norm_del2 = dot_vect_h(2, del, del, 6)
        norm_del = dsqrt(norm_del2)
        statev(12) = norm_del / parms(12)
     
    else
        statev(12) = 0.0
    endif

    return
endsubroutine

! ========================== Utilities =======================
! copy asv to qq_n
! changes intergranular strain from continuum to soil mechanics convention
subroutine move_asv_h(asv, nasv, qq_n)
! del has 6 components
    implicit none
    integer::nasv,i
    real(8)::asv(nasv), qq_n(nasv)
    do i = 1, nasv
        qq_n(i) = 0.0
    enddo
    ! intergranular strain tensor stored in qq_n(1:6)
    do i = 1, 6
        qq_n(i) = -asv(i)
    enddo
    ! void ratio stored in qq_n(7)
    qq_n(7) = asv(7)
    return
endsubroutine

! copy strain increment dstran into deps and computes volumetric strain increment
subroutine move_eps_h(dstran, ntens, deps, depsv)
! strains negative in compression
! deps has 6 components
    implicit none
    integer::ntens, i
    real(8)::deps(6), dstran(ntens), depsv
    do i = 1, ntens
        deps(i) = dstran(i)
    enddo
    depsv = deps(1) + deps(2) + deps(3)
    return
endsubroutine

! computes effective stress (sig) from total stress (stress) and pore pressure (pore)
subroutine move_sig_h(stress, ntens, pore, sig)
!   stress = total stress tensor (tension positive)
!   pore   = exc. pore pressure (undrained conds., compression positive)
!   sig    = effective stress (tension positive), sig has always 6 components
    implicit none
    integer::ntens, i
    real(8)::sig(6), stress(ntens), pore
    do i = 1, ntens
        if(i <= 3) then
            sig(i) = stress(i) + pore
        else
            sig(i) = stress(i)
        endif
    enddo
    return
endsubroutine

! evaluate norm of residual vector Res = ||y_hat - y_til||
subroutine norm_res_h(y_til, y_hat, ny, nasv, norm_R)
    implicit none
    integer ny,nasv,ng,k,i,testnan
    real(8)::y_til(ny),y_hat(ny),void_til,void_hat,del_void
    real(8)::err(ny),norm_R2,norm_R
    real(8)::norm_sig2,norm_q2,norm_sig,norm_q
    real(8)::sig_hat(6),sig_til(6),del_sig(6)
    real(8)::q_hat(nasv),q_til(nasv),del_q(nasv)
    real(8)::dot_vect_h
    
    ng = 6 * nasv
    k = 42 + nasv
    do i = 1, ny
        err(i) = 0.0
    enddo
    ! recover stress tensor and internal variables
    do i=1,6
        sig_hat(i)=y_hat(i)
        sig_til(i)=y_til(i)
        del_sig(i)=dabs(sig_hat(i)-sig_til(i))
    enddo
    do i=1,nasv-1
        q_hat(i)=y_hat(6+i)
        q_til(i)=y_til(6+i)
        del_q(i)=dabs(q_hat(i)-q_til(i))
    enddo

    void_hat=y_hat(6+nasv)
    void_til=y_til(6+nasv)
    del_void=dabs(void_hat-void_til)

    ! relative error norms
    norm_sig2=dot_vect_h(1,sig_hat,sig_hat,6)
    norm_q2=dot_vect_h(2,q_hat,q_hat,6)
    norm_sig=dsqrt(norm_sig2)
    norm_q=dsqrt(norm_q2)
    if(norm_sig>0.0) then
        do i=1,6
                err(i)=del_sig(i)/norm_sig
        enddo
    endif
    if(norm_q>0.0) then
        do i=1,nasv-1
            err(6+i)=del_q(i)/norm_q
        enddo
    endif
    err(6+nasv)=del_void/void_hat

    ! relative error norm
    norm_R2 = dot_vect_h(3,err,err,ny)
    norm_R = dsqrt(norm_R2)
    return
end

! ======================== math utilities ==========================

! calculate eigenvalues and eigenvectors
subroutine PrnSig_h(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
    Implicit real(8) (A-H,O-Z)
    Dimension S(*),xN1(*),xN2(*),xN3(*)
    if (iOpt == 1) then
        call Eig_3_h(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
    else
        call Eig_3a_h(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
    endif
    return
endsubroutine
      
! get Eigenvalues/Eigenvectors for 3*3 matrix
subroutine Eig_3_h(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
! stress vector St(): XX, YY, ZZ, XY, YZ, ZX
    implicit real(8) (A-H,O-Z)
    Dimension St(6),A(3,3),V(3,3), xN1(3),xN2(3),xN3(3)
    A(1,1) = St(1) ! xx
    A(1,2) = St(4) ! xy = yx
    A(1,3) = St(6) ! zx = xz
    A(2,1) = St(4) ! xy = yx
    A(2,2) = St(2) ! yy
    A(2,3) = St(5) ! zy = yz
    A(3,1) = St(6) ! zx = xz
    A(3,2) = St(5) ! zy = yz
    A(3,3) = St(3) ! zz

    ! Set V to unity matrix
    V(1,1) = 1
    V(2,1) = 0
    V(3,1) = 0
    V(1,2) = 0
    V(2,2) = 1
    V(3,2) = 0
    V(1,3) = 0
    V(2,3) = 0
    V(3,3) = 1

    abs_max_s = 0.0
    do i = 1, 3
        do j = 1, 3
            if (abs(a(i,j)) > abs_max_s) abs_max_s = abs(a(i,j))
        enddo
    enddo
    Tol = 1d-20 * abs_max_s
    it = 0
    itmax = 50
    Do While (it < itMax .and. abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) > Tol)
    it=it+1
    Do k=1,3
        If (k == 1) Then
        ip=1
        iq=2
        Else If (k ==2) Then
        ip=2
        iq=3
        Else
        ip=1
        iq=3
        endif
        If (abs(a(ip,iq)) > Tol) Then
        tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
        If (tau >=0.0) Then
            sign_tau=1.0
        Else
            sign_tau=-1.0
        endif
        t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
        c=1.0/sqrt(1.0+t*t)
        s=t*c
        a1p=c*a(1,ip)-s*a(1,iq)
        a2p=c*a(2,ip)-s*a(2,iq)
        a3p=c*a(3,ip)-s*a(3,iq)
        a(1,iq)=s*a(1,ip)+c*a(1,iq)
        a(2,iq)=s*a(2,ip)+c*a(2,iq)
        a(3,iq)=s*a(3,ip)+c*a(3,iq)
        a(1,ip)=a1p
        a(2,ip)=a2p
        a(3,ip)=a3p

        v1p=c*v(1,ip)-s*v(1,iq)
        v2p=c*v(2,ip)-s*v(2,iq)
        v3p=c*v(3,ip)-s*v(3,iq)
        v(1,iq)=s*v(1,ip)+c*v(1,iq)
        v(2,iq)=s*v(2,ip)+c*v(2,iq)
        v(3,iq)=s*v(3,ip)+c*v(3,iq)
        v(1,ip)=v1p
        v(2,ip)=v2p
        v(3,ip)=v3p

        ap1=c*a(ip,1)-s*a(iq,1)
        ap2=c*a(ip,2)-s*a(iq,2)
        ap3=c*a(ip,3)-s*a(iq,3)
        a(iq,1)=s*a(ip,1)+c*a(iq,1)
        a(iq,2)=s*a(ip,2)+c*a(iq,2)
        a(iq,3)=s*a(ip,3)+c*a(iq,3)
        a(ip,1)=ap1
        a(ip,2)=ap2
        a(ip,3)=ap3
        endif ! a(ip,iq)<>0
    enddo ! k
    enddo ! While
    ! principal values on diagonal of a
    S1 = a(1,1)
    S2 = a(2,2)
    S3 = a(3,3)
    ! Derived invariants
    P = (S1+S2+S3)/3
    Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

    ! Sort eigenvalues S1 <= S2 <= S3
    is1 = 1
    is2 = 2
    is3 = 3
    if (s1 > s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
    endif
    if (s2 > s3) then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
    endif
    if (s1 > s2) then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
    endif
    do i = 1, 3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
    enddo
    return
endsubroutine

! Get Eigenvalues (no Eigenvectors) for 3*3 matrix
Subroutine Eig_3a_h(iOpt,St,S1,S2,S3,P,Q)
! Stress vector XX, YY, ZZ, XY, YZ, ZX
    Implicit real(8) (A-H,O-Z)
    Dimension St(6),A(3,3)
    
    A(1,1) = St(1) ! xx
    A(1,2) = St(4) ! xy = yx
    A(1,3) = St(6) ! zx = xz

    A(2,1) = St(4) ! xy = yx
    A(2,2) = St(2) ! yy
    A(2,3) = St(5) ! zy = yz

    A(3,1) = St(6) ! zx = xz
    A(3,2) = St(5) ! zy = yz
    A(3,3) = St(3) ! zz

    abs_max_s=0.0
    Do i=1,3
    Do j=1,3
        if (abs(a(i,j)) > abs_max_s) abs_max_s=abs(a(i,j))
    enddo
    enddo
    Tol = 1d-20 * abs_max_s
    If (iOpt==1) Tol = 1d-50*abs_max_s
    it=0
    itmax = 50

    do while (it < itmax .and. &
        abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) > Tol)

    it=it+1
    do k = 1,3
        If (k == 1) then
            ip=1
            iq=2
        elseIf (k == 2) then
            ip=2
            iq=3
        else
            ip=1
            iq=3
        endif

        If (abs(a(ip,iq)) > Tol) then ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau >=0.0) Then
                sign_tau=1.0
            Else
                sign_tau=-1.0
            endif
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
        endif ! a(ip,iq)<>0
    enddo ! k
    enddo ! While
    ! principal values on diagonal of a
    S1 = a(1,1)
    S2 = a(2,2)
    S3 = a(3,3)
    ! Derived invariants
    P = (S1+S2+S3)/3
    Q = dsqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )
    
    if (s1>s2) Then
        t   = s2
        s2  = s1
        s1  = t
    endif
    if (s2>s3) Then
        t   = s3
        s3  = s2
        s2  = t
    endif
        if (s1>s2) Then
        t   = s2
        s2  = s1
        s1  = t
    endif
    return
endsubroutine
      
! integration of linear elasticity
subroutine calc_elasti_h(y, n, nasv, dtsub, err_tol, maxnint, DTmin, deps_np1, &
    parms, nparms, nfev, elprsw, dtime, DDtan, youngel, nuel, error)
    implicit none
    logical::elprsw
    integer::n,nasv,nparms,i,ksubst,kreject,nfev
    integer::maxnint,error,error_RKF,tension,j
    real(8)::y(n),parms(nparms),dtsub,err_tol,DTmin
    real(8)::deps_np1(6),y_k(n),y_2(n),y_3(n),y_til(n)
    real(8)::y_hat(n),DDtan(6,6)
    real(8)::T_k,DT_k,dtime,II(6,6),krondelta(6)
    real(8)::kRK_1(n),kRK_2(n),kRK_3(n)
    real(8)::norm_R,S_hull,youngel,nuel,F_sig(6)
    
    ! initialize y_k vector and other variables
    do i = 1, n
        y_k(i) = 0.0
    enddo
    
    ! fourth order identity tensors in Voigt notation
    do i = 1, 6
        do j = 1, 6
            II(i,j) = 0.0
        enddo
    enddo
    II(1,1) = 1.0
    II(2,2) = 1.0
    II(3,3) = 1.0
    II(4,4) = 0.5
    II(5,5) = 0.5
    II(6,6) = 0.5

    krondelta(1) = 1.0
    krondelta(2) = 1.0
    krondelta(3) = 1.0
    krondelta(4) = 0.0
    krondelta(5) = 0.0
    krondelta(6) = 0.0

    ! elastic stiffness tensor 
    do i = 1,6
        do j = 1, 6
            DDtan(i,j) = (youngel/(1.0+nuel)) * (II(i,j) + nuel/(1.0-2.0*nuel) * krondelta(i) * krondelta(j));
        enddo
    enddo
    call matmul_h(DDtan, deps_np1, F_sig, 6, 6, 1)
    do i = 1, 6
        y(i) = y(i) + F_sig(i)
    enddo
    
    return
endsubroutine

! copy a(n) to b(n)
subroutine push_h(a, b, n)
    implicit none
    integer::i, n
    real(8)::a(n), b(n)
    do i = 1, n
        b(i) = a(i)
    enddo
    return
endsubroutine

! dot product of a 2nd order tensor, stored in Voigt notation
real(8) function dot_vect_h(flag, a, b, n)
! flag meaning:
!   1 -> stress in Voigt notation
!   2 -> strain in Voigt notation
!   3 -> ordinary dot product between R^n vectors
    implicit none
    integer::i, n, flag
    real(8)::a(n), b(n), coeff

    if (flag == 1) then ! stress tensor
        coeff = 2.0
    elseif(flag==2) then ! strain tensor
        coeff = 0.5
    else ! standard vectors
        coeff = 1.0
    endif
      
    dot_vect_h = 0.0
    do i = 1, n
        if(i <= 3) then
            dot_vect_h = dot_vect_h + a(i) * b(i)
        else
            dot_vect_h = dot_vect_h + coeff * a(i) * b(i)
        endif
    enddo
    return
endfunction

! matrix multiplication
subroutine matmul_h(a, b, c, l, m, n)
    implicit none
    integer::i, j, k, l, m, n
    real(8)::a(l, m),b(m, n),c(l, n)
    do i = 1, l
        do j = 1, n
            do k = 1, m
                c(i, j) = c(i, j) + a(i, k) * b(k, j)
            enddo
        enddo
    enddo
    return
endsubroutine
