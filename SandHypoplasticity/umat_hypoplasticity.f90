! user subroutine for Abaqus
!	Author: D. Masin, based on RKF23 implementation by C. Tamagnini
recursive subroutine umat_hypoplasticity(stress, statev, ddsdde, sse, spd, &
    scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
    predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, nprops, coords, &
    drot, pnewdt, celent, dfgrd0, dfgrd1, noel, npt, layer, kspt, kstep, kinc)
    implicit None
    
    character *80 cmname
    integer ntens, ndi, nshr, nstatv, nprops, noel, npt, layer, kspt, kstep, kinc, inittension

    Double Precision stress(ntens), statev(nstatv), ddsdde(ntens, ntens), ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens), time(2), predef(1), dpred(1), props(nprops), coords(3), drot(3, 3), dfgrd0(3, 3), dfgrd1(3, 3)
    Double Precision sse, spd, scd, rpl, drpldt, dtime, temp, dtemp, pnewdt, celent

! ... 1. nasvdim    = maximum number of additional state variables
!     2. tolintT    = prescribed error tolerance for the adaptive
!                     substepping scheme
!     3. maxnint    = maximum number of time substeps allowed.
!                     If the limit is exceeded abaqus is forced to reduce
!                     the overall time step size (cut-back)
!     4. DTmin      = minimum substeps size allowed.
!                     If the limit is exceeded abaqus is forced to reduce
!                     the overall time step size (cut-back)
!     5. perturb    = perturbation parameter for numerical computation of Jacobian matrices
!     6. nfasv      = number of first additional state variable in statev field
!     7. prsw       = switch for printing information
!
! ... declaration of local variables
  Logical prsw, elprsw
  integer i, error, maxnint, nfev, testnan, maxninttest
  integer nparms, nasvdim, nfasv, nydim, nasv, nyact, testing
!
  Double Precision dot_vect_h
!
  Double Precision parms(nprops), theta, tolintt, dtsub, dtmin, perturb
  Double Precision sig_n(6), sig_np1(6), ddtan(6, 6), pore
  Double Precision deps_np1(6), depsv_np1, norm_deps, tolintttest
  Double Precision norm_deps2, pp, qq, cos3t, i1, i2, i3, norm_d2, norm_d
  Double Precision youngel, tdepel0, tdepel1, nuel
  Double Precision eyoung0, eyoung1, nu0, nu1

  Parameter (nasvdim=15)
  Parameter (nydim=6+nasvdim)
  Parameter (tolintt=1.0D-1)
  Parameter (tolintttest=1.0D-1)
  Parameter (maxnint=200)
  Parameter (maxninttest=1000)
  Parameter (dtmin=1.0D-17)
  Parameter (perturb=1.0D-5)
  Parameter (nfasv=1)
  Parameter (prsw=.False.)

! ... additional state variables
  Double Precision asv(nasvdim)

! ... solution vector (stresses, additional state variables)
  Double Precision y(nydim), y_n(nydim), dy(nydim)

! ... Error Management:
!     error =  0 ... no problem in time integration
!     error =  1 ... problems in evaluation of the time rate, (e.g. undefined
!                    stress state), reduce time integration substeps
!     error =  3 ... problems in time integration, reduce abaqus load increment
!                    (cut-back)
!     error = 10 ... severe error, terminate calculation
  error = 0

! ... check material parameters and move them to array parms(nparms)
  Call check_parms_h(props, nprops, parms, nparms, error)

  elprsw = .False. 
! ... define number of additional state variables
! not used                                       
  nasv = 7
  nyact = 6 + nasv

! ... suggested time substep size, and initial excess pore pressure
  dtsub = statev(13)
  pore = -statev(8)

! ... init void ratio
  If (statev(7)<0.001) Then
    statev(7) = props(16)
  End If

! ... additional state variables
  Do i = 1, nasv
    asv(i) = statev(i-1+nfasv)
  End Do

! ... compute volume strain increment and effective stress
  Do i = 1, 6
    sig_n(i) = 0.0
    deps_np1(i) = 0.0
  End Do
  Call move_sig_h(stress, ntens, pore, sig_n)
  Call move_eps_h(dstran, ntens, deps_np1, depsv_np1)

  norm_d2 = dot_vect_h(2, deps_np1, deps_np1, 6)
  norm_d = sqrt(norm_d2)

! ... Time integration
  Call iniy_h(y, nydim, nasv, ntens, sig_n, asv)
  Call push_h(y, y_n, nydim)

! ... check whether the initial state is not tensile
  inittension = 0
  Call check_rkf_h(inittension, y, nyact, nasv, parms, nparms)

! ... Switch for elasticity in the case tensile stress is reached
  youngel = 1.0

! ... local integration using adaptive RKF23 method, consistent Jacobian and error estimation
  If ((dtsub<=0.0D0) .Or. (dtsub>dtime)) Then
    dtsub = dtime
  End If

  testing = 0
! for ddsdde only
  If (norm_d==0.0) testing = 2

  nfev = & 
    0
! initialisation                                         
  If (inittension==0) Then

    If (testing==2) Then
! no need to update stress, for ddsdde only
      Do i = 1, nyact
        y(i) = y_n(i)
      End Do
! ... Normal RKF23 integration, main integration
    Else
! testing == 0                                         
      Call rkf23_update_h(y, nyact, nasv, dtsub, tolintt, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, dtime, error)
    End If

    If (error==3) Then
! keep stress the same, we are likely close to tensile region
      Do i = 1, nyact
        y(i) = y_n(i)
      End Do
!             write(6, *) 'tensile state'
    Else If (error==10) Then
      Call wrista_h(2, y, nydim, deps_np1, dtime, coords, statev, nstatv, parms, nparms, noel, npt, ndi, nshr, kstep, kinc)
      Call xit_h

    End If
! ... compute ddsdde
! end error                                             
    Call perturbate_h(y_n, y, nyact, nasv, dtsub, tolintt, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, theta, ntens, ddtan, dtime, error)

! in tensile state, calc elastic
  Else
! inittension == 1                                           
    youngel = 1000.0
    nuel = 0.3
    Call calc_elasti_h(y, nyact, nasv, dtsub, tolintt, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, dtime, ddtan, youngel, nuel, error)

  End If
! ... update dtsub and nfev
! end inittension                                           
  If (dtsub<=0.0D0) Then
    dtsub = 0.0D0
  Else If (dtsub>=dtime) Then
    dtsub = dtime
  End If
  statev(13) = dtsub
  statev(10) = dfloat(nfev)

! ... convert solution (stress + cons. tangent) to abaqus format
!     update pore pressure and compute total stresses
  Call solout_h(stress, ntens, asv, nasv, ddsdde, y, nydim, pore, depsv_np1, parms, nparms, ddtan)

! ... updated vector of additional state variables to abaqus statev vector
  Do i = 1, nasv
    statev(i-1+nfasv) = asv(i)
  End Do

! ... transfer additional information to statev vector
  Do i = 1, 6
    sig_np1(i) = y(i)
  End Do
  pp = -(sig_np1(1)+sig_np1(2)+sig_np1(3))/3.0

  statev(8) = -pore
  statev(9) = pp

! ... update statev variables
  If (inittension==0) Then
    Call calc_statev_h(sig_np1, statev, parms, nparms, nasv, nasvdim, deps_np1)
  End If

! ... complete time integration
  Return
End Subroutine umat_hypoplasticity

! ---------------------------- helper functions ------------------------------
!-----------------------------------------------------------------------------
recursive Subroutine check_parms_h(props, nprops, parms, nparms, error)
!-----------------------------------------------------------------------------
! checks input material parameters
!
! written 10/2004 (Tamagnini & Sellari)
!-----------------------------------------------------------------------------
  implicit None

  integer nprops, nparms, i, error
  Double Precision props(nprops), parms(nprops)
  Double Precision phi_deg, phi, hs, en, ed0, ec0, ei0, alpha, beta
  Double Precision m_r, m_t, r_uc, beta_r, chi, bulk_w, p_t

  Double Precision zero, one, four, pi, pi_deg
  Parameter (zero=0.0D0, one=1.0D0, four=4.0D0, pi_deg=180.0D0)

  pi = four*datan(one)

  nparms = nprops
  Do i = 1, nprops
    parms(i) = props(i)
  End Do

! ... recover material parameters
  phi_deg = parms(1)
  hs = parms(3)
  en = parms(4)
  ed0 = parms(5)
  ec0 = parms(6)
  ei0 = parms(7)
  alpha = parms(8)
  beta = parms(9)
  m_r = parms(10)
  m_t = parms(11)
  r_uc = parms(12)
  beta_r = parms(13)
  chi = parms(14)
  bulk_w = parms(15)
  p_t = parms(2)

  phi = phi_deg*pi/pi_deg
  parms(1) = phi

  If (phi<=zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'phi = ', phi
    error = 10
    Return
  End If

  If (m_r<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'm_R = ', m_r
    error = 10
    Return
  End If

  If (m_t<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'm_T = ', m_t
    error = 10
    Return
  End If

  If (r_uc<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'r_uc = ', r_uc
    error = 10
    Return
  End If

  If (beta_r<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'beta_r = ', beta_r
    error = 10
    Return
  End If

  If (chi<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'chi = ', chi
    error = 10
    Return
  End If

  If (bulk_w<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'bulk_w = ', bulk_w
    error = 10
    Return
  End If

  If (p_t<zero) Then
    Write (6, *) 'ERROR: subroutine CHECK_PARMS:'
    Write (6, *) 'p_t = ', p_t
    error = 10
    Return
  End If

  Return
End Subroutine check_parms_h

!------------------------------------------------------------------------------
recursive Double Precision Function dot_vect_h(flag, a, b, n)
!------------------------------------------------------------------------------
! dot product of a 2nd order tensor, stored in Voigt notation
! created 10/2004 (Tamagnini & Sellari)
!
! flag = 1 -> vectors are stresses in Voigt notation
! flag = 2 -> vectors are strains in Voigt notation
! flag = 3 -> ordinary dot product between R^n vectors
!------------------------------------------------------------------------------
  implicit None
  integer i, n, flag
  Double Precision a(n), b(n)
  Double Precision zero, half, one, two, coeff
  Parameter (zero=0.0D0, half=0.5D0, one=1.0D0, two=2.0D0)

  If (flag==1) Then
! ... stress tensor (or the like)
    coeff = two
  Else If (flag==2) Then
! ... strain tensor (or the like)
    coeff = half
  Else
! ... standard vectors
    coeff = one
  End If

  dot_vect_h = zero
  Do i = 1, n
    If (i<=3) Then
      dot_vect_h = dot_vect_h + a(i)*b(i)
    Else
      dot_vect_h = dot_vect_h + coeff*a(i)*b(i)
    End If
  End Do

  Return
End Function dot_vect_h

!-----------------------------------------------------------------------------
recursive Subroutine get_f_sig_q_h(sig, q, nasv, parms, nparms, deps, f_sig, f_q, error)
!-----------------------------------------------------------------------------
!
!  finds vectors F_sigma and F_q in F(y)
!
!  written 6/2005 (Tamagnini, Sellari & Miriano)
!-----------------------------------------------------------------------------
  implicit None
  Double Precision dot_vect_h
  integer nparms, nasv, ii
  Double Precision sig(6), q(nasv), parms(nparms), deps(6)
  Double Precision mm(6, 6), hh(nasv, 6), f_sig(6), f_q(nasv)
  Double Precision ll(6, 6), nn(6), norm_d, norm_d2
  integer istrain, error

! ... compute tangent operators
  If (parms(10)<=0.5) Then
    istrain = 0
  Else
    istrain = 1
  End If

  Call get_tan_h(deps, sig, q, nasv, parms, nparms, mm, hh, ll, nn, istrain, error)

! ... compute F_sig=MM*deps
  If (istrain==1) Then
    Call matmul_h(mm, deps, f_sig, 6, 6, 1)
  Else
    Call matmul_h(ll, deps, f_sig, 6, 6, 1)
    norm_d2 = dot_vect_h(2, deps, deps, 6)
    norm_d = sqrt(norm_d2)
    Do ii = 1, 6
      f_sig(ii) = f_sig(ii) + nn(ii)*norm_d
    End Do
  End If

! ... compute F_q=HH*deps
  Call matmul_h(hh, deps, f_q, nasv, 6, 1)

  Return
End Subroutine get_f_sig_q_h

!-----------------------------------------------------------------------------
recursive Subroutine get_tan_h(deps, sig, q, nasv, parms, nparms, mm, hh, ll, nn, istrain, error)
!-----------------------------------------------------------------------------
!  computes matrices M and H for Masin hypoplastic model for clays
!  version with intergranular strains
!  NOTE: stress and strain convention: tension and extension positive
!  written 6/2005 (Tamagnini & Sellari)
!-----------------------------------------------------------------------------
  implicit None

  integer nparms, nasv, i, j, error
  Double Precision dot_vect_h
  Double Precision sig(6), q(nasv), parms(nparms), deps(6)
  Double Precision eta(6), eta_dev(6), del(6), void, sig_star(6)
  Double Precision eta_del(6), eta_delta(6), eta_eps(6)
  Double Precision norm_del, norm_del2, norm_deps, norm_deps2, eta_dn2
  Double Precision pp, qq, cos3t, i1, i2, i3, tanpsi
  Double Precision a, a2, ff, fd, fs
  Double Precision num, den, af, fa2, eta_n2, norm_m, norm_m2
  Double Precision ii(6, 6), iu(6, 6)
  Double Precision mm(6, 6), hh(nasv, 6), ll(6, 6), nn(6), aa(6, 6), m(6)
  integer istrain
  Double Precision m_dir(6), m_dir1(6), leta(6), h_del(6, 6), h_e(6)
  Double Precision load, rho
  Double Precision temp1, temp2, temp3, temp4
  Double Precision phi, hs, en, ed0, ec0, ei0, alpha, beta, r_uc
  Double Precision m_r, m_t, beta_r, chi, bulk_w, p_t, sinphi, sinphi2
  Double Precision ec, ed, ei, bauer, fb, fe, sq2, sq3, sq6, az

  Double Precision zero, tiny, half, one, two, three, six, eight, nine
  Double Precision onethird, sqrt3, twosqrt2, sqrt2, oneeight, ln2m1
  Parameter (zero=0.0D0, one=1.0D0, two=2.0D0, three=3.0D0, six=6.0D0)
  Parameter (tiny=1.0D-17, half=0.5D0, eight=8.0D0, nine=9.0D0)
  Parameter (sq2=1.4142135623730951455D0, sq3=1.7320508075688771931D0, sq6=2.4494897427831778813D0)

  onethird = one/three
  sqrt3 = dsqrt(three)
  twosqrt2 = two*dsqrt(two)
  sqrt2 = dsqrt(two)
  oneeight = one/eight
  onethird = one/three
  ln2m1 = one/dlog(two)

  Do i = 1, 6
    Do j = 1, 6
      mm(i, j) = zero
      ll(i, j) = zero
      ii(i, j) = zero
      iu(i, j) = zero
      h_del(i, j) = zero
    End Do
    eta_del(i) = zero
    eta_delta(i) = zero
    eta_eps(i) = zero
  End Do
!
  Do i = 1, nasv
    Do j = 1, 6
      hh(i, j) = zero
    End Do
  End Do

! ... fourth order identity tensors in Voigt notation
  ii(1, 1) = one
  ii(2, 2) = one
  ii(3, 3) = one
  ii(4, 4) = half
  ii(5, 5) = half
  ii(6, 6) = half

  iu(1, 1) = one
  iu(2, 2) = one
  iu(3, 3) = one
  iu(4, 4) = one
  iu(5, 5) = one
  iu(6, 6) = one

! ... recover material parameters
  phi = parms(1)
  hs = parms(3)
  en = parms(4)
  ed0 = parms(5)
  ec0 = parms(6)
  ei0 = parms(7)
  alpha = parms(8)
  beta = parms(9)
  m_r = parms(10)
  m_t = parms(11)
  r_uc = parms(12)
  beta_r = parms(13)
  chi = parms(14)
  bulk_w = parms(15)
  p_t = parms(2)

  sinphi = dsin(phi)
  sinphi2 = sinphi*sinphi

! ... recover internal state variables
  del(1) = q(1)
  del(2) = q(2)
  del(3) = q(3)
  del(4) = q(4)
  del(5) = q(5)
  del(6) = q(6)
  void = q(7)

! ... axis translation due to cohesion (p_t>0)
  sig_star(1) = sig(1) - p_t
  sig_star(2) = sig(2) - p_t
  sig_star(3) = sig(3) - p_t
  sig_star(4) = sig(4)
  sig_star(5) = sig(5)
  sig_star(6) = sig(6)

! ... strain increment and intergranular strain directions
  norm_deps2 = dot_vect_h(2, deps, deps, 6)
  norm_del2 = dot_vect_h(2, del, del, 6)
  norm_deps = dsqrt(norm_deps2)
  norm_del = dsqrt(norm_del2)

  If (norm_del>=tiny) Then
    Do i = 1, 6
      eta_del(i) = del(i)/norm_del
    End Do
  End If

  eta_delta(1) = eta_del(1)
  eta_delta(2) = eta_del(2)
  eta_delta(3) = eta_del(3)
  eta_delta(4) = half*eta_del(4)
  eta_delta(5) = half*eta_del(5)
  eta_delta(6) = half*eta_del(6)

  If (norm_deps>=tiny) Then
    Do i = 1, 6
      eta_eps(i) = deps(i)/norm_deps
    End Do
  End If

! ... auxiliary stress tensors
  Call inv_sig_h(sig_star, pp, qq, cos3t, i1, i2, i3)
!     T_cap
  eta(1) = sig_star(1)/i1
  eta(2) = sig_star(2)/i1
  eta(3) = sig_star(3)/i1
  eta(4) = sig_star(4)/i1
  eta(5) = sig_star(5)/i1
  eta(6) = sig_star(6)/i1
!     T_cap_dot
  eta_dev(1) = eta(1) - onethird
  eta_dev(2) = eta(2) - onethird
  eta_dev(3) = eta(3) - onethird
  eta_dev(4) = eta(4)
  eta_dev(5) = eta(5)
  eta_dev(6) = eta(6)

! ... functions a and F
  eta_dn2 = dot_vect_h(1, eta_dev, eta_dev, 6)
  tanpsi = sqrt3*dsqrt(eta_dn2)
  temp1 = oneeight*tanpsi*tanpsi + (two-tanpsi*tanpsi)/(two+sqrt2*tanpsi*cos3t)
  temp2 = tanpsi/twosqrt2

  a = sqrt3*(three-sin(phi))/(twosqrt2*sin(phi))
  a2 = a*a
  ff = dsqrt(temp1) - temp2

! ... barotropy and pyknotropy functions
  bauer = dexp(-(-i1/hs)**en)
  ed = ed0*bauer
  ec = ec0*bauer
  ei = ei0*bauer

!     fb
  temp1 = three + a*a - a*sq3*((ei0-ed0)/(ec0-ed0))**alpha
  If (temp1<zero) Stop 'factor fb not defined'
  fb = hs/en/temp1*(one+ei)/ei*(ei0/ec0)**beta*(-i1/hs)**(one-en)

!     fe
  fe = (ec/void)**beta

  fs = fb*fe

!     fd
  If (void>=ed) Then
    fd = ((void-ed)/(ec-ed))**alpha
  Else
    fd = 0.0
  End If

! ... tensor L
  eta_n2 = dot_vect_h(1, eta, eta, 6)
  Do i = 1, 6
    Do j = 1, 6
      ll(i, j) = (ii(i,j)*ff*ff+a2*eta(i)*eta(j))/eta_n2
    End Do
  End Do

! ... tensor NN
  Do i = 1, 6
    nn(i) = ff*a*(eta(i)+eta_dev(i))/eta_n2
  End Do

! ... BEGIN intergranular STRAIN
  If (istrain==1) Then

! ... loading function
    load = dot_vect_h(2, eta_del, eta_eps, 6)
! ... intergranular strain--related tensors
    rho = norm_del/r_uc
    If (rho>one) Then
      rho = one
    End If

    Call matmul_h(ll, eta_del, leta, 6, 6, 1)

! ... tangent stiffness M(sig,q,eta_eps)
    temp1 = ((rho**chi)*m_t+(one-rho**chi)*m_r)*fs

    If (load>zero) & ! loading                            
      Then
      temp2 = (rho**chi)*(one-m_t)*fs
      temp3 = (rho**chi)*fs*fd
      Do i = 1, 6
        Do j = 1, 6
          aa(i, j) = temp2*leta(i)*eta_delta(j) + temp3*nn(i)*eta_delta(j)
          mm(i, j) = temp1*ll(i, j) + aa(i, j)
        End Do
      End Do
! reverse loading                                        
    Else
      temp4 = (rho**chi)*(m_r-m_t)*fs
      Do i = 1, 6
        Do j = 1, 6
          aa(i, j) = temp4*leta(i)*eta_delta(j)
          mm(i, j) = temp1*ll(i, j) + aa(i, j)
        End Do
      End Do
    End If

! ... intergranular strain evolution function
!     NOTE: H_del transforms a strain-like vector into a strain-like vector
!           eta_del(i) instead of eta_delta(i)
!           I = 6x6 unit matrix
    If (load>zero) Then
      Do i = 1, 6
        Do j = 1, 6
          h_del(i, j) = iu(i, j) - (rho**beta_r)*eta_del(i)*eta_delta(j)
        End Do
      End Do
    Else
      Do i = 1, 6
        h_del(i, i) = one
      End Do
    End If

! ... void ratio evolution function (tension positive)
    Do i = 1, 6
      If (i<=3) Then
        h_e(i) = one + void
      Else
        h_e(i) = zero
      End If
    End Do

! ... assemble hardening matrix
!     for clay hypoplasticity only, can be ignored for sand
    Do i = 1, nasv
      If (i<=6) Then
        Do j = 1, 6
          hh(i, j) = h_del(i, j)
        End Do
      Else
        Do j = 1, 6
          hh(i, j) = h_e(j)
        End Do
      End If
    End Do

! ... not accouting for intergranular strain
  Else If (istrain==0) Then

    Do i = 1, 6
      If (i<=3) Then
        h_e(i) = one + void
      Else
        h_e(i) = zero
      End If
    End Do

    Do i = 1, nasv
      If (i<=6) Then
        Do j = 1, 6
          hh(i, j) = 0.0
        End Do
      Else
        Do j = 1, 6
          hh(i, j) = h_e(j)
        End Do
      End If
    End Do

! ... end istrain/noistrain switch
  End If

  Do i = 1, 6

    Do j = 1, 6
      ll(i, j) = ll(i, j)*fs
    End Do

    nn(i) = nn(i)*fs*fd

  End Do

  Return
End Subroutine get_tan_h

!-----------------------------------------------------------------------------
recursive Subroutine iniy_h(y, nydim, nasv, ntens, sig, qq)
!-----------------------------------------------------------------------------
! initializes the vector of state variables
!-----------------------------------------------------------------------------
  implicit None
  integer i, nydim, nasv, ntens
  Double Precision y(nydim), qq(nasv), sig(ntens)

  Do i = 1, nydim
    y(i) = 0.0
  End Do

  Do i = 1, ntens
    y(i) = sig(i)
  End Do

!     additional state variables
  Do i = 1, nasv
    y(6+i) = qq(i)
  End Do

  Return
End Subroutine iniy_h

!------------------------------------------------------------------------------
recursive Subroutine inv_eps_h(eps, eps_v, eps_s, sin3t)
!------------------------------------------------------------------------------
! calculate invariants of strain tensor
!------------------------------------------------------------------------------
  implicit None
  integer i
  Double Precision eps(6), edev(6), edev2(6), ev3
  Double Precision tredev3, eps_v, eps_s, sin3t
  Double Precision norm2, numer, denom
  Double Precision zero, one, two, three, six
  Double Precision onethird, twothirds, sqrt6
  Data zero, one, two, three, six/0.0D0, 1.0D0, 2.0D0, 3.0D0, 6.0D0/
  onethird = one/three
  twothirds = two/three
  sqrt6 = dsqrt(six)

! ... volumetric strain
  eps_v = eps(1) + eps(2) + eps(3)
  ev3 = onethird*eps_v

! ... deviator strain
  edev(1) = eps(1) - ev3
  edev(2) = eps(2) - ev3
  edev(3) = eps(3) - ev3
  edev(4) = eps(4)/two
  edev(5) = eps(5)/two
  edev(6) = eps(6)/two

! ... second invariant
  norm2 = edev(1)*edev(1) + edev(2)*edev(2) + edev(3)*edev(3) + two*(edev(4)*edev(4)+edev(5)*edev(5)+edev(6)*edev(6))
  eps_s = dsqrt(twothirds*norm2)

! ... components of (edev_ij)(edev_jk)
  edev2(1) = edev(1)*edev(1) + edev(4)*edev(4) + edev(5)*edev(5)
  edev2(2) = edev(4)*edev(4) + edev(2)*edev(2) + edev(6)*edev(6)
  edev2(3) = edev(6)*edev(6) + edev(5)*edev(5) + edev(3)*edev(3)
  edev2(4) = two*(edev(1)*edev(4)+edev(4)*edev(2)+edev(6)*edev(5))
  edev2(5) = two*(edev(5)*edev(1)+edev(6)*edev(4)+edev(3)*edev(5))
  edev2(6) = two*(edev(4)*edev(5)+edev(2)*edev(6)+edev(6)*edev(3))

! ... Lode angle
  If (eps_s==zero) Then
    sin3t = -one
  Else
    tredev3 = zero
    Do i = 1, 6
      tredev3 = tredev3 + edev(i)*edev2(i)
    End Do
    numer = sqrt6*tredev3
    denom = (dsqrt(norm2))**3
    sin3t = numer/denom
    If (dabs(sin3t)>one) Then
      sin3t = sin3t/dabs(sin3t)
    End If
  End If

  Return
End Subroutine inv_eps_h

!------------------------------------------------------------------------------
recursive Subroutine inv_sig_h(sig, pp, qq, cos3t, i1, i2, i3)
!------------------------------------------------------------------------------
! calculate invariants of stress tensor
! NOTE: Voigt notation is used with the following index conversion
!       11 -> 1
!       22 -> 2
!       33 -> 3
!       12 -> 4
!       13 -> 5
!       23 -> 6
!------------------------------------------------------------------------------
  implicit None
!
  Double Precision sig(6), sdev(6)
  Double Precision eta(6), eta_d(6), eta_d2(6)
  Double Precision xmin1, xmin2, xmin3
  Double Precision tretadev3, pp, qq, cos3t, i1, i2, i3
  Double Precision norm2, norm2sig, norm2eta, numer, denom
!
  Double Precision half, one, two, three, six
  Double Precision onethird, threehalves, sqrt6, tiny
!
  Double Precision dot_vect_h
!
  Data half, one/0.5D0, 1.0D0/
  Data two, three, six/2.0D0, 3.0D0, 6.0D0/
  Data tiny/1.0D-18/
!
! ... some constants
!
  onethird = one/three
  threehalves = three/two
  sqrt6 = dsqrt(six)
!
! ... trace and mean stress
!
  i1 = sig(1) + sig(2) + sig(3)
  pp = onethird*i1
!
! ... deviator stress
!
  sdev(1) = sig(1) - pp
  sdev(2) = sig(2) - pp
  sdev(3) = sig(3) - pp
  sdev(4) = sig(4)
  sdev(5) = sig(5)
  sdev(6) = sig(6)
!
! ... normalized stress and dev. normalized stress
!

  If (i1/=0) Then
    eta(1) = sig(1)/i1
    eta(2) = sig(2)/i1
    eta(3) = sig(3)/i1
    eta(4) = sig(4)/i1
    eta(5) = sig(5)/i1
    eta(6) = sig(6)/i1
  Else
    eta(1) = sig(1)/tiny
    eta(2) = sig(2)/tiny
    eta(3) = sig(3)/tiny
    eta(4) = sig(4)/tiny
    eta(5) = sig(5)/tiny
    eta(6) = sig(6)/tiny
  End If
!
  eta_d(1) = eta(1) - onethird
  eta_d(2) = eta(2) - onethird
  eta_d(3) = eta(3) - onethird
  eta_d(4) = eta(4)
  eta_d(5) = eta(5)
  eta_d(6) = eta(6)
!
! ... second invariants
!
  norm2 = dot_vect_h(1, sdev, sdev, 6)
  norm2sig = dot_vect_h(1, sig, sig, 6)
  norm2eta = dot_vect_h(1, eta_d, eta_d, 6)
!
  qq = dsqrt(threehalves*norm2)
  i2 = half*(norm2sig-i1*i1)
!
! ... components of (eta_d_ij)(eta_d_jk)
!
  eta_d2(1) = eta_d(1)*eta_d(1) + eta_d(4)*eta_d(4) + eta_d(5)*eta_d(5)
  eta_d2(2) = eta_d(4)*eta_d(4) + eta_d(2)*eta_d(2) + eta_d(6)*eta_d(6)
  eta_d2(3) = eta_d(6)*eta_d(6) + eta_d(5)*eta_d(5) + eta_d(3)*eta_d(3)
  eta_d2(4) = eta_d(1)*eta_d(4) + eta_d(4)*eta_d(2) + eta_d(6)*eta_d(5)
  eta_d2(5) = eta_d(5)*eta_d(1) + eta_d(6)*eta_d(4) + eta_d(3)*eta_d(5)
  eta_d2(6) = eta_d(4)*eta_d(5) + eta_d(2)*eta_d(6) + eta_d(6)*eta_d(3)
!
! ... Lode angle
!
  If (norm2eta<tiny) Then
!
    cos3t = -one
!
  Else
!
    tretadev3 = dot_vect_h(1, eta_d, eta_d2, 6)
!
    numer = -sqrt6*tretadev3
    denom = (dsqrt(norm2eta))**3
    cos3t = numer/denom
    If (dabs(cos3t)>one) Then
      cos3t = cos3t/dabs(cos3t)
    End If
!
  End If
!
! ... determinant
!
  xmin1 = sig(2)*sig(3) - sig(6)*sig(6)
  xmin2 = sig(4)*sig(3) - sig(6)*sig(5)
  xmin3 = sig(4)*sig(6) - sig(5)*sig(2)
!
  i3 = sig(1)*xmin1 - sig(4)*xmin2 + sig(5)*xmin3

!
  Return
End Subroutine inv_sig_h

!------------------------------------------------------------------------------
recursive Subroutine matmul_h(a, b, c, l, m, n)
!------------------------------------------------------------------------------
! matrix multiplication
!------------------------------------------------------------------------------
  implicit None
  integer i, j, k, l, m, n
  Double Precision a(l, m), b(m, n), c(l, n)

  Do i = 1, l
    Do j = 1, n
      c(i, j) = 0.0D0
      Do k = 1, m
        c(i, j) = c(i, j) + a(i, k)*b(k, j)
      End Do
    End Do
  End Do

  Return
End Subroutine matmul_h

!-----------------------------------------------------------------------------
recursive Subroutine move_asv_h(asv, nasv, qq_n)
!-----------------------------------------------------------------------------
! move internal variables in vector qq_n and changes intergranular strain
! from continuum to soil mechanics convention
! NOTE: del has always 6 components
!
! written 6/2005 (Tamagnini, Sellari & Miriano)
!-----------------------------------------------------------------------------
  implicit None
  integer nasv, i
  Double Precision asv(nasv), qq_n(nasv), zero
  Parameter (zero=0.0D0)

  Do i = 1, nasv
    qq_n(i) = zero
  End Do

! ... intergranular strain tensor stored in qq_n(1:6)
  Do i = 1, 6
    qq_n(i) = -asv(i)
  End Do

! ... void ratio stored in qq_n(7)
  qq_n(7) = asv(7)

  Return
End Subroutine move_asv_h

!-----------------------------------------------------------------------------
recursive Subroutine move_eps_h(dstran, ntens, deps, depsv)
!-----------------------------------------------------------------------------
! Move strain increment dstran into deps and computes
! volumetric strain increment
!
! NOTE: all strains negative in compression; deps has always 6 components
!
! written 7/2005 (Tamagnini, Sellari & Miriano)
!-----------------------------------------------------------------------------
  implicit None
  integer ntens, i
  Double Precision deps(6), dstran(ntens), depsv
!
  Do i = 1, ntens
    deps(i) = dstran(i)
  End Do
!
  depsv = deps(1) + deps(2) + deps(3)
!
  Return
End Subroutine move_eps_h
!-----------------------------------------------------------------------------
recursive Subroutine move_sig_h(stress, ntens, pore, sig)
!-----------------------------------------------------------------------------
! computes effective stress from total stress (stress) and pore pressure (pore)
!
! NOTE: stress = total stress tensor (tension positive)
!         pore   = exc. pore pressure (undrained conds., compression positive)
!         sig    = effective stress (tension positive)
!
!       sig has always 6 components
!
! written 7/2005 (Tamagnini, Sellari & Miriano)
!-----------------------------------------------------------------------------
  implicit None
  integer ntens, i
  Double Precision sig(6), stress(ntens), pore, zero
!
  Parameter (zero=0.0D0)
!
  Do i = 1, 6
    sig(i) = zero
  End Do
!
  Do i = 1, ntens
    If (i<=3) Then
      sig(i) = stress(i) + pore
    Else
      sig(i) = stress(i)
    End If
  End Do
!
  Return
End Subroutine move_sig_h
!-----------------------------------------------------------------------------
recursive Subroutine norm_res_h(y_til, y_hat, ny, nasv, norm_r)
!-----------------------------------------------------------------------------
!  evaluate norm of residual vector Res=||y_hat-y_til||
!
!  written 6/2005 (Tamagnini, Sellari & Miriano)
!-----------------------------------------------------------------------------
  implicit None
!
  integer ny, nasv, ng, k, i, testnan
!
  Double Precision y_til(ny), y_hat(ny), void_til, void_hat, del_void
  Double Precision err(ny), norm_r2, norm_r
  Double Precision norm_sig2, norm_q2, norm_sig, norm_q
  Double Precision sig_hat(6), sig_til(6), del_sig(6)
  Double Precision q_hat(nasv), q_til(nasv), del_q(nasv)
  Double Precision dot_vect_h, zero
!
  Parameter (zero=0.0D0)
!
  ng = 6*nasv
  k = 42 + nasv
!
  Do i = 1, ny
    err(i) = zero
  End Do
!
! ... recover stress tensor and internal variables
!
  Do i = 1, 6
    sig_hat(i) = y_hat(i)
    sig_til(i) = y_til(i)
    del_sig(i) = dabs(sig_hat(i)-sig_til(i))
  End Do
!
  Do i = 1, nasv - 1
    q_hat(i) = y_hat(6+i)
    q_til(i) = y_til(6+i)
    del_q(i) = dabs(q_hat(i)-q_til(i))
  End Do
!
  void_hat = y_hat(6+nasv)
  void_til = y_til(6+nasv)
  del_void = dabs(void_hat-void_til)
!
! ... relative error norms
!
  norm_sig2 = dot_vect_h(1, sig_hat, sig_hat, 6)
  norm_q2 = dot_vect_h(2, q_hat, q_hat, 6)
  norm_sig = dsqrt(norm_sig2)
  norm_q = dsqrt(norm_q2)
!
  If (norm_sig>zero) Then
    Do i = 1, 6
      err(i) = del_sig(i)/norm_sig
    End Do
  End If
!
  If (norm_q>zero) Then
    Do i = 1, nasv - 1
      err(6+i) = del_q(i)/norm_q
    End Do
  End If
!
  err(6+nasv) = del_void/void_hat
!
! ... global relative error norm
!
  norm_r2 = dot_vect_h(3, err, err, ny)
  norm_r = dsqrt(norm_r2)
!
  testnan = 0
  Call umatisnan_h(norm_sig, testnan)
  Call umatisnan_h(norm_q, testnan)
  Call umatisnan_h(void_hat, testnan)
  If (testnan==1) Then
    norm_r = 1.D20
  End If

  Return
End Subroutine norm_res_h

!-----------------------------------------------------------------------------
recursive Subroutine perturbate_h(y_n, y_np1, n, nasv, dtsub, err_tol, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, theta, ntens, dd, dtime, error)
!-----------------------------------------------------------------------------
!  compute numerically consistent tangent stiffness
!  written 12/2005 (Tamagnini)
!-----------------------------------------------------------------------------
  implicit None

  Logical elprsw
  integer ntens, jj, kk, i
  integer n, nasv, nparms, nfev
  integer maxnint, error
  Double Precision y_n(n), y_np1(n), y_star(n), parms(nparms)
  Double Precision dtsub, err_tol, dtmin, dtime
  Double Precision theta, sig(6), q(nasv)
  Double Precision deps_np1(6), deps_star(6)
  Double Precision dsig(6), dd(6, 6), hhtmp(nasv, 6)
  Double Precision ll(6, 6), nn(6)
  integer istrain

  Double Precision zero
  Parameter (zero=0.0D0)

! ... initialize DD and y_star
  If (parms(10)<=0.5) Then
    istrain = 0
  Else
    istrain = 1
  End If

  Do kk = 1, 6
    Do jj = 1, 6
      dd(kk, jj) = zero
    End Do
  End Do

  Do i = 1, 6
    sig(i) = y_n(i)
  End Do

  Do i = 1, nasv
    q(i) = y_n(6+i)
  End Do

  Call push_h(y_n, y_star, n)

  If (error/=10) Then
    Call get_tan_h(deps_np1, sig, q, nasv, parms, nparms, dd, hhtmp, ll, nn, istrain, error)
  End If

  If (istrain==0) Then
    Do kk = 1, 6
      Do jj = 1, 6
        dd(kk, jj) = ll(kk, jj)
      End Do
    End Do
  Else
    Do kk = 1, 6
      Do jj = 1, 6
        dd(kk, jj) = parms(10)*ll(kk, jj)
      End Do
    End Do
  End If

  Return
End Subroutine perturbate_h

!-----------------------------------------------------------------------------
recursive Subroutine push_h(a, b, n)
!-----------------------------------------------------------------------------
! push vector a into vector b
!-----------------------------------------------------------------------------
  implicit None
  integer i, n
  Double Precision a(n), b(n)
  Do i = 1, n
    b(i) = a(i)
  End Do
  Return
End Subroutine push_h

!-----------------------------------------------------------------------------
recursive Subroutine rhs_h(y, ny, nasv, parms, nparms, deps, krk, nfev, error)
!-----------------------------------------------------------------------------
! calculate coefficient kRK from current state y and strain increment deps
! Masin hypoplastic model for clays with intergranular strains ???
! written 12/2005 (Tamagnini & Sellari)
!-----------------------------------------------------------------------------
  implicit None

  integer error, ny, nparms, nasv, i, nfev
  Double Precision y(ny), krk(ny), parms(nparms), deps(6)
  Double Precision sig(6), q(nasv)
  Double Precision f_sig(6), f_q(nasv)

  Double Precision zero, one, two, four
  Parameter (zero=0.0D0, one=1.0D0, two=2.0D0, four=4.0D0)

! ... update counter for the number of function f(y) evaluations
  nfev = nfev + 1

! ... initialize kRK
  Do i = 1, ny
    krk(i) = zero
  End Do

! ... recover current state variables (sig,q)
  Do i = 1, 6
    sig(i) = y(i)
  End Do

  Do i = 1, nasv
    q(i) = y(6+i)
  End Do

! ... build F_sig(6) and F_q(nasv) vectors and move them into kRK
  Call get_f_sig_q_h(sig, q, nasv, parms, nparms, deps, f_sig, f_q, error)

  If (error==10) Return

  Do i = 1, 6
    krk(i) = f_sig(i)
  End Do

  Do i = 1, nasv
    krk(6+i) = f_q(i)
  End Do

  Return
End Subroutine rhs_h

!-----------------------------------------------------------------------------
recursive Subroutine rkf23_update_h(y, n, nasv, dtsub, err_tol, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, dtime, error)
!-----------------------------------------------------------------------------
!
!  numerical solution of y'=f(y)
!  explicit, adapive RKF23 scheme with local time step extrapolation
!
!  Tamagnini, Sellari & Miriano 6/2005
!
!-----------------------------------------------------------------------------
  implicit None
!
  Logical elprsw
  integer n, nasv, nparms, i, ksubst, kreject, nfev
  integer maxnint, error, error_rkf
  Double Precision y(n), parms(nparms), dtsub, err_tol, dtmin
  Double Precision deps_np1(6), y_k(n), y_2(n), y_3(n), y_til(n)
  Double Precision y_hat(n)
  Double Precision t_k, dt_k, dtime
  Double Precision krk_1(n), krk_2(n), krk_3(n)
  Double Precision norm_r, s_hull
!
  Double Precision zero, half, one, two, three, four, six
  Double Precision ptnine, onesixth, onethird, twothirds, temp
  Parameter (zero=0.0D0, one=1.0D0, two=2.0D0)
  Parameter (three=3.0D0, four=4.0D0, six=6.0D0)
  Parameter (half=0.5D0, ptnine=0.9D0)
  onesixth = one/six
  onethird = one/three
  twothirds = two/three

! ... start of update process
  error_rkf = 0
  t_k = zero
  dt_k = dtsub/dtime
  ksubst = 0
  kreject = 0
  nfev = 0
  Do i = 1, n
    y_k(i) = y(i)
  End Do

! ... start substepping
  Do While (t_k<one)

    ksubst = ksubst + 1

! ... check for maximum number of substeps
    If (ksubst>maxnint) Then
!              write(6,*) 'number of substeps ', ksubst,
!     &                   ' is too big, step rejected'
      error = 3
      Return
    End If

! ... build RK functions
          !write(*, *) y_k(1), y_k(2), y_k(3), y_k(4), y_k(5), y_k(6)
    Call check_rkf_h(error_rkf, y_k, n, nasv, parms, nparms)

    If (error_rkf==1) Then
      error = 3
      Return
    Else
      Call rhs_h(y_k, n, nasv, parms, nparms, deps_np1, krk_1, nfev, error)
    End If

    If (error==10) Return

! ... find y_2
    temp = half*dt_k
    Do i = 1, n
      y_2(i) = y_k(i) + temp*krk_1(i)
    End Do

          !write(*, *) y_2(1), y_2(2), y_2(3), y_2(4), y_2(5), y_2(6)
    Call check_rkf_h(error_rkf, y_2, n, nasv, parms, nparms)

    If (error_rkf==1) Then
      error = 3
      Return
    Else
      Call rhs_h(y_2, n, nasv, parms, nparms, deps_np1, krk_2, nfev, error)
    End If

    If (error==10) Return

! ... find y_3
    Do i = 1, n
      y_3(i) = y_k(i) - dt_k*krk_1(i) + two*dt_k*krk_2(i)
    End Do

          !write(*, *) y_3(1), y_3(2), y_3(3), y_3(4), y_3(5), y_3(6)
    Call check_rkf_h(error_rkf, y_3, n, nasv, parms, nparms)

    If (error_rkf==1) Then
      error = 3
      Return
    Else
      Call rhs_h(y_3, n, nasv, parms, nparms, deps_np1, krk_3, nfev, error)
    End If

    If (error==10) Return

! ... approx. solutions of 2nd (y_til) and 3rd (y_hat) order
    Do i = 1, n
      y_til(i) = y_k(i) + dt_k*krk_2(i)
      y_hat(i) = y_k(i) + dt_k*(onesixth*krk_1(i)+twothirds*krk_2(i)+onesixth*krk_3(i))
    End Do

! ... local error estimate
    Call norm_res_h(y_til, y_hat, n, nasv, norm_r)
!     check if output y_hat can be used as an input into the next step
    Call check_rkf_h(error_rkf, y_hat, n, nasv, parms, nparms)

    If (error_rkf/=0) Then
      error = 3
      Return
    End If

! ... time step size estimator according to Hull
    If (norm_r/=0) Then
      s_hull = ptnine*dt_k*(err_tol/norm_r)**onethird
    Else
      s_hull = 1
    End If

    If (norm_r<err_tol) Then
! ... substep is accepted, update y_k and T_k and estimate new substep size DT_k

      Do i = 1, n
        y_k(i) = y_hat(i)
      End Do

      t_k = t_k + dt_k
      dt_k = min(four*dt_k, s_hull) ! adjust DT_k                   
      dtsub = dt_k*dtime
      dt_k = min(one-t_k, dt_k)

    Else
! ... substep is not accepted, recompute with new (smaller) substep size DT

      dt_k = max(dt_k/four, s_hull)

! ... check for minimum step size
      If (dt_k<dtmin) Then
!              write(6,*) 'substep size ', DT_k,
!     &                   ' is too small, step rejected'
        error = 3
        Return
      End If

    End If

! ... bottom of while loop
  End Do

! ... recover final state
  Do i = 1, n
    y(i) = y_k(i)
  End Do

  Return
End Subroutine rkf23_update_h

! Checks if RKF23 solout vector y is OK for hypoplasticity
recursive Subroutine check_rkf_h(error_rkf, y, ny, nasv, parms, nparms)
  implicit None
  integer error_rkf, ny, nasv, i, nparms, testnan, iopt
  Double Precision y(ny), parms(nparms)
  Double Precision sig(6), pmean, sig_star(6)
  Double Precision xn1(3), xn2(3), xn3(3), s(3), p, q, tmin
  Double Precision p_t, minstress

  integer intv(10)
  Double Precision realv(10)
  character *8 charv(10)

  minstress = 0.9
  p_t = parms(2)
  Do i = 1, 6
    sig(i) = y(i)
  End Do

  sig_star(1) = sig(1) - p_t
  sig_star(2) = sig(2) - p_t
  sig_star(3) = sig(3) - p_t
  sig_star(4) = sig(4)
  sig_star(5) = sig(5)
  sig_star(6) = sig(6)
  pmean = -(sig_star(1)+sig_star(2)+sig_star(3))/3

!     check for positive mean stress
  If (pmean<=minstress) Then
    error_rkf = 1
  End If

!     calculate minimum principal stress
  iopt = 0
  Call prnsig_h(iopt, sig_star, xn1, xn2, xn3, s(1), s(2), s(3), p, q)
  tmin = 1.0D+20
  Do i = 1, 3
    If (tmin>=-s(i)) Then
      tmin = -s(i)
    End If
  End Do

!     check for tension
  If (tmin<=minstress) Then
    error_rkf = 1
  End If

  Return
End Subroutine check_rkf_h

!-----------------------------------------------------------------------------
recursive Subroutine solout_h(stress, ntens, asv, nasv, ddsdde, y, nydim, pore, depsv_np1, parms, nparms, dd)
!-----------------------------------------------------------------------------
! copy the vector of state variables to umat output
! modified 7/2005 (Tamagnini, Sellari)
!
! NOTE: solid mechanics convention for stress and strain components
!       pore is always positive in compression
!-----------------------------------------------------------------------------
  implicit None
  integer nydim, nasv, nparms, ntens, i, j
  Double Precision y(nydim), asv(nasv), stress(ntens)
  Double Precision ddsdde(ntens, ntens), dd(6, 6)
  Double Precision parms(nparms), bulk_w, pore, depsv_np1

  bulk_w = parms(15)

! ... update excess pore pressure (if undrained conditions), compression positive
  pore = pore - bulk_w*depsv_np1

! updated total stresses (effective stresses stored in y(1:6))
  Do i = 1, ntens
    If (i<=3) Then
      stress(i) = y(i) - pore
    Else
      stress(i) = y(i)
    End If
  End Do

! additional state variables (first 6 components are intergranular strains)
  Do i = 1, nasv
    asv(i) = y(6+i)
  End Do

! consistent tangent stiffness
  Do j = 1, ntens
    Do i = 1, ntens
      ddsdde(i, j) = dd(i, j)
    End Do
  End Do

  Do j = 1, 3
    Do i = 1, 3
      ddsdde(i, j) = ddsdde(i, j) + bulk_w
    End Do
  End Do

  Return
End Subroutine solout_h


!-----------------------------------------------------------------------------
recursive Subroutine wrista_h(mode, y, nydim, deps_np1, dtime, coords, statev, nstatv, parms, nparms, noel, npt, ndi, nshr, kstep, kinc)
!-----------------------------------------------------------------------------
! ... subroutine for managing output messages
!
!     mode
!
!     all = writes:             kstep, kinc, noel, npt
!       2   = writes also:      error message,coords(3),parms(nparms),ndi,nshr,stress(nstress)
!                                               deps(nstress),dtime,statev(nstatv)
!     3   = writes also:        stress(nstress),deps(nstress),dtime,statev(nstatv)
!-----------------------------------------------------------------------------
  implicit None
  integer mode, nydim, nstatv, nparms, noel, npt, ndi, nshr, kstep, kinc, i
  Double Precision y(nydim), statev(nstatv), parms(nparms)
  Double Precision deps_np1(6), coords(3), dtime

! ... writes for mode = 2
  If (mode==2) Then
    Write (6, *) '==================================================='
    Write (6, *) 'ERROR: abaqus job failed during call of UMAT'
    Write (6, *) '==================================================='
    Write (6, *) 'state dump:'
    Write (6, *)
  End If

! ... writes for all mode values
  Write (1, 111) 'Step: ', kstep, 'increment: ', kinc, 'element: ', noel, 'Integration point: ', npt
  Write (6, *)

! ... writes for mode = 2
  If (mode==2) Then
    Write (6, *) 'Co-ordinates of material point:'
    Write (1, 104) 'x1 = ', coords(1), ' x2 = ', coords(2), ' x3 = ', coords(3)
    Write (6, *)
    Write (6, *) 'Material parameters:'
    Write (6, *)
    Do i = 1, nparms
      Write (1, 105) 'prop(', i, ') = ', parms(i)
    End Do
    Write (6, *)
    Write (1, 102) 'No. of mean components:  ', ndi
    Write (1, 102) 'No. of shear components: ', nshr
    Write (6, *)
  End If

! ... writes for mode = 2 or 3
  If ((mode==2) .Or. (mode==3)) Then
    Write (6, *) 'Stresses:'
    Write (6, *)
    Write (1, 101) 'sigma(1) = ', y(1)
    Write (1, 101) 'sigma(2) = ', y(2)
    Write (1, 101) 'sigma(3) = ', y(3)
    Write (1, 101) 'sigma(4) = ', y(4)
    Write (1, 101) 'sigma(5) = ', y(5)
    Write (1, 101) 'sigma(6) = ', y(6)
    Write (6, *)
    Write (6, *) 'Strain increment:'
    Write (6, *)
    Write (1, 101) 'deps_np1(1) = ', deps_np1(1)
    Write (1, 101) 'deps_np1(2) = ', deps_np1(2)
    Write (1, 101) 'deps_np1(3) = ', deps_np1(3)
    Write (1, 101) 'deps_np1(4) = ', deps_np1(4)
    Write (1, 101) 'deps_np1(5) = ', deps_np1(5)
    Write (1, 101) 'deps_np1(6) = ', deps_np1(6)
    Write (6, *)
    Write (6, *) 'Time increment:'
    Write (6, *)
    Write (1, 108) 'dtime = ', dtime
    Write (6, *)
    Write (6, *) 'Internal variables:'
    Write (6, *)
    Write (1, 109) 'del(1) = ', statev(1)
    Write (1, 109) 'del(2) = ', statev(2)
    Write (1, 109) 'del(3) = ', statev(3)
    Write (1, 109) 'del(4) = ', statev(4)
    Write (1, 109) 'del(5) = ', statev(5)
    Write (1, 109) 'del(6) = ', statev(6)
    Write (1, 109) 'void   = ', statev(7)
    Write (6, *)
    Write (6, *) '==================================================='

  End If
  Return

  101 Format (1X, A15, E11.4)
  102 Format (1X, A25, I1)
  103 Format (1X, A7, I5)
  104 Format (1X, 3(A5,F10.4,2X))
  105 Format (1X, A5, I2, A4, F20.3)
  106 Format (1X, 3(A9,F12.4,2X))
  107 Format (1X, 3(A10,F12.4,2X))
  108 Format (1X, A8, F12.4)
  109 Format (1X, A6, F10.4)
  110 Format (1X, A5, F10.4)
  111 Format (1X, A6, I4, 2X, A11, I4, 2X, A9, I10, 2X, A19, I4)
End Subroutine wrista_h

!-----------------------------------------------------------------------------
recursive Subroutine calc_statev_h(stress, statev, parms, nparms, nasv, nasvdim, deps)
!-----------------------------------------------------------------------------
!
!  computes additional state variables for postprocessing
!
!-----------------------------------------------------------------------------
  implicit None

  Logical elprsw
  integer ntens, jj, kk, i
  integer n, nasv, nparms, nfev, nasvdim
  integer maxnint, error

  Double Precision parms(nparms), dot_vect_h
  Double Precision stress(6), statev(nasvdim)
  Double Precision deps(6), tmax, tmin
  Double Precision mm(6, 6), hhtmp(nasv, 6)
  Double Precision ll(6, 6), nn(6)
  integer istrain
  Double Precision zero, two, four, iopt, three
  Double Precision i1, i2, i3, cos3t, pp, qq
  Double Precision sin2phi, sinphi, sig_star(6), p_t
  Double Precision norm_del, norm_del2, del(6)

  Parameter (zero=0.0D0, two=2.0D0, four=4.0D0, three=3.0D0)

! ... calc phimob (statev 11) from Matsuoka-Nakai YS
  p_t = parms(2)
  Do i = 1, 3
    sig_star(i) = stress(i) - p_t
  End Do
  Do i = 4, 6
    sig_star(i) = stress(i)
  End Do
  Call inv_sig_h(sig_star, pp, qq, cos3t, i1, i2, i3)

  If (i3/=0) Then
    sin2phi = (9.0D0+i1*i2/i3)/(1.0D0+i1*i2/i3)
  Else
    sin2phi = 0.0
  End If

  If (sin2phi<0) Then
    sin2phi = 0.0
  End If
  If (sin2phi>1) Then
    sin2phi = 1.0
  End If
  sinphi = sqrt(sin2phi)

  statev(11) = asin(sinphi)*180.0D0/3.141592D0

! ... cal normalized length of intergranular strain rho (statev 12)
  If (parms(10)<=0.5) Then
    istrain = 0
  Else
    istrain = 1
  End If

  If (istrain==1) Then

    Do i = 1, 6
      del(i) = statev(i)
    End Do

    norm_del2 = dot_vect_h(2, del, del, 6)
    norm_del = dsqrt(norm_del2)
    statev(12) = norm_del/parms(12)

  Else
    statev(12) = 0.0
  End If

  Return
End Subroutine calc_statev_h

! *** checks whether number is NaN ***
recursive Subroutine umatisnan_h(chcknum, testnan)
  Double Precision chcknum
  integer testnan
  If (.Not. (chcknum>=0. .Or. chcknum<0.)) testnan = 1
  If (chcknum>1.D30) testnan = 1
  If (chcknum<-1.D30) testnan = 1
  If (chcknum/=chcknum) testnan = 1
      !debug
  testnan = 0
  Return
End Subroutine umatisnan_h

recursive Subroutine xit_h
  Stop
  Return
End Subroutine xit_h

!***********************************************************************
recursive Subroutine prnsig_h(iopt, s, xn1, xn2, xn3, s1, s2, s3, p, q)
  implicit Double Precision (A-H, O-Z)
  Dimension s(*), xn1(*), xn2(*), xn3(*)

  If (iopt==1) Then
    Call eig_3_h(0, s, xn1, xn2, xn3, s1, s2, s3, p, q) ! with Eigenvectors  
  Else
    Call eig_3a_h(0, s, s1, s2, s3, p, q) ! no Eigenvectors               
  End If
  Return
End Subroutine prnsig_h

!***********************************************************************
recursive Subroutine eig_3_h(iopt, st, xn1, xn2, xn3, s1, s2, s3, p, q)
  implicit Double Precision (A-H, O-Z)
  Dimension st(6), a(3, 3), v(3, 3), xn1(3), xn2(3), xn3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
  a(1, 1) = st(1) ! xx                                               
  a(1, 2) = st(4) ! xy = yx                                          
  a(1, 3) = st(6) 
! zx = xz                                          
  a(2, 1) = st(4) ! xy = yx                                          
  a(2, 2) = st(2) ! yy                                               
  a(2, 3) = st(5) 
! zy = yz                                          
  a(3, 1) = st(6) ! zx = xz                                          
  a(3, 2) = st(5) ! zy = yz                                          
  a(3, 3) = st(3) 
      ! Set V to unity matrix
! zz                                               
  v(1, 1) = 1
  v(2, 1) = 0
  v(3, 1) = 0

  v(1, 2) = 0
  v(2, 2) = 1
  v(3, 2) = 0

  v(1, 3) = 0
  v(2, 3) = 0
  v(3, 3) = 1


  abs_max_s = 0.0
  Do i = 1, 3
    Do j = 1, 3
      If (abs(a(i,j))>abs_max_s) abs_max_s = abs(a(i,j))
    End Do
  End Do
  tol = 1D-20*abs_max_s
  it = 0
  itmax = 50
  Do While (it<itmax .And. abs(a(1,2))+abs(a(2,3))+abs(a(1,3))>tol)
    it = it + 1
    Do k = 1, 3
      If (k==1) Then
        ip = 1
        iq = 2
      Else If (k==2) Then
        ip = 2
        iq = 3
      Else
        ip = 1
        iq = 3
      End If
      If (abs(a(ip,iq))>tol) Then
        tau = (a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
        If (tau>=0.0) Then
          sign_tau = 1.0
        Else
          sign_tau = -1.0
        End If
        t = sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
        c = 1.0/sqrt(1.0+t*t)
        s = t*c
        a1p = c*a(1, ip) - s*a(1, iq)
        a2p = c*a(2, ip) - s*a(2, iq)
        a3p = c*a(3, ip) - s*a(3, iq)
        a(1, iq) = s*a(1, ip) + c*a(1, iq)
        a(2, iq) = s*a(2, ip) + c*a(2, iq)
        a(3, iq) = s*a(3, ip) + c*a(3, iq)
        a(1, ip) = a1p
        a(2, ip) = a2p
        a(3, ip) = a3p

        v1p = c*v(1, ip) - s*v(1, iq)
        v2p = c*v(2, ip) - s*v(2, iq)
        v3p = c*v(3, ip) - s*v(3, iq)
        v(1, iq) = s*v(1, ip) + c*v(1, iq)
        v(2, iq) = s*v(2, ip) + c*v(2, iq)
        v(3, iq) = s*v(3, ip) + c*v(3, iq)
        v(1, ip) = v1p
        v(2, ip) = v2p
        v(3, ip) = v3p

        ap1 = c*a(ip, 1) - s*a(iq, 1)
        ap2 = c*a(ip, 2) - s*a(iq, 2)
        ap3 = c*a(ip, 3) - s*a(iq, 3)
        a(iq, 1) = s*a(ip, 1) + c*a(iq, 1)
        a(iq, 2) = s*a(ip, 2) + c*a(iq, 2)
        a(iq, 3) = s*a(ip, 3) + c*a(iq, 3)
        a(ip, 1) = ap1
        a(ip, 2) = ap2
        a(ip, 3) = ap3
! a(ip,iq)<>0                                          
      End If
! k                                                      
    End Do
      ! principal values on diagonal of a
  End Do
! While                                                    
  s1 = a(1, 1)
  s2 = a(2, 2)
  s3 = a(3, 3)
      ! Derived invariants
  p = (s1+s2+s3)/3
  q = sqrt(((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/2)

      ! Sort eigenvalues S1 <= S2 <= S3
  is1 = 1
  is2 = 2
  is3 = 3
  If (s1>s2) Then
    t = s2
    s2 = s1
    s1 = t
    it = is2
    is2 = is1
    is1 = it
  End If
  If (s2>s3) Then
    t = s3
    s3 = s2
    s2 = t
    it = is3
    is3 = is2
    is2 = it
  End If
  If (s1>s2) Then
    t = s2
    s2 = s1
    s1 = t
    it = is2
    is2 = is1
    is1 = it
  End If
  Do i = 1, 3
    xn1(i) = v(i, is1) ! first  column                               
    xn2(i) = v(i, is2) ! second column                               
    xn3(i) = v(i, is3) ! third  column                               
  End Do
  Return

End Subroutine eig_3_h
! Eig_3                                                       
recursive Subroutine eig_3a_h(iopt, st, s1, s2, s3, p, q) ! xN1,xN2,xN3,          
  implicit Double Precision (A-H, O-Z)
  Dimension st(6), a(3, 3)       !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
!  V(3,3),xN1(3),xN2(3),xN3(3)           
  a(1, 1) = st(1) ! xx                                               
  a(1, 2) = st(4) ! xy = yx                                          
  a(1, 3) = st(6) 
! zx = xz                                          
  a(2, 1) = st(4) ! xy = yx                                          
  a(2, 2) = st(2) ! yy                                               
  a(2, 3) = st(5) 
! zy = yz                                          
  a(3, 1) = st(6) ! zx = xz                                          
  a(3, 2) = st(5) ! zy = yz                                          
  a(3, 3) = st(3) 
! zz                                               
  abs_max_s = 0.0
  Do i = 1, 3
    Do j = 1, 3
      If (abs(a(i,j))>abs_max_s) abs_max_s = abs(a(i,j))
    End Do
  End Do
  tol = 1D-20*abs_max_s
  If (iopt==1) tol = 1D-50*abs_max_s
  it = 0
  itmax = 50

  Do While (it<itmax .And. abs(a(1,2))+abs(a(2,3))+abs(a(1,3))>tol)

    it = it + 1
    Do k = 1, 3
      If (k==1) Then
        ip = 1
        iq = 2
      Else If (k==2) Then
        ip = 2
        iq = 3
      Else
        ip = 1
        iq = 3
      End If

      If (abs(a(ip,iq))>tol) & ! ongelijk nul ?     
        Then
        tau = (a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
        If (tau>=0.0) Then
          sign_tau = 1.0
        Else
          sign_tau = -1.0
        End If
        t = sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
        c = 1.0/sqrt(1.0+t*t)
        s = t*c
        a1p = c*a(1, ip) - s*a(1, iq)
        a2p = c*a(2, ip) - s*a(2, iq)
        a3p = c*a(3, ip) - s*a(3, iq)
        a(1, iq) = s*a(1, ip) + c*a(1, iq)
        a(2, iq) = s*a(2, ip) + c*a(2, iq)
        a(3, iq) = s*a(3, ip) + c*a(3, iq)
        a(1, ip) = a1p
        a(2, ip) = a2p
        a(3, ip) = a3p

        ap1 = c*a(ip, 1) - s*a(iq, 1)
        ap2 = c*a(ip, 2) - s*a(iq, 2)
        ap3 = c*a(ip, 3) - s*a(iq, 3)
        a(iq, 1) = s*a(ip, 1) + c*a(iq, 1)
        a(iq, 2) = s*a(ip, 2) + c*a(iq, 2)
        a(iq, 3) = s*a(ip, 3) + c*a(iq, 3)
        a(ip, 1) = ap1
        a(ip, 2) = ap2
        a(ip, 3) = ap3
! a(ip,iq)<>0                                          
      End If
! k                                                      
    End Do
      ! principal values on diagonal of a
  End Do
! While                                                    
  s1 = a(1, 1)
  s2 = a(2, 2)
  s3 = a(3, 3)
      ! Derived invariants
  p = (s1+s2+s3)/3
  q = sqrt(((s1-s2)**2+(s2-s3)**2+(s3-s1)**2)/2)

  If (s1>s2) Then
    t = s2
    s2 = s1
    s1 = t
  End If
  If (s2>s3) Then
    t = s3
    s3 = s2
    s2 = t
  End If
  If (s1>s2) Then
    t = s2
    s2 = s1
    s1 = t
  End If
  Return

End Subroutine eig_3a_h
!-----------------------------------------------------------------------------
! Eig_3a                                                      
recursive Subroutine calc_elasti_h(y, n, nasv, dtsub, err_tol, maxnint, dtmin, deps_np1, parms, nparms, nfev, elprsw, dtime, ddtan, youngel, nuel, error)
!-----------------------------------------------------------------------------
!  numerical solution of y'=f(y)
!  explicit, adapive RKF23 scheme with local time step extrapolation
!
!  Tamagnini, Sellari & Miriano 6/2005
!-----------------------------------------------------------------------------
  implicit None

  Logical elprsw
  integer n, nasv, nparms, i, ksubst, kreject, nfev
  integer maxnint, error, error_rkf, tension, j
!
  Double Precision y(n), parms(nparms), dtsub, err_tol, dtmin
  Double Precision deps_np1(6), y_k(n), y_2(n), y_3(n), y_til(n)
  Double Precision y_hat(n), ddtan(6, 6)
  Double Precision t_k, dt_k, dtime, ii(6, 6), krondelta(6)
  Double Precision krk_1(n), krk_2(n), krk_3(n)
  Double Precision norm_r, s_hull, youngel, nuel, f_sig(6)

  Double Precision zero, half, one, two, three, four, six
  Double Precision ptnine, onesixth, onethird, twothirds, temp
  Parameter (zero=0.0D0, one=1.0D0, two=2.0D0, three=3.0D0)
  Parameter (four=4.0D0, six=6.0D0, half=0.5D0, ptnine=0.9D0)

  onesixth = one/six
  onethird = one/three
  twothirds = two/three

! ... initialize y_k vector and other variables
  Do i = 1, n
    y_k(i) = zero
  End Do

! ... fourth order identity tensors in Voigt notation
  Do i = 1, 6
    Do j = 1, 6
      ii(i, j) = zero
    End Do
  End Do

  ii(1, 1) = one
  ii(2, 2) = one
  ii(3, 3) = one
  ii(4, 4) = half
  ii(5, 5) = half
  ii(6, 6) = half

  krondelta(1) = one
  krondelta(2) = one
  krondelta(3) = one
  krondelta(4) = zero
  krondelta(5) = zero
  krondelta(6) = zero

! ... Elastic stiffness tensor
  Do i = 1, 6
    Do j = 1, 6
      ddtan(i, j) = (youngel/(1+nuel))*(ii(i,j)+nuel/(1-2*nuel)*krondelta(i)*krondelta(j))
    End Do
  End Do

  Call matmul_h(ddtan, deps_np1, f_sig, 6, 6, 1)

  Do i = 1, 6
    y(i) = y(i) + f_sig(i)
  End Do

  Return
End Subroutine calc_elasti_h
