module isostasy_benchmarks

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use isostasy_defs, only : sp, dp, wp, pi

    implicit none

    type isos_analytical_elva_disk_load_class 
        integer  :: ifunc   = 0                      ! choice of function f(x) to integrate (integrand of analytical solution, in case there are several) 
        integer  :: n_evals = 0                      ! number of function evaluations
        integer  :: method  = 6                      ! integration method (#points in gaussian quadrature)
        procedure(func_1d),pointer :: fun => null()  ! function f(x) to be integrated (at given time t)
        procedure(gauss_func),pointer :: gg => null() ! the gauss quadrature formula to use
        real(wp) :: a      = 0._wp                   ! lower limit of integration
        real(wp) :: b      = 0._wp                   ! upper limit of integration (may be less than a)
        real(wp) :: tol    = 0._wp                   ! the requested relative error tolerance.
        real(wp) :: val    = 0._wp                   ! the value of x, only used for multiple integration to pass to the inner integrals
        real(wp) :: r      = 0._wp                   ! distance of disk center at which function f(x) is evaluated
        real(wp) :: time   = 0._wp                   ! time at which function is f(x) evaluated
        real(wp) :: r0     = 0._wp                   ! [m] 
        real(wp) :: h0     = 0._wp                   ! [m] 
        real(wp) :: D_lith = 0.e+24                  ! [N m]   
        real(wp) :: eta    = 0.e+21                  ! [Pa s]

        real(wp) :: rho_ice 
        real(wp) :: rho_uppermantle 
        real(wp) :: g 

        real(wp), allocatable :: kappa_mod(:)
        real(wp), allocatable :: rad(:) 
        integer,  allocatable :: lrad(:,:)
        
    contains

        private
        
        procedure :: dgauss_generic                                    ! core integration routine. refactored from
                                                                       ! SLATEC with selectable quadrature method
        procedure, public :: initialize => initialize_integration_class ! to set up the class
        procedure, public :: integrate  => integrate_1d                 ! to integrate function `fun`
        
    end type isos_analytical_elva_disk_load_class
    

    abstract interface
        
        function func_1d(me,x) result(f)
            ! 1d user function f(x)
            import :: wp,isos_analytical_elva_disk_load_class
            implicit none
            class(isos_analytical_elva_disk_load_class),intent(inout) :: me
            real(wp), intent(in)                      :: x
            real(wp)                                  :: f
        end function func_1d

        function gauss_func(me, x, h) result(f)
            ! gauss quadrature formula
            import :: wp,isos_analytical_elva_disk_load_class
            implicit none
            class(isos_analytical_elva_disk_load_class),intent(inout)  :: me
            real(wp), intent(in)                    :: x
            real(wp), intent(in)                    :: h
            real(wp)                                :: f
        end function gauss_func
       
    end interface

    private
    public :: isosbench_elva_disk

contains

    subroutine isosbench_elva_disk(z_bed,r0,h0,eta,dx,D_lith,rho_ice,rho_uppermantle,g,time)
        ! This will calculate the analytical solution 'elva_disk' for a specific time

        implicit none

        real(wp), intent(INOUT) :: z_bed(:,:)   ! [m]
        real(wp), intent(IN)    :: r0           ! Radius, r0=1000e3 m (1000 km) by default
        real(wp), intent(IN)    :: h0           ! Ice thickness, h0=1000 m by default
        real(wp), intent(IN)    :: eta          ! Viscosity, eta=1e21 Pa s by default
        real(wp), intent(IN)    :: dx 
        real(wp), intent(IN)    :: D_lith
        real(wp), intent(IN)    :: rho_ice
        real(wp), intent(IN)    :: rho_uppermantle
        real(wp), intent(IN)    :: g
        real(wp), intent(IN)    :: time 

        ! Local variables
        integer  :: nx, ny 
        
        type(isos_analytical_elva_disk_load_class) :: ana        ! object containing all info relative to function to be integrated for analytical solution

        ! Parameters for analytical ELVA disk solution - Alex - class apart?                              
        real(wp), parameter :: kappa_min = 0.  ! [] Minimum kappa for gaussian quadrature       
        real(wp), parameter :: kappa_max = 0.1 ! [] Maximum kappa for gaussian quadrature       
        real(wp), parameter :: dk = 1.e-3      ! [] Step in kappa for gaussian quadrature    

        ! Initialization =======

        nx = size(z_bed,1)
        ny = size(z_bed,2) 

        ! Store several parameters for the analytical integrand
        ana%r0         = r0
        ana%h0         = h0
        ana%D_lith     = D_lith
        ana%eta        = eta
        
        ana%rho_ice    = rho_ice
        ana%rho_uppermantle      = rho_uppermantle
        ana%g          = g
        
        ! Calculate some arrays that are also needed 
        call calc_analytical_viscous_disk_params(ana%kappa_mod,ana%rad,ana%lrad, &
                                                            kappa_min,kappa_max,dk,nx,ny,dx)                                    

        ! Updating ===========

        call calc_analytical_viscous_disk(ana,z_bed,dx,time)  

        write(*,*) "isostasy benchmark elva_disk: ", time, minval(z_bed), maxval(z_bed)

        return

    end subroutine isosbench_elva_disk


!=========================================================
!
! Routines for analytical ELVA disk
!
!=========================================================

    subroutine calc_analytical_viscous_disk_params(kappa_mod,r,lr,kappa_min,kappa_max,dk,nx,ny,dx)  
       
        real(wp), allocatable, intent(OUT)  :: kappa_mod(:)
        real(wp), allocatable, intent(OUT)  :: r(:) 
        integer,  allocatable, intent(OUT)  :: lr(:,:)

        real(wp), intent(IN)                :: kappa_min
        real(wp), intent(IN)                :: kappa_max
        real(wp), intent(IN)                :: dk
        integer,  intent(IN)                :: nx
        integer,  intent(IN)                :: ny
        real(wp), intent(IN)                :: dx
        
        ! Local variables
        integer, allocatable                :: n(:,:)
        real(wp)                            :: xd, yd
        integer                             :: i, j, ip, iq, ic, jc, k, nk, l, nl
        real(wp), allocatable               :: dist2c(:,:)
        
        nk = int((kappa_max-kappa_min)/dk)
        nl = int(max(nx,ny)*sqrt(2.)/2) + 2

        allocate(kappa_mod(nk+1))
        allocate(dist2c(nx,ny))
        allocate(r(nl))
        allocate(lr(nx,ny))
        allocate(n(nx,ny))
      
        do k = 1, nk+1
            kappa_mod(k) = kappa_min + dk * (k-1)
        end do

      
        ic = (nx-1)/2 + 1
        jc = (ny-1)/2 + 1

        do i = 1, nx
        do j = 1, ny
            xd = dx*(i-ic)
            yd = dx*(j-jc)
            dist2c(i,j) = sqrt(xd**2 + yd**2)
        end do
        end do

        ! Remap to (nx * xy) grid (w = w(r,t) so all points with same r have same w; 
        ! this reduces computational time      

        do l = 1, nl 
            r(l) =  dx * (l-1)
        end do
     
        do i = 1, nx 
        do j = 1, ny
            n(i,j) = 0       
            do l = 1, nl-1
                if (dist2c(i,j).ge.r(l) .and. dist2c(i,j) .lt. r(l+1) ) then
                    lr(i,j) = l 
                    n(i,j) = n(i,j) + 1
                end if
            end do
            if (n(i,j).ne.1) then
                write(*,*) "==> error in radial distance allocation", n(i,j), i, j
                stop
            end if
        end do
        end do

        deallocate(n)
         
        return
      
    end subroutine calc_analytical_viscous_disk_params
    
    subroutine calc_analytical_viscous_disk(me,w,dx,t) 

        ! Calculate analytical solution for displacement for the asthenosphere viscous disk 
        ! u(r,t) as in Bueler et al (2007), eq. 17
        ! remap into original (nx * ny) grid

        implicit none

        class(isos_analytical_elva_disk_load_class), intent(INOUT) :: me
        real(wp), intent(OUT)         :: w(:,:)
        real(wp), intent(IN)          :: dx
        real(wp), intent(IN)          :: t
        
        ! Local variables
        real(wp), allocatable         :: wr(:)
        integer(kind=4), allocatable  :: n(:,:)
        
        real(wp)                      :: ans
        real(wp)                      :: err
        integer                       :: method    ! quadrature method to use for x
        integer                       :: ierr      ! error code
        
        integer(kind=4)               :: i, j, k, l, nx, ny, nk, nl

        real(wp), parameter :: tol = 1.e-6 !  error tolerance
        
        nx = size(w,1)
        ny = size(w,2)
        nk = size(me%kappa_mod)-1
        nl = size(me%rad)

        allocate(wr(nl)) 
        allocate(n(nx,ny))

        method       = me%method  ! 6-point gaussian quadrature     
        me%time      = t          ! [hours] 
        me%n_evals   = 0

        w =  0.0
        wr = 0.0
      
        do l = 1, nl ! number of points neccessary to cover whole domain radially (ca. diagonal)
         
            me%r = me%rad(l)

            do k = 1, nk

                ! sets parameters but most notably integration limits      

                call initialize_integration_class(me, fx=analytical_integrand, &
                            xl=me%kappa_mod(k), xu=me%kappa_mod(k+1), tolx=tol, methodx=method)

                ! reset number of function evaluations:

                me%n_evals = 0

                ! integrate the function over specified interval

                call integrate_1d (me, ans, ierr, err) 

                wr(l)  = wr(l) + ans

            end do ! nk 
        end do ! nl

        do i = 1, nx
        do j = 1, ny
            w(i,j) = wr(me%lrad(i,j))  
        end do
        end do
      
        w  = w - 0.25*(w(1,1)+w(nx,ny)+w(1,ny)+w(nx,1)) 

        deallocate(n)

        return
         
    end subroutine calc_analytical_viscous_disk

    subroutine initialize_integration_class(me,fx,xl,xu,tolx,methodx)

        implicit none

        class(isos_analytical_elva_disk_load_class),intent(inout)  :: me

        procedure(func_1d)    :: fx       !! 1d function: f(x)
        real(wp),intent(in)   :: xl       !! x integration lower bound
        real(wp),intent(in)   :: xu       !! x integration upper bound
        real(wp),intent(in)   :: tolx     !! error tolerance for dx integration
        integer,intent(in)    :: methodx  !! quadrature method to use for x

        ! select quadrature rule (only g6 implemented; for others cut and paste)
        select case (methodx)
        case(6);  me%gg => g6
        !    case(8);  me%gg => g8  
        !    case(10); me%gg => g10
        !    case(12); me%gg => g12
        !    case(14); me%gg => g14
        case default
        error stop 'invalid quadrature method in initialize_integration_class'
        end select

        me%fun  => fx         !function f(x) to integrate
        me%tol  = tolx        !tolerance
        me%a    = xl          !lower bound
        me%b    = xu          !upper bound
        me%method = methodx   !method

    return
    
    end subroutine initialize_integration_class

    subroutine initialize_analytical_integrand(me,r0,h0,D_lith,eta,rho_ice,rho_uppermantle,g) 

        implicit none

        class(isos_analytical_elva_disk_load_class),intent(out)  :: me
        real(wp), intent(in)   ::  h0, r0, D_lith, eta, rho_ice, rho_uppermantle, g 

        me%r0         = r0
        me%h0         = h0
        me%D_lith     = D_lith
        me%eta        = eta
        
        me%rho_ice    = rho_ice
        me%rho_uppermantle      = rho_uppermantle
        me%g          = g
        
        return

    end subroutine initialize_analytical_integrand

  function analytical_integrand(me,x) result(f)

        implicit none

        class(isos_analytical_elva_disk_load_class),intent(inout)  :: me
        real(wp), intent(in)  :: x ! (kappa)
        real(wp)              :: f
        real(wp) :: beta, h0, r0, eta, D_lith, t_sec, r_eq_zero
        real(wp) :: rho_ice, rho_uppermantle, g 

        r0     = me%r0
        h0     = me%h0
        D_lith = me%D_lith
        eta    = me%eta
        t_sec  = me%time*365*24*3600

        rho_ice = me%rho_ice 
        rho_uppermantle   = me%rho_uppermantle 
        g       = me%g 

        beta = rho_uppermantle*g + D_lith*(x**4)

        r_eq_zero = 0.

        f   =  rho_ice*g*h0*r0*( exp(-beta*t_sec/(2*eta*x)) - 1.) * bessel_j1(x*r0) *  bessel_j0(x*me%r) / beta

        ! asymptotic solution (t --> infty ) 
        !    f   =  rho_ice*g*h0*r0*( - 1.) * bessel_j1(x*r0) *  bessel_j0(x*me%r) / beta

        !number of function evaluations
        me%n_evals = me%n_evals + 1          

    end function analytical_integrand


    subroutine integrate_1d (me, ans, ierr, err)
        ! 1d integral by Jacob Williams.

        implicit none

        class(isos_analytical_elva_disk_load_class),intent(inout)  :: me
        real(wp),intent(out)  :: ans
        integer,intent(out)   :: ierr
        real(wp),intent(out)  :: err


        !call the low-level routine 
        call me%dgauss_generic(me%a, me%b, me%tol, ans, ierr, err)

        return
    
    end subroutine integrate_1d

!===========================================================
!
!  Routines for Gaussian quadrature
!
! ==========================================================
  
!  Integrate a real function of one variable over a finite
!  interval using the specified adaptive algorithm.
!  Intended primarily for high accuracy
!  integration or integration of smooth functions.
!
!### License
!  * SLATEC is public domain software: http://www.netlib.org/slatec/guide
!
!### See also
!  * Original sourcecode from: http://www.netlib.org/slatec/src/dgaus8.f
!
!### Author
!  * Jones, R. E., (SNLA) -- Original SLATEC code.
!  * Jacob Williams : 1/20/2020 : refactored to modern Fortran and generalized.
!
!@note This function is recursive.
!      [It can call itself indirectly during double integration]


  recursive subroutine dgauss_generic (me, lb, ub, error_tol, ans, ierr, err)

    implicit none

    class(isos_analytical_elva_disk_load_class),intent(inout)  :: me

    real(wp),intent(in)   :: lb         !! lower bound of the integration
    real(wp),intent(in)   :: ub         !! upper bound of the integration
    real(wp),intent(in)   :: error_tol  !! is a requested pseudorelative error tolerance.  normally
                                        !! pick a value of abs(error_tol) so that
                                        !! dtol < abs(error_tol) <= 1.0e-3 where dtol is the larger
                                        !! of 1.0e-18 and the real(wp) unit roundoff d1mach(4).
                                        !! ans will normally have no more error than abs(error_tol)
                                        !! times the integral of the absolute value of fun(x).  usually,
                                        !! smaller values of error_tol yield more accuracy and require
                                        !! more function evaluations.
    real(wp),intent(out)  :: ans        !! computed value of integral
    integer,intent(out)   :: ierr       !! status code:
                                        !!
                                        !!  * normal codes:
                                        !!    * 1 : `ans` most likely meets requested error tolerance,
                                        !!      or `lb=ub`.
                                        !!    * -1 : `lb` and `ub` are too nearly equal to allow normal
                                        !!      integration. `ans` is set to zero.
                                        !!  * abnormal code:
                                        !!    * 2 : `ans` probably does not meet requested error tolerance.
    real(wp),intent(out)  :: err        !! an estimate of the absolute error in `ans`.
                                        !! the estimated error is solely for information to the user and
                                        !! should not be used as a correction to the computed integral.

    real(wp),parameter  :: sq2      = sqrt(2._wp)
    real(wp),parameter  :: ln2      = log(2._wp)
    integer,parameter   :: nlmn     = 1                    !! ??
    integer,parameter   :: kmx      = 5000                 !! ??
    integer,parameter   :: kml      = 6                    !! ??
    real(wp),parameter  :: magic    = 0.30102000_wp        !! ??
    integer,parameter   :: iwork    = 60                   !! size of the work arrays. ?? Why 60 ??
    real(wp),parameter  :: bb       = radix(1_wp)          !! machine constant
    real(wp),parameter  :: d1mach4  = bb**(1-digits(1_wp)) !! machine constant
    real(wp),parameter  :: d1mach5  = log10(bb)            !! machine constant

    integer                   :: k,l,lmn,lmx,mxl,nbits,nib,nlmx
    real(wp)                  :: ae,anib,area,c,ee,ef,eps,est,gl,glr,tol
    real(wp),dimension(iwork) :: aa,hh,vl,gr
    integer,dimension(iwork)  :: lr

    ans = 0_wp
    ierr = 1
    err = 0_wp
    if (lb == ub) return
    aa = 0_wp
    hh = 0_wp
    vl = 0_wp
    gr = 0_wp
    lr = 0
    k = digits(1._wp)
    anib = d1mach5*k/magic
    nbits = anib
    nlmx = min(60,(nbits*5)/8)         ! ... is this the same 60 as iwork???
    lmx = nlmx
    lmn = nlmn
    if (ub /= 0_wp) then
        if (sign(1._wp,ub)*lb > 0_wp) then
            c = abs(1._wp-lb/ub)
            if (c <= 0.1_wp) then
                if (c <= 0_wp) return
                anib = 0.5_wp - log(c)/ln2
                nib = anib
                lmx = min(nlmx,nbits-nib-7)
                if (lmx < 1) then
                    ! lb and ub are too nearly equal to allow
                    ! normal integration [ans is set to zero]
                    ierr = -1
                    return
                end if
                lmn = min(lmn,lmx)
            end if
        end if
    end if
    tol = max(abs(error_tol),2._wp**(5-nbits))/2._wp
    if (error_tol == 0_wp) tol = sqrt(d1mach4)
    eps = tol
    hh(1) = (ub-lb)/4. 
    aa(1) = lb
    lr(1) = 1
    l = 1
    est = me%gg(aa(l)+2._wp*hh(l),2._wp*hh(l))
    k = 8
    area = abs(est)
    ef = 0.5 !one_half
    mxl = 0

    !compute refined estimates, estimate the error, etc.
    main : do

        gl = me%gg(aa(l)+hh(l),hh(l))
        gr(l) = me%gg(aa(l)+3._wp*hh(l),hh(l))
        k = k + 16
        area = area + (abs(gl)+abs(gr(l))-abs(est))
        glr = gl + gr(l)
        ee = abs(est-glr)*ef
        ae = max(eps*area,tol*abs(glr))
        if (ee-ae > 0_wp) then
            !consider the left half of this level
            if (k > kmx) lmx = kml
            if (l >= lmx) then
                mxl = 1
            else
                l = l + 1
                eps = eps*0.5_wp
                ef = ef/sq2
                hh(l) = hh(l-1)*0.5_wp
                lr(l) = -1
                aa(l) = aa(l-1)
                est = gl
                cycle main
            end if
        end if

        err = err + (est-glr)
        if (lr(l) > 0) then
            !return one level
            ans = glr
            do
                if (l <= 1) exit main ! finished
                l = l - 1
                eps = eps*2._wp
                ef = ef*sq2
                if (lr(l) <= 0) then
                    vl(l) = vl(l+1) + ans
                    est = gr(l-1)
                    lr(l) = 1
                    aa(l) = aa(l) + 4._wp*hh(l)
                    cycle main
                end if
                ans = vl(l+1) + ans
            end do
        else
            !proceed to right half at this level
            vl(l) = glr
            est = gr(l-1)
            lr(l) = 1
            aa(l) = aa(l) + 4._wp*hh(l)
            cycle main
        end if

    end do main

    if ((mxl/=0) .and. (abs(err)>2._wp*tol*area)) ierr = 2 ! ans is probably insufficiently accurate

    return
    
    end subroutine dgauss_generic

!
!  6-point method.
!
!### See also
!  * Coefficients from:
!    http://processingjs.nihongoresources.com/bezierinfo/legendre-gauss-values.php  

    function g6(me,x,h) result(f)

      implicit none

      class(isos_analytical_elva_disk_load_class),intent(inout)  :: me
      real(wp), intent(in)                    :: x
      real(wp), intent(in)                    :: h
      real(wp)                                :: f

      !> abscissae:
      real(wp),dimension(3),parameter ::  a = [   0.6612093864662645136613995950199053470064485643&
                                                &951700708145267058521834966071431009442864037464&
                                                &614564298883716392751466795573467722253804381723&
                                                &198010093367423918538864300079016299442625145884&
                                                &9024557188219703863032236201173523213570221879361&
                                                &8906974301231555871064213101639896769013566165126&
                                                &1150514997832_wp,&
                                                &0.2386191860831969086305017216807119354186106301&
                                                &400213501813951645742749342756398422492244272573&
                                                &491316090722230970106872029554530350772051352628&
                                                &872175189982985139866216812636229030578298770859&
                                                &440976999298617585739469216136216592222334626416&
                                                &400139367778945327871453246721518889993399000945&
                                                &408150514997832_wp,&
                                                &0.9324695142031520278123015544939946091347657377&
                                                &122898248725496165266135008442001962762887399219&
                                                &259850478636797265728341065879713795116384041921&
                                                &786180750210169211578452038930846310372961174632&
                                                &524612619760497437974074226320896716211721783852&
                                                &305051047442772222093863676553669179038880252326&
                                                &771150514997832_wp ]
    !> weights:
    real(wp),dimension(3),parameter ::  w = [   0.36076157304813860756983351383771611166152189274&
                                                &674548228973924023714003783726171832096220198881&
                                                &934794311720914037079858987989027836432107077678&
                                                &721140858189221145027225257577711260007323688285&
                                                &916316028951118005174081368554707448247248610118&
                                                &325993144981721640242558677752676819993095031068&
                                                &73150514997832_wp,&
                                                0.46791393457269104738987034398955099481165560576&
                                                &921053531162531996391420162039812703111009258479&
                                                &198230476626878975479710092836255417350295459356&
                                                &355927338665933648259263825590180302812735635025&
                                                &362417046193182590009975698709590053347408007463&
                                                &437682443180817320636917410341626176534629278889&
                                                &17150514997832_wp,&
                                                0.17132449237917034504029614217273289352682250148&
                                                &404398239863543979894576054234015464792770542638&
                                                &866975211652206987440430919174716746217597462964&
                                                &922931803144845206713510916832108437179940676688&
                                                &721266924855699404815942932735702498405343382418&
                                                &236324411837461039120523911904421970357029774978&
                                                &12150514997832_wp ]

    f = h * ( w(1)*( me%fun(x-a(1)*h) + me%fun(x+a(1)*h) ) + &
              w(2)*( me%fun(x-a(2)*h) + me%fun(x+a(2)*h) ) + &
              w(3)*( me%fun(x-a(3)*h) + me%fun(x+a(3)*h) ) )

    return
    
  end function g6

end module isostasy_benchmarks





subroutine r4mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R4MAT_PRINT prints an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = * )  title

  call r4mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine r4mat_print

subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R4MAT_PRINT_SOME prints some of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 4 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 4 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r4mat_print_some




subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! R4MAT_PRINT prints an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of rows in A.
!
!    Input, integer ( kind = 4 ) N, the number of columns in A.
!
!    Input, real ( kind = 4 ) A(M,N), the matrix.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = * )  title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end subroutine r8mat_print

subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R4MAT_PRINT_SOME prints some of an R4MAT.
!
!  Discussion:
!
!    An R4MAT is an array of R4 values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
!
!    Input, real ( kind = 4 ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
!
!    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer ( kind = 4 ), parameter :: incx = 5
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i2hi
  integer ( kind = 4 ) i2lo
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) j2hi
  integer ( kind = 4 ) j2lo
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( a(i,j) == real ( int ( a(i,j) ), kind = 4 ) ) then
          write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
        else
          write ( ctemp(j2), '(g14.6)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end subroutine r8mat_print_some

