
module isostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    ! Note: currently routine `isos_par_load` has dependency on the nml.f90 module 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 

    real(wp), parameter :: g  = 9.81                ! [m/s^2]
    real(wp), parameter :: pi = 3.14159265359

    real(wp), parameter :: rho_ice  = 910.0         ! [kg/m^3]
    real(wp), parameter :: rho_sw   = 1028.0        ! [kg/m^3] 
    real(wp), parameter :: rho_w    = 1000.0        ! [kg/m^3] 
    real(wp), parameter :: rho_a    = 3300.0        ! [kg/m^3] 3370 used by Coulon et al (2021)
!mmr----------------------------------------------------------------
    real(wp), parameter :: nu       = 0.5   !0.25   ! [-]   Poisson's ratio: 0.25 for Coulon et al (2021) | 0.5 for Bueler et al (2007) 
    real(wp), parameter :: E        = 66.0  !100.0  ! [GPa] Young's modulus | 100 for Coulon et al (2021) | 66  for Bueler et al (2007) 
    real(wp), parameter :: eta      = 1.e21         ! [Pa s] Asthenosphere's viscosity, Bueler et al (2007)
!mmr----------------------------------------------------------------    
    real(wp), parameter :: r_earth  = 6.378e6       ! [m]  Earth's radius, Coulon et al (2021) 
    real(wp), parameter :: m_earth  = 5.972e24      ! [kg] Earth's mass,   Coulon et al (2021) 
    
    type isos_param_class 
        integer            :: method            ! Type of isostasy to use
        real(wp)           :: dt_lith           ! [yr] Timestep to recalculate equilibrium lithospheric displacement
        real(wp)           :: dt_step           ! [yr] Timestep to recalculate bedrock uplift and rate
        real(wp)           :: He_lith           ! [km] Effective elastic thickness of lithosphere
        real(wp)           :: tau               ! [yr] Asthenospheric relaxation constant
        
        ! Internal parameters 
        real(wp) :: L_w                       ! [m] Lithosphere flexural length scale (for method=2)
        integer  :: nr                        ! [-] Radius of neighborhood for convolution, in number of grid points        
        real(wp) :: time_lith                 ! [yr] Current model time of last update of equilibrium lithospheric displacement
        real(wp) :: time_step                 ! [yr] Current model time of last update of bedrock uplift
!mmr----------------------------------------------------------------
        real(wp) :: mu                        ! [1/m] 2pi/L                                         
        real(wp) :: dx                        ! [m] Horizontal resolution                           

        ! Parameters for analytical ELVA disk solution - Alex - class apart?                              
        real(wp) :: kappa_min = 0.  ! [] Minimum kappa for gaussian quadrature       
        real(wp) :: kappa_max = 0.1 ! [] Maximum kappa for gaussian quadrature       
        real(wp) :: dk = 1.e-3      ! [] Step in kappa for gaussian quadrature       
!mmr----------------------------------------------------------------
    end type 

    type isos_state_class 
        
        real(wp), allocatable :: tau(:,:)           ! [yr] Asthenospheric relaxation timescale field
        real(wp), allocatable :: He_lith(:,:)       ! [m]  Effective elastic thickness of the lithosphere
        real(wp), allocatable :: D_lith(:,:)        ! [N-m] Lithosphere flexural rigidity
        
        real(wp), allocatable :: kei(:,:)           ! Kelvin function filter values 
        real(wp), allocatable :: G0(:,:)            ! Green's function values
       
        real(wp), allocatable :: z_bed(:,:)         ! Bedrock elevation         [m]
        real(wp), allocatable :: dzbdt(:,:)         ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: z_bed_ref(:,:)     ! Reference (unweighted) bedrock 
        real(wp), allocatable :: q0(:,:)            ! Reference load
        real(wp), allocatable :: w0(:,:)            ! Reference equilibrium displacement
        real(wp), allocatable :: q1(:,:)            ! Current load          
        real(wp), allocatable :: w1(:,:)            ! Current equilibrium displacement
!mmr----------------------------------------------------------------
!mmr Alex -  class apart?   (this is for the ELVA numerical solution)             
        real(wp), allocatable :: kappa(:,:)         ! sqrt(p^2 + q^2) as in Bueler et al 2007                      
        real(wp), allocatable :: beta(:,:)          ! beta as in Bueler et al 2007                                 
        real(wp), allocatable :: w2(:,:)            ! Current viscous equilibrium displacement
        
! mmr Alex  - class apart? (this is for tge ELVA disk analytical solution)      
        real(wp), allocatable    :: w2_ana(:,:)        ! Current analytical viscous equilibrium displacement                               
        real(wp), allocatable    :: kappa_mod(:)       ! module of kappa (for analytical disk-load solution)                               
        real(wp), allocatable    :: dist2c(:,:)        ! distance of individual points to disk center (for analytical disk-load solution)  
        real(wp), allocatable    :: r(:)               ! radial distance to disk center (for analytical disk-load solution)                
        integer(wp), allocatable :: lr(:,:)         ! index for remapping (for analytical disk-load solution)                           
!mmr----------------------------------------------------------------
    end type isos_state_class

    type isos_class
        type(isos_param_class) :: par
        type(isos_state_class) :: now      
    end type

    type isos_analytical_elva_disk_load_class 
       integer  :: ifunc   = 0                      ! choice of function f(x) to integrate (integrand of analytical solution, in case there are several) 
       integer  :: n_evals = 0                      ! number of function evaluations
       integer  :: method  = 6                      ! integration method (#points in gaussian quadrature)
       procedure(func_1d),pointer :: fun => null()  ! function f(x) to be integrated (at given time t)
       procedure(gauss_func),pointer :: g => null() ! the gauss quadrature formula to use
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


        contains

        private

        procedure :: dgauss_generic                                    ! core integration routine. refactored from
                                                                       ! SLATEC with selectable quadrature method
        procedure,public :: initialize => initialize_integration_class ! to set up the class
        procedure,public :: integrate  => integrate_1d                 ! to integrate function `fun`
        
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

     type(isos_analytical_elva_disk_load_class) :: ana        ! object containing all info relative to function to be integrated for analytical solution
       
    private
    public :: isos_class 
    public :: isos_init
    public :: isos_init_state 
    public :: isos_update 
    public :: isos_end  
    public :: isos_set_field

    public :: isos_analytical_elva_disk_load_class 

   
  contains
    
!mmr----------------------------------------------------------------    
!mmr     subroutine isos_init(isos,filename,nx,ny,dx)
    subroutine isos_init(isos,filename,nx,ny,dx,calc_analytical)
!mmr----------------------------------------------------------------    

        implicit none 

        type(isos_class), intent(OUT) :: isos 
        character(len=*), intent(IN)  :: filename 
        integer,  intent(IN) :: nx, ny 
        real(wp), intent(IN) :: dx
!mmr----------------------------------------------------------------    
        logical, intent(IN) :: calc_analytical 
!mmr----------------------------------------------------------------    
        ! Local variables
        real(wp) :: radius_fac
        real(wp) :: filter_scaling
        real(wp) :: D_lith_const
!mmr----------------------------------------------------------------    
        real(wp) :: r0      
        real(wp) :: h0      
        real(wp) :: eta    
!mmr----------------------------------------------------------------    

        
        ! First, load parameters from parameter file `filename`
        call isos_par_load(isos%par,filename)
        

        select case(isos%par%method)
            case(2)
                ! ELRA method is being used, which requires a constant
                ! value of L_w, D_Lith and thus He_lith everywhere 
                
                ! Calculate the flexural rigidity based on the 
                ! effective elastic thickness of the lithosphere (He_lith),
                ! Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
                D_lith_const = (E*1e9) * (isos%par%He_lith*1e3)**3 / (12.0*(1.0-nu**2))

                ! Calculate the flexural length scale
                ! (Coulon et al, 2021, Eq. in text after Eq. 3)
                ! Note: will be on the order of 100km
                isos%par%L_w = (D_lith_const / (rho_a*g))**0.25 

                ! Modify the relative radius to use for the regional filter
                ! depending on whether we want to include the forbulge at
                ! large radii. Large radius makes the code run slower though too. 
                ! Set filter_scaling to < 1.0 to adjust for values of 
                ! radius_fac <~ 5.0 

                ! Larger radius, no filter scaling needed
                ! radius_fac      = 6.0 
                ! filter_scaling  = 1.0 

                ! Smaller radius, filter scaling to improve accuracy
                radius_fac      = 4.0 
                filter_scaling  = 0.91 

                ! Calculate radius of grid points to use for regional elastic plate filter
                ! See Greve and Blatter (2009), Chpt 8, page 192 for methodology 
                ! and Le Muer and Huybrechts (1996). It seems that this value
                ! should be 5-6x radius of relative stiffness to capture the forebuldge
                ! further out from the depression near the center. 
                ! Note: previous implementation stopped at 400km, hard coded. 
                isos%par%nr = int(radius_fac*isos%par%L_w/dx)+1
                

                ! Now, initialize isos variables 
                call isos_allocate(isos%now,nx,ny,nr=isos%par%nr)
            

                ! Calculate the Kelvin function filter 
                call calc_kei_filter_2D(isos%now%kei,L_w=isos%par%L_w,dx=dx,dy=dx)

                ! Apply scaling to adjust magnitude of filter for radius of filter cutoff
                isos%now%kei = filter_scaling*isos%now%kei

                ! Calculate the reference Green's function values
                call calc_greens_function_scaling(isos%now%G0,isos%now%kei, &
                                                isos%par%L_w,D_lith_const,dx=dx,dy=dx)

                ! Populate 2D D_lith field for aesethics too
                isos%now%D_lith = D_lith_const 

                write(*,*) "isos_init:: summary"
                write(*,*) "    He_lith (km): ", isos%par%He_lith 
                write(*,*) "    D_lith (N m): ", D_lith_const
                write(*,*) "    L_w (m):      ", isos%par%L_w 
                write(*,*) "    nr:           ", isos%par%nr
                
                write(*,*) "    range(kei): ", minval(isos%now%kei),    maxval(isos%now%kei)
                write(*,*) "    range(G0):  ", minval(isos%now%G0),     maxval(isos%now%G0)

!mmr----------------------------------------------------------------
            case(4) 
                ! ELRA/ELVA method is being used, which require a constant 
                ! value of L_w, D_Lith and thus He_lith everywhere         
                
                ! Calculate the flexural rigidity based on the 
                ! effective elastic thickness of the lithosphere (He_lith),
                ! Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
                D_lith_const = (E*1e9) * (isos%par%He_lith*1e3)**3 / (12.0*(1.0-nu**2))

               ! Calculate the flexural length scale
                ! (Coulon et al, 2021, Eq. in text after Eq. 3)
                ! Note: will be on the order of 100km
                isos%par%L_w = (D_lith_const / (rho_a*g))**0.25 

                ! Modify the relative radius to use for the regional filter
                ! depending on whether we want to include the forbulge at
                ! large radii. Large radius makes the code run slower though too. 
                ! Set filter_scaling to < 1.0 to adjust for values of 
                ! radius_fac <~ 5.0 

                ! Larger radius, no filter scaling needed
                ! radius_fac      = 6.0 
                ! filter_scaling  = 1.0 

                ! Smaller radius, filter scaling to improve accuracy
                radius_fac      = 4.0 
                filter_scaling  = 0.91 

                ! Calculate radius of grid points to use for regional elastic plate filter
                ! See Greve and Blatter (2009), Chpt 8, page 192 for methodology 
                ! and Le Muer and Huybrechts (1996). It seems that this value
                ! should be 5-6x radius of relative stiffness to capture the forebuldge
                ! further out from the depression near the center. 
                ! Note: previous implementation stopped at 400km, hard coded. 
                isos%par%nr = int(radius_fac*isos%par%L_w/dx)+1
                

                ! Now, initialize isos variables 
                call isos_allocate(isos%now,nx,ny,nr=isos%par%nr)
            

                ! Calculate the Kelvin function filter 
                call calc_kei_filter_2D(isos%now%kei,L_w=isos%par%L_w,dx=dx,dy=dx)

                ! Apply scaling to adjust magnitude of filter for radius of filter cutoff
                isos%now%kei = filter_scaling*isos%now%kei
                
!mmr recheck - Elastic solution is different in Bueler et al. 2007, need to recalculate Green's coefficients
                
                ! Calculate the reference Green's function values
                call calc_greens_function_scaling(isos%now%G0,isos%now%kei, &
                                                isos%par%L_w,D_lith_const,dx=dx,dy=dx)

                ! Populate 2D D_lith field for aesethics too
                isos%now%D_lith = D_lith_const

                ! Store horizonal resolution
                isos%par%dx = dx
                
                ! Calculate parameters needed for elastic lithosphere viscous asthenosphere (ELVA)                                        
                ! solution as in Bueler et al 2007 (eq 11)                                                                                
!
                call calc_asthenosphere_viscous_params(nx,ny,isos%par%dx,isos%now%D_lith,isos%par%mu,isos%now%kappa,isos%now%beta)        

                if (calc_analytical) then   
                   call calc_analytical_asthenosphere_viscous_disk_params(nx,ny,isos%par%dx,isos%par%kappa_min,isos%par%kappa_max,isos%par%dk,&
                        isos%now%kappa_mod,isos%now%dist2c,isos%now%r,isos%now%lr)                                    

                   r0     = 1000.0e3 ! [m] recheck - include into routine?
                   h0     = 1000.0   ! [m] 
                   eta    = 1.e+21   ! [Pa s]

                endif

                write(*,*) "isos_init:: summary"
                write(*,*) "    He_lith (km): ", isos%par%He_lith 
                write(*,*) "    D_lith (N m): ", D_lith_const
                write(*,*) "    L_w (m):      ", isos%par%L_w 
                write(*,*) "    nr:           ", isos%par%nr
                
                write(*,*) "    range(kei): ", minval(isos%now%kei),   maxval(isos%now%kei)
                write(*,*) "    range(G0):  ", minval(isos%now%G0),  maxval(isos%now%G0)

                write(*,*) "    mu (1/m):  ", isos%par%mu                                                       
                write(*,*) "    dx (m):    ", isos%par%dx                                                       
                write(*,*) "    kappa_min: ", isos%par%kappa_min                                                
                write(*,*) "    kappa_max: ", isos%par%kappa_max                                                
                write(*,*) "    dk ( ):    ", isos%par%dk                                                       

                write(*,*) "    range(kappa): ", minval(isos%now%kappa),  maxval(isos%now%kappa)                
                write(*,*) "    range(beta): ", minval(isos%now%beta),  maxval(isos%now%beta)                   

                if (calc_analytical) then                
                   write(*,*) "    range(kappa_mod): ", minval(isos%now%kappa_mod),  maxval(isos%now%kappa_mod)    
                   call initialize_analytical_integrand(ana,r0,h0,D_lith_const,eta)
                endif           
!mmr----------------------------------------------------------------
            
             case DEFAULT 

                ! Set elastic length scale to zero (not used)
                isos%par%L_w = 0.0 

                ! Set a small number of points for radius to avoid allocating
                ! any large, unused arrays
                isos%par%nr = 1 

                ! Now, initialize isos variables 
                call isos_allocate(isos%now,nx,ny,nr=isos%par%nr)
                
                ! Set regional filter fields to zero (not used)
                isos%now%kei = 0.0 
                isos%now%G0  = 0.0
                
!mmr----------------------------------------------------------------
!mmr - Alex, these should go within case = 4 but isos is initialized just after                
                ! Set kappa to zero (only used for ELVA)                
                isos%now%kappa = 0.0                                    
                ! Set mu to zero (only used for ELVA)                   
                isos%par%mu = 0.0                                       
                ! Initializw (only used for ELVA)                   
                isos%par%dx = dx
!mmr----------------------------------------------------------------

             end select

         ! Set He_lith and tau to constant fields initially.
        ! If these will be spatially varying fields, they should
        ! be defined after calling `isos_init` and before calling
        ! `isos_init_state`. 
        isos%now%He_lith    = isos%par%He_lith      ! [m]
        isos%now%tau        = isos%par%tau          ! [yr]

        
        ! Intially ensure all variables are zero 
        isos%now%z_bed_ref  = 0.0
        isos%now%z_bed      = 0.0 
        isos%now%dzbdt      = 0.0 

        isos%now%w0         = 0.0   
        isos%now%w1         = 0.0
!mmr----------------------------------------------------------------        
        isos%now%w2         = 0.0 
        isos%now%w2_ana     = 0.0 
!mmr----------------------------------------------------------------        
        ! Set time to very large value in the future 
        isos%par%time_lith = 1e10 
        isos%par%time_step = 1e10
        
        write(*,*) "isos_init:: complete." 
        
        return 

    end subroutine isos_init

    subroutine isos_init_state(isos,z_bed,H_ice,z_sl,z_bed_ref,H_ice_ref,z_sl_ref,time,calc_analytical)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: z_bed(:,:)            ! [m] Current bedrock elevation 
        real(wp), intent(IN) :: H_ice(:,:)            ! [m] Current ice thickness  
        real(wp), intent(IN) :: z_sl(:,:)             ! [m] Current sea level 
        real(wp), intent(IN) :: z_bed_ref(:,:)        ! [m] Reference bedrock elevation (with known load)
        real(wp), intent(IN) :: H_ice_ref(:,:)        ! [m] Reference ice thickness (associated with reference z_bed)
        real(wp), intent(IN) :: z_sl_ref(:,:)         ! [m] Reference sea level (associated with reference z_bed)
        real(wp), intent(IN) :: time                  ! [a] Initial time 
!mmr----------------------------------------------------------------
        logical :: calc_analytical
!mmr----------------------------------------------------------------        
        ! Store initial bedrock field 
        isos%now%z_bed = z_bed 
        
        ! Store reference bedrock field
        isos%now%z_bed_ref = z_bed_ref 
        isos%now%dzbdt     = 0.0 

        
        select case(isos%par%method)

            case(0,1) 
                ! 0: Steady-state lithospheric depression 
                ! 1: LLRA - Local lithosphere, relaxing asthenosphere

                ! Local lithosphere (LL)
                call calc_litho_local(isos%now%w0,isos%now%q0,z_bed_ref,H_ice_ref,z_sl_ref)
                
            case(2)
                ! 2: ELRA - Elastic lithosphere, relaxing asthenosphere

                ! Elastic lithosphere
                call calc_litho_regional(isos%now%w0,isos%now%q0,z_bed_ref,H_ice_ref,z_sl_ref,isos%now%G0)

            case(3)
                ! 3: EGIA - 2D Elastic lithosphere, relaxing asthenosphere - to do!

                ! Re-calculate the flexural rigidity based on the spatially variable
                ! effective elastic thickness of the lithosphere (He_lith),
                ! Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
                ! Note: this 2D field is only needed for the EGIA model (method=3) below. 
                isos%now%D_lith = (E*1e9) * (isos%now%He_lith*1e3)**3 / (12.0*(1.0-nu**2))


                ! Elastic lithosphere
                !call calc_litho_regional(isos%now%w0,isos%now%q0,z_bed_ref,H_ice_ref,z_sl_ref,isos%now%G0)
                
                write(*,*) "isos_init_state:: Error: method=3 not implemented yet."
                stop
                
!mmr----------------------------------------------------------------                
            case(4) 
                ! ELVA - viscous half-space asthenosphere overlain by                                       
                ! elastic plate lithosphere with uniform constants                                          

               ! Elastic lithosphere (EL)

!mmr recheck - Alex - ignored for the moment, need to recalculate Green's coefficients (solution is different in Bueler et al 2007)
               
!mmr               call calc_litho_regional(isos%now%w0,isos%now%q0,z_bed_ref,H_ice_ref,z_sl_ref,isos%now%G0)
!mmr               call calc_litho_local(isos%now%w0,isos%now%q0,z_bed_ref,H_ice_ref,z_sl_ref)  !mmr
                                             
                ! Viscous (half-space) asthenosphere (VA)                                                   

!mmr recheck - calculate plans here forth and back?

!mmr               call make_fft_plans(isos%now%q0,plan_fwd,plan_bck) 
!mmr----------------------------------------------------------------               
               
         end select 

        ! Define initial time of isostasy model 
        ! (set time_lith earlier, so that it is definitely updated on the first timestep)
        isos%par%time_lith = time - isos%par%dt_lith
        isos%par%time_step = time 

        ! Call isos_update to diagnose rate of change
        ! (no change to z_bed will be applied since isos%par%time==time)
!mmr----------------------------------------------------------------
!mmr         call isos_update(isos,H_ice,z_sl,time)
        call isos_update(isos,H_ice,z_sl,time,calc_analytical) 
!mmr----------------------------------------------------------------               
        write(*,*) "isos_init_state:: "
        write(*,*) "  Initial time:   ", isos%par%time_step 
        write(*,*) "  range(He_lith): ", minval(isos%now%He_lith), maxval(isos%now%He_lith)
        write(*,*) "  range(tau):     ", minval(isos%now%tau),     maxval(isos%now%tau)
        write(*,*) "  range(w0):      ", minval(isos%now%w0),      maxval(isos%now%w0)
        write(*,*) "  range(w1):      ", minval(isos%now%w1),      maxval(isos%now%w1)
        write(*,*) "  range(z_bed):   ", minval(isos%now%z_bed),   maxval(isos%now%z_bed)
!mmr----------------------------------------------------------------         
        write(*,*) "  range(w2):      ", minval(isos%now%w2),      maxval(isos%now%w2)      
        write(*,*) "  range(w2_ana):  ", minval(isos%now%w2_ana),  maxval(isos%now%w2_ana)  
!mmr----------------------------------------------------------------

        ! Make sure tau does not contain zero values, if so
        ! output an error for diagnostics. Don't kill the program
        ! so that external program can still write output as needed. 
        if (minval(isos%now%tau) .le. 0.0) then 
            write(error_unit,*) "isos_init_state:: Error: tau initialized with zero values present. &
            &This will lead to the model crashing."
            !stop 
        end if 

        return 

    end subroutine isos_init_state

!mmr----------------------------------------------------------------
!mmr     subroutine isos_update(isos,H_ice,z_sl,time)
    subroutine isos_update(isos,H_ice,z_sl,time,calc_analytical) 
!mmr----------------------------------------------------------------
        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: H_ice(:,:)        ! [m] Current ice thickness 
        real(wp), intent(IN) :: z_sl(:,:)         ! [m] Current sea level 
        real(wp), intent(IN) :: time              ! [a] Current time  
!mmr----------------------------------------------------------------
        logical, intent(IN)   :: calc_analytical
!mmr----------------------------------------------------------------
        ! Local variables 
        real(wp) :: dt, dt_now  
        integer  :: n, nstep 
        logical  :: update_equil
!mmr----------------------------------------------------------------
! mmr Alex - introduced update analytical solution to reduce computational time
        logical  :: update_analytical 
!mmr----------------------------------------------------------------       
        
        ! Step 0: determine current timestep and number of iterations
        dt = time - isos%par%time_step 

        ! Get maximum number of iterations needed to reach desired time
        nstep = ceiling( (time - isos%par%time_step) / isos%par%dt_step )
        nstep = max(nstep,1)

        ! Loop over iterations until maximum time is reached
        do n = 1, nstep 

            ! Get current dt (either total time or maximum allowed timestep)
            dt_now = min(time-isos%par%time_step, isos%par%dt_step)

            ! Only update equilibrium displacement height (w1) if 
            ! enough time has passed (to save on computations)
            if ( (isos%par%time_step+dt_now) - isos%par%time_lith .ge. isos%par%dt_lith) then
                update_equil = .TRUE. 
            else 
                update_equil = .FALSE. 
            end if 

!mmr----------------------------------------------------------------
           ! Only update analytical displacement height (w_ana) if 
           ! enough time has passed (to save on computations)
             if ( (isos%par%time_step+dt_now) .ge. time) then 
                 update_analytical = .TRUE.                   
             else                                             
                 update_analytical = .FALSE.                  
             end if 
!mmr----------------------------------------------------------------
            
            ! Step 1: diagnose equilibrium displacement and rate of bedrock uplift
            select case(isos%par%method)

                case(0)
                    ! Steady-state lithosphere

                    isos%now%w1    = isos%now%w0
                    isos%now%dzbdt = 0.0 

                case(1)
                    ! Local lithosphere, relaxing asthenosphere (LLRA)

                    ! Local lithosphere (LL)
                    ! (update every time because it is cheap)
                    call calc_litho_local(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl)

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt,isos%now%z_bed,isos%now%z_bed_ref, &
                                                    w_b=isos%now%w1-isos%now%w0,tau=isos%now%tau)

                case(2)
                    ! Elastic lithosphere, relaxing asthenosphere (ELRA)
                    
                    ! Elastic lithosphere (EL)
                    if (update_equil) then 
                        call calc_litho_regional(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%now%G0)
                    end if 

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt,isos%now%z_bed,isos%now%z_bed_ref, &
                                                    w_b=isos%now%w1-isos%now%w0,tau=isos%now%tau)

                case(3) 
                    ! Elementary GIA model (spatially varying ELRA with geoid - to do!)

                    ! 2D elastic lithosphere (2DEL)
                    if (update_equil) then
                       !call calc_litho_regional(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%now%G0)
 
                        write(*,*) "isos_update:: Error: method=3 not implemented yet."
                        stop 

                    end if 

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt,isos%now%z_bed,isos%now%z_bed_ref, &
                                                    w_b=isos%now%w1-isos%now%w0,tau=isos%now%tau)

!mmr----------------------------------------------------------------
                 case(4) 
                    ! ELVA - viscous half-space asthenosphere overlain by                                                   
                    ! elastic plate lithosphere with uniform constants                                                      
                    
                    ! Elastic lithosphere (EL)                                                                              
                    if (update_equil) then                                                                              

!mmr recheck - Alex - ignored for the moment, need to recalculate Green's coefficients (solution is different in Bueler et al 2007)
                       
!mmr                    call calc_litho_regional(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%now%G0)              
!mmr                    call calc_litho_local(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl)                         

                    end if
                    
                    ! Viscous (half-space) asthenosphere                                                                             

                    isos%now%q1 = rho_ice*g*H_ice !mmr recheck - imposing load here because there's no EL component

                    call calc_asthenosphere_viscous(isos%now%dzbdt,isos%now%q1,isos%now%w2,    &                            
                          dt_now,isos%par%mu,isos%now%kappa,isos%now%beta)                                                  

!mmr Alex  - introduced to reduce computational time in analytical solution
                    if (calc_analytical .and. update_analytical) then                                                       
                       call calc_analytical_asthenosphere_viscous_disk(ana,isos%par%dx,   &                                 
                            isos%par%time_step + dt_now, &                                   !mmr  Alex - this is not very clean (from me)
                            isos%now%kappa_mod,isos%now%dist2c,isos%now%r,isos%now%lr,isos%now%w2_ana)  
                     endif
!mmr----------------------------------------------------------------
                 end select

            if (update_equil) then 
                ! Update current time_lith value 
                isos%par%time_lith = isos%par%time_step + dt_now 
            end if 

            ! Step 2: update bedrock elevation and current model time
            if (dt_now .gt. 0.0) then

               isos%now%z_bed = isos%now%z_bed + isos%now%dzbdt*dt_now !mmr recheck -  contains only viscous component 
               
               isos%par%time_step = isos%par%time_step + dt_now
               
            end if

            ! Write a summary to check timestepping
            !write(*,*) "isos: ", n, time, dt_now, isos%par%time_step, isos%par%time_lith, update_analytical

            if ( abs(time-isos%par%time_step) .lt. 1e-5) then 
                ! Desired time has reached, exit the loop 
                isos%par%time_step = time 
                exit 
            end if 

         end do

        return 

    end subroutine isos_update

    subroutine isos_end(isos)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 

        ! Deallocate isos object 
        call isos_deallocate(isos%now)

        return 

    end subroutine isos_end

    subroutine isos_par_load(par,filename)

        use nml 

        implicit none

        type(isos_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename 

        call nml_read(filename,"isostasy","method",         par%method)
        call nml_read(filename,"isostasy","dt_lith",        par%dt_lith)
        call nml_read(filename,"isostasy","dt_step",        par%dt_step)
        call nml_read(filename,"isostasy","He_lith",        par%He_lith)
        call nml_read(filename,"isostasy","tau",            par%tau)
        
        return

    end subroutine isos_par_load

    subroutine isos_allocate(now,nx,ny,nr)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nr  

        ! Local variables
        integer :: nfilt 

        nfilt = 2*nr+1 

        ! First ensure arrays are not allocated
        call isos_deallocate(now)

        ! Allocate arrays
        allocate(now%tau(nx,ny))
        allocate(now%He_lith(nx,ny))
        allocate(now%D_lith(nx,ny))

        allocate(now%kei(nfilt,nfilt))
        allocate(now%G0(nfilt,nfilt))
!mmr----------------------------------------------------------------        
        allocate(now%kappa(nx,ny)) 
        allocate(now%beta(nx,ny))  
!mmr----------------------------------------------------------------                
        allocate(now%z_bed(nx,ny))
        allocate(now%dzbdt(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        allocate(now%q0(nx,ny))
        allocate(now%w0(nx,ny))
        allocate(now%q1(nx,ny))
        allocate(now%w1(nx,ny))
!mmr----------------------------------------------------------------        
        allocate(now%w2(nx,ny))     
        allocate(now%w2_ana(nx,ny)) 
!mmr----------------------------------------------------------------                 
        return 

    end subroutine isos_allocate

    subroutine isos_deallocate(now)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 

        if (allocated(now%tau))         deallocate(now%tau)
        if (allocated(now%He_lith))     deallocate(now%He_lith)
        if (allocated(now%D_lith))      deallocate(now%D_lith)
        
        if (allocated(now%kei))         deallocate(now%kei)
        if (allocated(now%G0))          deallocate(now%G0)
!mmr----------------------------------------------------------------        
        if (allocated(now%kappa))       deallocate(now%kappa) 
        if (allocated(now%beta))        deallocate(now%beta)  
!mmr----------------------------------------------------------------                
        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%dzbdt))       deallocate(now%dzbdt)
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)
        if (allocated(now%q0))          deallocate(now%q0)
        if (allocated(now%w0))          deallocate(now%w0)
        if (allocated(now%q1))          deallocate(now%q1)
        if (allocated(now%w1))          deallocate(now%w1)
!mmr----------------------------------------------------------------        
        if (allocated(now%w2))          deallocate(now%w2)     
        if (allocated(now%w2_ana))      deallocate(now%w2_ana) 
!mmr----------------------------------------------------------------                
        return 

    end subroutine isos_deallocate


    subroutine isos_set_field(var,var_values,mask_values,mask,dx,sigma)
        ! Impose multiple var values according to the corresponding
        ! locations provided in a mask. 
        ! Additionally impose Gaussian smoothing via the
        ! smoothing radius sigma.

        implicit none 

        real(wp), intent(OUT) :: var(:,:) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:,:) 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: sigma 

        ! Local variables
        integer :: j, n 

        ! Safety check 
        if (sigma .le. dx) then 
            write(*,*) "isos_set_field:: Error: sigma must be larger than dx."
            write(*,*) "dx    = ", dx 
            write(*,*) "sigma = ", sigma 
            stop 
        end if 

        ! Determine how many values should be assigned
        n = size(var_values,1)

        ! Initially set var=0 everywhere
        var = 0.0_wp 

        ! Loop over unique var values and assign them
        ! to correct regions of domain.
        do j = 1, n 

            where(mask .eq. mask_values(j)) var = var_values(j)

        end do

        ! Apply Gaussian smoothing as desired
        call smooth_gauss_2D(var,dx=dx,sigma=sigma)
        
        return

    end subroutine isos_set_field

    subroutine smooth_gauss_2D(var,dx,sigma,mask_apply,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2, k 
        integer  :: imx, ipx, jmx, jpx 
        real(wp), allocatable :: filter0(:,:), filter(:,:) 
        real(wp), allocatable :: var_old(:,:) 
        logical,  allocatable :: mask_apply_local(:,:) 
        logical,  allocatable :: mask_use_local(:,:) 

        real(wp), allocatable :: var_ext(:,:), var_ref_ext(:,:) 

        nx    = size(var,1)
        ny    = size(var,2)

        ! Determine half-width of filter as 2-sigma
        n2 = ceiling(2.0_wp*sigma / dx)

        ! Get total number of points for filter window in each direction
        n = 2*n2+1
        
        allocate(var_old(nx,ny))
        allocate(mask_apply_local(nx,ny))
        allocate(mask_use_local(nx,ny))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        allocate(var_ext(-n2:nx+n2,-n2:ny+n2))
        allocate(var_ref_ext(-n2:nx+n2,-n2:ny+n2))
        
        ! Check whether mask_apply is available 
        if (present(mask_apply)) then 
            ! use mask_use to define neighborhood points
            
            mask_apply_local = mask_apply 

        else
            ! Assume that everywhere should be smoothed

            mask_apply_local = .TRUE.
        
        end if

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply_local
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = var 

        var_ref_ext = -9999.0 

        var_ref_ext(1:nx,1:ny) = var 
        do i = 0, -n2, -1
            k = -i+1
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 
        do i = nx+1, nx+n2 
            k = nx + ((nx+1)-i)
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 

        do j = 0, -n2, -1
            k = -j+1
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 
        do j = ny+1, ny+n2 
            k = ny + ((ny+1)-j)
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 

        if (count(var_ref_ext .eq. -9999.0) .gt. 0) then 
            write(*,*) "Missing points!"
            stop 
        end if 

        do j = 1, ny
        do i = 1, nx


            !if (mask_apply_local(i,j)) then 
                ! Apply smoothing to this point 

                ! Limit filter input to neighbors of interest
                filter = filter0 
                !where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2) ) filter = 0.0

                ! If neighbors are available, normalize and perform smoothing  
                if (sum(filter) .gt. 0.0) then 
                    filter = filter/sum(filter)
                    var_ext(i,j) = sum(var_ref_ext(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            !end if 

        end do 
        end do 

        ! Get variable on normal grid 
        var = var_ext(1:nx,1:ny)

        return 

    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,  intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
        integer  :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return 

    end function gauss_values

    subroutine calc_greens_function_scaling(G0,kei2D,L_w,D_lith,dx,dy)
        ! The Green's function (Eq. 3 of Coulon et al, 2021)
        ! gives displacement G in [m] as a function of the distance
        ! r from the point load P_b [Pa]. 

        ! Here G0 is calculated, which is G without including the point load.
        ! G0 has units of [m N-1]. 
        ! This can then be multiplied with the actual magnitude of the
        ! point load to obtain G.
        ! G = G0 * P_b = [m N-1] * [Pa] = [m]. 

        ! Note that L_w contains information about rho_a. 

        implicit none

        real(wp), intent(OUT) :: G0(:,:) 
        real(wp), intent(IN)  :: kei2D(:,:) 
        real(wp), intent(IN)  :: L_w 
        real(wp), intent(IN)  :: D_lith 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        
        G0 = -L_w**2 / (2.0*pi*D_lith) * kei2D * (dx*dy)

        return

    end subroutine calc_greens_function_scaling

    ! === isos physics routines ======================================

    elemental subroutine calc_asthenosphere_relax(dzbdt,z_bed,z_bed_ref,w_b,tau)
        ! Calculate rate of change of vertical bedrock height
        ! from a relaxing asthenosphere.

        implicit none

        real(wp), intent(OUT) :: dzbdt 
        real(wp), intent(IN)  :: z_bed 
        real(wp), intent(IN)  :: z_bed_ref
        real(wp), intent(IN)  :: w_b        ! w_b = w1-w0
        real(wp), intent(IN)  :: tau

        dzbdt = -((z_bed-z_bed_ref) + w_b) / tau
        
        return

    end subroutine calc_asthenosphere_relax

    elemental subroutine calc_litho_local(w,q,z_bed,H_ice,z_sl)
        ! Calculate the local lithospheric loading from ice or ocean weight 
        ! in units of [Pa] and local equilibrium displacement w [m].

        implicit none 

        real(wp), intent(OUT) :: w
        real(wp), intent(OUT) :: q
        real(wp), intent(IN)  :: z_bed, H_ice, z_sl 

        if (rho_ice*H_ice.ge.rho_sw*(z_sl-z_bed)) then
            ! Ice or land

            q = rho_ice*g*H_ice

        else
            ! Ocean

            q = rho_sw*g*(z_sl-z_bed)

        end if

        ! Scale to get local displacement given the load q
        w = q / (rho_a*g) 

        return 

    end subroutine calc_litho_local

    subroutine calc_litho_regional(w1,q1,z_bed,H_ice,z_sl,G0)
        ! Calculate the load on the lithosphere as
        ! distributed on an elastic plate. 

        implicit none

        real(wp), intent(INOUT) :: w1(:,:)      ! [m] Lithospheric displacement
        real(wp), intent(INOUT) :: q1(:,:)      ! [Pa] Lithospheric load
        real(wp), intent(IN)    :: z_bed(:,:)   ! [m] Bed elevation
        real(wp), intent(IN)    :: H_ice(:,:)   ! [m] Ice thickness 
        real(wp), intent(IN)    :: z_sl(:,:)    ! [m] Sea level 
        real(wp), intent(IN)    :: G0(:,:)      ! Regional filter function
        
        ! Calculate local lithospheric load and displacement first
        call calc_litho_local(w1,q1,z_bed,H_ice,z_sl)

        ! Convolve the estimated point load with the regional
        ! filter to obtain the distributed load w1. 
        call convolve_load_elastic_plate(w1,q1,G0)

        return

    end subroutine calc_litho_regional

    subroutine convolve_load_elastic_plate(w1,q1,G0)
        ! Spread the load q1 [Pa] from each point in the grid
        ! via the regional Green's function scaling G0 [m N-1]

        implicit none

        real(wp), intent(OUT) :: w1(:,:)        ! [m] Lithospheric displacement
        real(wp), intent(IN)  :: q1(:,:)        ! [Pa] Lithospheric loading
        real(wp), intent(IN)  :: G0(:,:)        ! Regional scaling filter

        ! Local variables
        !integer :: ip, jp, lpx, lpy
        real(wp), allocatable :: q1_ext(:,:)
        real(wp), allocatable :: w_reg(:,:)

        integer :: i, j, nx ,ny, nr

        nx = size(w1,1)
        ny = size(w1,2)

        ! Size of regional neighborhood 
        nr = (size(G0,1)-1)/2 

        ! Populate load on extended grid
        allocate(q1_ext(1-nr:nx+nr,1-nr:ny+nr))

        ! First fill in main grid points with current point load
        q1_ext(1:nx,1:ny) = q1 

        ! Populate the extended grid points
        do i = 1, nx
            q1_ext(i,1-nr:0)=q1_ext(i,1)
            q1_ext(i,ny+1:ny+nr)=q1_ext(i,ny)
        end do
        do j = 1, ny
            q1_ext(1-nr:0,j)=q1_ext(1,j)
            q1_ext(NX+1:NX+nr,j)=q1_ext(nx,j)
        end do
        
        ! Populate the extended grid corner points     
        q1_ext(1-nr:0,1-nr:0)         = q1_ext(1,1)
        q1_ext(1-nr:0,ny+1:ny+nr)     = q1_ext(1,ny)
        q1_ext(nx+1:nx+nr,1-nr:0)     = q1_ext(nx,1)
        q1_ext(nx+1:nx+nr,ny+1:ny+nr) = q1_ext(nx,ny)

        ! ----- allocation de w_reg  et de croix -----------

        allocate(w_reg(-nr:nr,-nr:nr))

        do j = 1, ny
        do i = 1, nx

            ! Apply the neighborhood scaling to the deflection in the neighborhood
            w_reg = G0 * q1_ext(i-nr:i+nr,j-nr:j+nr)

            ! Sum to get total deflection at current point due to all neighbors
            w1(i,j) = sum(w_reg)

        end do
        end do
        
        return

    end subroutine convolve_load_elastic_plate

! ========================================================
!
! Functions related to the zero-order Kelvin function
!
! ========================================================

    subroutine calc_kei_filter_2D(filt,L_w,dx,dy)
        ! Calculate 2D Kelvin function (kei) kernel

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
        real(wp), intent(IN)  :: L_w 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy   

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp), allocatable :: rn_vals(:) 
        real(wp), allocatable :: kei_vals(:) 
        
        real(wp) :: kei_test_0
        real(wp) :: kei_test_1 

        ! Get size of filter array and half-width
        n  = size(filt,1) 
        n2 = (n-1)/2 

        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_kei_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_kei_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop  
        end if 

        ! Use tabulated values saved directly into module
        ! to avoid external dependency on input file
        call load_kei_values(rn_vals,kei_vals)

        ! Load tabulated kei values from file
        ! filename="input/kelvin.res"
        ! call load_kei_values_file(rn_vals,kei_vals,filename)

        ! TESTING to compare inline kei value calculation 
        ! to tabulated values.
        ! So far, the calc_kei_value does not work.
if (.FALSE.) then 
        do i = 1, size(rn_vals,1)
            r = rn_vals(i)*L_w
            kei_test_0 = get_kei_value(r,L_w,rn_vals,kei_vals)
            kei_test_1 = calc_kei_value(r,L_w)
            write(*,*) "kei: ", r, rn_vals(i), kei_test_0, kei_test_1 
        end do 
        stop 
end if 

        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center, normalized by L_w
        ! and impose correct Kelvin function value. 

        do j = -n2, n2 
        do i = -n2, n2

            x  = i*dx 
            y  = j*dy 
            r  = sqrt(x**2+y**2)

            ! Get actual index of array
            i1 = i+1+n2 
            j1 = j+1+n2 

            ! Get correct kei value for this point
            filt(i1,j1) = get_kei_value(r,L_w,rn_vals,kei_vals)

        end do 
        end do 
        
        return 

    end subroutine calc_kei_filter_2D

    function get_kei_value(r,L_w,rn_vals,kei_vals) result(kei)

        implicit none

        real(wp), intent(IN) :: r           ! [m] Radius from point load 
        real(wp), intent(IN) :: L_w         ! [m] Flexural length scale
        real(wp), intent(IN) :: rn_vals(:)  ! [-] Tabulated normalised radius values (r/L_w)
        real(wp), intent(IN) :: kei_vals(:) ! [-] Tabulated normalised Kelvin function values
        real(wp) :: kei

        ! Local variables 
        integer :: k, n 
        real(wp) :: rn_now 
        real(wp) :: wt 

        n = size(rn_vals,1) 

        ! Get current normalized radius from point load
        rn_now = r / L_w 

        if (rn_now .gt. rn_vals(n)) then 

            kei = kei_vals(n)

        else 

            do k = 1, n-1
                if (rn_now .ge. rn_vals(k) .and. rn_now .lt. rn_vals(k+1)) exit
            end do 

            ! Linear interpolation to get current kei value
            kei = kei_vals(k) &
                + (rn_now-rn_vals(k))/(rn_vals(k+1)-rn_vals(k))*(kei_vals(k+1)-kei_vals(k))   

        end if 

        ! Diagnostics (write only works if function is changed to a subroutine!)
        !write(*,*) "get_kei_value: ", L_w, r, rn_now, k, rn_vals(k), kei_vals(k) 

        return

    end function get_kei_value

    function calc_kei_value(r,L_w) result(kei)
        ! This function is based on Greve and Blatter (2009), Eq. 8.34.
        ! Combined with 8.11, it should be possible to obtain
        ! an analytical expression for the Kelvin function (kei). 
        ! This could eventually be used as a replacement for loading 
        ! the tabulated values from a file. 

        ! So far, this approach is not giving the right answer. 
        
        implicit none 

        real(wp), intent(IN) :: r       ! [m] Radius from point load 
        real(wp), intent(IN) :: L_w     ! [m] Flexural length scale
        real(wp) :: kei 

        ! Local variables
        real(wp) :: alpha 
        real(wp) :: f_now 
        real(wp) :: fac 

        alpha = r / (sqrt(2.0)*L_w)
        f_now = exp(-alpha)*(cos(alpha)+sin(alpha))

        fac = sqrt(2.0)**3 * L_w * (pi/4.0)
        
        kei = fac * f_now 

        ! Note: this doesn't give the right values of kei!!!

        return

    end function calc_kei_value

    subroutine load_kei_values(rn,kei)
        ! Originally obtained from the IMSL special math library function AKEIO
        ! Table columns: x         kei(x)

        real(wp), allocatable, intent(OUT) :: rn(:) 
        real(wp), allocatable, intent(OUT) :: kei(:) 
        
        ! Local variables
        integer :: i 

        integer, parameter :: n = 1001
        real(wp) :: tbl(n*2)

        if (allocated(rn))  deallocate(rn)
        if (allocated(kei)) deallocate(kei) 

        allocate(rn(n))
        allocate(kei(n)) 

        tbl = [  0.0000,  -0.78540E+00, &
                 0.0100,  -0.78526E+00, &
                 0.0200,  -0.78490E+00, &
                 0.0300,  -0.78436E+00, &
                 0.0400,  -0.78366E+00, &
                 0.0500,  -0.78283E+00, &
                 0.0600,  -0.78186E+00, &
                 0.0700,  -0.78077E+00, &
                 0.0800,  -0.77957E+00, &
                 0.0900,  -0.77826E+00, &
                 0.1000,  -0.77685E+00, &
                 0.1100,  -0.77534E+00, &
                 0.1200,  -0.77375E+00, &
                 0.1300,  -0.77206E+00, &
                 0.1400,  -0.77029E+00, &
                 0.1500,  -0.76844E+00, &
                 0.1600,  -0.76652E+00, &
                 0.1700,  -0.76452E+00, &
                 0.1800,  -0.76246E+00, &
                 0.1900,  -0.76032E+00, &
                 0.2000,  -0.75812E+00, &
                 0.2100,  -0.75587E+00, &
                 0.2200,  -0.75355E+00, &
                 0.2300,  -0.75117E+00, &
                 0.2400,  -0.74874E+00, &
                 0.2500,  -0.74625E+00, &
                 0.2600,  -0.74372E+00, &
                 0.2700,  -0.74113E+00, &
                 0.2800,  -0.73850E+00, &
                 0.2900,  -0.73582E+00, &
                 0.3000,  -0.73310E+00, &
                 0.3100,  -0.73034E+00, &
                 0.3200,  -0.72753E+00, &
                 0.3300,  -0.72469E+00, &
                 0.3400,  -0.72181E+00, &
                 0.3500,  -0.71889E+00, &
                 0.3600,  -0.71594E+00, &
                 0.3700,  -0.71295E+00, &
                 0.3800,  -0.70993E+00, &
                 0.3900,  -0.70688E+00, &
                 0.4000,  -0.70380E+00, &
                 0.4100,  -0.70069E+00, &
                 0.4200,  -0.69755E+00, &
                 0.4300,  -0.69439E+00, &
                 0.4400,  -0.69120E+00, &
                 0.4500,  -0.68799E+00, &
                 0.4600,  -0.68475E+00, &
                 0.4700,  -0.68149E+00, &
                 0.4800,  -0.67821E+00, &
                 0.4900,  -0.67490E+00, &
                 0.5000,  -0.67158E+00, &
                 0.5100,  -0.66824E+00, &
                 0.5200,  -0.66488E+00, &
                 0.5300,  -0.66150E+00, &
                 0.5400,  -0.65811E+00, &
                 0.5500,  -0.65470E+00, &
                 0.5600,  -0.65128E+00, &
                 0.5700,  -0.64784E+00, &
                 0.5800,  -0.64439E+00, &
                 0.5900,  -0.64093E+00, &
                 0.6000,  -0.63745E+00, &
                 0.6100,  -0.63396E+00, &
                 0.6200,  -0.63046E+00, &
                 0.6300,  -0.62696E+00, &
                 0.6400,  -0.62344E+00, &
                 0.6500,  -0.61991E+00, &
                 0.6600,  -0.61638E+00, &
                 0.6700,  -0.61284E+00, &
                 0.6800,  -0.60929E+00, &
                 0.6900,  -0.60574E+00, &
                 0.7000,  -0.60218E+00, &
                 0.7100,  -0.59861E+00, &
                 0.7200,  -0.59504E+00, &
                 0.7300,  -0.59146E+00, &
                 0.7400,  -0.58789E+00, &
                 0.7500,  -0.58431E+00, &
                 0.7600,  -0.58072E+00, &
                 0.7700,  -0.57714E+00, &
                 0.7800,  -0.57355E+00, &
                 0.7900,  -0.56996E+00, &
                 0.8000,  -0.56637E+00, &
                 0.8100,  -0.56278E+00, &
                 0.8200,  -0.55919E+00, &
                 0.8300,  -0.55560E+00, &
                 0.8400,  -0.55201E+00, &
                 0.8500,  -0.54842E+00, &
                 0.8600,  -0.54483E+00, &
                 0.8700,  -0.54125E+00, &
                 0.8800,  -0.53767E+00, &
                 0.8900,  -0.53409E+00, &
                 0.9000,  -0.53051E+00, &
                 0.9100,  -0.52694E+00, &
                 0.9200,  -0.52337E+00, &
                 0.9300,  -0.51980E+00, &
                 0.9400,  -0.51624E+00, &
                 0.9500,  -0.51269E+00, &
                 0.9600,  -0.50914E+00, &
                 0.9700,  -0.50559E+00, &
                 0.9800,  -0.50205E+00, &
                 0.9900,  -0.49852E+00, &
                 1.0000,  -0.49499E+00, &
                 1.0100,  -0.49147E+00, &
                 1.0200,  -0.48796E+00, &
                 1.0300,  -0.48445E+00, &
                 1.0400,  -0.48096E+00, &
                 1.0500,  -0.47746E+00, &
                 1.0600,  -0.47398E+00, &
                 1.0700,  -0.47050E+00, &
                 1.0800,  -0.46704E+00, &
                 1.0900,  -0.46358E+00, &
                 1.1000,  -0.46013E+00, &
                 1.1100,  -0.45669E+00, &
                 1.1200,  -0.45326E+00, &
                 1.1300,  -0.44984E+00, &
                 1.1400,  -0.44642E+00, &
                 1.1500,  -0.44302E+00, &
                 1.1600,  -0.43963E+00, &
                 1.1700,  -0.43625E+00, &
                 1.1800,  -0.43287E+00, &
                 1.1900,  -0.42951E+00, &
                 1.2000,  -0.42616E+00, &
                 1.2100,  -0.42282E+00, &
                 1.2200,  -0.41950E+00, &
                 1.2300,  -0.41618E+00, &
                 1.2400,  -0.41287E+00, &
                 1.2500,  -0.40958E+00, &
                 1.2600,  -0.40630E+00, &
                 1.2700,  -0.40303E+00, &
                 1.2800,  -0.39977E+00, &
                 1.2900,  -0.39653E+00, &
                 1.3000,  -0.39329E+00, &
                 1.3100,  -0.39007E+00, &
                 1.3200,  -0.38686E+00, &
                 1.3300,  -0.38367E+00, &
                 1.3400,  -0.38048E+00, &
                 1.3500,  -0.37731E+00, &
                 1.3600,  -0.37416E+00, &
                 1.3700,  -0.37101E+00, &
                 1.3800,  -0.36788E+00, &
                 1.3900,  -0.36477E+00, &
                 1.4000,  -0.36166E+00, &
                 1.4100,  -0.35858E+00, &
                 1.4200,  -0.35550E+00, &
                 1.4300,  -0.35244E+00, &
                 1.4400,  -0.34939E+00, &
                 1.4500,  -0.34635E+00, &
                 1.4600,  -0.34333E+00, &
                 1.4700,  -0.34033E+00, &
                 1.4800,  -0.33734E+00, &
                 1.4900,  -0.33436E+00, &
                 1.5000,  -0.33140E+00, &
                 1.5100,  -0.32845E+00, &
                 1.5200,  -0.32551E+00, &
                 1.5300,  -0.32259E+00, &
                 1.5400,  -0.31969E+00, &
                 1.5500,  -0.31680E+00, &
                 1.5600,  -0.31392E+00, &
                 1.5700,  -0.31106E+00, &
                 1.5800,  -0.30821E+00, &
                 1.5900,  -0.30538E+00, &
                 1.6000,  -0.30257E+00, &
                 1.6100,  -0.29976E+00, &
                 1.6200,  -0.29698E+00, &
                 1.6300,  -0.29421E+00, &
                 1.6400,  -0.29145E+00, &
                 1.6500,  -0.28871E+00, &
                 1.6600,  -0.28598E+00, &
                 1.6700,  -0.28327E+00, &
                 1.6800,  -0.28057E+00, &
                 1.6900,  -0.27789E+00, &
                 1.7000,  -0.27523E+00, &
                 1.7100,  -0.27258E+00, &
                 1.7200,  -0.26994E+00, &
                 1.7300,  -0.26732E+00, &
                 1.7400,  -0.26472E+00, &
                 1.7500,  -0.26213E+00, &
                 1.7600,  -0.25956E+00, &
                 1.7700,  -0.25700E+00, &
                 1.7800,  -0.25446E+00, &
                 1.7900,  -0.25193E+00, &
                 1.8000,  -0.24942E+00, &
                 1.8100,  -0.24692E+00, &
                 1.8200,  -0.24444E+00, &
                 1.8300,  -0.24197E+00, &
                 1.8400,  -0.23952E+00, &
                 1.8500,  -0.23709E+00, &
                 1.8600,  -0.23467E+00, &
                 1.8700,  -0.23226E+00, &
                 1.8800,  -0.22987E+00, &
                 1.8900,  -0.22750E+00, &
                 1.9000,  -0.22514E+00, &
                 1.9100,  -0.22280E+00, &
                 1.9200,  -0.22047E+00, &
                 1.9300,  -0.21816E+00, &
                 1.9400,  -0.21586E+00, &
                 1.9500,  -0.21358E+00, &
                 1.9600,  -0.21131E+00, &
                 1.9700,  -0.20906E+00, &
                 1.9800,  -0.20683E+00, &
                 1.9900,  -0.20461E+00, &
                 2.0000,  -0.20240E+00, &
                 2.0100,  -0.20021E+00, &
                 2.0200,  -0.19803E+00, &
                 2.0300,  -0.19587E+00, &
                 2.0400,  -0.19373E+00, &
                 2.0500,  -0.19160E+00, &
                 2.0600,  -0.18948E+00, &
                 2.0700,  -0.18738E+00, &
                 2.0800,  -0.18530E+00, &
                 2.0900,  -0.18323E+00, &
                 2.1000,  -0.18117E+00, &
                 2.1100,  -0.17913E+00, &
                 2.1200,  -0.17711E+00, &
                 2.1300,  -0.17510E+00, &
                 2.1400,  -0.17310E+00, &
                 2.1500,  -0.17112E+00, &
                 2.1600,  -0.16915E+00, &
                 2.1700,  -0.16720E+00, &
                 2.1800,  -0.16526E+00, &
                 2.1900,  -0.16334E+00, &
                 2.2000,  -0.16143E+00, &
                 2.2100,  -0.15954E+00, &
                 2.2200,  -0.15766E+00, &
                 2.2300,  -0.15579E+00, &
                 2.2400,  -0.15394E+00, &
                 2.2500,  -0.15211E+00, &
                 2.2600,  -0.15028E+00, &
                 2.2700,  -0.14847E+00, &
                 2.2800,  -0.14668E+00, &
                 2.2900,  -0.14490E+00, &
                 2.3000,  -0.14314E+00, &
                 2.3100,  -0.14138E+00, &
                 2.3200,  -0.13965E+00, &
                 2.3300,  -0.13792E+00, &
                 2.3400,  -0.13621E+00, &
                 2.3500,  -0.13452E+00, &
                 2.3600,  -0.13283E+00, &
                 2.3700,  -0.13117E+00, &
                 2.3800,  -0.12951E+00, &
                 2.3900,  -0.12787E+00, &
                 2.4000,  -0.12624E+00, &
                 2.4100,  -0.12463E+00, &
                 2.4200,  -0.12303E+00, &
                 2.4300,  -0.12144E+00, &
                 2.4400,  -0.11986E+00, &
                 2.4500,  -0.11830E+00, &
                 2.4600,  -0.11676E+00, &
                 2.4700,  -0.11522E+00, &
                 2.4800,  -0.11370E+00, &
                 2.4900,  -0.11219E+00, &
                 2.5000,  -0.11070E+00, &
                 2.5100,  -0.10921E+00, &
                 2.5200,  -0.10774E+00, &
                 2.5300,  -0.10629E+00, &
                 2.5400,  -0.10484E+00, &
                 2.5500,  -0.10341E+00, &
                 2.5600,  -0.10199E+00, &
                 2.5700,  -0.10059E+00, &
                 2.5800,  -0.99193E-01, &
                 2.5900,  -0.97812E-01, &
                 2.6000,  -0.96443E-01, &
                 2.6100,  -0.95086E-01, &
                 2.6200,  -0.93742E-01, &
                 2.6300,  -0.92410E-01, &
                 2.6400,  -0.91090E-01, &
                 2.6500,  -0.89782E-01, &
                 2.6600,  -0.88486E-01, &
                 2.6700,  -0.87202E-01, &
                 2.6800,  -0.85930E-01, &
                 2.6900,  -0.84670E-01, &
                 2.7000,  -0.83422E-01, &
                 2.7100,  -0.82185E-01, &
                 2.7200,  -0.80960E-01, &
                 2.7300,  -0.79747E-01, &
                 2.7400,  -0.78545E-01, &
                 2.7500,  -0.77354E-01, &
                 2.7600,  -0.76175E-01, &
                 2.7700,  -0.75007E-01, &
                 2.7800,  -0.73850E-01, &
                 2.7900,  -0.72705E-01, &
                 2.8000,  -0.71571E-01, &
                 2.8100,  -0.70447E-01, &
                 2.8200,  -0.69335E-01, &
                 2.8300,  -0.68234E-01, &
                 2.8400,  -0.67143E-01, &
                 2.8500,  -0.66064E-01, &
                 2.8600,  -0.64995E-01, &
                 2.8700,  -0.63937E-01, &
                 2.8800,  -0.62889E-01, &
                 2.8900,  -0.61852E-01, &
                 2.9000,  -0.60825E-01, &
                 2.9100,  -0.59809E-01, &
                 2.9200,  -0.58803E-01, &
                 2.9300,  -0.57808E-01, &
                 2.9400,  -0.56823E-01, &
                 2.9500,  -0.55848E-01, &
                 2.9600,  -0.54882E-01, &
                 2.9700,  -0.53927E-01, &
                 2.9800,  -0.52982E-01, &
                 2.9900,  -0.52047E-01, &
                 3.0000,  -0.51122E-01, &
                 3.0100,  -0.50206E-01, &
                 3.0200,  -0.49300E-01, &
                 3.0300,  -0.48404E-01, &
                 3.0400,  -0.47518E-01, &
                 3.0500,  -0.46641E-01, &
                 3.0600,  -0.45773E-01, &
                 3.0700,  -0.44915E-01, &
                 3.0800,  -0.44066E-01, &
                 3.0900,  -0.43226E-01, &
                 3.1000,  -0.42395E-01, &
                 3.1100,  -0.41574E-01, &
                 3.1200,  -0.40762E-01, &
                 3.1300,  -0.39958E-01, &
                 3.1400,  -0.39164E-01, &
                 3.1500,  -0.38379E-01, &
                 3.1600,  -0.37602E-01, &
                 3.1700,  -0.36834E-01, &
                 3.1800,  -0.36075E-01, &
                 3.1900,  -0.35324E-01, &
                 3.2000,  -0.34582E-01, &
                 3.2100,  -0.33849E-01, &
                 3.2200,  -0.33124E-01, &
                 3.2300,  -0.32407E-01, &
                 3.2400,  -0.31699E-01, &
                 3.2500,  -0.30999E-01, &
                 3.2600,  -0.30307E-01, &
                 3.2700,  -0.29623E-01, &
                 3.2800,  -0.28947E-01, &
                 3.2900,  -0.28279E-01, &
                 3.3000,  -0.27620E-01, &
                 3.3100,  -0.26968E-01, &
                 3.3200,  -0.26324E-01, &
                 3.3300,  -0.25688E-01, &
                 3.3400,  -0.25059E-01, &
                 3.3500,  -0.24438E-01, &
                 3.3600,  -0.23825E-01, &
                 3.3700,  -0.23219E-01, &
                 3.3800,  -0.22621E-01, &
                 3.3900,  -0.22030E-01, &
                 3.4000,  -0.21446E-01, &
                 3.4100,  -0.20870E-01, &
                 3.4200,  -0.20301E-01, &
                 3.4300,  -0.19739E-01, &
                 3.4400,  -0.19185E-01, &
                 3.4500,  -0.18637E-01, &
                 3.4600,  -0.18096E-01, &
                 3.4700,  -0.17563E-01, &
                 3.4800,  -0.17036E-01, &
                 3.4900,  -0.16516E-01, &
                 3.5000,  -0.16003E-01, &
                 3.5100,  -0.15496E-01, &
                 3.5200,  -0.14996E-01, &
                 3.5300,  -0.14503E-01, &
                 3.5400,  -0.14016E-01, &
                 3.5500,  -0.13536E-01, &
                 3.5600,  -0.13063E-01, &
                 3.5700,  -0.12595E-01, &
                 3.5800,  -0.12134E-01, &
                 3.5900,  -0.11680E-01, &
                 3.6000,  -0.11231E-01, &
                 3.6100,  -0.10789E-01, &
                 3.6200,  -0.10353E-01, &
                 3.6300,  -0.99224E-02, &
                 3.6400,  -0.94983E-02, &
                 3.6500,  -0.90801E-02, &
                 3.6600,  -0.86678E-02, &
                 3.6700,  -0.82614E-02, &
                 3.6800,  -0.78608E-02, &
                 3.6900,  -0.74659E-02, &
                 3.7000,  -0.70767E-02, &
                 3.7100,  -0.66932E-02, &
                 3.7200,  -0.63152E-02, &
                 3.7300,  -0.59428E-02, &
                 3.7400,  -0.55758E-02, &
                 3.7500,  -0.52143E-02, &
                 3.7600,  -0.48582E-02, &
                 3.7700,  -0.45075E-02, &
                 3.7800,  -0.41620E-02, &
                 3.7900,  -0.38217E-02, &
                 3.8000,  -0.34867E-02, &
                 3.8100,  -0.31567E-02, &
                 3.8200,  -0.28319E-02, &
                 3.8300,  -0.25121E-02, &
                 3.8400,  -0.21973E-02, &
                 3.8500,  -0.18875E-02, &
                 3.8600,  -0.15826E-02, &
                 3.8700,  -0.12825E-02, &
                 3.8800,  -0.98717E-03, &
                 3.8900,  -0.69663E-03, &
                 3.9000,  -0.41081E-03, &
                 3.9100,  -0.12965E-03, &
                 3.9200,   0.14689E-03, &
                 3.9300,   0.41886E-03, &
                 3.9400,   0.68630E-03, &
                 3.9500,   0.94926E-03, &
                 3.9600,   0.12078E-02, &
                 3.9700,   0.14619E-02, &
                 3.9800,   0.17117E-02, &
                 3.9900,   0.19572E-02, &
                 4.0000,   0.21984E-02, &
                 4.0100,   0.24354E-02, &
                 4.0200,   0.26682E-02, &
                 4.0300,   0.28969E-02, &
                 4.0400,   0.31216E-02, &
                 4.0500,   0.33421E-02, &
                 4.0600,   0.35587E-02, &
                 4.0700,   0.37713E-02, &
                 4.0800,   0.39800E-02, &
                 4.0900,   0.41848E-02, &
                 4.1000,   0.43858E-02, &
                 4.1100,   0.45830E-02, &
                 4.1200,   0.47764E-02, &
                 4.1300,   0.49662E-02, &
                 4.1400,   0.51522E-02, &
                 4.1500,   0.53346E-02, &
                 4.1600,   0.55135E-02, &
                 4.1700,   0.56887E-02, &
                 4.1800,   0.58605E-02, &
                 4.1900,   0.60288E-02, &
                 4.2000,   0.61936E-02, &
                 4.2100,   0.63551E-02, &
                 4.2200,   0.65131E-02, &
                 4.2300,   0.66679E-02, &
                 4.2400,   0.68194E-02, &
                 4.2500,   0.69676E-02, &
                 4.2600,   0.71126E-02, &
                 4.2700,   0.72544E-02, &
                 4.2800,   0.73931E-02, &
                 4.2900,   0.75287E-02, &
                 4.3000,   0.76613E-02, &
                 4.3100,   0.77908E-02, &
                 4.3200,   0.79173E-02, &
                 4.3300,   0.80408E-02, &
                 4.3400,   0.81614E-02, &
                 4.3500,   0.82792E-02, &
                 4.3600,   0.83940E-02, &
                 4.3700,   0.85061E-02, &
                 4.3800,   0.86153E-02, &
                 4.3900,   0.87218E-02, &
                 4.4000,   0.88256E-02, &
                 4.4100,   0.89267E-02, &
                 4.4200,   0.90252E-02, &
                 4.4300,   0.91210E-02, &
                 4.4400,   0.92142E-02, &
                 4.4500,   0.93049E-02, &
                 4.4600,   0.93930E-02, &
                 4.4700,   0.94786E-02, &
                 4.4800,   0.95618E-02, &
                 4.4900,   0.96426E-02, &
                 4.5000,   0.97209E-02, &
                 4.5100,   0.97969E-02, &
                 4.5200,   0.98705E-02, &
                 4.5300,   0.99418E-02, &
                 4.5400,   0.10011E-01, &
                 4.5500,   0.10078E-01, &
                 4.5600,   0.10142E-01, &
                 4.5700,   0.10205E-01, &
                 4.5800,   0.10265E-01, &
                 4.5900,   0.10323E-01, &
                 4.6000,   0.10379E-01, &
                 4.6100,   0.10433E-01, &
                 4.6200,   0.10485E-01, &
                 4.6300,   0.10534E-01, &
                 4.6400,   0.10582E-01, &
                 4.6500,   0.10628E-01, &
                 4.6600,   0.10672E-01, &
                 4.6700,   0.10714E-01, &
                 4.6800,   0.10754E-01, &
                 4.6900,   0.10792E-01, &
                 4.7000,   0.10829E-01, &
                 4.7100,   0.10863E-01, &
                 4.7200,   0.10896E-01, &
                 4.7300,   0.10927E-01, &
                 4.7400,   0.10957E-01, &
                 4.7500,   0.10984E-01, &
                 4.7600,   0.11010E-01, &
                 4.7700,   0.11034E-01, &
                 4.7800,   0.11057E-01, &
                 4.7900,   0.11078E-01, &
                 4.8000,   0.11097E-01, &
                 4.8100,   0.11115E-01, &
                 4.8200,   0.11132E-01, &
                 4.8300,   0.11146E-01, &
                 4.8400,   0.11160E-01, &
                 4.8500,   0.11172E-01, &
                 4.8600,   0.11182E-01, &
                 4.8700,   0.11191E-01, &
                 4.8800,   0.11199E-01, &
                 4.8900,   0.11205E-01, &
                 4.9000,   0.11210E-01, &
                 4.9100,   0.11213E-01, &
                 4.9200,   0.11215E-01, &
                 4.9300,   0.11216E-01, &
                 4.9400,   0.11216E-01, &
                 4.9500,   0.11214E-01, &
                 4.9600,   0.11211E-01, &
                 4.9700,   0.11207E-01, &
                 4.9800,   0.11202E-01, &
                 4.9900,   0.11195E-01, &
                 5.0000,   0.11188E-01, &
                 5.0100,   0.11179E-01, &
                 5.0200,   0.11169E-01, &
                 5.0300,   0.11158E-01, &
                 5.0400,   0.11146E-01, &
                 5.0500,   0.11133E-01, &
                 5.0600,   0.11119E-01, &
                 5.0700,   0.11103E-01, &
                 5.0800,   0.11087E-01, &
                 5.0900,   0.11070E-01, &
                 5.1000,   0.11052E-01, &
                 5.1100,   0.11033E-01, &
                 5.1200,   0.11013E-01, &
                 5.1300,   0.10992E-01, &
                 5.1400,   0.10970E-01, &
                 5.1500,   0.10947E-01, &
                 5.1600,   0.10924E-01, &
                 5.1700,   0.10899E-01, &
                 5.1800,   0.10874E-01, &
                 5.1900,   0.10848E-01, &
                 5.2000,   0.10821E-01, &
                 5.2100,   0.10794E-01, &
                 5.2200,   0.10765E-01, &
                 5.2300,   0.10736E-01, &
                 5.2400,   0.10706E-01, &
                 5.2500,   0.10676E-01, &
                 5.2600,   0.10644E-01, &
                 5.2700,   0.10612E-01, &
                 5.2800,   0.10579E-01, &
                 5.2900,   0.10546E-01, &
                 5.3000,   0.10512E-01, &
                 5.3100,   0.10477E-01, &
                 5.3200,   0.10442E-01, &
                 5.3300,   0.10406E-01, &
                 5.3400,   0.10370E-01, &
                 5.3500,   0.10333E-01, &
                 5.3600,   0.10295E-01, &
                 5.3700,   0.10257E-01, &
                 5.3800,   0.10218E-01, &
                 5.3900,   0.10179E-01, &
                 5.4000,   0.10139E-01, &
                 5.4100,   0.10099E-01, &
                 5.4200,   0.10058E-01, &
                 5.4300,   0.10017E-01, &
                 5.4400,   0.99754E-02, &
                 5.4500,   0.99333E-02, &
                 5.4600,   0.98907E-02, &
                 5.4700,   0.98477E-02, &
                 5.4800,   0.98043E-02, &
                 5.4900,   0.97605E-02, &
                 5.5000,   0.97163E-02, &
                 5.5100,   0.96717E-02, &
                 5.5200,   0.96268E-02, &
                 5.5300,   0.95814E-02, &
                 5.5400,   0.95358E-02, &
                 5.5500,   0.94897E-02, &
                 5.5600,   0.94434E-02, &
                 5.5700,   0.93967E-02, &
                 5.5800,   0.93498E-02, &
                 5.5900,   0.93025E-02, &
                 5.6000,   0.92550E-02, &
                 5.6100,   0.92071E-02, &
                 5.6200,   0.91590E-02, &
                 5.6300,   0.91107E-02, &
                 5.6400,   0.90621E-02, &
                 5.6500,   0.90132E-02, &
                 5.6600,   0.89641E-02, &
                 5.6700,   0.89149E-02, &
                 5.6800,   0.88653E-02, &
                 5.6900,   0.88156E-02, &
                 5.7000,   0.87657E-02, &
                 5.7100,   0.87156E-02, &
                 5.7200,   0.86653E-02, &
                 5.7300,   0.86149E-02, &
                 5.7400,   0.85643E-02, &
                 5.7500,   0.85135E-02, &
                 5.7600,   0.84626E-02, &
                 5.7700,   0.84116E-02, &
                 5.7800,   0.83604E-02, &
                 5.7900,   0.83091E-02, &
                 5.8000,   0.82577E-02, &
                 5.8100,   0.82062E-02, &
                 5.8200,   0.81546E-02, &
                 5.8300,   0.81029E-02, &
                 5.8400,   0.80511E-02, &
                 5.8500,   0.79993E-02, &
                 5.8600,   0.79473E-02, &
                 5.8700,   0.78953E-02, &
                 5.8800,   0.78433E-02, &
                 5.8900,   0.77912E-02, &
                 5.9000,   0.77390E-02, &
                 5.9100,   0.76868E-02, &
                 5.9200,   0.76346E-02, &
                 5.9300,   0.75824E-02, &
                 5.9400,   0.75301E-02, &
                 5.9500,   0.74778E-02, &
                 5.9600,   0.74256E-02, &
                 5.9700,   0.73733E-02, &
                 5.9800,   0.73210E-02, &
                 5.9900,   0.72687E-02, &
                 6.0000,   0.72165E-02, &
                 6.0100,   0.71643E-02, &
                 6.0200,   0.71121E-02, &
                 6.0300,   0.70599E-02, &
                 6.0400,   0.70078E-02, &
                 6.0500,   0.69557E-02, &
                 6.0600,   0.69036E-02, &
                 6.0700,   0.68517E-02, &
                 6.0800,   0.67997E-02, &
                 6.0900,   0.67479E-02, &
                 6.1000,   0.66961E-02, &
                 6.1100,   0.66443E-02, &
                 6.1200,   0.65927E-02, &
                 6.1300,   0.65411E-02, &
                 6.1400,   0.64896E-02, &
                 6.1500,   0.64382E-02, &
                 6.1600,   0.63869E-02, &
                 6.1700,   0.63357E-02, &
                 6.1800,   0.62846E-02, &
                 6.1900,   0.62336E-02, &
                 6.2000,   0.61827E-02, &
                 6.2100,   0.61320E-02, &
                 6.2200,   0.60813E-02, &
                 6.2300,   0.60308E-02, &
                 6.2400,   0.59804E-02, &
                 6.2500,   0.59301E-02, &
                 6.2600,   0.58800E-02, &
                 6.2700,   0.58300E-02, &
                 6.2800,   0.57801E-02, &
                 6.2900,   0.57304E-02, &
                 6.3000,   0.56808E-02, &
                 6.3100,   0.56313E-02, &
                 6.3200,   0.55820E-02, &
                 6.3300,   0.55329E-02, &
                 6.3400,   0.54839E-02, &
                 6.3500,   0.54351E-02, &
                 6.3600,   0.53865E-02, &
                 6.3700,   0.53380E-02, &
                 6.3800,   0.52897E-02, &
                 6.3900,   0.52415E-02, &
                 6.4000,   0.51936E-02, &
                 6.4100,   0.51458E-02, &
                 6.4200,   0.50982E-02, &
                 6.4300,   0.50507E-02, &
                 6.4400,   0.50035E-02, &
                 6.4500,   0.49564E-02, &
                 6.4600,   0.49096E-02, &
                 6.4700,   0.48629E-02, &
                 6.4800,   0.48164E-02, &
                 6.4900,   0.47701E-02, &
                 6.5000,   0.47240E-02, &
                 6.5100,   0.46781E-02, &
                 6.5200,   0.46324E-02, &
                 6.5300,   0.45869E-02, &
                 6.5400,   0.45416E-02, &
                 6.5500,   0.44965E-02, &
                 6.5600,   0.44516E-02, &
                 6.5700,   0.44070E-02, &
                 6.5800,   0.43625E-02, &
                 6.5900,   0.43183E-02, &
                 6.6000,   0.42742E-02, &
                 6.6100,   0.42304E-02, &
                 6.6200,   0.41868E-02, &
                 6.6300,   0.41434E-02, &
                 6.6400,   0.41003E-02, &
                 6.6500,   0.40573E-02, &
                 6.6600,   0.40146E-02, &
                 6.6700,   0.39721E-02, &
                 6.6800,   0.39298E-02, &
                 6.6900,   0.38878E-02, &
                 6.7000,   0.38459E-02, &
                 6.7100,   0.38044E-02, &
                 6.7200,   0.37630E-02, &
                 6.7300,   0.37219E-02, &
                 6.7400,   0.36809E-02, &
                 6.7500,   0.36403E-02, &
                 6.7600,   0.35998E-02, &
                 6.7700,   0.35596E-02, &
                 6.7800,   0.35196E-02, &
                 6.7900,   0.34799E-02, &
                 6.8000,   0.34404E-02, &
                 6.8100,   0.34011E-02, &
                 6.8200,   0.33621E-02, &
                 6.8300,   0.33233E-02, &
                 6.8400,   0.32847E-02, &
                 6.8500,   0.32464E-02, &
                 6.8600,   0.32083E-02, &
                 6.8700,   0.31705E-02, &
                 6.8800,   0.31329E-02, &
                 6.8900,   0.30955E-02, &
                 6.9000,   0.30584E-02, &
                 6.9100,   0.30215E-02, &
                 6.9200,   0.29849E-02, &
                 6.9300,   0.29484E-02, &
                 6.9400,   0.29123E-02, &
                 6.9500,   0.28764E-02, &
                 6.9600,   0.28407E-02, &
                 6.9700,   0.28052E-02, &
                 6.9800,   0.27700E-02, &
                 6.9900,   0.27351E-02, &
                 7.0000,   0.27004E-02, &
                 7.0100,   0.26659E-02, &
                 7.0200,   0.26317E-02, &
                 7.0300,   0.25977E-02, &
                 7.0400,   0.25639E-02, &
                 7.0500,   0.25304E-02, &
                 7.0600,   0.24971E-02, &
                 7.0700,   0.24641E-02, &
                 7.0800,   0.24313E-02, &
                 7.0900,   0.23988E-02, &
                 7.1000,   0.23665E-02, &
                 7.1100,   0.23344E-02, &
                 7.1200,   0.23026E-02, &
                 7.1300,   0.22710E-02, &
                 7.1400,   0.22397E-02, &
                 7.1500,   0.22086E-02, &
                 7.1600,   0.21777E-02, &
                 7.1700,   0.21471E-02, &
                 7.1800,   0.21167E-02, &
                 7.1900,   0.20865E-02, &
                 7.2000,   0.20566E-02, &
                 7.2100,   0.20270E-02, &
                 7.2200,   0.19975E-02, &
                 7.2300,   0.19683E-02, &
                 7.2400,   0.19393E-02, &
                 7.2500,   0.19106E-02, &
                 7.2600,   0.18821E-02, &
                 7.2700,   0.18538E-02, &
                 7.2800,   0.18258E-02, &
                 7.2900,   0.17980E-02, &
                 7.3000,   0.17705E-02, &
                 7.3100,   0.17431E-02, &
                 7.3200,   0.17160E-02, &
                 7.3300,   0.16891E-02, &
                 7.3400,   0.16625E-02, &
                 7.3500,   0.16361E-02, &
                 7.3600,   0.16099E-02, &
                 7.3700,   0.15839E-02, &
                 7.3800,   0.15582E-02, &
                 7.3900,   0.15327E-02, &
                 7.4000,   0.15074E-02, &
                 7.4100,   0.14824E-02, &
                 7.4200,   0.14575E-02, &
                 7.4300,   0.14329E-02, &
                 7.4400,   0.14086E-02, &
                 7.4500,   0.13844E-02, &
                 7.4600,   0.13604E-02, &
                 7.4700,   0.13367E-02, &
                 7.4800,   0.13132E-02, &
                 7.4900,   0.12899E-02, &
                 7.5000,   0.12669E-02, &
                 7.5100,   0.12440E-02, &
                 7.5200,   0.12214E-02, &
                 7.5300,   0.11990E-02, &
                 7.5400,   0.11768E-02, &
                 7.5500,   0.11548E-02, &
                 7.5600,   0.11330E-02, &
                 7.5700,   0.11114E-02, &
                 7.5800,   0.10901E-02, &
                 7.5900,   0.10689E-02, &
                 7.6000,   0.10480E-02, &
                 7.6100,   0.10272E-02, &
                 7.6200,   0.10067E-02, &
                 7.6300,   0.98637E-03, &
                 7.6400,   0.96626E-03, &
                 7.6500,   0.94634E-03, &
                 7.6600,   0.92663E-03, &
                 7.6700,   0.90712E-03, &
                 7.6800,   0.88781E-03, &
                 7.6900,   0.86870E-03, &
                 7.7000,   0.84979E-03, &
                 7.7100,   0.83108E-03, &
                 7.7200,   0.81256E-03, &
                 7.7300,   0.79424E-03, &
                 7.7400,   0.77611E-03, &
                 7.7500,   0.75818E-03, &
                 7.7600,   0.74043E-03, &
                 7.7700,   0.72288E-03, &
                 7.7800,   0.70553E-03, &
                 7.7900,   0.68835E-03, &
                 7.8000,   0.67137E-03, &
                 7.8100,   0.65458E-03, &
                 7.8200,   0.63797E-03, &
                 7.8300,   0.62154E-03, &
                 7.8400,   0.60530E-03, &
                 7.8500,   0.58925E-03, &
                 7.8600,   0.57337E-03, &
                 7.8700,   0.55768E-03, &
                 7.8800,   0.54216E-03, &
                 7.8900,   0.52682E-03, &
                 7.9000,   0.51166E-03, &
                 7.9100,   0.49668E-03, &
                 7.9200,   0.48187E-03, &
                 7.9300,   0.46724E-03, &
                 7.9400,   0.45278E-03, &
                 7.9500,   0.43849E-03, &
                 7.9600,   0.42437E-03, &
                 7.9700,   0.41042E-03, &
                 7.9800,   0.39664E-03, &
                 7.9900,   0.38303E-03, &
                 8.0000,   0.36958E-03, &
                 8.0100,   0.35630E-03, &
                 8.0200,   0.34319E-03, &
                 8.0300,   0.33023E-03, &
                 8.0400,   0.31744E-03, &
                 8.0500,   0.30481E-03, &
                 8.0600,   0.29234E-03, &
                 8.0700,   0.28003E-03, &
                 8.0800,   0.26787E-03, &
                 8.0900,   0.25588E-03, &
                 8.1000,   0.24403E-03, &
                 8.1100,   0.23234E-03, &
                 8.1200,   0.22081E-03, &
                 8.1300,   0.20942E-03, &
                 8.1400,   0.19819E-03, &
                 8.1500,   0.18711E-03, &
                 8.1600,   0.17617E-03, &
                 8.1700,   0.16539E-03, &
                 8.1800,   0.15475E-03, &
                 8.1900,   0.14425E-03, &
                 8.2000,   0.13390E-03, &
                 8.2100,   0.12369E-03, &
                 8.2200,   0.11363E-03, &
                 8.2300,   0.10370E-03, &
                 8.2400,   0.93920E-04, &
                 8.2500,   0.84274E-04, &
                 8.2600,   0.74766E-04, &
                 8.2700,   0.65395E-04, &
                 8.2800,   0.56159E-04, &
                 8.2900,   0.47057E-04, &
                 8.3000,   0.38089E-04, &
                 8.3100,   0.29254E-04, &
                 8.3200,   0.20549E-04, &
                 8.3300,   0.11975E-04, &
                 8.3400,   0.35296E-05, &
                 8.3500,  -0.47876E-05, &
                 8.3600,  -0.12978E-04, &
                 8.3700,  -0.21042E-04, &
                 8.3800,  -0.28982E-04, &
                 8.3900,  -0.36798E-04, &
                 8.4000,  -0.44491E-04, &
                 8.4100,  -0.52062E-04, &
                 8.4200,  -0.59513E-04, &
                 8.4300,  -0.66844E-04, &
                 8.4400,  -0.74057E-04, &
                 8.4500,  -0.81152E-04, &
                 8.4600,  -0.88131E-04, &
                 8.4700,  -0.94994E-04, &
                 8.4800,  -0.10174E-03, &
                 8.4900,  -0.10838E-03, &
                 8.5000,  -0.11490E-03, &
                 8.5100,  -0.12131E-03, &
                 8.5200,  -0.12762E-03, &
                 8.5300,  -0.13381E-03, &
                 8.5400,  -0.13989E-03, &
                 8.5500,  -0.14587E-03, &
                 8.5600,  -0.15174E-03, &
                 8.5700,  -0.15750E-03, &
                 8.5800,  -0.16316E-03, &
                 8.5900,  -0.16872E-03, &
                 8.6000,  -0.17418E-03, &
                 8.6100,  -0.17953E-03, &
                 8.6200,  -0.18478E-03, &
                 8.6300,  -0.18994E-03, &
                 8.6400,  -0.19499E-03, &
                 8.6500,  -0.19995E-03, &
                 8.6600,  -0.20481E-03, &
                 8.6700,  -0.20957E-03, &
                 8.6800,  -0.21425E-03, &
                 8.6900,  -0.21882E-03, &
                 8.7000,  -0.22331E-03, &
                 8.7100,  -0.22770E-03, &
                 8.7200,  -0.23200E-03, &
                 8.7300,  -0.23621E-03, &
                 8.7400,  -0.24033E-03, &
                 8.7500,  -0.24437E-03, &
                 8.7600,  -0.24832E-03, &
                 8.7700,  -0.25218E-03, &
                 8.7800,  -0.25595E-03, &
                 8.7900,  -0.25964E-03, &
                 8.8000,  -0.26325E-03, &
                 8.8100,  -0.26677E-03, &
                 8.8200,  -0.27022E-03, &
                 8.8300,  -0.27358E-03, &
                 8.8400,  -0.27686E-03, &
                 8.8500,  -0.28006E-03, &
                 8.8600,  -0.28318E-03, &
                 8.8700,  -0.28623E-03, &
                 8.8800,  -0.28920E-03, &
                 8.8900,  -0.29209E-03, &
                 8.9000,  -0.29491E-03, &
                 8.9100,  -0.29765E-03, &
                 8.9200,  -0.30033E-03, &
                 8.9300,  -0.30292E-03, &
                 8.9400,  -0.30545E-03, &
                 8.9500,  -0.30791E-03, &
                 8.9600,  -0.31029E-03, &
                 8.9700,  -0.31261E-03, &
                 8.9800,  -0.31486E-03, &
                 8.9900,  -0.31704E-03, &
                 9.0000,  -0.31915E-03, &
                 9.0100,  -0.32120E-03, &
                 9.0200,  -0.32318E-03, &
                 9.0300,  -0.32510E-03, &
                 9.0400,  -0.32696E-03, &
                 9.0500,  -0.32875E-03, &
                 9.0600,  -0.33048E-03, &
                 9.0700,  -0.33215E-03, &
                 9.0800,  -0.33375E-03, &
                 9.0900,  -0.33530E-03, &
                 9.1000,  -0.33679E-03, &
                 9.1100,  -0.33822E-03, &
                 9.1200,  -0.33959E-03, &
                 9.1300,  -0.34091E-03, &
                 9.1400,  -0.34217E-03, &
                 9.1500,  -0.34337E-03, &
                 9.1600,  -0.34452E-03, &
                 9.1700,  -0.34561E-03, &
                 9.1800,  -0.34665E-03, &
                 9.1900,  -0.34764E-03, &
                 9.2000,  -0.34858E-03, &
                 9.2100,  -0.34946E-03, &
                 9.2200,  -0.35030E-03, &
                 9.2300,  -0.35108E-03, &
                 9.2400,  -0.35182E-03, &
                 9.2500,  -0.35251E-03, &
                 9.2600,  -0.35315E-03, &
                 9.2700,  -0.35374E-03, &
                 9.2800,  -0.35428E-03, &
                 9.2900,  -0.35478E-03, &
                 9.3000,  -0.35524E-03, &
                 9.3100,  -0.35565E-03, &
                 9.3200,  -0.35601E-03, &
                 9.3300,  -0.35633E-03, &
                 9.3400,  -0.35661E-03, &
                 9.3500,  -0.35685E-03, &
                 9.3600,  -0.35704E-03, &
                 9.3700,  -0.35720E-03, &
                 9.3800,  -0.35731E-03, &
                 9.3900,  -0.35739E-03, &
                 9.4000,  -0.35742E-03, &
                 9.4100,  -0.35742E-03, &
                 9.4200,  -0.35737E-03, &
                 9.4300,  -0.35730E-03, &
                 9.4400,  -0.35718E-03, &
                 9.4500,  -0.35703E-03, &
                 9.4600,  -0.35684E-03, &
                 9.4700,  -0.35662E-03, &
                 9.4800,  -0.35636E-03, &
                 9.4900,  -0.35607E-03, &
                 9.5000,  -0.35574E-03, &
                 9.5100,  -0.35539E-03, &
                 9.5200,  -0.35500E-03, &
                 9.5300,  -0.35457E-03, &
                 9.5400,  -0.35412E-03, &
                 9.5500,  -0.35363E-03, &
                 9.5600,  -0.35312E-03, &
                 9.5700,  -0.35258E-03, &
                 9.5800,  -0.35200E-03, &
                 9.5900,  -0.35140E-03, &
                 9.6000,  -0.35077E-03, &
                 9.6100,  -0.35011E-03, &
                 9.6200,  -0.34942E-03, &
                 9.6300,  -0.34871E-03, &
                 9.6400,  -0.34797E-03, &
                 9.6500,  -0.34721E-03, &
                 9.6600,  -0.34642E-03, &
                 9.6700,  -0.34560E-03, &
                 9.6800,  -0.34476E-03, &
                 9.6900,  -0.34390E-03, &
                 9.7000,  -0.34301E-03, &
                 9.7100,  -0.34210E-03, &
                 9.7200,  -0.34117E-03, &
                 9.7300,  -0.34021E-03, &
                 9.7400,  -0.33923E-03, &
                 9.7500,  -0.33824E-03, &
                 9.7600,  -0.33722E-03, &
                 9.7700,  -0.33618E-03, &
                 9.7800,  -0.33512E-03, &
                 9.7900,  -0.33404E-03, &
                 9.8000,  -0.33294E-03, &
                 9.8100,  -0.33182E-03, &
                 9.8200,  -0.33069E-03, &
                 9.8300,  -0.32953E-03, &
                 9.8400,  -0.32836E-03, &
                 9.8500,  -0.32717E-03, &
                 9.8600,  -0.32597E-03, &
                 9.8700,  -0.32474E-03, &
                 9.8800,  -0.32351E-03, &
                 9.8900,  -0.32225E-03, &
                 9.9000,  -0.32098E-03, &
                 9.9100,  -0.31970E-03, &
                 9.9200,  -0.31840E-03, &
                 9.9300,  -0.31709E-03, &
                 9.9400,  -0.31576E-03, &
                 9.9500,  -0.31442E-03, &
                 9.9600,  -0.31307E-03, &
                 9.9700,  -0.31170E-03, &
                 9.9800,  -0.31032E-03, &
                 9.9900,  -0.30893E-03, &
                10.0000,  -0.30752E-03  ]

        do i = 1, n
            rn(i)  = tbl((i-1)*2+1)
            kei(i) = tbl((i-1)*2+2)
        end do 

        return

    end subroutine load_kei_values

    subroutine load_kei_values_file(rn,kei,filename)

        implicit none

        real(wp), allocatable,  intent(OUT):: rn(:)
        real(wp), allocatable,  intent(OUT):: kei(:)
        character(len=*),       intent(IN) :: filename

        ! Local variables 
        integer :: k 
        integer, parameter :: ntot    = 1001
        integer, parameter :: filenum = 177

        if (allocated(rn))  deallocate(rn)
        if (allocated(kei)) deallocate(kei) 

        allocate(rn(ntot))
        allocate(kei(ntot)) 

        ! fonction de kelvin
        ! lecture de la table kei qui est tous les 0.01 entre 0 et 10
        ! STEPK=100=1/ecart 
        !cdc modification du chemin maintenant fonction de dir_inp
        ! trim(dir_inp)//'kelvin.res'
        open(filenum,file=trim(filename))
        read(filenum,*)  ! Skip first line
        do k = 1, ntot
            read(filenum,*) rn(k),kei(k)
        end do
        close(filenum)

        return

    end subroutine load_kei_values_file

!mmr----------------------------------------------------------------
!mmr all below is new 


! mmr ========================================================
!
! mmr Functions related to the ELVA
!
! mmr ========================================================

 
    subroutine calc_asthenosphere_viscous_params(nx,ny,dx,D_lith,mu,kappa,beta) 
       
      integer(kind=4), intent(IN)              :: nx, ny
      real(wp), intent(IN)                     :: dx
      real(wp), intent(IN)                     :: D_lith(:,:)

      real(wp), intent(OUT)                     :: mu      
      real(wp), allocatable, intent(OUT)        :: kappa(:,:), beta(:,:) 
      
      real(wp)                                  :: xd, yd
      integer(kind=4)                           :: i, j, ip, iq, ic, jc 
      
      allocate(kappa(nx,ny))
      allocate(beta(nx,ny))

      mu = 2.*pi/((nx-1)*dx)
      
      kappa = 0.0
      
      ic = (nx-1)/2 + 1  !mmr Alex - I have already calculated these in test_isostasy, I don't like to recalculate it here
      jc = (ny-1)/2 + 1

      do i = 1, nx
            if (i.le.ic) then 
            ip = i-1
         else
            ip = nx-i+1
         endif
         do j = 1, ny
               if (j.le.jc) then
               iq = j-1  
            else
               iq = ny-j+1
            endif
            kappa(i,j)  = (ip*ip + iq*iq)**0.5
            beta(i,j)   = rho_a*g + D_lith(i,j)*(mu**4)*kappa(i,j)**4
         enddo
      enddo

         return
      
    end subroutine calc_asthenosphere_viscous_params

 
    subroutine calc_asthenosphere_viscous(dzbdt,q,w,dt,mu,kappa,beta)
        ! Calculate rate of change of vertical bedrock height
        ! from a viscous half-space asthenosphere.
        ! Contains viscous component only.

      use, intrinsic :: iso_c_binding
      implicit none
      include 'fftw3.f03'

 
      integer(kind=4), parameter :: nd = 2
      
      real(wp), intent(OUT)   :: dzbdt(:,:)
      real(wp), intent(IN)    :: q(:,:)
      real(wp), intent(INOUT) :: w(:,:)
      real(wp), intent(IN)    :: dt
      real(wp), intent(IN)    :: mu
      real(wp), intent(IN)    :: kappa(:,:)
      real(wp), intent(IN)    :: beta(:,:)     

      
      real(wp), allocatable   :: w0(:,:)
      real(wp), allocatable   :: q_hat(:,:)
      real(wp), allocatable   :: w_hat(:,:)

      real(c_double), allocatable :: data_in(:,:), data_inbis(:,:)
      complex(c_double_complex), allocatable :: data_out(:,:)

      complex(dp), allocatable   :: q_hat_c(:,:)
      complex(dp), allocatable   :: w_hat_c(:,:)

      real(wp), allocatable   :: w_hat_c_re(:,:), w_hat_c_im(:,:)      

      real(wp)                :: dt_sec
      integer (kind = 4)      :: l, m, i, j

      logical  :: fft_r2r, fft_c2c


      dt_sec = dt * 3600*24*365 ! [s] 

     
      l = size(q,1)
      m = size(q,2)

      
      allocate(w0(l,m))
      allocate(q_hat(l,m))
      allocate(w_hat(l,m))
      allocate(q_hat_c(l,m/2+1)) 
      allocate(w_hat_c(l,m/2+1))

      allocate(w_hat_c_re(l,m/2+1))
      allocate(w_hat_c_im(l,m/2+1))


      
      
      !  Initialize 
      
      w0 = w
      
      q_hat    = 0.0
      w_hat    = 0.0
      w_hat_c  = 0.0
      dzbdt    = 0.0
      
      fft_r2r = .true.
      fft_c2c = .false.
      
      if (fft_r2r) then

         ! fft of load
         ! 

         call calc_fft_forward_r2r(q,q_hat)

         ! fft of displacement

         call calc_fft_forward_r2r(w,w_hat)

         ! displacement fft at timestep n+1    

         if (dt.gt.0)  w_hat = ( ( 2.*eta*mu*kappa - (dt_sec/2.)*beta)*w_hat + dt_sec*q_hat)/(2*eta*mu*kappa + (dt_sec/2.)*beta)

         ! Inverse fft to obtain displacement at n+1

         call calc_fft_backward_r2r(w_hat,w)

      else if (fft_c2c) then

         ! fft of load

         call calc_fft_forward_c2c(q,q_hat)

         ! fft of displacement

         call calc_fft_forward_c2c(w,w_hat)

         ! displacement fft at timestep n+1    

         if (dt.gt.0)  w_hat = ( ( 2.*eta*mu*kappa - (dt_sec/2.)*beta)*w_hat + dt_sec*q_hat)/(2*eta*mu*kappa + (dt_sec/2.)*beta)

         ! Inverse fft to obtain displacement at n+1

         call calc_fft_backward_c2c(w_hat,w)

      else  

         write(*,*) 'Error, you have to choose an option to calculate the FFTs'
         stop
      endif
         
    ! Impose boundary conditions 

      w  = w - 0.25*(w(1,1)+w(l,m)+w(1,m)+w(l,1)) 

      ! Rate of viscous asthenosphere uplift

      if (dt.gt.0.)  dzbdt = -(w-w0)/dt


      deallocate(w0)
      deallocate(q_hat)
      deallocate(w_hat)
      deallocate(q_hat_c)
      deallocate(w_hat_c)

      return

    end subroutine calc_asthenosphere_viscous

    subroutine make_fft_plans(in,plan_fwd,plan_bck) 
 
      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03' 
      
      real(wp), intent(IN)       :: in(:,:)
      type(c_ptr), intent(OUT)   :: plan_fwd, plan_bck
      complex (dp), allocatable  :: in_aux(:,:)
      complex (dp), allocatable  :: out_aux(:,:)
      integer (kind = 4)         :: l, m

      
      l    = size(in,1)
      m    = size(in,2)

      if(l.ne.m) then
         print*,'Dimensions do not match, stopping now'
         stop
      endif

      in_aux = in

! r2r      
      plan_fwd = fftw_plan_dft_2d(l,m,in_aux,out_aux,-1,FFTW_ESTIMATE)
      plan_bck = fftw_plan_dft_2d(l,m,out_aux,in_aux,+1,FFTW_ESTIMATE)

      return

    end subroutine make_fft_plans


   subroutine calc_fft_forward_r2r(in,out)

      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03'  

      real(wp), intent(IN)       :: in(:,:)
      real(wp), intent(OUT)      :: out(:,:)
      
      real(wp), allocatable      :: rec(:,:)
      real(dp), allocatable      :: in_aux(:,:)
      real(dp), allocatable      :: out_aux(:,:) 
      real(dp), allocatable      :: rec_aux(:,:)
      real(dp)                   :: dx, cc
      type(c_ptr)                :: plan
      integer(kind=4)            :: m,n
      logical                    :: print_check
      

      m = size(in,1)
      n = size(in,2)

      print_check = .false.

      allocate(in_aux(m,n))
      allocate(out_aux(m,n))
      allocate(rec(m,n))
      allocate(rec_aux(m,n))


! http://www.fftw.org/fftw3_doc/The-Discrete-Hartley-Transform.html
      
!      The discrete Hartley transform (DHT) is an invertible linear
!      transform closely related to the DFT. In the DFT, one
!      multiplies each input by cos - i * sin (a complex exponential),
!      whereas in the DHT each input is multiplied by simply cos +
!      sin. Thus, the DHT transforms n real numbers to n real numbers,
!      and has the convenient property of being its own inverse. In
!      FFTW, a DHT (of any positive n) can be specified by an r2r kind
!      of FFTW_DHT.

!      Like the DFT, in FFTW the DHT is unnormalized, so computing a
!      DHT of size n followed by another DHT of the same size will
!      result in the original array multiplied by n.

!      The DHT was originally proposed as a more efficient alternative
!      to the DFT for real data, but it was subsequently shown that a
!      specialized DFT (such as FFTW’s r2hc or r2c transforms) could
!      be just as fast. In FFTW, the DHT is actually computed by
!      post-processing an r2hc transform, so there is ordinarily no
!      reason to prefer it from a performance perspective.5 However,
!      we have heard rumors that the DHT might be the most appropriate
!      transform in its own right for certain applications, and we
!      would be very interested to hear from anyone who finds it
!      useful.


      in_aux = in 
      plan = fftw_plan_r2r_2d(m,n,in_aux,out_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)

      call fftw_execute_r2r(plan, in_aux, out_aux)
!mmr      call fftw_destroy_plan(plan)
      out = real(out_aux/sqrt(m*n*1.))

      if (print_check) then
         call r4mat_print_some ( m, n, in, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the original data:' )
         plan =  fftw_plan_r2r_2d(m,n,out_aux,rec_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)
         call fftw_execute_r2r(plan, out_aux, rec_aux)
         rec_aux = rec_aux/(m*n)
         rec = real(rec_aux,wp) 
         call fftw_destroy_plan(plan)
         call r4mat_print_some (m, n, rec, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the recovered data:' )
      endif

      deallocate(in_aux)
      deallocate(out_aux)
      deallocate(rec)
      deallocate(rec_aux)
      
      return
      
    end subroutine calc_fft_forward_r2r

    subroutine calc_fft_backward_r2r(in,out)

      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03'  

      real(wp), intent(IN)       :: in(:,:)
      real(wp), intent(OUT)      :: out(:,:)

      real(wp), allocatable      :: rec(:,:)
      real(dp), allocatable      :: in_aux(:,:)
      real(dp), allocatable      :: out_aux(:,:) 
      real(dp), allocatable      :: rec_aux(:,:)
      real(dp)                   :: dx, cc
      type(c_ptr)                :: plan
      integer(kind=4)            :: m,n
      logical                    :: print_check

      m = size(in,1)
      n = size(in,2)


      allocate(in_aux(m,n))
      allocate(out_aux(m,n))
      allocate(rec(m,n))
      allocate(rec_aux(m,n))


      in_aux = in
      plan = fftw_plan_r2r_2d(m,n,in_aux,out_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)

      call fftw_execute_r2r(plan, in_aux, out_aux)
!mmr      call fftw_destroy_plan(plan)
      out = real(out_aux/sqrt(m*n*1.)) 

      if (print_check) then
         call r4mat_print_some ( m, n, in, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the original data:' )
         plan =  fftw_plan_r2r_2d(m,n,out_aux,rec_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)
         call fftw_execute_r2r(plan, out_aux, rec_aux)
         rec_aux = rec_aux/(m*n)
         rec = real(rec_aux,wp) 
         call fftw_destroy_plan(plan)
         call r4mat_print_some ( m, n, rec, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the recovered data:' )
      endif

      deallocate(in_aux)
      deallocate(out_aux)
      deallocate(rec)
      deallocate(rec_aux)

      return
      
    end subroutine calc_fft_backward_r2r

    
    ! http://www.fftw.org/fftw3_doc/Real_002ddata-DFTs.html
    ! FFTW computes an unnormalized transform: computing an r2c
    ! followed by a c2r transform (or vice versa) will result in the
    ! original data multiplied by the size of the transform (the
    ! product of the logical dimensions). An r2c transform produces
    ! the same output as a FFTW_FORWARD complex DFT of the same input,
    ! and a c2r transform is correspondingly equivalent to
    ! FFTW_BACKWARD.
    
 
    subroutine calc_fft_forward_c2c(in,out)

      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03'  

      real(wp), intent(IN)       :: in(:,:)
      real(wp), intent(OUT)      :: out(:,:)

      real(wp), allocatable      :: rec(:,:)
      complex(dp), allocatable   :: in_aux(:,:)   
      complex(dp), allocatable   :: out_aux(:,:)     
      real(dp), allocatable      :: rec_aux(:,:)
      real(dp)                   :: dx, cc
      type(c_ptr)                :: plan
      integer(kind=4)            :: m, n

      m = size(in,1)
      n = size(in,2)

      
      allocate(in_aux(m,n))
      allocate(out_aux(m,n))
      allocate(rec(m,n))
      allocate(rec_aux(m,n))


      in_aux = in
      plan = fftw_plan_dft_2d(m,n,in_aux,out_aux,-1,FFTW_ESTIMATE)

      call fftw_execute_dft(plan, in_aux, out_aux)                 
      call fftw_destroy_plan(plan)
      out = real(out_aux/sqrt(m*n*1.)) 


      deallocate(in_aux)
      deallocate(out_aux)
      deallocate(rec)
      deallocate(rec_aux)

      return
      
    end subroutine calc_fft_forward_c2c

    subroutine calc_fft_backward_c2c(in,out)

      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03'  

      real(wp), intent(IN)       :: in(:,:)
      real(wp), intent(OUT)      :: out(:,:)

      real(wp), allocatable      :: rec(:,:)
      complex(dp), allocatable   :: in_aux(:,:)   
      complex(dp), allocatable   :: out_aux(:,:)     
      real(dp), allocatable      :: rec_aux(:,:)
      real(dp)                   :: dx, cc
      type(c_ptr)                :: plan
      integer(kind=4)            :: m, n

      m = size(in,1)
      n = size(in,2)


      allocate(in_aux(m,n))
      allocate(out_aux(m,n))
      allocate(rec(m,n))
      allocate(rec_aux(m,n))


      in_aux = in
      plan = fftw_plan_dft_2d(m,n,in_aux,out_aux,1,FFTW_ESTIMATE)

      call fftw_execute_dft(plan, in_aux, out_aux)                 
      call fftw_destroy_plan(plan)
      out = real(out_aux/sqrt(m*n*1.)) 


      deallocate(in_aux)
      deallocate(out_aux)
      deallocate(rec)
      deallocate(rec_aux)

      return

    end subroutine calc_fft_backward_c2c


!=========================================================
!
! Routines for analytical ELVA disk
!
!=========================================================

    
    subroutine calc_analytical_asthenosphere_viscous_disk_params(nx,ny,dx,kappa_min,kappa_max,dk,kappa_mod,dist2c,r,lr)  
       
      integer(kind=4), intent(IN)              :: nx, ny
      real(wp), intent(IN)                     :: dx
      real(wp), intent(IN)                     :: kappa_min, kappa_max, dk

      real(wp), allocatable, intent(OUT)        :: kappa_mod(:), dist2c(:,:), r(:) 
      integer(kind=4), allocatable, intent(OUT) :: lr(:,:)
      
      integer(kind=4), allocatable              :: n(:,:)
      real(wp)                                  :: xd, yd
      integer(kind=4)                           :: i, j, ip, iq, ic, jc, k, nk, l, nl


      nk = int((kappa_max-kappa_min)/dk)
      
      nl = int(nx*sqrt(2.)/2) + 2
      
      
      allocate(kappa_mod(nk+1))
      allocate(dist2c(nx,ny))
      allocate(r(nl))
      allocate(lr(nx,ny))
      allocate(n(nx,ny))
      
      do k = 1, nk+1
         kappa_mod(k) = kappa_min + dk * (k-1)
      enddo

      
      ic = (nx-1)/2 + 1  !mmr Alex - I have already calculated these in test_isostasy, I don't like to recalculate it here
      jc = (ny-1)/2 + 1

      do i = 1, nx
         do j = 1, ny
            xd = dx*(i-ic)
            yd = dx*(j-jc)
            dist2c(i,j) = sqrt(xd**2 + yd**2)  !mmr Alex - I have already calculated it in test_isostasy, I don't like to recalculate it here
         enddo
      enddo

               
      ! remap to (nx * xy) grid (w = w(r,t) so all points with same r have same w; this reduces computational time)
      

      do l = 1, nl 
         r(l) =  dx * (l-1)
      enddo
     
         do i = 1, nx 
            do j = 1, ny
               n(i,j) = 0       
               do l = 1, nl-1
                  if (dist2c(i,j).ge.r(l) .and. dist2c(i,j) .lt. r(l+1) ) then
                     lr(i,j) = l 
                     n(i,j) = n(i,j) + 1
                  endif
               enddo
               if (n(i,j).ne.1) then
                  print*,'==> error in radial distance allocation', n(i,j), i, j
                  stop
               endif
            enddo
         enddo

         deallocate(n)
         
         return
      
       end subroutine calc_analytical_asthenosphere_viscous_disk_params
       
    subroutine calc_analytical_asthenosphere_viscous_disk(me,dx,t,kappa_mod,dist2c,r,lr,w) 

      ! Calculate analytical solution for displacement for the asthenosphere viscous disk 
      ! u(r,t) as in Bueler et al (2007), eq. 17
      ! remap into original (nx * ny) grid
      
      implicit none

      class(isos_analytical_elva_disk_load_class), intent(INOUT) :: me

      real(wp),parameter :: tol       = 1.e-3  ! error tolerance
      
      real(wp), intent(IN)          :: dx, t
      real(wp), intent(IN)          :: kappa_mod(:), dist2c(:,:), r(:)
      integer(kind=4), intent(IN)   :: lr(:,:)
      real(wp), intent(OUT)         :: w(:,:)

      real(wp), allocatable         :: wr(:)
      integer(kind=4), allocatable  :: n(:,:)
      
      real(wp)                      :: ans
      real(wp)                      :: err
      integer                       :: method    ! quadrature method to use for x
      integer                       :: ierr      ! error code
      
      integer(kind=4)               :: i, j, k, l, nx, ny, nk, nl
      
      nx = size(w,1)
      ny = size(w,2)
      nk = size(kappa_mod)-1
      nl = size(r)

      allocate(wr(nl)) 
      allocate(n(nx,ny))

      method       = me%method  ! 6-point gaussian quadrature     
      me%time      = t          ! [hours] 
      me%n_evals   = 0
      
      w =  0.0
      wr = 0.0
      
      do l = 1, nl ! number of points neccessary to cover whole domain radially (ca. diagonal)
         
         me%r = r(l)

         do k = 1, nk

      ! sets parameters but most notably integration limits      

            call initialize_integration_class(me, fx=analytical_integrand, xl=kappa_mod(k), xu=kappa_mod(k+1), tolx=tol, methodx=method)

            !reset number of function evaluations:

            me%n_evals = 0

            !integrate the function over specified interval

            call integrate_1d (me, ans, ierr, err) 

            wr(l)  = wr(l) + ans

         enddo ! nk 
      enddo ! nl

      do i = 1, nx
         do j = 1, ny
            w(i,j) = wr(lr(i,j))  
         enddo
      enddo
      
      w  = w - 0.25*(w(1,1)+w(nx,ny)+w(1,ny)+w(nx,1)) 

      deallocate(n)

         return
         
  end subroutine calc_analytical_asthenosphere_viscous_disk

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
    case(6);  me%g => g6
       !    case(8);  me%g => g8  
       !    case(10); me%g => g10
       !    case(12); me%g => g12
       !    case(14); me%g => g14
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

    subroutine initialize_analytical_integrand(me,r0,h0,D_lith,eta) 

    implicit none

    class(isos_analytical_elva_disk_load_class),intent(out)  :: me
    real(wp), intent(in)   ::  h0, r0, D_lith, eta 

    me%r0         = r0
    me%h0         = h0
    me%D_lith     = D_lith
    me%eta        = eta
    
    return

    end subroutine initialize_analytical_integrand

  function analytical_integrand(me,x) result(f)

    implicit none

    class(isos_analytical_elva_disk_load_class),intent(inout)  :: me
    real(wp), intent(in)  :: x ! (kappa)
    real(wp)              :: f
    real(wp) beta, h0, r0, eta, D_lith, t_sec, r_eq_zero

    r0     = me%r0
    h0     = me%h0
    D_lith = me%D_lith
    eta    = me%eta
    t_sec  = me%time*365*24*3600
    
    beta = rho_a*g + D_lith*(x**4)

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
    est = me%g(aa(l)+2._wp*hh(l),2._wp*hh(l))
    k = 8
    area = abs(est)
    ef = 0.5 !one_half
    mxl = 0

    !compute refined estimates, estimate the error, etc.
    main : do

        gl = me%g(aa(l)+hh(l),hh(l))
        gr(l) = me%g(aa(l)+3._wp*hh(l),hh(l))
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

!mmr----------------------------------------------------------------
    
end module isostasy




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

