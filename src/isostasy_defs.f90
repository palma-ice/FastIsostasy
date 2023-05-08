module isostasy_defs

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = dp ! sp dp mmr2 

    real(wp), parameter :: pi = 3.14159265359

    
    type isos_param_class 
        integer            :: method            ! Type of isostasy to use
        real(wp)           :: dt_lith           ! [yr] Timestep to recalculate equilibrium lithospheric displacement
        real(wp)           :: dt_step           ! [yr] Timestep to recalculate bedrock uplift and rate
        real(wp)           :: He_lith           ! [km] Effective elastic thickness of lithosphere
        real(wp)           :: visc              ! [Pa s] Asthenosphere viscosity (constant)
        real(wp)           :: tau               ! [yr] Asthenospheric relaxation constant
        
        ! Physical constants
        real(wp) :: E 
        real(wp) :: nu 
        real(wp) :: rho_ice 
        real(wp) :: rho_sw 
        real(wp) :: rho_w
        real(wp) :: rho_a 
        real(wp) :: g 
        real(wp) :: r_earth 
        real(wp) :: m_earth 

        ! Internal parameters 
        integer  :: nx
        integer  :: ny
        real(wp) :: dx                        ! [m] Horizontal resolution                           
        real(wp) :: L_w                       ! [m] Lithosphere flexural length scale (for method=2)
        integer  :: nr                        ! [-] Radius of neighborhood for convolution, in number of grid points        
        real(wp) :: mu                        ! [1/m] 2pi/L                                         
        
        real(wp) :: time_lith                 ! [yr] Current model time of last update of equilibrium lithospheric displacement
        real(wp) :: time_step                 ! [yr] Current model time of last update of bedrock uplift
        
    end type 

    type isos_state_class 
        
        real(wp), allocatable :: He_lith(:,:)       ! [m]  Effective elastic thickness of the lithosphere
        real(wp), allocatable :: D_lith(:,:)        ! [N-m] Lithosphere flexural rigidity
        real(wp), allocatable :: eta_eff(:,:)       ! [Pa-s] Effective asthenosphere viscosity
        real(wp), allocatable :: tau(:,:)           ! [yr] Asthenospheric relaxation timescale field
        
        real(wp), allocatable :: kei(:,:)           ! Kelvin function filter values 
        real(wp), allocatable :: G0(:,:)            ! Green's function values
       
        real(wp), allocatable :: z_bed(:,:)         ! Bedrock elevation         [m]
        real(wp), allocatable :: dzbdt(:,:)         ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: z_bed_ref(:,:)     ! Reference (unweighted) bedrock 
        real(wp), allocatable :: q0(:,:)            ! Reference load
        real(wp), allocatable :: w0(:,:)            ! Reference equilibrium displacement
        real(wp), allocatable :: q1(:,:)            ! Current load          
        real(wp), allocatable :: w1(:,:)            ! Current equilibrium displacement
        real(wp), allocatable :: w2(:,:)            ! Current viscous equilibrium displacement (ELVA)
        
    end type isos_state_class

    type isos_class
        type(isos_param_class) :: par
        type(isos_state_class) :: now      
    end type

    public :: isos_param_class
    public :: isos_state_class
    public :: isos_class 

contains


end module isostasy_defs
