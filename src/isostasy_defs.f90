module isostasy_defs

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = dp

    real(wp), parameter :: pi = 3.14159265359

    
    type isos_param_class
        logical            :: interactive_sealevel
        logical            :: correct_distortion
        integer            :: method                ! Type of isostasy to use
        real(wp)           :: dt_prognostics        ! [yr] Timestep to recalculate equilibrium lithospheric displacement
        real(wp)           :: dt_diagnostics        ! [yr] Timestep to recalculate bedrock uplift and rate
        
        real(wp)           :: He_lith               ! [km] Effective elastic thickness of lithosphere
        real(wp)           :: visc                  ! [Pa s] Asthenosphere viscosity (constant)
        real(wp)           :: tau                   ! [yr] Asthenospheric relaxation constant

        ! Physical constants
        real(wp) :: rho_water
        real(wp) :: rho_ice 
        real(wp) :: rho_seawater 
        real(wp) :: rho_uppermantle
        real(wp) :: rho_litho 
        
        real(wp) :: E 
        real(wp) :: nu
        real(wp) :: compressibility_correction

        real(wp) :: g 
        real(wp) :: r_earth 
        real(wp) :: m_earth
        logical            :: static_load       ! [-] Load static / transient
        character(len=56)  :: visc_method       ! [-] Method use to prescribe asthenosphere's viscosity field
        real(wp)           :: visc_c            ! [Pa s]  Viscosity in channel between elastic lithosphere and viscous asthenosphere (for LV-ELVA only)
        real(wp)           :: thck_c            ! [km]    Thickness of channel between elastic lithosphere and viscous asthenosphere (for LV-ELVA only)
        integer            :: n_lev             ! [-]     Number of layers within viscous asthenosphere (for LV-ELVA only)
        character(len=56)  :: rigidity_method   ! [-] Method use to prescribe lithosphere's rigidity field
        
        real(wp) :: L_w                         ! [m] Lithosphere flexural length scale (for method=2)
        integer  :: nr                          ! [-] Radius of neighborhood for convolution, in number of grid points        
        real(wp) :: mu                          ! [1/m] 2pi/L
        
        real(wp) :: time_diagnostics            ! [yr] Current model time of last update of equilibrium lithospheric displacement
        real(wp) :: time_prognostics            ! [yr] Current model time of last update of bedrock uplift
        
    end type 

    type isos_state_class 
        integer               :: count_updates      ! [1] Number of sea-level updates since beginnning of simulation
        real(wp)              :: t                  ! [yr] Time
        real(wp)              :: bsl                ! [m] Barystatic sea level
        
        real(wp), allocatable :: He_lith(:,:)       ! [m]  Effective elastic thickness of the lithosphere
        real(wp), allocatable :: D_lith(:,:)        ! [N-m] Lithosphere flexural rigidity
        real(wp), allocatable :: eta(:,:,:)           ! [Pa-s]3D mantle viscosity
        real(wp), allocatable :: eta_eff(:,:)       ! [Pa-s] Effective asthenosphere viscosity
        real(wp), allocatable :: tau(:,:)           ! [yr] Asthenospheric relaxation timescale field
       
        real(wp), allocatable :: z_bed(:,:)         ! Bedrock elevation         [m]
        real(wp), allocatable :: dzbdt(:,:)         ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: q(:,:)             ! [Pa] Load
        real(wp), allocatable :: w(:,:)             ! Current viscous displacement
        real(wp), allocatable :: w_equilibrium(:,:) ! Current equilibrium displacement (XLRA)

        real(wp), allocatable :: ssh_perturb(:,:)  ! [m] sea-surface height perturbation
        real(wp), allocatable :: Haf(:,:)           ! [m] Ice thickness above floatation
        real(wp), allocatable :: Hice(:,:)          ! [m] Thickness of ice column
        real(wp), allocatable :: Hsw(:,:)           ! [m] Thickness of seawater column
        real(wp), allocatable :: u(:,:)             ! [m] Viscous displacement
        real(wp), allocatable :: ue(:,:)            ! [m] Elastic displacement

        real(wp), allocatable :: seasurfaceheight(:,:)  ! [m] sea-surface height
        real(wp), allocatable :: canom_load(:,:)    ! [kg m^-2] Load column anomaly
        real(wp), allocatable :: canom_full(:,:)    ! [kg m^-2] Full column anomaly
        real(wp), allocatable :: mass_anom(:,:)     ! [kg] Mass anomaly
        real(wp), allocatable :: maskocean(:,:)     ! [1] Ocean mask
        real(wp), allocatable :: maskactive(:,:)    ! [1] Active mask
        real(wp), allocatable :: maskgrounded(:,:)  ! [1] Grounded mask
    end type isos_state_class

    type isos_domain_class
        integer                 :: i1
        integer                 :: i2
        integer                 :: j1
        integer                 :: j2
        integer                 :: convo_offset
        integer                 :: nx
        integer                 :: ny
        real(wp)                :: dx
        real(wp)                :: dy
        real(wp)                :: mu

        real(wp), allocatable   :: dx_matrix(:,:)   ! [m] Length difference in x between neighboring cells
        real(wp), allocatable   :: dy_matrix(:,:)   ! [m] Length difference in y between neighboring cells
        real(wp), allocatable   :: A(:,:)           ! [m^2] Cell area
        real(wp), allocatable   :: K(:,:)           ! [1] Distortion matrix
        real(wp), allocatable   :: kappa(:,:)       ! Pseudodifferential operator

        ! Auxiliary arrays for fft computation
        real(wp), allocatable       :: real_in_aux(:,:)
        real(wp), allocatable       :: real_out_aux(:,:)
        complex(wp), allocatable    :: cplx_in_aux(:,:)
        complex(wp), allocatable    :: cplx_out_aux(:,:)

        type(c_ptr)             :: forward_fftplan_r2r
        type(c_ptr)             :: backward_fftplan_r2r
        type(c_ptr)             :: forward_dftplan_r2c
        type(c_ptr)             :: backward_dftplan_c2r

        real(wp), allocatable :: kei(:,:)           ! Kelvin function filter values
        real(wp), allocatable :: G0(:,:)            ! Green's function values
        real(wp), allocatable :: GF(:,:)            ! Green's function values (Farrell 1972)
        real(wp), allocatable :: GN(:,:)            ! Green's function for geoid values
    end type isos_domain_class

    type isos_class
        type(isos_param_class)  :: par
        type(isos_domain_class) :: domain
        type(isos_state_class)  :: now
        type(isos_state_class)  :: ref
    end type

    public :: isos_param_class
    public :: isos_domain_class
    public :: isos_state_class
    public :: isos_class

    contains ! nothing but definition of derived types

end module isostasy_defs
