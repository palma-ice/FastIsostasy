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
        integer            :: method                ! Computation method for viscous displacement
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
        real(wp) :: Vden_factor
        
        real(wp) :: E 
        real(wp) :: nu
        real(wp) :: compressibility_correction

        real(wp) :: sec_per_year
        real(wp) :: g
        real(wp) :: r_earth
        real(wp) :: m_earth
        real(wp) :: A_ocean_pd

        character(len=56)  :: visc_method       ! [-] Method use to prescribe asthenosphere's viscosity field
        real(wp)           :: visc_c            ! [Pa s] Channel viscosity (for LV-ELVA only)
        real(wp)           :: thck_c            ! [km] Channel thickness (for LV-ELVA only)
        integer            :: n_lev             ! [-] Number of layers for LV-ELVA
        character(len=56)  :: rigidity_method   ! [-] Method use to prescribe lithosphere's rigidity field
        
        real(wp) :: L_w                         ! [m] Lithosphere flexural length scale (for method=2)
        real(wp) :: time_diagnostics            ! [yr] Current model time of last diagnostic update
        real(wp) :: time_prognostics            ! [yr] Current model time of last prognostic update

    end type

    type isos_domain_class
        integer                 :: i1
        integer                 :: i2
        integer                 :: j1
        integer                 :: j2
        integer                 :: offset
        integer                 :: nx
        integer                 :: ny
        integer                 :: nsq
        real(wp)                :: dx
        real(wp)                :: dy
        real(wp)                :: mu

        real(wp), allocatable   :: dx_matrix(:, :)   ! [m] K * dx
        real(wp), allocatable   :: dy_matrix(:, :)   ! [m] K * dy
        real(wp), allocatable   :: A(:, :)           ! [m^2] Cell area
        real(wp), allocatable   :: K(:, :)           ! [1] Distortion matrix
        real(wp), allocatable   :: kappa(:, :)       ! [1] Pseudodifferential operator
        logical,  allocatable   :: maskactive(:, :)  ! [1] Active mask

        real(wp), allocatable   :: He_lith(:, :)       ! [m] Elastic thickness of the lithosphere
        real(wp), allocatable   :: D_lith(:, :)        ! [N-m] Lithosphere flexural rigidity
        real(wp), allocatable   :: eta(:, :,:)         ! [Pa-s] 3D mantle viscosity
        real(wp), allocatable   :: eta_eff(:, :)       ! [Pa-s] Effective mantle viscosity
        real(wp), allocatable   :: tau(:, :)           ! [yr] Asthenospheric relaxation timescale field

        type(c_ptr)             :: forward_fftplan_r2r
        type(c_ptr)             :: backward_fftplan_r2r
        type(c_ptr)             :: forward_dftplan_r2c
        type(c_ptr)             :: backward_dftplan_c2r

        real(wp), allocatable   :: kei(:, :)   ! Kelvin function filter values
        real(wp), allocatable   :: GV(:, :)    ! Green's function values
        real(wp), allocatable   :: GE(:, :)    ! Green's function for elastic displacement (Farrell 1972)
        real(wp), allocatable   :: GN(:, :)    ! Green's function for ssh_perturb

        complex(wp), allocatable :: FGV(:, :)    ! FFT of GV
        complex(wp), allocatable :: FGE(:, :)    ! FFT of GE
        complex(wp), allocatable :: FGN(:, :)    ! FFT of GN

    end type isos_domain_class

    type isos_state_class 
        real(wp)              :: t                  ! [yr] Time
        real(wp)              :: bsl                ! [m] Barystatic sea level
       
        real(wp), allocatable       :: z_bed(:, :)         ! Bedrock elevation         [m]
        real(wp), allocatable       :: dzbdt(:, :)         ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable       :: q(:, :)             ! [Pa] Load
        real(wp), allocatable       :: w(:, :)             ! Current viscous displacement
        real(wp), allocatable       :: w_equilibrium(:, :) ! Current viscous equilibrium displacement (XLRA)
        real(wp), allocatable       :: we(:, :)            ! [m] Elastic displacement
        complex(wp), allocatable    :: cplx_out_aux(:, :)

        real(wp), allocatable :: Haf(:, :)           ! [m] Ice thickness above floatation
        real(wp), allocatable :: Hice(:, :)          ! [m] Thickness of ice column
        real(wp), allocatable :: Hseawater(:, :)     ! [m] Thickness of seawater column
        ! real(wp), allocatable :: Hsediment(:, :)     ! [m] Thickness of sediment column

        real(wp), allocatable :: ssh(:, :)           ! [m] sea-surface height
        real(wp), allocatable :: ssh_perturb(:, :)   ! [m] sea-surface height perturbation
        real(wp), allocatable :: canom_load(:, :)    ! [kg m^-2] Load column anomaly
        real(wp), allocatable :: canom_full(:, :)    ! [kg m^-2] Full column anomaly
        real(wp), allocatable :: mass_anom(:, :)     ! [kg] Mass anomaly

        logical, allocatable :: maskocean(:, :)      ! [1] Ocean mask
        logical, allocatable :: maskgrounded(:, :)   ! [1] Grounded mask
        logical, allocatable :: maskcontinent(:, :)  ! [1] Continent mask
    end type isos_state_class

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
