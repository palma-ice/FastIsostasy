module isostasy_defs

    use, intrinsic :: iso_c_binding
    use precision
    implicit none
    include 'fftw3.f03'

    real(wp), parameter :: pi = 3.14159265359

    type isos_param_class
        logical     :: heterogeneous_ssh    ! Sea level varies spatially?
        logical     :: interactive_sealevel ! Sea level interacts with solid Earth deformation?
        logical     :: include_elastic      ! Include elastic deformation?
        logical     :: correct_distortion   ! Account for distortion of projection?
        integer     :: method               ! Computation method for viscous displacement
        real(wp)    :: dt_prognostics       ! [yr] Timestep to recalculate prognostics
        real(wp)    :: dt_diagnostics       ! [yr] Timestep to recalculate diagnostics
        real(wp)    :: min_pad              ! [km] Padding around domain

        character(len=56)       :: mantle       ! [-] Method to prescribe viscosity field
        character(len=56)       :: lithosphere  ! [-] Method to prescribe lithospheric thickness field
        character(len=64)       :: layering
        character(len=64)       :: lumping
        character(len=56)       :: viscosity_scaling_method
        real(wp)                :: viscosity_scaling
        real(wp)                :: rheo_smoothing_radius ! [km] Smoothing radius for rheology
        
        real(wp), allocatable   :: zl(:)            ! [km] Layer boundaries
        real(wp), allocatable   :: viscosities(:)   ! [Pa s] Layer viscosities
        integer                 :: nl               ! [-] Number of layers for LV-ELVA
        real(wp)                :: dl
        real(wp)                :: eta_ref          ! [Pa s] Reference viscosity for scaling
        real(wp)                :: tau              ! [yr] Asthenospheric relaxation time

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

        character(len=256)       :: restart
        character(len=256)       :: mask_file
        character(len=256)       :: rheology_file
        character(len=64)        :: litho_thickness_varname ! Variable name in netcdf file
        character(len=64)        :: log10_viscosity_varname ! Variable name in netcdf file
        character(len=64)        :: stddev_log10_viscosity_varname ! Variable name in netcdf file
        logical                  :: use_restart

        real(wp) :: L_w                         ! [m] Lithosphere flexural length scale (for method=2)
        real(wp) :: time_diagnostics            ! [yr] Current model time of last diagnostic update
        real(wp) :: time_prognostics            ! [yr] Current model time of last prognostic update

        ! Internal parameter
        logical :: ref_was_set
        
    end type

    type isos_domain_class
        integer                 :: nx_ice, ny_ice, n_pad_x, n_pad_y, n_pad_xy
        integer                 :: i1, i2, j1, j2
        integer                 :: icrop1, icrop2, jcrop1, jcrop2
        integer                 :: offset
        integer                 :: nx, ny
        real(wp)                :: dx, dy

        real(wp), allocatable :: xc_ice(:)        ! [m] X coordinates of ice domain
        real(wp), allocatable :: yc_ice(:)        ! [m] Y coordinates of ice domain
        real(wp), allocatable :: x_ice(:, :)    ! [m] X coordinates of ice domain
        real(wp), allocatable :: y_ice(:, :)    ! [m] Y coordinates of ice domain
        real(wp), allocatable   :: xc(:)            ! [m]
        real(wp), allocatable   :: yc(:)            ! [m]
        real(wp), allocatable   :: x(:, :)          ! [m]
        real(wp), allocatable   :: y(:, :)          ! [m]

        real(wp), allocatable   :: dx_matrix(:, :)  ! [m] K * dx
        real(wp), allocatable   :: dy_matrix(:, :)  ! [m] K * dy
        real(wp), allocatable   :: A(:, :)          ! [m^2] Cell area
        real(wp), allocatable   :: K(:, :)          ! [1] Distortion matrix
        real(wp), allocatable   :: kappa(:, :)      ! [1] Pseudodifferential operator
        real(wp), allocatable   :: R(:, :)          ! [1] Scaling operator for viscosity variations
        logical,  allocatable   :: maskactive(:, :) ! [1] Active mask

        real(wp), allocatable   :: boundaries(:, :, :)
        real(wp), allocatable   :: He_lith(:, :)    ! [m] Elastic thickness of the lithosphere
        real(wp), allocatable   :: D_lith(:, :)     ! [N-m] Lithosphere flexural rigidity
        real(wp), allocatable   :: eta(:, :, :)     ! [Pa-s] 3D mantle viscosity
        real(wp), allocatable   :: eta_eff(:, :)    ! [Pa-s] Effective mantle viscosity
        real(wp), allocatable   :: tau(:, :)        ! [yr] Asthenospheric relaxation timescale field

        type(c_ptr)             :: forward_fftplan_r2r
        type(c_ptr)             :: backward_fftplan_r2r
        type(c_ptr)             :: forward_dftplan_r2c
        type(c_ptr)             :: backward_dftplan_c2r

        real(wp), allocatable   :: kei(:, :)   ! Kelvin function filter values
        real(wp), allocatable   :: GV(:, :)    ! Green's function values
        real(wp), allocatable   :: GE(:, :)    ! Green's function for elastic displacement (Farrell 1972)
        real(wp), allocatable   :: GN(:, :)    ! Green's function for dz_ss

        complex(dp), allocatable :: FGV(:, :)    ! FFT of GV
        complex(dp), allocatable :: FGE(:, :)    ! FFT of GE
        complex(dp), allocatable :: FGN(:, :)    ! FFT of GN

    end type isos_domain_class

    type isos_state_class
        real(wp)              :: t              ! [yr] Time
        real(wp)              :: bsl            ! [m] BSL
        real(wp)              :: V_af           ! [m^3] Volume contribution of ice above floatation
        real(wp)              :: V_pov          ! [m^3] Potential ocean volume contribution
        real(wp)              :: V_den          ! [m^3] Density-driven ocean volume contribution
        real(wp)              :: deltaV_bsl     ! [m^3] Volume contrib. from domain wrt previous time step
        real(wp)              :: dbsl_total     ! [m^3] Cumulative BSL contrib. from domain

        real(wp), allocatable :: z_bed(:, :)          ! Bedrock elevation         [m]
        real(wp), allocatable :: dwdt(:, :)           ! Rate of bedrock uplift    [m/a]
        real(wp), allocatable :: w(:, :)              ! Current viscous displacement
        real(wp), allocatable :: w_equilibrium(:, :)  ! Current viscous equilibrium displacement (XLRA)
        real(wp), allocatable :: we(:, :)             ! [m] Elastic displacement

        real(wp), allocatable :: Haf(:, :)              ! [m] Ice thickness above floatation
        real(wp), allocatable :: Hice(:, :)             ! [m] Thickness of ice column
        ! real(wp), allocatable :: Hsediment(:, :)      ! [m] Thickness of sediment column

        real(wp), allocatable :: z_ss(:, :)           ! [m] sea-surface height
        real(wp), allocatable :: dz_ss(:, :)   ! [m] sea-surface height perturbation
        real(wp), allocatable :: canom_ice(:, :)
        real(wp), allocatable :: canom_sediment(:, :)
        real(wp), allocatable :: canom_seawater(:, :)
        real(wp), allocatable :: canom_litho(:, :)
        real(wp), allocatable :: canom_mantle(:, :)
        real(wp), allocatable :: canom_load(:, :)    ! [kg m^-2] Load column anomaly
        real(wp), allocatable :: canom_full(:, :)    ! [kg m^-2] Full column anomaly
        real(wp), allocatable :: mass_anom(:, :)     ! [kg] Mass anomaly

        logical, allocatable :: maskocean(:, :)      ! [1] Ocean mask
        logical, allocatable :: maskgrounded(:, :)   ! [1] Grounded mask
        logical, allocatable :: maskcontinent(:, :)  ! [1] Continent mask
    end type isos_state_class

    type isos_out_class
        real(wp), allocatable       :: Hice(:, :)
        real(wp), allocatable       :: canom_full(:, :)
        real(wp), allocatable       :: dwdt(:, :)
        real(wp), allocatable       :: w(:, :)
        real(wp), allocatable       :: we(:, :)
        real(wp), allocatable       :: w_equilibrium(:, :)

        real(wp), allocatable       :: z_ss(:, :)
        real(wp), allocatable       :: dz_ss(:, :)
        real(wp), allocatable       :: z_bed(:, :)

        logical, allocatable        :: maskocean(:, :)
        logical, allocatable        :: maskgrounded(:, :)
        logical, allocatable        :: maskcontinent(:, :)
    end type isos_out_class

    type isos_class
        type(isos_param_class)  :: par
        type(isos_domain_class) :: domain
        type(isos_state_class)  :: now
        type(isos_state_class)  :: ref
        type(isos_out_class)    :: out
    end type

    public :: sp, dp, wp
    public :: isos_param_class
    public :: isos_domain_class
    public :: isos_state_class
    public :: isos_class

    contains ! nothing but definition of derived types

end module isostasy_defs
