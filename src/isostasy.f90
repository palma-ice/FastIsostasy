!# TODO: We should rename the module to something more specific. For instance:
! RegionalGIA
module isostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    ! Note: currently routine `isos_par_load` has dependency on the nml.f90 module 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use, intrinsic :: iso_c_binding

    use isostasy_defs, only : wp, pi, isos_param_class, isos_state_class, &
        isos_domain_class, isos_output_class, isos_class
    use green_functions
    use kelvin_function
    use isos_utils
    use lv_xlra
    use convolutions
    use lv_elva
    use sealevel
    use ncio

    implicit none
    include 'fftw3.f03'

    private
    public :: isos_param_class
    public :: isos_state_class
    public :: isos_domain_class
    public :: isos_class
    public :: isos_init
    public :: isos_init_state
    public :: isos_update
    public :: isos_end

    contains
    
    subroutine isos_init(isos, filename, group, nx, ny, dx, dy, E, nu, &
        rho_water, rho_ice, rho_seawater, rho_uppermantle, rho_litho, g, r_earth, m_earth, &
        A_ocean_pd, visc_method, rigidity_method, boundaries, viscosities, &
        interactive_sealevel, correct_distortion, method, dt_prognostics, dt_diagnostics, K)

        implicit none

        type(isos_class), intent(INOUT) :: isos
        character(len=*), intent(IN)    :: filename
        character(len=*), intent(IN)    :: group
        integer,  intent(IN)  :: nx
        integer,  intent(IN)  :: ny
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        real(wp), intent(IN), optional  :: E
        real(wp), intent(IN), optional  :: nu
        real(wp), intent(IN), optional  :: rho_water
        real(wp), intent(IN), optional  :: rho_ice
        real(wp), intent(IN), optional  :: rho_seawater
        real(wp), intent(IN), optional  :: rho_uppermantle
        real(wp), intent(IN), optional  :: rho_litho
        real(wp), intent(IN), optional  :: g
        real(wp), intent(IN), optional  :: r_earth
        real(wp), intent(IN), optional  :: m_earth
        real(wp), intent(IN), optional  :: A_ocean_pd

        character(len=*), intent(IN), optional  :: visc_method
        character(len=*), intent(IN), optional  :: rigidity_method
        real(wp), intent(IN), optional  :: boundaries(:)
        real(wp), intent(IN), optional  :: viscosities(:)

        real(wp), intent(IN), optional  :: dt_prognostics
        real(wp), intent(IN), optional  :: dt_diagnostics
        logical, intent(IN), optional   :: interactive_sealevel
        logical, intent(IN), optional   :: correct_distortion
        integer, intent(IN), optional   :: method
        real(wp), intent(IN), optional  :: K(:, :)


        ! Local variables
        integer                 :: n, nbsl, i, j, l, ncz    ! ncz = helper to load 3D fields
        real(wp), allocatable   :: z(:)
        real(wp), allocatable   :: helper_convo(:, :)
        real(wp), allocatable   :: eta_raw(:, :, :)
        real(wp)                :: D_lith_const
        character*256           :: filename_laty
        character*256           :: ocean_surface_file

        ! First, load parameters from parameter file `filename`
        write(*,*) "Defining params..."
        call isos_par_load(isos%par, filename, group)

        if (present(interactive_sealevel))  isos%par%interactive_sealevel   = interactive_sealevel
        if (present(correct_distortion))    isos%par%correct_distortion     = correct_distortion
        if (present(method))                isos%par%method                 = method
        if (present(dt_prognostics))        isos%par%dt_prognostics         = dt_prognostics
        if (present(dt_diagnostics))        isos%par%dt_diagnostics         = dt_diagnostics

        ! Overwrite physical constants with arguments if available
        if (present(E))                 isos%par%E                  = E
        if (present(nu))                isos%par%nu                 = nu
        if (present(rho_ice))           isos%par%rho_ice            = rho_ice
        if (present(rho_seawater))      isos%par%rho_seawater       = rho_seawater
        if (present(rho_water))         isos%par%rho_water          = rho_water
        if (present(rho_uppermantle))   isos%par%rho_uppermantle    = rho_uppermantle
        if (present(rho_litho))         isos%par%rho_litho          = rho_litho
        if (present(g))                 isos%par%g                  = g
        if (present(r_earth))           isos%par%r_earth            = r_earth
        if (present(m_earth))           isos%par%m_earth            = m_earth
        if (present(A_ocean_pd))        isos%par%A_ocean_pd         = A_ocean_pd

        if (present(rigidity_method))   isos%par%rigidity_method    = rigidity_method
        if (present(visc_method))       isos%par%visc_method        = visc_method
        if (present(boundaries))        isos%par%boundaries         = boundaries
        if (present(viscosities))       isos%par%viscosities        = viscosities

        write(*,*) "Allocating fields..."
        n = max(nx, ny)
        ocean_surface_file = "isostasy_data/tools/ocean_surface/OceanSurfaceETOPO2022.nc"
        ! nbsl = nc_size(ocean_surface_file, "z")
        call allocate_isos(isos, n, n, isos%par%nl, nbsl)

        ! Init scalar fields of domain
        isos%domain%nx = n
        isos%domain%ny = n
        isos%domain%dx = dx
        isos%domain%dy = dy

        ! Init ocean surface interpolator
        ! call nc_read(ocean_surface_file, "z", isos%domain%bsl_vec, start=[1], count=[nbsl])
        ! call nc_read(ocean_surface_file, "A", isos%domain%A_ocean_vec, start=[1], &
        !     count=[nbsl])

        call calc_cropindices(isos%domain%icrop1, isos%domain%icrop2, &
            isos%domain%jcrop1, isos%domain%jcrop2, nx, ny)

        ! Init plans
        write(*,*) "Computing FFT plans..."
        allocate(helper_convo(2*n-1, 2*n-1))
        helper_convo = 0.0
        isos%domain%forward_fftplan_r2r = fftw_plan_r2r_2d(n, n, isos%now%w, &
            isos%now%w, FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)
        isos%domain%backward_fftplan_r2r = fftw_plan_r2r_2d(n, n, isos%now%w, &
            isos%now%w, FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)
        isos%domain%forward_dftplan_r2c = fftw_plan_dft_r2c_2d(2*n-1, 2*n-1, &
            helper_convo, isos%domain%FGE, 1)
        isos%domain%backward_dftplan_c2r = fftw_plan_dft_c2r_2d(2*n-1, 2*n-1, &
            isos%domain%FGE, helper_convo, 1)

        !# TODO: provide square fields for params even if yelmo grid is rectangular.
        ! Init domain
        write(*,*) "Complementing domain informations..."
        if (isos%par%correct_distortion) then
            if (present(K)) then
                isos%domain%K(:, :) = 1
            else
                print*, 'Provide distortion matrix for it to be accounted for.'
                stop
            endif
        else
            isos%domain%K(:, :) = 1
        endif

        isos%domain%dx_matrix = isos%domain%dx * isos%domain%K
        isos%domain%dy_matrix = isos%domain%dy * isos%domain%K
        isos%domain%A = isos%domain%dx_matrix * isos%domain%dy_matrix
        call convenient_calc_convolution_indices(isos%domain)

        ! Init parameters
        isos%par%sec_per_year = 3600.0 * 24.0 * 365.25           ! [s/a]
        isos%par%compressibility_correction = 1.5 / (1 + isos%par%nu)
        call calc_density_correction_factor(isos%par)
        call calc_homogeneous_rigidity(D_lith_const, isos%par%E, &
            isos%par%boundaries(1), isos%par%nu)
        call calc_flexural_lengthscale(isos%par%L_w, D_lith_const, &
            isos%par%rho_uppermantle, isos%par%g)

        isos%domain%D_lith = D_lith_const

        call calc_elastic_green(isos%domain%GE, dx=isos%domain%dx, dy=isos%domain%dx)
        call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GE, &
        isos%domain%FGE, isos%domain%nx, isos%domain%ny)
        if (isos%par%interactive_sealevel) then
            call calc_ssh_green(isos%domain%GN, isos%par%m_earth, &
                isos%par%r_earth, dx=isos%domain%dx, dy=isos%domain%dx)
            call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GN, &
                isos%domain%FGN, isos%domain%nx, isos%domain%ny)
        endif

        select case(isos%par%method)

            ! LV-ELRA method is being used, which allows heterogeneous value of D_Lith(x, y)
            ! and tau(x, y). ELRA is a particular case of LV-ELRA with constant values.
            case(2)
                write(*,*) "Using (laterally-variable) ELRA..."

                ! Calculate the Kelvin function filter
                call calc_kei_filter_2D(isos%domain%kei, L_w=isos%par%L_w, &
                    dx=isos%domain%dx, dy=isos%domain%dx)

                ! Calculate the reference Green's function values
                call calc_viscous_green(isos%domain%GV, isos%domain%kei, &
                    isos%par%L_w, D_lith_const, dx=isos%domain%dx, dy=isos%domain%dx)

                write(*,*) "Initialising viscous Green kernel..."
                call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GV, &
                    isos%domain%FGV, isos%domain%nx, isos%domain%ny)

                !# TODO: allow heterogeneous tau
                isos%domain%tau        = isos%par%tau          ! [yr]

            ! LV-ELVA method is being used, which allows heterogeneous value of D_Lith(x, y)
            ! and eta_eff(x, y). ELVA is a particular case of LV-ELVA with constant values.
            case(3)
                write(*,*) "Using (laterally-variable) ELVA..."
                call convenient_calc_kappa(isos%domain)
                isos%domain%maskactive  = .true.

                write(*,*) "Choosing rigidity field..."
                select case(trim(isos%par%rigidity_method))
                    !# TODO: this part of the code should be in test_isostasy.f90
                    case("uniform")
                        isos%domain%D_lith   = D_lith_const
                        isos%domain%He_lith  = isos%par%boundaries(1)

                    case("gaussian_plus")
                        call calc_gaussian_rigidity(isos%domain%He_lith, &
                            1.5*isos%par%boundaries(1), &
                            isos%par%boundaries(1), sign=1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                    case("gaussian_minus")
                        call calc_gaussian_rigidity(isos%domain%He_lith, &
                            1.5*isos%par%boundaries(1), &
                            isos%par%boundaries(1), sign=1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                    case("laty")
                        filename_laty = "isostasy_data/earth_structure/ANT-32KM_Latychev.nc"
                        call nc_read(filename_laty, "lithos_thck", isos%domain%He_lith, &
                            start=[1, 1], count=[isos%domain%nx, isos%domain%ny])

                        ! calc_heterogeneous_rigidity assumes km, not m. Therefore conversion.
                        isos%domain%He_lith = isos%domain%He_lith*1.e-3
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                    case DEFAULT
                        ! do nothing, eta was set above
                        print*,'what is the default?' 
                        stop
                end select

                call calc_layerboundaries(isos%domain%boundaries, isos%domain%He_lith, &
                    isos%par%boundaries)

                write(*,*) "Choosing viscosity field..."
                call layered_viscosity(isos%domain%eta, isos%par%viscosities)
                select case(trim(isos%par%visc_method))
                !# TODO: this part of the code should be in test_isostasy.f90
                    case("uniform")
                        call calc_effective_viscosity(isos%domain%eta_eff, isos%domain%eta, &
                            isos%domain%dx, isos%domain%dy, isos%domain%boundaries)

                    case("gaussian_plus")
                        call calc_gaussian_viscosity(isos%domain%eta_eff, &
                            isos%par%viscosities(1), +1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)

                    case("gaussian_minus")
                        call calc_gaussian_viscosity(isos%domain%eta_eff, &
                            isos%par%viscosities(1), -1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)

                    case("laty")
                        !# TODO: loading fields should be more general (with 2D interpolation)
                        filename_laty = "isostasy_data/earth_structure/ANT-32KM_Latychev.nc"
                        ncz = nc_size(filename_laty,"zc")
                        allocate(z(ncz))
                        allocate(eta_raw(n, n, ncz))
                        call nc_read(filename_laty, "zc", z, start=[1], count=[ncz])
                        call nc_read(filename_laty, "log10_mantle_visc", eta_raw, &
                            start=[1, 1, 1], count=[n, n, ncz])
                        where (eta_raw < 16) eta_raw = 1e16

                        do l = 1, isos%par%nl
                            do i = 1, nx
                            do j = 1, ny
                                call interp_0d(z, eta_raw(i, j, :), &
                                    isos%domain%boundaries(i, j, l), isos%domain%eta(i, j, l))
                            end do
                            end do
                        end do

                        isos%domain%eta = 10. ** (isos%domain%eta)
                        call calc_effective_viscosity(isos%domain%eta_eff, isos%domain%eta, &
                            isos%domain%dx, isos%domain%dy, isos%domain%boundaries)
                        
                        ! use w_equilibrium as helper to load mask since it is not used.
                        filename_laty = "isostasy_data/tools/masks/ANT-32KM_activemask.nc"
                        call nc_read(filename_laty, "M", &
                            isos%now%w_equilibrium, start=[1, 1], &
                            count=[isos%domain%nx, isos%domain%ny])
                        where(isos%now%w_equilibrium .lt. 0.5) isos%domain%maskactive = .false.

                    case DEFAULT
                        ! do nothing, eta was set above
                        print*,'what is the default?' 
                        stop
                end select

                isos%domain%eta_eff = isos%domain%eta_eff * isos%par%compressibility_correction

                write(*,*) "isos_init:: summary"
                write(*,*) "    E               : ",    isos%par%E 
                write(*,*) "    nu              : ",    isos%par%nu
                write(*,*) "    rho_ice         : ",    isos%par%rho_ice
                write(*,*) "    rho_seawater    : ",    isos%par%rho_seawater
                write(*,*) "    rho_water       : ",    isos%par%rho_water
                write(*,*) "    rho_uppermantle : ",    isos%par%rho_uppermantle
                write(*,*) "    rho_litho       : ",    isos%par%rho_litho
                write(*,*) "    g               : ",    isos%par%g
                write(*,*) "    r_earth         : ",    isos%par%r_earth
                write(*,*) "    m_earth         : ",    isos%par%m_earth
                write(*,*) "    A_ocean_pd      : ",    isos%par%A_ocean_pd
                
                write(*,*) "    range(eta_eff):  ", minval(isos%domain%eta_eff),    maxval(isos%domain%eta_eff) 
                write(*,*) "    range(He_lith):  ", minval(isos%domain%He_lith),    maxval(isos%domain%He_lith)
                
                write(*,*) "    L_w (m):      ", isos%par%L_w
                write(*,*) "    nx:           ", isos%domain%nx
                write(*,*) "    ny:           ", isos%domain%ny
                write(*,*) "    dx (m):       ", isos%domain%dx
                write(*,*) "    dy (m):       ", isos%domain%dy

                write(*,*) "    range(kei): ", minval(isos%domain%kei),     maxval(isos%domain%kei)
                write(*,*) "    range(GV):  ", minval(isos%domain%GV),      maxval(isos%domain%GV)
                write(*,*) "    range(GE):  ", minval(isos%domain%GE),      maxval(isos%domain%GE)
                write(*,*) "    range(GN):  ", minval(isos%domain%GN),      maxval(isos%domain%GN)

            case DEFAULT 

                ! Set elastic length scale to zero (not used)
                isos%par%L_w = 0.0
                
                ! Set regional filter fields to zero (not used)
                isos%domain%kei = 0.0 
                isos%domain%GV  = 0.0
                isos%domain%GE  = 0.0
                isos%domain%GN  = 0.0
                
                isos%domain%tau         = isos%par%tau          ! [yr]
                isos%domain%maskactive  = .true.

            end select

        !# TODO: the initialisation below should be more generic and coherent with isos_init_state.
        ! isos%now%bsl            = -1e10
        isos%now%bsl            = 0.0

        ! Set time to very large value in the future 
        isos%par%time_diagnostics = 1e10 
        isos%par%time_prognostics = 1e10

        ! This is then changed by isos_init_state
        isos%now%z_bed          = 0.0
        isos%now%dwdt           = 0.0
        isos%now%w              = 0.0
        isos%now%w_equilibrium  = 0.0
        isos%now%we             = 0.0

        isos%now%Haf            = 0.0
        isos%now%Hice           = 0.0

        isos%now%rsl            = -1e6_wp
        isos%now%ssh            = 0.0
        isos%now%ssh_perturb    = 0.0
        isos%now%canom_load     = 0.0
        isos%now%canom_full     = 0.0
        isos%now%mass_anom      = 0.0

        isos%now%maskocean      = .false.
        isos%now%maskgrounded   = .false.
        isos%now%maskcontinent  = .false.

        call cropdomain2output(isos%output, isos%domain)
        write(*,*) "isos_init:: complete." 

        return
    end subroutine isos_init

    subroutine isos_init_state(isos, z_bed, H_ice, ssh, time)

        implicit none

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: z_bed(:, :)         ! [m] Current bedrock elevation
        real(wp), intent(IN) :: H_ice(:, :)         ! [m] Current ice thickness
        real(wp), intent(IN) :: ssh(:, :)           ! [m] Current sea-surface height
        real(wp), intent(IN) :: time                ! [a] Initial time 

        ! Store initial bedrock field
        call extendice2isostasy(isos%now, z_bed, H_ice, ssh, &
            isos%domain%icrop1, isos%domain%icrop2, &
            isos%domain%jcrop1, isos%domain%jcrop2)
        ! call copy_state(isos%ref, isos%now)

        ! Define initial time of isostasy model
        ! (set time_diagnostics earlier, so that it is definitely updated on the first timestep)
        isos%par%time_diagnostics = time - isos%par%dt_diagnostics
        isos%par%time_prognostics = time

        write(*,*) "Calling first update..."
        ! Call isos_update to diagnose rate of change
        ! (no change to z_bed will be applied since isos%par%time==time)
        ! isos%now%ssh = isos%now%bsl + isos%ref%ssh + isos%now%ssh_perturb
        isos%now%rsl = isos%now%ssh - isos%now%z_bed
        call calc_masks(isos)
        call calc_columnanoms_load(isos)
        call calc_columnanoms_solidearth(isos)
        call copy_sparsestate(isos%ref, isos%now)
        call isos_update(isos, H_ice, time)
        call copy_state(isos%ref, isos%now)

        write(*,*) "isos_init_state: "
        write(*,*) "    Initial time:   ",  isos%par%time_prognostics
        write(*,*) "    range(He_lith): ",  minval(isos%domain%He_lith), maxval(isos%domain%He_lith)
        write(*,*) "    range(tau):     ",  minval(isos%domain%tau), maxval(isos%domain%tau)
        write(*,*) "    range(w):       ",  minval(isos%now%w), maxval(isos%now%w)
        write(*,*) "    range(w_eq):    ",  minval(isos%now%w_equilibrium), maxval(isos%now%w_equilibrium)
        write(*,*) "  range(z_bed):     ",  minval(isos%now%z_bed), maxval(isos%now%z_bed)

        if ((minval(isos%domain%tau) .le. 0.0) .and. (isos%par%method .le. 2)) then
            write(error_unit,*) "isos_init_state:: Error: tau initialized with zero values present. &
            &This will lead to the model crashing."
        end if
        
        return

    end subroutine isos_init_state

    subroutine isos_update(isos, H_ice, time, dwdt_corr) 

        implicit none

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN)            :: H_ice(:, :)      ! [m] Current ice thickness
        real(wp), intent(IN)            :: time             ! [a] Current time
        real(wp), intent(IN), optional  :: dwdt_corr(:, :) ! [m/yr] Basal topography displacement rate (ie, to relax from low resolution to high resolution) 

        ! Local variables
        real(wp) :: dt, dt_now
        integer  :: n, nstep
        logical  :: update_diagnostics
 
        isos%now%Hice = 0.0
        isos%now%Hice(isos%domain%icrop1:isos%domain%icrop2, &
            isos%domain%jcrop1:isos%domain%jcrop2) = H_ice

        ! Step 0: determine current timestep and number of iterations
        dt = time - isos%par%time_prognostics

        ! Get maximum number of iterations needed to reach desired time
        nstep = ceiling( (time - isos%par%time_prognostics) / isos%par%dt_prognostics )
        nstep = max(nstep, 1)

        ! Loop over iterations until maximum time is reached
        do n = 1, nstep

            ! Get current dt (either total time or maximum allowed timestep)
            dt_now = min(time-isos%par%time_prognostics, isos%par%dt_prognostics)

            ! Only update diagnostics if enough time has passed (to save on computations)
            if ( (isos%par%time_prognostics+dt_now) - isos%par%time_diagnostics .ge. &
                isos%par%dt_diagnostics) then
                update_diagnostics = .TRUE.
            else
                update_diagnostics = .FALSE.
            end if

            call calc_columnanoms_load(isos)

            ! update elastic response
            if (update_diagnostics) then
                call precomputed_fftconvolution(isos%now%we, isos%domain%FGE, &
                    isos%now%canom_load * isos%par%g * isos%domain%K ** 2.0, &
                    isos%domain%i1, isos%domain%i2, &
                    isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                    isos%domain%nx, isos%domain%ny, &
                    isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)
            endif

            call calc_columnanoms_solidearth(isos)

            ! write(*,*) "Updating the sea-level..."
            if (update_diagnostics .and. isos%par%interactive_sealevel) then
                ! write(*,*) "Updating the ssh perturbation..."
                call calc_mass_anom(isos)
                call precomputed_fftconvolution(isos%now%ssh_perturb, isos%domain%FGN, &
                    isos%now%mass_anom, isos%domain%i1, isos%domain%i2, &
                    isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                    isos%domain%nx, isos%domain%ny, &
                    isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)
                
                ! write(*,*) "Updating sea-level contributions..."
                ! call calc_sl_contribution(isos)
                isos%now%ssh = isos%now%bsl + isos%ref%ssh + isos%now%ssh_perturb
                isos%now%rsl = isos%now%ssh - isos%now%z_bed

                ! write(*,*) "Updating masks..."
                call calc_masks(isos)
            endif

            ! Need to re-compute the load after updateing the ssh
            call calc_columnanoms_load(isos)
            call calc_columnanoms_solidearth(isos)

            ! Step 1: diagnose equilibrium displacement and rate of bedrock uplift
            select case(isos%par%method)

                ! Steady-state lithosphere
                case(0)
                    call calc_llra_equilibrium(isos%now%w_equilibrium, &
                        isos%now%canom_load, isos%par%rho_uppermantle)
                        isos%now%w    = isos%now%w_equilibrium
                        isos%now%dwdt = 0.0

                ! LLRA
                case(1)
                    call calc_llra_equilibrium(isos%now%w_equilibrium, &
                        isos%now%canom_load, isos%par%rho_uppermantle)
                    call calc_asthenosphere_relax(isos%now%dwdt, isos%now%w, &
                        isos%now%w_equilibrium, isos%domain%tau)

                ! LV-ELRA (Coulon et al. 2021).
                ! Gives ELRA (LeMeur and Huybrechts 1996) if tau = const.
                case(2)
                    
                    if (update_diagnostics) then
                        call precomputed_fftconvolution(isos%now%w_equilibrium, isos%domain%FGV, &
                            isos%now%canom_load * isos%par%g * isos%domain%K ** 2.0, &
                            isos%domain%i1, isos%domain%i2, &
                            isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                            isos%domain%nx, isos%domain%ny, &
                            isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)
                        call apply_zerobc_at_corners(isos%now%w_equilibrium, isos%domain%nx, &
                            isos%domain%ny)
                    end if
                    call calc_asthenosphere_relax(isos%now%dwdt, isos%now%w, &
                        isos%now%w_equilibrium, isos%domain%tau)

                ! LV-ELVA (Swierczek-jereczek et al. 2024).
                ! Gives ELVA (Bueler et al. 2007) if eta = const.
                case(3)
                    call calc_lvelva(isos%now%dwdt, isos%now%w, isos%now%canom_full, &
                        isos%domain%maskactive, isos%par%g, isos%par%nu, isos%domain%D_lith, &
                        isos%domain%eta_eff, isos%domain%kappa, isos%domain%nx, isos%domain%ny, &
                        isos%domain%dx_matrix, isos%domain%dy_matrix, isos%par%sec_per_year, &
                        isos%domain%forward_fftplan_r2r, isos%domain%backward_fftplan_r2r)
                end select

            if (update_diagnostics) then 
                ! Update current time_diagnostics value
                isos%par%time_diagnostics = isos%par%time_prognostics + dt_now
            end if

            ! Step 2: update bedrock elevation and current model time
            if (dt_now .gt. 0.0) then

                !# TODO: do we need this? If yes, lines below should be adapted adequatly
                ! Additionally apply bedrock adjustment field
                if (present(dwdt_corr)) then 
                    isos%now%z_bed = isos%now%z_bed + dwdt_corr*dt_now
                end if

                isos%now%w = isos%now%w + isos%now%dwdt*dt_now
                call apply_zerobc_at_corners(isos%now%w, isos%domain%nx, isos%domain%ny)
                isos%now%z_bed = isos%ref%z_bed + isos%now%w + isos%now%we
                isos%par%time_prognostics = isos%par%time_prognostics + dt_now

            end if

            if ( abs(time-isos%par%time_prognostics) .lt. 1e-5) then 
                ! Final time has been reached, exit the loop 
                isos%par%time_prognostics = time 
                exit
            end if
        end do

        call cropstate2output(isos%output, isos%now, isos%domain%icrop1, isos%domain%icrop2, &
            isos%domain%jcrop1, isos%domain%jcrop2)

        ! rsl = isos%output%ssh - isos%output%z_bed

        return
    end subroutine isos_update

    subroutine isos_end(isos)

        implicit none 
        type(isos_class), intent(INOUT) :: isos 

        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call deallocate_isos_domain(isos%domain)
        call deallocate_isos_output(isos%output)

        return
    end subroutine isos_end

    subroutine isos_par_load(par, filename, group)

        use nml

        implicit none

        type(isos_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename 
        character(len=*),       intent(IN)  :: group 

        call nml_read(filename,group,"interactive_sealevel",par%interactive_sealevel)
        call nml_read(filename,group,"correct_distortion",  par%correct_distortion)
        call nml_read(filename,group,"method",              par%method)
        call nml_read(filename,group,"dt_diagnostics",      par%dt_diagnostics)
        call nml_read(filename,group,"dt_prognostics",      par%dt_prognostics)
        call nml_read(filename,group,"visc_method",         par%visc_method) 
        call nml_read(filename,group,"rigidity_method",     par%rigidity_method)

        call nml_read(filename,group,"nl",                  par%nl)
        allocate(par%boundaries(par%nl))
        allocate(par%viscosities(par%nl))

        call nml_read(filename,group,"boundaries",          par%boundaries)
        call nml_read(filename,group,"viscosities",         par%viscosities)
        call nml_read(filename,group,"tau",                 par%tau)

        ! Physical constants
        call nml_read(filename,group,"E",               par%E)
        call nml_read(filename,group,"nu",              par%nu)
        call nml_read(filename,group,"rho_water",       par%rho_water)
        call nml_read(filename,group,"rho_ice",         par%rho_ice)
        call nml_read(filename,group,"rho_seawater",    par%rho_seawater)
        call nml_read(filename,group,"rho_uppermantle", par%rho_uppermantle)
        call nml_read(filename,group,"rho_litho",       par%rho_litho)
        call nml_read(filename,group,"g",               par%g)
        call nml_read(filename,group,"r_earth",         par%r_earth)
        call nml_read(filename,group,"m_earth",         par%m_earth)
        call nml_read(filename,group,"A_ocean_pd",      par%A_ocean_pd)

        return
    end subroutine isos_par_load

    subroutine allocate_isos(isos, nx, ny, nl, nbsl)
        implicit none
        type(isos_class), intent(INOUT) :: isos
        integer, intent(IN)             :: nx, ny, nl, nbsl

        ! First ensure arrays are not allocated
        call deallocate_isos_domain(isos%domain)
        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call deallocate_isos_output(isos%output)
        call allocate_isos_domain(isos%domain, nx, ny, nl, nbsl)
        call allocate_isos_state(isos%now, nx, ny)
        call allocate_isos_state(isos%ref, nx, ny)
        call allocate_isos_output(isos%output, nx, ny)

        return
    end subroutine allocate_isos

    subroutine deallocate_isos_domain(domain)
        implicit none 
        type(isos_domain_class), intent(INOUT) :: domain

        if (allocated(domain%bsl_vec))      deallocate(domain%bsl_vec)
        if (allocated(domain%A_ocean_vec))  deallocate(domain%A_ocean_vec)

        if (allocated(domain%dx_matrix))    deallocate(domain%dx_matrix)
        if (allocated(domain%dy_matrix))    deallocate(domain%dy_matrix)
        if (allocated(domain%A))            deallocate(domain%A)
        if (allocated(domain%K))            deallocate(domain%K)
        if (allocated(domain%kappa))        deallocate(domain%kappa)
        if (allocated(domain%maskactive))   deallocate(domain%maskactive)

        if (allocated(domain%He_lith))           deallocate(domain%He_lith)
        if (allocated(domain%D_lith))            deallocate(domain%D_lith)
        if (allocated(domain%eta))               deallocate(domain%eta)
        if (allocated(domain%eta_eff))           deallocate(domain%eta_eff)
        if (allocated(domain%tau))               deallocate(domain%tau)

        if (allocated(domain%kei))       deallocate(domain%kei)
        if (allocated(domain%GV))        deallocate(domain%GV)
        if (allocated(domain%GE))        deallocate(domain%GE)
        if (allocated(domain%GN))        deallocate(domain%GN)

        if (allocated(domain%FGV))        deallocate(domain%FGV)
        if (allocated(domain%FGE))        deallocate(domain%FGE)
        if (allocated(domain%FGN))        deallocate(domain%FGN)
        return
    end subroutine deallocate_isos_domain

    subroutine deallocate_isos_state(state)
        implicit none 
        type(isos_state_class), intent(INOUT) :: state

        if (allocated(state%z_bed))             deallocate(state%z_bed)
        if (allocated(state%dwdt))              deallocate(state%dwdt)
        if (allocated(state%w))                 deallocate(state%w)
        if (allocated(state%w_equilibrium))     deallocate(state%w_equilibrium)
        if (allocated(state%we))                deallocate(state%we)

        if (allocated(state%ssh_perturb))       deallocate(state%ssh_perturb)
        if (allocated(state%Haf))               deallocate(state%Haf)
        if (allocated(state%Hice))              deallocate(state%Hice)

        if (allocated(state%rsl))               deallocate(state%rsl)
        if (allocated(state%ssh))               deallocate(state%ssh)
        if (allocated(state%canom_load))        deallocate(state%canom_load)
        if (allocated(state%canom_full))        deallocate(state%canom_full)
        if (allocated(state%mass_anom))         deallocate(state%mass_anom)

        if (allocated(state%maskocean))         deallocate(state%maskocean)
        if (allocated(state%maskgrounded))      deallocate(state%maskgrounded)
        if (allocated(state%maskcontinent))     deallocate(state%maskcontinent)

        return 
    end subroutine deallocate_isos_state

    subroutine deallocate_isos_output(output)
        implicit none 
        type(isos_output_class), intent(INOUT) :: output

        if (allocated(output%He_lith))          deallocate(output%He_lith)
        if (allocated(output%D_lith))           deallocate(output%D_lith)
        if (allocated(output%eta_eff))          deallocate(output%eta_eff)
        if (allocated(output%tau))              deallocate(output%tau)
        if (allocated(output%kappa))            deallocate(output%kappa)
        if (allocated(output%kei))              deallocate(output%kei)
        if (allocated(output%GE))               deallocate(output%GE)
        if (allocated(output%GV))               deallocate(output%GV)
        if (allocated(output%GN))               deallocate(output%GN)

        if (allocated(output%Hice))             deallocate(output%Hice)
        if (allocated(output%canom_full))       deallocate(output%canom_full)
        if (allocated(output%dwdt))            deallocate(output%dwdt)

        if (allocated(output%w))                deallocate(output%w)
        if (allocated(output%we))               deallocate(output%we)
        if (allocated(output%w_equilibrium))    deallocate(output%w_equilibrium)

        if (allocated(output%rsl))              deallocate(output%rsl)
        if (allocated(output%ssh))              deallocate(output%ssh)
        if (allocated(output%ssh_perturb))      deallocate(output%ssh_perturb)
        if (allocated(output%z_bed))            deallocate(output%z_bed)

        if (allocated(output%maskocean))        deallocate(output%maskocean)
        if (allocated(output%maskgrounded))     deallocate(output%maskgrounded)
        if (allocated(output%maskcontinent))    deallocate(output%maskcontinent)
        return
    end subroutine deallocate_isos_output

    subroutine allocate_isos_domain(domain, nx, ny, nl, nbsl)
        implicit none
        integer, intent(IN) :: nx, ny, nl, nbsl
        type(isos_domain_class), intent(INOUT) :: domain

        allocate(domain%bsl_vec(nbsl))
        allocate(domain%A_ocean_vec(nbsl))

        allocate(domain%dx_matrix(nx, ny))
        allocate(domain%dy_matrix(nx, ny))
        allocate(domain%A(nx, ny))
        allocate(domain%K(nx, ny))
        allocate(domain%kappa(nx, ny))
        allocate(domain%maskactive(nx, ny))

        allocate(domain%He_lith(nx, ny))
        allocate(domain%D_lith(nx, ny))
        allocate(domain%eta_eff(nx, ny))
        allocate(domain%tau(nx, ny))
        allocate(domain%eta(nx, ny, nl))
        allocate(domain%boundaries(nx, ny, nl))

        allocate(domain%kei(nx, ny))
        allocate(domain%GV(nx, ny))
        allocate(domain%GE(nx, ny))
        allocate(domain%GN(nx, ny))

        allocate(domain%FGV(2*nx-1, 2*ny-1))
        allocate(domain%FGE(2*nx-1, 2*ny-1))
        allocate(domain%FGN(2*nx-1, 2*ny-1))

        return
    end subroutine allocate_isos_domain

    subroutine allocate_isos_state(state, nx, ny)
        implicit none
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny
        type(isos_state_class), intent(INOUT) :: state

        allocate(state%z_bed(nx, ny))
        allocate(state%dwdt(nx, ny))
        allocate(state%w(nx, ny))
        allocate(state%w_equilibrium(nx, ny))
        allocate(state%we(nx, ny))

        allocate(state%ssh_perturb(nx, ny))
        allocate(state%Haf(nx, ny))
        allocate(state%Hice(nx, ny))

        allocate(state%rsl(nx, ny))
        allocate(state%ssh(nx, ny))
        allocate(state%canom_load(nx, ny))
        allocate(state%canom_full(nx, ny))
        allocate(state%mass_anom(nx, ny))

        allocate(state%maskocean(nx, ny))
        allocate(state%maskgrounded(nx, ny))
        allocate(state%maskcontinent(nx, ny))

        return
    end subroutine allocate_isos_state

    subroutine allocate_isos_output(output, nx, ny)
        implicit none 
        type(isos_output_class), intent(INOUT) :: output
        integer, intent(IN) :: nx, ny

        allocate(output%He_lith(nx, ny))
        allocate(output%D_lith(nx, ny))
        allocate(output%eta_eff(nx, ny))
        allocate(output%tau(nx, ny))
        allocate(output%kappa(nx, ny))
        allocate(output%kei(nx, ny))
        allocate(output%GE(nx, ny))
        allocate(output%GV(nx, ny))
        allocate(output%GN(nx, ny))

        allocate(output%Hice(nx, ny))
        allocate(output%canom_full(nx, ny))
        allocate(output%dwdt(nx, ny))

        allocate(output%w(nx, ny))
        allocate(output%we(nx, ny))
        allocate(output%w_equilibrium(nx, ny))

        allocate(output%rsl(nx, ny))
        allocate(output%ssh(nx, ny))
        allocate(output%ssh_perturb(nx, ny))
        allocate(output%z_bed(nx, ny))

        allocate(output%maskocean(nx, ny))
        allocate(output%maskgrounded(nx, ny))
        allocate(output%maskcontinent(nx, ny))
        return
    end subroutine allocate_isos_output

end module isostasy