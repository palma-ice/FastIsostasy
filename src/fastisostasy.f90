module fastisostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    ! The details of the present implementation are described in
    ! Swierczek-Jereczek et al. (2024)

    ! Note: `isos_par_load` has dependency on the nml.f90 module 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
    use, intrinsic :: iso_c_binding
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    
    use isostasy_defs, only : sp, dp, wp, pi, isos_param_class, isos_state_class, &
        isos_domain_class, isos_out_class, isos_class
    use green_functions
    use kelvin_function
    use isos_utils
    use lv_xlra
    use convolutions
    use ncio
    use isostasy_io
    use lv_elva
    use sealevel

    implicit none
    include 'fftw3.f03'

    private
    public :: isos_param_class
    public :: isos_state_class
    public :: isos_domain_class
    public :: isos_class
    public :: isos_init
    public :: isos_init_ref
    public :: isos_init_state
    public :: isos_update
    public :: isos_end
    public :: isos_restart_write
    public :: isos_restart_read

    contains
    
    subroutine isos_init(isos, filename, group, nx_ice, ny_ice, dx, dy, rho_ice, K)

        implicit none

        type(isos_class), intent(INOUT) :: isos
        character(len=*), intent(IN)    :: filename
        character(len=*), intent(IN)    :: group
        integer,  intent(IN)  :: nx_ice
        integer,  intent(IN)  :: ny_ice
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        real(wp), intent(IN), optional  :: rho_ice
        real(wp), intent(IN), optional  :: K(:, :)

        ! Local variables
        logical                 :: correct_distortion
        integer                 :: nbsl, i, j, l, ncz    ! ncz = helper to load 3D fields
        real(wp), allocatable   :: z(:), depth(:)
        real(dp), allocatable   :: buffer_2n(:, :)
        real(dp), allocatable   :: buffer_n(:, :)
        real(wp), allocatable   :: eta_raw(:, :, :)
        real(wp)                :: D_lith_const
        real(wp)                :: viscosity_scaling
        logical, allocatable    :: mask_inner(:, :), mask_outer(:, :), mask_all(:, :)

        ! write(*,*) "Defining params..."
        call isos_par_load(isos%par, filename, group)

        if (present(K))         isos%par%correct_distortion     = .true.
        if (present(rho_ice))   isos%par%rho_ice                = rho_ice
        
        ! write(*,*) "Padding domain..."
        call init_domain_size(isos%domain, nx_ice, ny_ice, dx, dy, isos%par%min_pad)
        write(*,*) "icrop1: ", isos%domain%icrop1
        write(*,*) "icrop2: ", isos%domain%icrop2
        write(*,*) "jcrop1: ", isos%domain%jcrop1
        write(*,*) "jcrop2: ", isos%domain%jcrop2
        call convenient_calc_convolution_indices(isos%domain)

        if (isos%par%variable_ocean_surface) then
            nbsl = nc_size(isos%par%ocean_surface_file, "z")
            write(*,*) "Variable ocean surface not supported yet in the fortran version!"
            stop
        else
            nbsl = 1
        end if

        call allocate_isos(isos, nx_ice, ny_ice, isos%par%nl, nbsl)
        call init_dims(isos%domain%xc, isos%domain%x, isos%domain%nx, dx)
        call init_dims(isos%domain%yc, isos%domain%y, isos%domain%ny, dy)
        isos%domain%K(:, :) = 1
        if (present(K))         isos%domain%K                   = K

        allocate(mask_inner(isos%domain%nx, isos%domain%ny))
        allocate(mask_outer(isos%domain%nx, isos%domain%ny))
        allocate(mask_all(isos%domain%nx, isos%domain%ny))
        mask_inner = .false.
        mask_outer = .false.
        mask_all = .true.
        mask_inner(isos%domain%icrop1:isos%domain%icrop2, &
            isos%domain%jcrop1:isos%domain%jcrop2) = .true.
        mask_outer = .not. mask_inner

        if (isos%par%variable_ocean_surface) then
            ! write(*,*) "Reading ocean surface file..."
            call nc_read(isos%par%ocean_surface_file, "z", isos%domain%bsl_vec, &
                start=[1], count=[nbsl])
            call nc_read(isos%par%ocean_surface_file, "A", isos%domain%A_ocean_vec, &
                start=[1], count=[nbsl])
        end if

        ! write(*,*) "Initializing FFT plans..."
        allocate(buffer_2n(2*isos%domain%nx-1, 2*isos%domain%ny-1))
        buffer_2n = 0.0
        allocate(buffer_n(isos%domain%nx, isos%domain%ny))
        buffer_n = 0.0

        isos%domain%forward_fftplan_r2r = fftw_plan_r2r_2d(isos%domain%nx, isos%domain%ny, &
            buffer_n, buffer_n, FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)

        isos%domain%backward_fftplan_r2r = fftw_plan_r2r_2d(isos%domain%nx, isos%domain%ny, &
            buffer_n, buffer_n, FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)

        isos%domain%forward_dftplan_r2c = fftw_plan_dft_r2c_2d(2*isos%domain%nx-1, &
            2*isos%domain%ny-1, buffer_2n, isos%domain%FGE, 1)

        isos%domain%backward_dftplan_c2r = fftw_plan_dft_c2r_2d(2*isos%domain%nx-1, &
            2*isos%domain%ny-1, isos%domain%FGE, buffer_2n, 1)

        ! write(*,*) "Initializing distorted domain..."
        isos%domain%dx_matrix = isos%domain%dx * isos%domain%K
        isos%domain%dy_matrix = isos%domain%dy * isos%domain%K
        isos%domain%A = isos%domain%dx_matrix * isos%domain%dy_matrix

        ! write(*,*) "Initializing physical constants..."
        isos%par%sec_per_year = 3600.0 * 24.0 * 365.25           ! [s/a]
        isos%par%compressibility_correction = 1.5 / (1 + isos%par%nu)
        call calc_density_correction_factor(isos%par)

        ! write(*,*) "Initializing lithosphere..."
        isos%domain%He_lith = isos%par%zl(1)
        call calc_homogeneous_rigidity(D_lith_const, isos%par%E, &
            isos%par%zl(1), isos%par%nu)
        call calc_flexural_lengthscale(isos%par%L_w, D_lith_const, &
            isos%par%rho_uppermantle, isos%par%g)
        isos%domain%D_lith = D_lith_const

        ! write(*,*) "Initializing elastic Green's function..."
        call calc_elastic_green(isos%domain%GE, dx=isos%domain%dx, dy=isos%domain%dx)
        call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GE, &
            isos%domain%FGE, isos%domain%nx, isos%domain%ny)

        ! write(*,*) "Initializing gravitational Green's function..."
        if (isos%par%interactive_sealevel) then
            call calc_z_ss_green(isos%domain%GN, isos%par%m_earth, &
                isos%par%r_earth, dx=isos%domain%dx, dy=isos%domain%dx)
            call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GN, &
                isos%domain%FGN, isos%domain%nx, isos%domain%ny)
        endif

        ! write(*,*) "Initializing maskactive with the help of buffer..."
        buffer_n = 0.0
        isos%domain%maskactive(:, :) = .false.
        if (trim(isos%par%mask_file) .eq. "None" .or. &
            trim(isos%par%mask_file) .eq. "none" .or. &
            trim(isos%par%mask_file) .eq. "no") then
            
            if (isos%par%interactive_sealevel .eqv. .false.) then
                isos%domain%maskactive = .true.
            else if (isos%par%min_pad > 100.e3) then
                isos%domain%maskactive(isos%domain%icrop1:isos%domain%icrop2, &
                    isos%domain%jcrop1:isos%domain%jcrop2) = .true.
            else
                write(*,*) "interactive_sealevel=.true. requires a mask file or at least 100km of padding."
                stop
            end if
        else
            call nc_read(isos%par%mask_file, "M", &
                buffer_n(isos%domain%icrop1:isos%domain%icrop2, &
                isos%domain%jcrop1:isos%domain%jcrop2), start=[1, 1], &
                count=[nx_ice, ny_ice])
            where(buffer_n .gt. 0.5) isos%domain%maskactive = .true.
        end if

        select case(isos%par%method)

            ! LV-ELRA method is being used, which allows heterogeneous value of D_Lith(x, y)
            ! and tau(x, y). ELRA is a particular case of LV-ELRA with constant values.
            case(2)
                ! write(*,*) "Using (laterally-variable) ELRA..."

                ! Calculate the Kelvin function filter
                call calc_kei_filter_2D(isos%domain%kei, L_w=isos%par%L_w, &
                    dx=isos%domain%dx, dy=isos%domain%dx)

                ! Calculate the reference Green's function values
                call calc_viscous_green(isos%domain%GV, isos%domain%kei, &
                    isos%par%L_w, D_lith_const, dx=isos%domain%dx, dy=isos%domain%dx)

                ! write(*,*) "Initialising viscous Green kernel..."
                call precompute_kernel(isos%domain%forward_dftplan_r2c, isos%domain%GV, &
                    isos%domain%FGV, isos%domain%nx, isos%domain%ny)

                !# TODO: allow heterogeneous tau
                isos%domain%tau        = isos%par%tau          ! [yr]

            ! LV-ELVA method is being used, which allows heterogeneous value of D_Lith(x, y)
            ! and eta_eff(x, y). ELVA is a particular case of LV-ELVA with constant values.
            case(3)
                ! write(*,*) "Using (laterally-variable) ELVA..."
                call convenient_calc_kappa(isos%domain)

                write(*,*) "Choosing rigidity field..."
                select case(trim(isos%par%lithosphere))
                    !# TODO: this part of the code should be in test_isostasy.f90
                    case("uniform")
                        isos%domain%D_lith   = D_lith_const
                        isos%domain%He_lith  = isos%par%zl(1)

                    case("gaussian_plus")
                        call calc_gaussian_rigidity(isos%domain%He_lith, &
                            1.5*isos%par%zl(1), &
                            isos%par%zl(1), sign=1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                    case("gaussian_minus")
                        call calc_gaussian_rigidity(isos%domain%He_lith, &
                            1.5*isos%par%zl(1), &
                            isos%par%zl(1), sign=1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                    case("rheology_file")
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!!!! This gives an error when the input file does not have padding !!!!
                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        isos%domain%He_lith = isos%par%zl(1) * 1e-3     ! [km]
                        call nc_read(isos%par%rheology_file, "litho_thickness", &
                            isos%domain%He_lith(isos%domain%icrop1:isos%domain%icrop2, &
                            isos%domain%jcrop1:isos%domain%jcrop2), start=[1, 1], &
                            count=[nx_ice, ny_ice])

                        call smooth_gauss_2D(isos%domain%He_lith, isos%domain%dx, &
                            4 * isos%domain%dx, mask_outer, mask_all)
                        
                        ! rigidity computation assumes He_lith is in km...
                        call calc_heterogeneous_rigidity(isos%domain%D_lith, isos%par%E, &
                            isos%domain%He_lith, isos%par%nu)

                        ! ... but the remaining routines assumes m!
                        isos%domain%He_lith = isos%domain%He_lith * 1e3

                    case DEFAULT
                        ! do nothing, eta was set above
                        print*,'what is the default?'
                        stop
                end select

                if (isos%par%layering .eq. "parallel") then
                    call calc_parallel_layerboundaries(isos%domain%boundaries, &
                        isos%domain%He_lith, isos%par%nl, isos%par%dl)
                else if (isos%par%layering .eq. "equalized") then
                    call calc_equalized_layerboundaries(isos%domain%boundaries, &
                        isos%domain%He_lith, isos%par%nl, isos%par%zl)
                else if (isos%par%layering .eq. "folded") then
                    call calc_folded_layerboundaries(isos%domain%boundaries, &
                        isos%domain%He_lith, isos%par%nl, isos%par%zl)
                else
                    print*, "Unknown layering scheme: ", isos%par%layering
                    stop
                end if

                write(*,*) "Choosing viscosity field..."
                call layered_viscosity(isos%domain%eta, isos%par%viscosities)
                select case(trim(isos%par%mantle))
                !# TODO: this part of the code should be in test_isostasy.f90
                    case("uniform")
                        call calc_effective_viscosity(isos%domain%eta_eff, isos%domain%eta, &
                            isos%domain%xc, isos%domain%yc, isos%domain%boundaries)

                    case("gaussian_plus")
                        call calc_gaussian_viscosity(isos%domain%eta_eff, &
                            isos%par%viscosities(1), +1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)

                    case("gaussian_minus")
                        call calc_gaussian_viscosity(isos%domain%eta_eff, &
                            isos%par%viscosities(1), -1._wp, &
                            dx=isos%domain%dx, dy=isos%domain%dx)

                    case("rheology_file")
                        ncz = nc_size(isos%par%rheology_file, "zc")
                        
                        allocate(z(ncz))
                        allocate(depth(ncz))
                        allocate(eta_raw(nx_ice, ny_ice, ncz))
                        call nc_read(isos%par%rheology_file, "zc", z, start=[1], count=[ncz])
                        do i = 1, ncz
                            depth(i) = z(1) - z(i)
                        end do

                        call nc_read(isos%par%rheology_file, "log10_eta", eta_raw, &
                            start=[1, 1, 1], count=[nx_ice, ny_ice, ncz])
                        where (eta_raw < 16) eta_raw = 16
                        where (eta_raw > 23) eta_raw = 23

                        do l = 1, isos%par%nl
                            do i = 1, isos%domain%nx
                            do j = 1, isos%domain%ny
                                call interp_0d(depth, eta_raw(i, j, :), &
                                    isos%domain%boundaries(i, j, l), isos%domain%eta(i, j, l))
                            end do
                            end do
                        end do

                        isos%domain%eta = 10. ** (isos%domain%eta)

                        call calc_effective_viscosity(isos%domain%eta_eff, isos%domain%eta, &
                            isos%domain%xc, isos%domain%yc, isos%domain%boundaries)

                        call smooth_gauss_2D(isos%domain%eta_eff, isos%domain%dx, &
                            4 * isos%domain%dx, mask_outer, mask_all)
                        
                        if (isos%par%viscosity_scaling .eq. -1) then
                            isos%domain%eta_eff(:, :) = minval(isos%domain%eta_eff)
                        elseif (isos%par%viscosity_scaling .eq. 0) then
                            isos%domain%eta_eff = isos%domain%eta_eff
                        elseif (isos%par%viscosity_scaling .eq. 1) then
                            isos%domain%eta_eff = maxval(isos%domain%eta_eff)
                        else
                            isos%domain%eta_eff = viscosity_scaling * isos%domain%eta_eff
                        end if

                    case DEFAULT
                        ! do nothing, eta was set above
                        print*,'what is the default?' 
                        stop
                end select

                isos%domain%eta_eff = isos%domain%eta_eff * isos%par%compressibility_correction

                write(*,*) "isos_init:: summary"
                ! write(*,*) "    E               : ",    isos%par%E 
                ! write(*,*) "    nu              : ",    isos%par%nu
                ! write(*,*) "    rho_ice         : ",    isos%par%rho_ice
                ! write(*,*) "    rho_seawater    : ",    isos%par%rho_seawater
                ! write(*,*) "    rho_water       : ",    isos%par%rho_water
                ! write(*,*) "    rho_uppermantle : ",    isos%par%rho_uppermantle
                ! write(*,*) "    rho_litho       : ",    isos%par%rho_litho
                ! write(*,*) "    g               : ",    isos%par%g
                ! write(*,*) "    r_earth         : ",    isos%par%r_earth
                ! write(*,*) "    m_earth         : ",    isos%par%m_earth
                ! write(*,*) "    A_ocean_pd      : ",    isos%par%A_ocean_pd
                
                ! write(*,*) "    range(eta_eff):  ", minval(isos%domain%eta_eff),    maxval(isos%domain%eta_eff) 
                ! write(*,*) "    range(He_lith):  ", minval(isos%domain%He_lith),    maxval(isos%domain%He_lith)
                
                ! write(*,*) "    L_w (m):      ", isos%par%L_w
                ! write(*,*) "    nx:           ", isos%domain%nx
                ! write(*,*) "    ny:           ", isos%domain%ny
                ! write(*,*) "    dx (m):       ", isos%domain%dx
                ! write(*,*) "    dy (m):       ", isos%domain%dy

                ! write(*,*) "    range(kei): ", minval(isos%domain%kei),     maxval(isos%domain%kei)
                ! write(*,*) "    range(GV):  ", minval(isos%domain%GV),      maxval(isos%domain%GV)
                ! write(*,*) "    range(GE):  ", minval(isos%domain%GE),      maxval(isos%domain%GE)
                ! write(*,*) "    range(GN):  ", minval(isos%domain%GN),      maxval(isos%domain%GN)

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

        ! Set time to very large value in the future 
        isos%par%time_diagnostics = 1e10
        isos%par%time_prognostics = 1e10

        ! Initialize all values to zero initially for safety
        call zero_init_state(isos%now, z_bed_background=-1e6_wp)
        call zero_init_state(isos%ref, z_bed_background=-1e6_wp)

        ! Initially the reference state has not been defined
        isos%par%ref_was_set = .FALSE.

        ! write(*,*) "out%He_lith", size(isos%out%He_lith, 1), size(isos%out%He_lith, 2)
        call cropdomain2out(isos%out, isos%domain)
        ! write(*,*) "out%He_lith", size(isos%out%He_lith, 1), size(isos%out%He_lith, 2)

        ! write(*,*) "isos_init:: complete." 

        return
    end subroutine isos_init

    subroutine zero_init_state(state, z_bed_background)
        implicit none

        type(isos_state_class), intent(INOUT) :: state
        real(wp), intent(IN), optional :: z_bed_background

        state%bsl            = 0.0
        state%V_af           = 0.0
        state%V_den          = 0.0
        state%V_pov          = 0.0

        ! This is then changed by isos_init_state
        if (present(z_bed_background)) then
            state%z_bed = z_bed_background
        else
            state%z_bed = 0.0
        end if
        state%dwdt              = 0.0
        state%w                 = 0.0
        state%w_equilibrium     = 0.0
        state%we                = 0.0

        state%H_above_bsl       = 0.0
        state%Haf               = 0.0
        state%Hice              = 0.0

        state%rsl               = 0.0 ! -1e6_wp
        state%z_ss              = 0.0
        state%dz_ss             = 0.0
        state%canom_load        = 0.0
        state%canom_full        = 0.0
        state%mass_anom         = 0.0

        state%maskocean         = .false.
        state%maskgrounded      = .false.
        state%maskcontinent     = .false.

    end subroutine zero_init_state

    subroutine isos_init_ref(isos, z_bed, H_ice, bsl, dz_ss)
        ! Set reference state from external fields with [nxo,nyo] dimensions

        implicit none

        type(isos_class), intent(INOUT) :: isos
        real(wp), intent(IN) :: z_bed(:, :)                 ! [m] Bedrock elevation
        real(wp), intent(IN) :: H_ice(:, :)                 ! [m] Ice thickness
        real(wp), intent(IN), optional :: bsl               ! [a] Barystatic sea level
        real(wp), intent(IN), optional :: dz_ss(:, :)       ! [m] Sea surface perturbation
        
        isos%ref%z_bed = 0.0
        isos%ref%Hice  = 0.0

        call out2in(isos%ref%z_bed, z_bed, isos%domain)
        call out2in(isos%ref%Hice,  H_ice, isos%domain)

        if (present(bsl))   isos%ref%bsl   = bsl
        if (present(dz_ss)) isos%ref%dz_ss = dz_ss
        ! No need to set w, we, dz_ss and bsl to 0 because `zero_init_state(isos%ref)` in isos_init()
        
        ! Now the reference state has been defined
        isos%par%ref_was_set = .TRUE. 

        return

    end subroutine isos_init_ref

    subroutine isos_init_state(isos, z_bed, H_ice, time, w, we, bsl, dz_ss)

        implicit none
        
        type(isos_class), intent(INOUT) :: isos
        real(wp), intent(IN) :: time                        ! [a] Initial time
        real(wp), intent(IN) :: z_bed(:, :)                 ! [m] Current bedrock elevation
        real(wp), intent(IN) :: H_ice(:, :)                 ! [m] Current ice thickness

        real(wp), intent(IN), optional :: bsl               ! [a] Barystatic sea level
        real(wp), intent(IN), optional :: dz_ss(:, :)       ! [m] Sea surface perturbation
        real(wp), intent(IN), optional :: w(:, :)           ! [m] Viscous displacement
        real(wp), intent(IN), optional :: we(:, :)          ! [m] Elastic displacement

        ! write(*,*) "BSL ref: ", isos%ref%bsl

        if (isos%par%use_restart) then
            ! write(*,*) "Reading restart file..."
            call isos_restart_read(isos,isos%par%restart,time)

            ! Reference state has been set via the restart file
            isos%par%ref_was_set = .TRUE. 

        else
            call out2in(isos%now%z_bed, z_bed, isos%domain)
            call out2in(isos%now%Hice,  H_ice, isos%domain)

            if (present(bsl)) then
                isos%now%bsl = bsl
            else
                isos%now%bsl   = 0.0
            end if
            if (present(dz_ss)) then
                call out2in(isos%now%dz_ss, dz_ss, isos%domain)
            else
                isos%now%dz_ss = 0.0
            end if

            if (.not. isos%par%ref_was_set) then
                call isos_init_ref(isos, z_bed, H_ice, isos%now%bsl, isos%now%dz_ss)
                write(*,*) "isos_init_state:: reference state was set to initial state."
            end if

        end if

        call calc_H_above_bsl(isos%ref, isos%par)
        call calc_rsl(isos%ref)

        ! Store initial bedrock field
        call extendice2isostasy(isos%now, z_bed, H_ice, isos%domain)
        ! call copy_state(isos%ref, isos%now)

        ! Define initial time of isostasy model
        ! (set time_diagnostics earlier, so that it is definitely updated on the first timestep)
        isos%par%time_prognostics = time
        isos%par%time_diagnostics = time

        ! write(*,*) "Calling first mask update..."
        ! call calc_sl_contributions(isos)
        call calc_masks(isos)
        ! write(*,*) "Calling first update..."
        ! call copy_sparsestate(isos%ref, isos%now)
        call isos_update(isos, H_ice, time)

        ! write(*,*) "isos_init_state: "
        ! write(*,*) "    Initial time:   ",  isos%par%time_prognostics
        ! write(*,*) "    range(He_lith): ",  minval(isos%domain%He_lith), maxval(isos%domain%He_lith)
        ! write(*,*) "    range(tau):     ",  minval(isos%domain%tau), maxval(isos%domain%tau)
        ! write(*,*) "    range(w):       ",  minval(isos%now%w), maxval(isos%now%w)
        ! write(*,*) "    range(w_eq):    ",  minval(isos%now%w_equilibrium), maxval(isos%now%w_equilibrium)
        ! write(*,*) "  range(z_bed):     ",  minval(isos%now%z_bed), maxval(isos%now%z_bed)

        if ((minval(isos%domain%tau) .le. 0.0) .and. (isos%par%method .le. 2)) then
            write(error_unit,*) "isos_init_state:: Error: tau initialized with zero values present. &
            &This will lead to the model crashing."
        end if

        ! write(*,*) "isos_init_state:: complete."
        
        return

    end subroutine isos_init_state

    subroutine isos_update(isos, H_ice, time, bsl, dwdt_corr) 

        implicit none

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN)            :: H_ice(:, :)      ! [m] Current ice thickness
        real(wp), intent(IN)            :: time             ! [a] Current time
        real(wp), intent(IN), optional  :: bsl              ! [m] Barystatic sea level forcing
        real(wp), intent(IN), optional  :: dwdt_corr(:, :) ! [m/yr] Basal topography displacement rate (ie, to relax from low resolution to high resolution)

        ! Local variables
        real(wp) :: dt, dt_now
        integer  :: n, nstep
        logical  :: update_diagnostics

        real(wp), allocatable :: dwdt_corr_ext(:,:)

        ! write(*,*) "isos_update:: updating ice thickness..."
        isos%now%Hice = 0.0
        call out2in(isos%now%Hice, H_ice, isos%domain)
        ! write(*,*) "Extrema of H_ice: ", minval(H_ice), maxval(H_ice)

        ! write(*,*) "isos_update:: updating correction..."
        allocate(dwdt_corr_ext(isos%domain%nx,isos%domain%ny))
        dwdt_corr_ext = 1.0
        if (present(dwdt_corr)) then
            call out2in(dwdt_corr_ext, dwdt_corr, isos%domain)
        end if
        ! write(*,*) "Extrema of dwdt_corr: ", minval(dwdt_corr_ext), maxval(dwdt_corr_ext)

        ! Step 0: determine current timestep and number of iterations
        dt = time - isos%par%time_prognostics

        ! Get maximum number of iterations needed to reach desired time
        nstep = ceiling( dt / isos%par%dt_prognostics )
        nstep = max(nstep, 1)

        ! Loop over iterations until maximum time is reached
        do n = 1, nstep

            ! Get current dt (either total time or maximum allowed timestep)
            dt_now = min(dt, isos%par%dt_prognostics)
            ! write(*,*) time

            ! write(*,*) "isos_update:: updating diagnostic bool..."
            ! Only update diagnostics if enough time has passed (to save on computations)
            if ( (isos%par%time_prognostics+dt_now) - isos%par%time_diagnostics .ge. 1e-3) then
                update_diagnostics = .TRUE.
                isos%par%time_diagnostics = isos%par%time_diagnostics + isos%par%dt_diagnostics
            else
                update_diagnostics = .FALSE.
            end if

            ! write (*,*) "time_prog, time_diag, update_diag: ", &
            !     isos%par%time_prognostics, isos%par%time_diagnostics, update_diagnostics

            call nan_check(isos,1)
            
            ! write(*,*) "isos_update:: updating load anomalies..."
            call calc_columnanoms_load(isos)
            ! write(*,*) "Extrema of canom_load: ", minval(isos%now%canom_load), maxval(isos%now%canom_load)

            call nan_check(isos,2)

            if (update_diagnostics) then
                ! write(*,*) "Updating the elastic response..."
                call precomputed_fftconvolution(isos%now%we, isos%domain%FGE, &
                    isos%now%canom_load * isos%par%g * isos%domain%K ** 2.0, &
                    isos%domain%i1, isos%domain%i2, &
                    isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                    isos%domain%nx, isos%domain%ny, &
                    isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)
                ! write(*,*) "Extrema of we: ", minval(isos%now%we), maxval(isos%now%we)
            endif

            call nan_check(isos,3)

            call calc_columnanoms_solidearth(isos)

            if (update_diagnostics .and. isos%par%interactive_sealevel) then
                ! write(*,*) "Updating the SSH perturbation..."
                call calc_mass_anom(isos)
                call precomputed_fftconvolution(isos%now%dz_ss, isos%domain%FGN, &
                    isos%now%mass_anom, isos%domain%i1, isos%domain%i2, &
                    isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                    isos%domain%nx, isos%domain%ny, &
                    isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)
                ! write(*,*) "Extrema of dz_ss: ", minval(isos%now%dz_ss), maxval(isos%now%dz_ss)

                call nan_check(isos,4)

                ! write(*,*) "Updating the SSH..."
                call calc_z_ss(isos%now%z_ss, isos%now%bsl, isos%ref%z_ss, isos%now%dz_ss)
                ! write(*,*) "Extrema of z_ss: ", minval(isos%now%z_ss), maxval(isos%now%z_ss)

                ! write(*,*) "Updating the BSL..."
                if (present(bsl)) then
                    isos%now%bsl = bsl
                else
                    call calc_H_above_bsl(isos%now, isos%par)
                    call calc_sl_contributions(isos)

                    if (isos%par%variable_ocean_surface) then
                        call calc_bsl_variable_Aocean(isos)
                    else
                        call calc_bsl_constant_Aocean(isos)
                    end if
                    ! write(*,*) "BSL: ", isos%now%bsl
                end if
            end if

            ! write(*,*) "Updating RSL..."
            call calc_rsl(isos%now)
            ! write(*,*) "Extrema of rsl: ", minval(isos%now%rsl), maxval(isos%now%rsl)

            call nan_check(isos,5)

            ! write(*,*) "Updating the masks..."
            call calc_masks(isos)

            ! Need to re-compute the load after updating the elastic response, sea level and masks
            call calc_columnanoms_load(isos)
            ! write(*,*) "Extrema of canom_load: ", minval(isos%now%canom_load), maxval(isos%now%canom_load)
            call calc_columnanoms_solidearth(isos)
            ! write(*,*) "Extrema of canom_full: ", minval(isos%now%canom_full), maxval(isos%now%canom_full)

            ! Step 1: diagnose equilibrium displacement and rate of bedrock uplift
            select case(isos%par%method)

                ! Case 0 turns off vertical displacement
                case(0)

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
                        ! call add_columnanom(isos%par%rho_litho, isos%now%we, &
                        !     isos%ref%we, isos%now%canom_load, isos%domain%maskactive)

                        call precomputed_fftconvolution(isos%now%w_equilibrium, isos%domain%FGV, &
                            -isos%now%canom_load * isos%par%g * isos%domain%K ** 2.0, &
                            isos%domain%i1, isos%domain%i2, &
                            isos%domain%j1, isos%domain%j2, isos%domain%offset, &
                            isos%domain%nx, isos%domain%ny, &
                            isos%domain%forward_dftplan_r2c, isos%domain%backward_dftplan_c2r)

                        ! write(*,*) "Checksums: w_eq, FGV, canom_load", &
                        !     sum(isos%now%w_equilibrium), sum(isos%domain%FGV), sum(isos%now%canom_load)

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
                ! write(*,*) "Extrema of dwdt: ", minval(isos%now%dwdt), maxval(isos%now%dwdt)

            ! if (update_diagnostics) then 
            !     ! Update current time_diagnostics value
            !     isos%par%time_diagnostics = isos%par%time_prognostics + dt_now
            ! end if

            ! Step 2: update bedrock elevation and current model time
            ! write(*,*) dt_now

            call nan_check(isos,6)

            if (dt_now .gt. 0.0) then

                !# TODO: do we need this? If yes, lines below should be adapted adequatly
                ! Additionally apply bedrock adjustment field (zero if dwdt_corr not provided as argument)
                ! isos%now%z_bed = isos%now%z_bed + dwdt_corr_ext*dt_now

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

        call cropstate2out(isos%out, isos%now, isos%domain)
        ! write(*,*) "extrema of w: ", minval(isos%now%w), maxval(isos%now%w)
        ! write(*,*) "extrema of we: ", minval(isos%now%we), maxval(isos%now%we)

        ! ajr, diagnostics
        if (ieee_is_nan(maxval(isos%now%z_bed)) .or. ieee_is_nan(maxval(isos%out%z_bed))) then
            write(error_unit,*) "isos_update:: Error: NaNs detected."
            write(error_unit,*) "now%z_bed: ", minval(isos%now%z_bed), maxval(isos%now%z_bed)
            write(error_unit,*) "out%z_bed: ", minval(isos%out%z_bed), maxval(isos%out%z_bed)
            stop
        end if 

        return
    end subroutine isos_update

    subroutine nan_check(isos,step)
        implicit none
        type(isos_class), intent(IN) :: isos
        integer, intent(IN) :: step
        write(error_unit,"(a12,i4,10g10.3)") "nan-check: ", step, maxval(isos%now%z_bed), &
            maxval(isos%now%we),maxval(isos%now%canom_load * isos%par%g * isos%domain%K ** 2.0), &
            maxval(isos%now%dz_ss), maxval(isos%now%mass_anom), maxval(isos%now%dwdt)
        return
    end subroutine nan_check

    subroutine isos_end(isos)

        implicit none 
        type(isos_class), intent(INOUT) :: isos 

        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call deallocate_isos_domain(isos%domain)
        call deallocate_isos_out(isos%out)

        return
    end subroutine isos_end

    subroutine isos_par_load(par, filename, group)

        use nml

        implicit none

        type(isos_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename 
        character(len=*),       intent(IN)  :: group 

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

        call nml_read(filename,group,"interactive_sealevel",par%interactive_sealevel)
        call nml_read(filename,group,"correct_distortion",  par%correct_distortion)
        call nml_read(filename,group,"method",              par%method)
        call nml_read(filename,group,"dt_diagnostics",      par%dt_diagnostics)
        call nml_read(filename,group,"dt_prognostics",      par%dt_prognostics)
        call nml_read(filename,group,"pad",                 par%min_pad)

        call nml_read(filename,group,"nl",          par%nl)
        allocate(par%viscosities(par%nl))
        allocate(par%zl(par%nl))
        call nml_read(filename,group,"zl",          par%zl)

        if (par%method .lt. 3) then     ! ELRA and LLRA
            call nml_read(filename,group,"tau",             par%tau)

        else                            ! ELVA
            call nml_read(filename,group,"mantle",          par%mantle) 
            call nml_read(filename,group,"lithosphere",     par%lithosphere)
            call nml_read(filename,group,"layering",        par%layering)
            call nml_read(filename,group,"viscosities",     par%viscosities)
            
            if (par%layering .eq. "parallel") then
                call nml_read(filename,group,"dl",          par%dl)
            end if

            if ((par%mantle .eq. "rheology_file") .or. (par%lithosphere .eq. "rheology_file")) then
                call nml_read(filename,group,"rheology_file", par%rheology_file)
            end if
        end if

        call nml_read(filename, group, "mask_file",           par%mask_file)

        call nml_read(filename,group,"ocean_surface_file",    par%ocean_surface_file)
        if (trim(par%ocean_surface_file) .eq. "None" .or. &
            trim(par%ocean_surface_file) .eq. "none" .or. &
            trim(par%ocean_surface_file) .eq. "no") then 
            par%variable_ocean_surface = .false. 
        else 
            par%variable_ocean_surface = .true.
        end if

        call nml_read(filename,group,"restart",         par%restart)
        if (trim(par%restart) .eq. "None" .or. &
            trim(par%restart) .eq. "none" .or. &
            trim(par%restart) .eq. "no") then 
            par%use_restart = .FALSE. 
        else 
            par%use_restart = .TRUE.
        end if

        return
    end subroutine isos_par_load

    subroutine allocate_isos(isos, nx_out, ny_out, nl, nbsl)
        implicit none
        type(isos_class), intent(INOUT) :: isos
        integer, intent(IN)             :: nx_out, ny_out, nl, nbsl

        ! First ensure arrays are not allocated
        call deallocate_isos_domain(isos%domain)
        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call deallocate_isos_out(isos%out)
        call allocate_isos_domain(isos%domain, isos%domain%nx, isos%domain%ny, nl, nbsl)
        call allocate_isos_state(isos%now, isos%domain%nx, isos%domain%ny)
        call allocate_isos_state(isos%ref, isos%domain%nx, isos%domain%ny)
        call allocate_isos_out(isos%out, nx_out, ny_out)

        return
    end subroutine allocate_isos

    subroutine deallocate_isos_domain(domain)
        implicit none 
        type(isos_domain_class), intent(INOUT) :: domain

        if (allocated(domain%bsl_vec))      deallocate(domain%bsl_vec)
        if (allocated(domain%A_ocean_vec))  deallocate(domain%A_ocean_vec)

        if (allocated(domain%xc))           deallocate(domain%xc)
        if (allocated(domain%yc))           deallocate(domain%yc)
        if (allocated(domain%x))            deallocate(domain%x)
        if (allocated(domain%y))            deallocate(domain%y)

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

        if (allocated(state%dz_ss))       deallocate(state%dz_ss)

        if (allocated(state%H_above_bsl))       deallocate(state%H_above_bsl)
        if (allocated(state%Haf))               deallocate(state%Haf)
        if (allocated(state%Hice))              deallocate(state%Hice)

        if (allocated(state%rsl))               deallocate(state%rsl)
        if (allocated(state%z_ss))               deallocate(state%z_ss)
        if (allocated(state%canom_load))        deallocate(state%canom_load)
        if (allocated(state%canom_full))        deallocate(state%canom_full)
        if (allocated(state%mass_anom))         deallocate(state%mass_anom)

        if (allocated(state%maskocean))         deallocate(state%maskocean)
        if (allocated(state%maskgrounded))      deallocate(state%maskgrounded)
        if (allocated(state%maskcontinent))     deallocate(state%maskcontinent)

        return 
    end subroutine deallocate_isos_state

    subroutine deallocate_isos_out(out)
        implicit none 
        type(isos_out_class), intent(INOUT) :: out

        if (allocated(out%He_lith))          deallocate(out%He_lith)
        if (allocated(out%D_lith))           deallocate(out%D_lith)
        if (allocated(out%eta_eff))          deallocate(out%eta_eff)
        if (allocated(out%tau))              deallocate(out%tau)
        if (allocated(out%kappa))            deallocate(out%kappa)
        if (allocated(out%kei))              deallocate(out%kei)
        if (allocated(out%GE))               deallocate(out%GE)
        if (allocated(out%GV))               deallocate(out%GV)
        if (allocated(out%GN))               deallocate(out%GN)

        if (allocated(out%Hice))             deallocate(out%Hice)
        if (allocated(out%canom_full))       deallocate(out%canom_full)
        if (allocated(out%dwdt))            deallocate(out%dwdt)

        if (allocated(out%w))                deallocate(out%w)
        if (allocated(out%we))               deallocate(out%we)
        if (allocated(out%w_equilibrium))    deallocate(out%w_equilibrium)

        if (allocated(out%rsl))              deallocate(out%rsl)
        if (allocated(out%z_ss))              deallocate(out%z_ss)
        if (allocated(out%dz_ss))      deallocate(out%dz_ss)
        if (allocated(out%z_bed))            deallocate(out%z_bed)

        if (allocated(out%maskocean))        deallocate(out%maskocean)
        if (allocated(out%maskgrounded))     deallocate(out%maskgrounded)
        if (allocated(out%maskcontinent))    deallocate(out%maskcontinent)
        return
    end subroutine deallocate_isos_out

    subroutine allocate_isos_domain(domain, nx, ny, nl, nbsl)
        implicit none
        integer, intent(IN) :: nx, ny, nl, nbsl
        type(isos_domain_class), intent(INOUT) :: domain

        allocate(domain%bsl_vec(nbsl))
        allocate(domain%A_ocean_vec(nbsl))

        allocate(domain%xc(nx))
        allocate(domain%yc(ny))
        allocate(domain%x(nx, ny))
        allocate(domain%y(nx, ny))

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

        allocate(state%H_above_bsl(nx, ny))
        allocate(state%Haf(nx, ny))
        allocate(state%Hice(nx, ny))

        allocate(state%rsl(nx, ny))
        allocate(state%z_ss(nx, ny))
        allocate(state%dz_ss(nx, ny))
        allocate(state%canom_load(nx, ny))
        allocate(state%canom_full(nx, ny))
        allocate(state%mass_anom(nx, ny))

        allocate(state%maskocean(nx, ny))
        allocate(state%maskgrounded(nx, ny))
        allocate(state%maskcontinent(nx, ny))

        return
    end subroutine allocate_isos_state

    subroutine allocate_isos_out(out, nx, ny)
        implicit none 
        type(isos_out_class), intent(INOUT) :: out
        integer, intent(IN) :: nx, ny

        allocate(out%He_lith(nx, ny))
        allocate(out%D_lith(nx, ny))
        allocate(out%eta_eff(nx, ny))
        allocate(out%tau(nx, ny))
        allocate(out%kappa(nx, ny))
        allocate(out%kei(nx, ny))
        allocate(out%GE(nx, ny))
        allocate(out%GV(nx, ny))
        allocate(out%GN(nx, ny))

        allocate(out%Hice(nx, ny))
        allocate(out%canom_full(nx, ny))
        allocate(out%dwdt(nx, ny))

        allocate(out%w(nx, ny))
        allocate(out%we(nx, ny))
        allocate(out%w_equilibrium(nx, ny))

        allocate(out%rsl(nx, ny))
        allocate(out%z_ss(nx, ny))
        allocate(out%dz_ss(nx, ny))
        allocate(out%z_bed(nx, ny))

        allocate(out%maskocean(nx, ny))
        allocate(out%maskgrounded(nx, ny))
        allocate(out%maskcontinent(nx, ny))
        return
    end subroutine allocate_isos_out

end module fastisostasy