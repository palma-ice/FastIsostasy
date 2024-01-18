
module isostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    ! Note: currently routine `isos_par_load` has dependency on the nml.f90 module 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use isostasy_defs, only : wp, pi, isos_param_class, isos_state_class, isos_domain_class, isos_class
    use green_functions
    use kelvin_function
    use isos_utils

    use solver_lv_elva

    use ncio

    implicit none 

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
        A_ocean_pd, visc_method, visc_c, thck_c, n_lev, rigidity_method)

        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        type(isos_class), intent(OUT) :: isos
        character(len=*), intent(IN)  :: filename
        character(len=*), intent(IN)  :: group
        integer,  intent(IN)  :: nx
        integer,  intent(IN)  :: ny
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        real(wp), intent(IN), optional :: E
        real(wp), intent(IN), optional :: nu
        real(wp), intent(IN), optional :: rho_water
        real(wp), intent(IN), optional :: rho_ice
        real(wp), intent(IN), optional :: rho_seawater
        real(wp), intent(IN), optional :: rho_uppermantle
        real(wp), intent(IN), optional :: rho_litho
        real(wp), intent(IN), optional :: g 
        real(wp), intent(IN), optional :: r_earth 
        real(wp), intent(IN), optional :: m_earth
        real(wp), intent(IN), optional :: A_ocean_pd

        ! Viscous asthenosphere-related parameters
        character(len=*), intent(IN), optional :: visc_method
        real(wp), intent(IN), optional :: visc_c
        real(wp), intent(IN), optional :: thck_c
        integer, intent(IN), optional  :: n_lev
        
        !Lithosphere effective thickness (and therefore ridigidity)
        character(len=*), intent(IN), optional :: rigidity_method

        ! Local variables
        integer         :: nsq
        real(wp)        :: radius_fac
        real(wp)        :: filter_scaling
        real(wp)        :: D_lith_const
        character*256   :: filename_laty

        ! By default one level in asthenosphere
        isos%par%n_lev = 1.        
        
        ! First, load parameters from parameter file `filename`
        call isos_par_load(isos%par, filename, group)

        if (present(interactive_sealevel))  isos%par%interactive_sealevel   = interactive_sealevel
        if (present(correct_distortion))    isos%par%correct_distortion     = correct_distortion
        if (present(method))                isos%par%method                 = method
        if (present(dt_prognostics))        isos%par%dt_prognostics         = dt_prognostics
        if (present(dt_diagnostics))        isos%par%dt_diagnostics         = dt_diagnostics

        ! Overwrite physical constants with arguments if available
        if (present(E))                 isos%par%E                  = E
        if (present(nu))                isos%par%nu                 = nu
        if (present(rho_water))         isos%par%rho_water          = rho_water
        if (present(rho_ice))           isos%par%rho_ice            = rho_ice
        if (present(rho_seawater))      isos%par%rho_seawater       = rho_seawater
        if (present(rho_uppermantle))   isos%par%rho_uppermantle    = rho_uppermantle
        if (present(rho_litho))         isos%par%rho_litho          = rho_litho
        if (present(g))                 isos%par%g                  = g
        if (present(r_earth))           isos%par%r_earth            = r_earth
        if (present(m_earth))           isos%par%m_earth            = m_earth
        if (present(A_ocean_pd))        isos%par%A_ocean_pd         = A_ocean_pd

        if (present(rigidity_method))   isos%par%rigidity_method    = rigidity_method
        if (present(visc_method))       isos%par%visc_method        = visc_method

        if (present(visc_c))            isos%par%visc_c             = visc_c   
        if (present(thck_c))            isos%par%thck_c             = thck_c
        if (present(n_lev))             isos%par%n_lev              = n_lev

        ! Init scalar fields of domain
        ! nsq = max(nx, ny) 
        ! isos%domain%nsq = nsq
        isos%domain%nx = nx
        isos%domain%ny = ny
        isos%domain%dx = dx
        isos%domain%dy = dy
        isos%domain%mu = 2. * pi / ((nx-1) * dx)    ! 2 * pi / L
        
        call allocate_isos(isos)

        ! Init plans
        isos%domain%forward_fftplan_r2r = fftw_plan_r2r_2d(nx, ny, in_aux, out_aux, &
            FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)
        isos%domain%backward_fftplan_r2r = fftw_plan_r2r_2d(nx, ny, in_aux, out_aux, &
            FFTW_DHT, FFTW_DHT, FFTW_ESTIMATE)

        isos%domain%forward_dftplan_r2c = fftw_plan_dft_r2c_2d(m, n, in_aux, out_aux, 1) ! recheck - shouldnt this be -1 (forward)
        isos%domain%backward_dftplan_c2r = fftw_plan_dft_c2r_2d(m, n, in_aux, out_aux, 1)

        ! Init domain
        if (correct_distortion) then
            ! TODO: implement correct expression for K
            isos%domain%K(:, :) = 1
        else
            isos%domain%K(:, :) = 1
        endif

        isos%domain%dx_matrix = dx * isos%domain%K
        isos%domain%dy_matrix = dy * isos%domain%K
        isos%domain%A = dx_matrix * dy_matrix
        call convenient_calc_convolution_indices(isos%domain)
        call convenient_calc_kappa(isos%domain)

        ! Init parameters
        isos%par%sec_per_year = 3600.0*24.0*365.0               ! [s]
        isos%par%compressibility_correction = 1.5 / (1 + nu)    ! [1]

        select case(isos%par%method)


            case(2) ! 
                ! ELRA is being used, which require a constant
                ! value of L_w, D_Lith and thus He_lith everywhere 
                
                ! Calculate the flexural rigidity based on the 
                ! effective elastic thickness of the lithosphere (He_lith),
                ! Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
                D_lith_const = (isos%par%E*1e9) * (isos%par%He_lith*1e3)**3 / (12.0*(1.0-isos%par%nu**2))

                ! Calculate the flexural length scale
                ! (Coulon et al, 2021, Eq. in text after Eq. 3)
                ! Note: will be on the order of 100km
                isos%par%L_w = (D_lith_const / (isos%par%rho_uppermantle*isos%par%g))**0.25 

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
                isos%par%nr =   int(radius_fac*isos%par%L_w/isos%par%dx)+1
                                
                ! Now, initialize isos variables 
                call isos_allocate(isos%now,isos%par%nx,isos%par%ny,isos%par%n_lev,nr=isos%par%nr)

                ! Calculate the Kelvin function filter 
                call calc_kei_filter_2D(isos%domain%kei,L_w=isos%par%L_w,dx=isos%par%dx,dy=isos%par%dx)

                ! Apply scaling to adjust magnitude of filter for radius of filter cutoff
                isos%domain%kei = filter_scaling*isos%domain%kei

                ! Calculate the reference Green's function values
                call calc_greens_function_scaling(isos%domain%G0,isos%domain%kei,isos%par%L_w, &
                                                    D_lith_const,dx=isos%par%dx,dy=isos%par%dx)

                ! Populate 2D D_lith field to have it available
                isos%now%D_lith = D_lith_const
                isos%now%eta_eff    = isos%par%visc                                     
                
            case(4)

                ! LV-ELVA method is being used, which allows for a non-constant
                ! value of L_w, D_Lith and thus He_lith.
                ! ELVA is a particular case of LV-ELVA with constant values
                ! value of L_w, D_Lith and thus He_lith everywhere.
                ! Calculate the flexural rigidity based on the 
                ! effective elastic thickness of the lithosphere (He_lith),
                ! Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
                D_lith_const = (isos%par%E*1e9) * (isos%par%He_lith*1e3)**3 / (12.0*(1.0-isos%par%nu**2))

                ! Calculate the flexural length scale
                ! (Coulon et al, 2021, Eq. in text after Eq. 3)
                ! Note: will be on the order of 100km
                isos%par%L_w = (D_lith_const / (isos%par%rho_uppermantle*isos%par%g))**0.25 

                ! Modify the relative radius to use for the regional filter
                ! depending on whether we want to include the forbulge at
                ! large radii. Large radius makes the code run slower though too. 
                ! Set filter_scaling to < 1.0 to adjust for values of 
                ! radius_fac <~ 5.0 

                ! Larger radius, no filter scaling needed
                ! radius_fac      = 6.0 
                ! filter_scaling  = 1.0 

                ! Smaller radius, filter scaling to improve accuracy
                radius_fac      = 10.   !mmr recheck (larger so nr is also larger, higher resolution for elastic Green function) - 4.0 
                filter_scaling  = 0.91  !mmr recheck

                ! Calculate radius of grid points to use for regional elastic plate filter
                ! See Farrell(1972) for methodology. It seems that this value
                ! should be 5-6x radius of relative stiffness to capture the forebuldge
                ! further out from the depression near the center. 
                ! Note: previous implementation stopped at 400km, hard coded.

                isos%par%nr =  100 !100 !23 ! spare int(radius_fac*isos%par%L_w/isos%par%dx)+1 !mmr spare 100
                
                call calc_ge_filter_2D(isos%domain%GF, dx=isos%par%dx, dy=isos%par%dx) 
             
                if (isos%par%interactive_sealevel) then
                    ! Calculate the geoid's Green function (Coulon et al. 2021) to determine geoid displacement
                    call calc_gn_filter_2D(isos%domain%GN, isos%par%m_earth, &
                        isos%par%r_earth, dx=isos%par%dx, dy=isos%par%dx)
                endif
              
                ! Populate 2D D_lith field to have it available
                
                isos%now%D_lith = D_lith_const

                ! Select asthenosphere viscosity method
                
                select case(trim(isos%par%visc_method))

                case("uniform")

                   isos%now%eta_eff    = isos%par%visc
                   
                case("gaussian_plus")

                   call calc_gaussian_viscosity(isos%now%eta_eff,isos%par%visc,+1._wp,dx=isos%par%dx,dy=isos%par%dx)
                   
                case("gaussian_minus")

                   call calc_gaussian_viscosity(isos%now%eta_eff,isos%par%visc,-1._wp,dx=isos%par%dx,dy=isos%par%dx)
                                      
                case("viscous_channel")

                   isos%now%eta_eff    = isos%par%visc

                   call calc_effective_viscosity_3layer_channel(isos%now%eta_eff,isos%par%visc_c,isos%par%thck_c,isos%par%He_lith,&
                        isos%par%n_lev,isos%par%nu,isos%par%dx,isos%par%dx)

                case("laty")


                   filename_laty = "/Users/montoya/work/ice_data/Antarctica/ANT-32KM/ANT-32KM_latyparams.nc"

                   call nc_read(filename_laty,"log10_mantle_visc",isos%now%eta,start=[1,1,1],count=[isos%par%nx,isos%par%ny,isos%par%n_lev])

                   isos%now%eta = 10.**(isos%now%eta)
                   
                   call calc_effective_viscosity_3d(isos%now%eta_eff,isos%now%eta,isos%par%nu,isos%par%dx,isos%par%dx)
                                   
                case DEFAULT

                   ! do nothing, eta was set above
                   print*,'what is the default?' 
                   stop

                end select

                isos%now%eta_eff    = isos%now%eta_eff * isos%params%compressibility_correction

                !Select lithosphere rigidity method
                
                select case(trim(isos%par%rigidity_method))

                case("uniform")

                   isos%now%D_lith   = D_lith_const         ! [Pa s]
                   isos%now%He_lith  = isos%par%He_lith     ! [km]

                case("gaussian_plus")

                   isos%par%He_lith = 100.0   !km

                   call calc_gaussian_rigidity(isos%now%He_lith,1.5*isos%par%He_lith, isos%par%He_lith, sign=1._wp,dx=isos%par%dx,dy=isos%par%dx) 

                   isos%now%D_lith = (isos%par%E*1.e9) * (isos%now%He_lith*1.e3)**3 / (12.0*(1.0-isos%par%nu**2))

                case("gaussian_minus")

                   isos%par%He_lith = 100.0   !km

                   call calc_gaussian_rigidity(isos%now%He_lith,1.5*isos%par%He_lith, isos%par%He_lith, sign=-1._wp,dx=isos%par%dx,dy=isos%par%dx) 

                   isos%now%D_lith = (isos%par%E*1.e9) * (isos%now%He_lith*1.e3)**3 / (12.0*(1.0-isos%par%nu**2))

                case("laty")

                   filename_laty = "/Users/montoya/work/ice_data/Antarctica/ANT-32KM/ANT-32KM_latyparams.nc"

                   call nc_read(filename_laty,"lithos_thck",isos%now%He_lith,start=[1,1],count=[isos%par%nx,isos%par%ny])

                   isos%now%He_lith = isos%now%He_lith*1.e-3  ! * 0.1  ! recheck for stability!!!
                   isos%now%D_lith = (isos%par%E*1e9) * (isos%now%He_lith*1e3)**3 / (12.0*(1.0-isos%par%nu**2))
   
                case DEFAULT

                   ! do nothing, eta was set above
                   print*,'what is the default?' 
                   stop

                end select

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
                
                write(*,*) "    range(eta_eff):  ", minval(isos%now%eta_eff),     maxval(isos%now%eta_eff) 
                write(*,*) "    range(He_lith):  ", minval(isos%now%He_lith),     maxval(isos%now%He_lith)
                
                write(*,*) "    L_w (m):      ", isos%par%L_w 
                write(*,*) "    nr:           ", isos%par%nr
                write(*,*) "    nx:           ", isos%par%nx
                write(*,*) "    ny:           ", isos%par%ny
                write(*,*) "    dx (m):       ", isos%par%dx 

                write(*,*) "    range(kei): ", minval(isos%domain%kei),    maxval(isos%domain%kei)
                write(*,*) "    range(G0):  ", minval(isos%domain%G0),     maxval(isos%domain%G0)
                write(*,*) "    range(GF):  ", minval(isos%domain%GF),     maxval(isos%domain%GF)
                write(*,*) "    range(GN):  ", minval(isos%domain%GN),     maxval(isos%domain%GN)

            case DEFAULT 

                ! Set elastic length scale to zero (not used)
                isos%par%L_w = 0.0 

                ! Set a small number of points for radius to avoid allocating
                ! any large, unused arrays
                isos%par%nr = 1 
                
                ! Set regional filter fields to zero (not used)
                isos%domain%kei = 0.0 
                isos%domain%G0  = 0.0
                isos%domain%GF  = 0.0
                isos%domain%GN  = 0.0
                
            end select

        ! Set He_lith and tau to constant fields initially.
        ! If these will be spatially varying fields, they should
        ! be defined after calling `isos_init` and before calling
        ! `isos_init_state`.
        !mmr moved above        isos%now%He_lith    = isos%par%He_lith      ! [m]
             
        isos%now%tau        = isos%par%tau          ! [yr]

        ! Set effective viscosity to a constant field too
        ! Later this can be overwritten.
        !mmr moved above        isos%now%eta_eff    = isos%par%visc         ! [Pa s]

        ! TODO: this should 
        ! Intially ensure all variables are zero 
        isos%ref%z_bed      = 0.0
        isos%now%z_bed      = 0.0 
        isos%now%dzbdt      = 0.0 

        isos%now%w              = 0.0   
        isos%now%w_equilibrium  = 0.0      
        isos%now%ssh_perturb    = 0.0 

        ! Set time to very large value in the future 
        isos%par%time_diagnostics = 1e10 
        isos%par%time_prognostics = 1e10

        
        write(*,*) "isos_init:: complete." 

        
        return 

    end subroutine isos_init

    subroutine isos_init_state(isos, z_bed, H_ice, z_sl, z_bed_ref, H_ice_ref, z_sl_ref, time)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: z_bed(:,:)            ! [m] Current bedrock elevation 
        real(wp), intent(IN) :: H_ice(:,:)            ! [m] Current ice thickness  
        real(wp), intent(IN) :: z_sl(:,:)             ! [m] Current sea level 
        real(wp), intent(IN) :: z_bed_ref(:,:)        ! [m] Reference bedrock elevation (with known load)
        real(wp), intent(IN) :: H_ice_ref(:,:)        ! [m] Reference ice thickness (associated with reference z_bed)
        real(wp), intent(IN) :: z_sl_ref(:,:)         ! [m] Reference sea level (associated with reference z_bed)
        real(wp), intent(IN) :: time                  ! [a] Initial time 

        integer :: nx, ny

!        character*256 filename
        
        ! Store initial bedrock field 
        isos%now%z_bed = z_bed 
        
        ! Store reference bedrock field
        isos%ref%z_bed = z_bed_ref 
        isos%now%dzbdt     = 0.0

        nx = size(H_ice,1)
        ny = size(H_ice,2)

        ! Define initial time of isostasy model 
        ! (set time_diagnostics earlier, so that it is definitely updated on the first timestep)
        isos%par%time_diagnostics = time - isos%par%dt_diagnostics
        isos%par%time_prognostics = time 

        ! Call isos_update to diagnose rate of change
        ! (no change to z_bed will be applied since isos%par%time==time)

        call isos_update(isos,H_ice,z_sl,time) 
            
        write(*,*) "isos_init_state:: "
        write(*,*) "  Initial time:   ", isos%par%time_prognostics 
        write(*,*) "  range(He_lith): ", minval(isos%now%He_lith), maxval(isos%now%He_lith)
        write(*,*) "  range(tau):     ", minval(isos%now%tau),     maxval(isos%now%tau)
        write(*,*) "  range(w):      ", minval(isos%now%w),      maxval(isos%now%w)
        write(*,*) "  range(w_equilibrium):      ", minval(isos%now%w_equilibrium),      maxval(isos%now%w_equilibrium)
        write(*,*) "  range(z_bed):   ", minval(isos%now%z_bed),   maxval(isos%now%z_bed)

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

    subroutine isos_update(isos, H_ice, time, dzbdt_corr) 

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: H_ice(:,:)                  ! [m] Current ice thickness 
        real(wp), intent(IN) :: time                        ! [a] Current time  
        real(wp), intent(IN), optional :: dzbdt_corr(:,:)   ! [m/yr] Basal topography adjustment rate (ie, to relax from low resolution to high resolution) 

        ! Local variables
        
        real(wp) :: dt, dt_now  
        integer  :: n, nstep 
        logical  :: update_equil
 
        ! Step 0: determine current timestep and number of iterations
        dt = time - isos%par%time_prognostics 

        ! Get maximum number of iterations needed to reach desired time
        nstep = ceiling( (time - isos%par%time_prognostics) / isos%par%dt_prognostics )
        nstep = max(nstep,1)

        ! Loop over iterations until maximum time is reached
        do n = 1, nstep 

            ! Get current dt (either total time or maximum allowed timestep)
            dt_now = min(time-isos%par%time_prognostics, isos%par%dt_prognostics)

            ! Only update equilibrium displacement height if 
            ! enough time has passed (to save on computations)
            if ( (isos%par%time_prognostics+dt_now) - isos%par%time_diagnostics .ge. isos%par%dt_diagnostics) then
                update_equil = .TRUE. 
            else 
                update_equil = .FALSE. 
            end if 

            ! Step 1: diagnose equilibrium displacement and rate of bedrock uplift
            select case(isos%par%method)

                case(0)
                    ! Steady-state lithosphere

                    isos%now%w    = isos%now%w_equilibrium
                    isos%now%dzbdt = 0.0

                case(1)
                    ! Local lithosphere, relaxing asthenosphere (LLRA)

                    ! Local lithosphere (LL)
                    ! (update every time because it is cheap)
                    call calc_litho_local(isos%now%w, isos%now%q1, isos%now%z_bed, H_ice, z_sl, &
                        isos%par%rho_ice, isos%par%rho_seawater, isos%par%rho_uppermantle, isos%par%g)

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt, isos%now%z_bed, isos%ref%z_bed, &
                        w_b=isos%now%w-isos%ref%w_equilibrium, tau=isos%now%tau)

                 case(2)
                    ! Elastic lithosphere, relaxing asthenosphere (ELRA)
                    
                    ! Elastic lithosphere (EL)
                    if (update_equil) then 
                        call calc_litho_regional(isos%now%w, isos%now%q1, isos%now%z_bed, H_ice, &
                            z_sl,isos%domain%G0, isos%par%rho_ice, isos%par%rho_seawater, &
                            isos%par%rho_uppermantle,isos%par%rho_litho,isos%par%g)
                    end if 

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt, isos%now%z_bed, isos%ref%z_bed, &
                        w_b=isos%now%w-isos%ref%w_equilibrium, tau=isos%now%tau)
                    
                case(3) 
                    ! Elementary GIA model (spatially varying ELRA with geoid - to do!)

                    ! 2D elastic lithosphere (2DEL)
                    if (update_equil) then
                       !call calc_litho_regional(isos%now%w,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%domain%G0)
 
                        write(*,*) "isos_update:: Error: method=3 not implemented yet."
                        stop 

                    end if 

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt,isos%now%z_bed,isos%ref%z_bed, &
                                                    w_b=isos%now%w-isos%ref%w_equilibrium,tau=isos%now%tau)

                 case(4) 
                    ! ELVA - viscous half-space asthenosphere
                    ! or  LV-ELVA - laterally-variable viscosity asthenosphere
                    ! overlain by elastic plate lithosphere with uniform constants                                                      
                    
                    ! TODO update load here

                    if (isos%par%interactive_sealevel) then
                        call calc_geoid_displacement(isos%now%ssh_perturb, isos%now%mass_anom, isos%domain%GN)
                    else
                        isos%now%ssh_perturb = 0.
                    endif

                    call calc_lvelva_viscous_square(isos%now%dzbdt, isos%now%w, &
                        isos%now%w,isos%now%q1,isos%par%nu,isos%now%D_lith, isos%now%eta_eff, &
                        isos%par%rho_uppermantle,isos%par%rho_litho,isos%par%g,isos%par%dx,dt_now)

                 end select

            if (update_equil) then 
                ! Update current time_diagnostics value 
                isos%par%time_diagnostics = isos%par%time_prognostics + dt_now 
             end if

            ! Step 2: update bedrock elevation and current model time
            if (dt_now .gt. 0.0) then

               isos%now%w = isos%now%w + isos%now%dzbdt*dt_now
               isos%now%z_bed = isos%now%z_bed + isos%now%dzbdt*dt_now

                ! Additionally apply bedrock adjustment field
                if (present(dzbdt_corr)) then 
                    isos%now%z_bed = isos%now%z_bed + dzbdt_corr*dt_now
                end if 
                
                isos%par%time_prognostics = isos%par%time_prognostics + dt_now
               
            end if

            if ( abs(time-isos%par%time_prognostics) .lt. 1e-5) then 
                ! Desired time has reached, exit the loop 
                isos%par%time_prognostics = time 
                exit 
            end if 

         end do

        return 

    end subroutine isos_update

    subroutine isos_end(isos)

        implicit none 

        type(isos_class), intent(INOUT) :: isos 

        ! Deallocate isos object 
        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call deallocate_isos_domain(isos%domain)

        return 

    end subroutine isos_end

    subroutine isos_par_load(par,filename,group)

        use nml 

        implicit none

        type(isos_param_class), intent(OUT) :: par
        character(len=*),       intent(IN)  :: filename 
        character(len=*),       intent(IN)  :: group 

        call nml_read(filename,group,"interactive_sealevel",          par%interactive_sealevel)
        call nml_read(filename,group,"correct_distortion",          par%correct_distortion)
        call nml_read(filename,group,"method",          par%method)
        call nml_read(filename,group,"dt_diagnostics",  par%dt_diagnostics)
        call nml_read(filename,group,"dt_prognostics",  par%dt_prognostics)
        call nml_read(filename,group,"He_lith",         par%He_lith)
        call nml_read(filename,group,"visc",            par%visc)
        call nml_read(filename,group,"tau",             par%tau)

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

        ! Viscous asthenosphere-related parameters
        
        call nml_read(filename,group,"visc_method",     par%visc_method) 
        call nml_read(filename,group,"visc_c",          par%visc_c)              
        call nml_read(filename,group,"thck_c",          par%thck_c)
        call nml_read(filename,group,"n_lev",           par%n_lev)
        
        ! Lithosphere effective thickness (and therefore ridigidity)
        call nml_read(filename,group,"rigidity_method",  par%rigidity_method)

        return


    end subroutine isos_par_load

    subroutine allocate_isos(isos)
        implicit none 
        type(isos_class), intent(INOUT) :: isos

        ! First ensure arrays are not allocated
        call deallocate_isos_domain(isos%domain)
        call deallocate_isos_state(isos%now)
        call deallocate_isos_state(isos%ref)
        call allocate_isos_domain(isos%domain, isos%domain%nx, isos%domain%ny)
        call allocate_isos_state(isos%now, isos%domain%nx, isos%domain%ny)
        call allocate_isos_state(isos%ref, isos%domain%nx, isos%domain%ny)

        return

    end subroutine allocate_isos

    subroutine deallocate_isos_domain(domain)
        implicit none 
        type(isos_domain_class), intent(INOUT) :: domain

        if (allocated(domain%dx_matrix))    deallocate(domain%dx_matrix)
        if (allocated(domain%dy_matrix))    deallocate(domain%dy_matrix)
        if (allocated(domain%A))            deallocate(domain%A)
        if (allocated(domain%K))            deallocate(domain%K)
        if (allocated(domain%kappa))        deallocate(domain%kappa)

        if (allocated(domain%forward_fftplan_r2r))       deallocate(domain%forward_fftplan_r2r)
        if (allocated(domain%backward_fftplan_r2r))        deallocate(domain%backward_fftplan_r2r)
        if (allocated(domain%forward_dftplan_r2c))        deallocate(domain%forward_dftplan_r2c)
        if (allocated(domain%backward_dftplan_c2r))        deallocate(domain%backward_dftplan_c2r)

        if (allocated(domain%kei))       deallocate(domain%kei)
        if (allocated(domain%G0))        deallocate(domain%G0)
        if (allocated(domain%GF))        deallocate(domain%GF)
        if (allocated(domain%GN))        deallocate(domain%GN)

        return 

    end subroutine deallocate_isos_domain

    subroutine deallocate_isos_state(state)
        implicit none 
        type(isos_state_class), intent(INOUT) :: state

        if (allocated(state%He_lith))           deallocate(state%He_lith)
        if (allocated(state%D_lith))            deallocate(state%D_lith)
        if (allocated(state%eta))               deallocate(state%eta)
        if (allocated(state%eta_eff))               deallocate(state%eta_eff)
        if (allocated(state%tau))               deallocate(state%tau)

        if (allocated(state%z_bed))             deallocate(state%z_bed)
        if (allocated(state%dzbdt))             deallocate(state%dzbdt)
        if (allocated(state%q))                 deallocate(state%q)
        if (allocated(state%w))                 deallocate(state%w)
        if (allocated(state%w_equilibrium))     deallocate(state%w_equilibrium)
        if (allocated(state%we))                deallocate(state%we)

        if (allocated(state%ssh_perturb))       deallocate(state%ssh_perturb)
        if (allocated(state%Haf))               deallocate(state%Haf)
        if (allocated(state%Hice))              deallocate(state%Hice)
        if (allocated(state%Hsw))               deallocate(state%Hsw)

        if (allocated(state%ssh))               deallocate(state%ssh)
        if (allocated(state%canom_load))        deallocate(state%canom_load)
        if (allocated(state%canom_full))        deallocate(state%canom_full)
        if (allocated(state%mass_anom))         deallocate(state%mass_anom)
        if (allocated(state%maskocean))         deallocate(state%maskocean)
        if (allocated(state%maskactive))        deallocate(state%maskactive)
        if (allocated(state%maskgrounded))      deallocate(state%maskgrounded)

        return 

    end subroutine deallocate_isos_state

    subroutine allocate_isos_domain(domain, nx, ny)
        implicit none
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny
        type(isos_state_class), intent(INOUT) :: domain

        allocate(domain%dx_matrix(nx,ny))
        allocate(domain%dy_matrix(nx,ny))
        allocate(domain%A(nx,ny))
        allocate(domain%K(nx,ny))
        allocate(domain%kappa(nx,ny))

        allocate(domain%kei(nx,ny))
        allocate(domain%G0(nx,ny))
        allocate(domain%GF(nx,ny))
        allocate(domain%GN(nx,ny))

        return
    end subroutine allocate_isos_domain

    subroutine allocate_isos_state(state, nx, ny)
        implicit none
        integer, intent(IN) :: nx
        integer, intent(IN) :: ny
        type(isos_state_class), intent(INOUT) :: state

        integer :: nfilt
        integer :: nz
        nfilt = 2 ! dummy number
        nz = 1

        ! Allocate arrays
        allocate(state%He_lith(nx,ny))
        allocate(state%D_lith(nx,ny))
        
        if(nz.gt.1) then
           allocate(state%eta(nx,ny,nz))
        else
           allocate(state%eta(nx,ny,1))
        endif
        
        allocate(state%eta_eff(nx,ny))
        allocate(state%tau(nx,ny))

        allocate(state%z_bed(nx,ny))
        allocate(state%dzbdt(nx,ny))
        allocate(state%q(nx,ny))
        allocate(state%w(nx,ny))
        allocate(state%w_equilibrium(nx,ny))
        allocate(state%we(nx,ny))

        allocate(state%ssh_perturb(nx,ny))
        allocate(state%Haf(nx,ny))
        allocate(state%Hice(nx,ny))
        allocate(state%Hsw(nx,ny))

        allocate(state%ssh(nx,ny))
        allocate(state%canom_load(nx,ny))
        allocate(state%canom_full(nx,ny))
        allocate(state%mass_anom(nx,ny))
        allocate(state%maskocean(nx,ny))
        allocate(state%maskactive(nx,ny))
        allocate(state%maskgrounded(nx,ny))

        return

    end subroutine allocate_isos_state

end module isostasy