
module isostasy 
    ! Isostasy (Greek ísos "equal", stásis "standstill") is the state of 
    ! gravitational equilibrium between Earth's crust and mantle such that 
    ! the crust "floats" at an elevation that depends on its thickness and density.
    ! -- https://en.wikipedia.org/wiki/Isostasy

    ! Note: currently routine `isos_par_load` has dependency on the nml.f90 module 
    
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use isostasy_defs, only : wp, pi, isos_param_class, isos_state_class, isos_domain_class, isos_class
    
    use solver_lv_elva

    use kelvin_function

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
    public :: isos_set_field

    contains
    
    subroutine isos_init(isos, filename, group, nx, ny, dx, static_load, E, nu, &
        rho_water, rho_ice, rho_seawater, rho_uppermantle, rho_litho, g, r_earth, m_earth, &
        A_ocean_pd, visc_method, visc_c, thck_c, n_lev, rigidity_method, calc_geoid)

        implicit none 

        type(isos_class), intent(OUT) :: isos
        character(len=*), intent(IN)  :: filename
        character(len=*), intent(IN)  :: group
        integer,  intent(IN)  :: nx
        integer,  intent(IN)  :: ny
        real(wp), intent(IN)  :: dx

        logical, intent(IN), optional :: static_load

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

        !Geoid
        logical, intent(IN), optional :: calc_geoid
        
        ! Local variables

        real(wp) :: radius_fac
        real(wp) :: filter_scaling
        real(wp) :: D_lith_const

        character*256 :: filename_laty

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

        if (present(visc_method))       isos%par%visc_method        = visc_method

        if (present(visc_c))            isos%par%visc_c             = visc_c   
        if (present(thck_c))            isos%par%thck_c             = thck_c
        if (present(n_lev))             isos%par%n_lev              = n_lev

        if (present(rigidity_method))   isos%par%rigidity_method    = rigidity_method

        ! Allocate
        isos%domain%nx = nx
        isos%domain%ny = ny
        isos%domain%dx = dx
        isos%domain%dy = dy
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

        if ( mod(nx,2).eq.0 ) then
            isos%domain%i1 = nx/2
        else
            isos%domain%i1 = int(nx/2) + 1
        endif
        isos%domain%i2 = 2*nx - 1 - int(nx/2)
 
        if ( mod(ny,2).eq.0 ) then
            isos%domain%j1 = ny/2
        else
            isos%domain%j1 = int(ny/2) + 1
        endif
        isos%domain%j2 = 2*ny - 1 - int(ny/2)

        if ( mod(ny - nx, 2).eq.0 ) then
            isos%domain%convo_offset = (ny - nx) / 2
        else
            isos%domain%convo_offset = (ny - nx - 1) / 2
        endif

        isos%domain%kappa(:, :) = 0.0
        ic = (nx-1)/2 + 1
        jc = (ny-1)/2 + 1

        do i = 1, nx
            if (i.le.ic) then 
                ip = i-1
            else
                ip = nx-i+1
            end if
            do j = 1, ny
                if (j.le.jc) then
                    iq = j-1  
                else
                    iq = ny-j+1
                end if
                isos%domain%kappa(i,j)  = (ip*ip + iq*iq)**0.5
            end do
        end do
        isos%domain%kappa = 2 * pi / L * isos%domain%kappa
        isos%domain%kappa(1,1) = (isos%domain%kappa(1,2) + isos%domain%kappa(2,1)) / 2.0

        ! Init parameters
        isos%params%compressibility_correction = 1.5 / (1 + nu)

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
                call calc_kei_filter_2D(isos%now%kei,L_w=isos%par%L_w,dx=isos%par%dx,dy=isos%par%dx)

                ! Apply scaling to adjust magnitude of filter for radius of filter cutoff
                isos%now%kei = filter_scaling*isos%now%kei

                ! Calculate the reference Green's function values
                call calc_greens_function_scaling(isos%now%G0,isos%now%kei,isos%par%L_w, &
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
                
                call calc_ge_filter_2D(isos%now%GF,dx=isos%par%dx,dy=isos%par%dx) 
             
                if (isos%par%interactive_sealevel) then
                    ! Calculate the geoid's Green function (Coulon et al. 2021) to determine geoid displacement
                    call calc_gn_filter_2D(isos%now%GN,isos%par%m_earth,isos%par%r_earth, dx=isos%par%dx,dy=isos%par%dx)
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

                write(*,*) "    range(kei): ", minval(isos%now%kei),    maxval(isos%now%kei)
                write(*,*) "    range(G0):  ", minval(isos%now%G0),     maxval(isos%now%G0)
                write(*,*) "    range(GF):  ", minval(isos%now%GF),     maxval(isos%now%GF)
                write(*,*) "    range(GN):  ", minval(isos%now%GN),     maxval(isos%now%GN)

            case DEFAULT 

                ! Set elastic length scale to zero (not used)
                isos%par%L_w = 0.0 

                ! Set a small number of points for radius to avoid allocating
                ! any large, unused arrays
                isos%par%nr = 1 

                ! Now, initialize isos variables 
                call isos_allocate(isos%now,isos%par%nx,isos%par%ny,isos%par%n_lev,nr=isos%par%nr)
                
                ! Set regional filter fields to zero (not used)
                isos%now%kei = 0.0 
                isos%now%G0  = 0.0
                isos%now%GF  = 0.0
                isos%now%GN  = 0.0
                
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

        
        ! Intially ensure all variables are zero 
        isos%now%z_bed_ref  = 0.0
        isos%now%z_bed      = 0.0 
        isos%now%dzbdt      = 0.0 

        isos%now%w0         = 0.0   
        isos%now%w1         = 0.0      
        isos%now%w2         = 0.0
        isos%now%ssh_perturb         = 0.0 

        ! Set time to very large value in the future 
        isos%par%time_diagnostics = 1e10 
        isos%par%time_prognostics = 1e10

        
        write(*,*) "isos_init:: complete." 

        
        return 

    end subroutine isos_init

    subroutine isos_init_state(isos,z_bed,H_ice,z_sl,z_bed_ref,H_ice_ref,z_sl_ref,time)

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
        isos%now%z_bed_ref = z_bed_ref 
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
        write(*,*) "  range(w0):      ", minval(isos%now%w0),      maxval(isos%now%w0)
        write(*,*) "  range(w1):      ", minval(isos%now%w1),      maxval(isos%now%w1)
        write(*,*) "  range(w2):      ", minval(isos%now%w2),      maxval(isos%now%w2)
        write(*,*) "  range(z_bed):   ", minval(isos%now%z_bed),   maxval(isos%now%z_bed)

        ! Make sure tau does not contain zero values, if so
        ! output an error for diagnostics. Don't kill the program
        ! so that external program can still write output as needed. 
        if (minval(isos%now%tau) .le. 0.0) then 
            write(error_unit,*) "isos_init_state:: Error: tau initialized with zero values present. &
            &This will lead to the model crashing."
            !stop 
        end if 

!        print*,'hola'
!        stop
        
        return 

    end subroutine isos_init_state

    subroutine isos_update(isos,H_ice,z_sl,time,dzbdt_corr) 

        implicit none 

        type(isos_class), intent(INOUT) :: isos 
        real(wp), intent(IN) :: H_ice(:,:)                  ! [m] Current ice thickness 
        real(wp), intent(IN) :: z_sl(:,:)                   ! [m] Current sea level
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

            ! Only update equilibrium displacement height (w1) if 
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

                    isos%now%w1    = isos%now%w0
                    isos%now%dzbdt = 0.0 

                case(1)
                    ! Local lithosphere, relaxing asthenosphere (LLRA)

                    ! Local lithosphere (LL)
                    ! (update every time because it is cheap)
                    call calc_litho_local(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl, &
                                                    isos%par%rho_ice,isos%par%rho_seawater,isos%par%rho_uppermantle,isos%par%g)

                    ! Relaxing asthenosphere (RA)
                    call calc_asthenosphere_relax(isos%now%dzbdt,isos%now%z_bed,isos%now%z_bed_ref, &
                                                    w_b=isos%now%w1-isos%now%w0,tau=isos%now%tau)

                 case(2)
                    ! Elastic lithosphere, relaxing asthenosphere (ELRA)
                    
                    ! Elastic lithosphere (EL)
                    if (update_equil) then 
                        call calc_litho_regional(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%now%G0, &
                             isos%par%rho_ice,isos%par%rho_seawater,isos%par%rho_uppermantle,isos%par%rho_litho,isos%par%g)
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

                 case(4) 
                    ! ELVA - viscous half-space asthenosphere
                    ! or  LV-ELVA - laterally-variable viscosity asthenosphere
                    ! overlain by                                                   
                    ! elastic plate lithosphere with uniform constants                                                      
                    
                    ! Local lithosphere (LL)
                    ! (update every time because it is cheap)
                    ! Note: calculate the local lithospheric load here because 
                    ! the EL component is contained in the ELVA model solution
                    ! q1 will be used (the load), while w1 will not be used further.

                    call calc_litho_regional(isos%now%w1,isos%now%q1,isos%now%z_bed,H_ice,z_sl,isos%now%GF,
                        isos%par%rho_ice,isos%par%rho_seawater,isos%par%rho_uppermantle,isos%par%rho_litho,
                        isos%par%g)

                    call calc_lv_asthenosphere_viscous_square(isos%now%dzbdt,isos%now%w2,isos%now%w1,isos%now%q1,isos%par%nu,isos%now%D_lith, &
                         isos%now%eta_eff,isos%par%rho_uppermantle,isos%par%rho_litho,isos%par%g,isos%par%dx,dt_now)

                    if (isos%par%interactive_sealevel) then
                       call calc_geoid_displacement(isos%now%ssh_perturb,-isos%now%q1/isos%par%g - isos%now%w2*isos%par%rho_uppermantle - isos%now%w1*isos%par%rho_litho, &
                            isos%now%GN)
                    else
                       isos%now%ssh_perturb = 0.
                    endif

                 end select

            if (update_equil) then 
                ! Update current time_diagnostics value 
                isos%par%time_diagnostics = isos%par%time_prognostics + dt_now 
             end if

            ! Step 2: update bedrock elevation and current model time
            if (dt_now .gt. 0.0) then


               isos%now%w2 = isos%now%w2 + isos%now%dzbdt*dt_now
               
               isos%now%z_bed = isos%now%z_bed + isos%now%dzbdt*dt_now


!               print*,'hola dzbdt',  isos%now%dzbdt(60,60)

!                ! Additionally apply bedrock adjustment field
!               if (present(dzbdt_corr)) then 
!                    isos%now%z_bed = isos%now%z_bed + dzbdt_corr*dt_now
!                end if 


!                isos%now%w2 = isos%now%z_bed 
                
                isos%par%time_prognostics = isos%par%time_prognostics + dt_now
               
            end if

            ! Write a summary to check timestepping
!            write(*,*) "isos: ", n, time, dt_now, isos%par%time_prognostics, isos%par%time_diagnostics, isos%now%z_bed(60,60)


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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TODO: Remove this section
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine isos_allocate(now,nx,ny,nz,nr)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 
        integer, intent(IN) :: nx, ny, nz, nr  

        ! Local variables
        integer :: nfilt 

        nfilt = 2*nr+1 

        ! First ensure arrays are not allocated
        call isos_deallocate(now)

        ! Allocate arrays
        allocate(now%He_lith(nx,ny))
        allocate(now%D_lith(nx,ny))
        
        if(nz.gt.1) then
           allocate(now%eta(nx,ny,nz))
        else
           allocate(now%eta(nx,ny,1))
        endif
        
        allocate(now%eta_eff(nx,ny))
        allocate(now%tau(nx,ny))

        
        allocate(now%kei(nfilt,nfilt))
        allocate(now%G0(nfilt,nfilt))
! recheck convol 
!        allocate(now%GF(nfilt,nfilt))
        allocate(now%GF(nx,ny))
! recheck convol 
!        allocate(now%GN(nfilt,nfilt)) 
        allocate(now%GN(nx,ny))

        allocate(now%z_bed(nx,ny))
        allocate(now%dzbdt(nx,ny))
        allocate(now%z_bed_ref(nx,ny))
        allocate(now%q0(nx,ny))
        allocate(now%w0(nx,ny))
        allocate(now%q1(nx,ny))
        allocate(now%w1(nx,ny))     
        allocate(now%w2(nx,ny))
        allocate(now%ssh_perturb(nx,ny))     

        return 

    end subroutine isos_allocate

    subroutine isos_deallocate(now)

        implicit none 

        type(isos_state_class), intent(INOUT) :: now 

        if (allocated(now%He_lith))     deallocate(now%He_lith)
        if (allocated(now%D_lith))      deallocate(now%D_lith)
        if (allocated(now%eta_eff))     deallocate(now%eta_eff)
        if (allocated(now%eta))         deallocate(now%eta)
        if (allocated(now%tau))         deallocate(now%tau)
        
        if (allocated(now%kei))         deallocate(now%kei)
        if (allocated(now%G0))          deallocate(now%G0)

        if (allocated(now%GF))          deallocate(now%GF)

        if (allocated(now%GN))          deallocate(now%GN)

        if (allocated(now%z_bed))       deallocate(now%z_bed)
        if (allocated(now%dzbdt))       deallocate(now%dzbdt)
        if (allocated(now%z_bed_ref))   deallocate(now%z_bed_ref)
        if (allocated(now%q0))          deallocate(now%q0)
        if (allocated(now%w0))          deallocate(now%w0)
        if (allocated(now%q1))          deallocate(now%q1)
        if (allocated(now%w1))          deallocate(now%w1)    
        if (allocated(now%w2))          deallocate(now%w2)     

        return 

    end subroutine isos_deallocate
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
        return 

    end subroutine deallocate_isos_domain

    subroutine deallocate_isos_state(state)
        implicit none 
        type(isos_state_class), intent(INOUT) :: state

        if (allocated(state%He_lith))   deallocate(state%He_lith)
        if (allocated(state%D_lith))    deallocate(state%D_lith)
        if (allocated(state%eta))       deallocate(state%eta)
        if (allocated(state%tau))       deallocate(state%tau)

        if (allocated(state%kei))       deallocate(state%kei)
        if (allocated(state%G0))        deallocate(state%G0)
        if (allocated(state%GF))        deallocate(state%GF)
        if (allocated(state%GN))        deallocate(state%GN)

        if (allocated(state%z_bed))     deallocate(state%z_bed)
        if (allocated(state%dzbdt))     deallocate(state%dzbdt)
        if (allocated(state%q))         deallocate(state%q)
        if (allocated(state%w))         deallocate(state%w)
        if (allocated(state%w_equilibrium))    deallocate(state%w_equilibrium)

        if (allocated(state%ssh_perturb))    deallocate(state%ssh_perturb)
        if (allocated(state%Haf))       deallocate(state%Haf)
        if (allocated(state%Hice))      deallocate(state%Hice)
        if (allocated(state%Hsw))       deallocate(state%Hsw)
        if (allocated(state%u))         deallocate(state%u)
        if (allocated(state%ue))        deallocate(state%ue)

        if (allocated(state%seasurfaceheight))  deallocate(state%seasurfaceheight)
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

        
        allocate(state%kei(nfilt,nfilt))
        allocate(state%G0(nfilt,nfilt))
! recheck convol 
!        allocate(now%GF(nfilt,nfilt))
        allocate(state%GF(nx,ny))
! recheck convol 
!        allocate(state%GN(nfilt,nfilt)) 
        allocate(state%GN(nx,ny))

        allocate(state%z_bed(nx,ny))
        allocate(state%dzbdt(nx,ny))
        allocate(state%z_bed_ref(nx,ny))
        allocate(state%q(nx,ny))
        allocate(state%w(nx,ny))
        allocate(state%w_equilibrium(nx,ny))

        allocate(state%ssh_perturb(nx,ny))
        allocate(state%Haf(nx,ny))
        allocate(state%Hice(nx,ny))
        allocate(state%Hsw(nx,ny))
        allocate(state%u(nx,ny))
        allocate(state%ue(nx,ny))

        allocate(state%seasurfaceheight(nx,ny))
        allocate(state%canom_load(nx,ny))
        allocate(state%canom_full(nx,ny))
        allocate(state%mass_anom(nx,ny))
        allocate(state%maskocean(nx,ny))
        allocate(state%maskactive(nx,ny))
        allocate(state%maskgrounded(nx,ny))

        return

    end subroutine allocate_isos_state



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


    subroutine calc_greens_function_scaling(G0,kei2D,L_w,D_lith,dx,dy)
        ! The Green's function (Eq. 3 of Coulon et al, 2021)
        ! gives displacement G in [m] as a function of the distance
        ! r from the point load P_b [Pa]. 

        ! Here G0 is calculated, which is G without including the point load.
        ! G0 has units of [m N-1]. 
        ! This can then be multiplied with the actual magnitude of the
        ! point load to obtain G.
        ! G = G0 * P_b = [m N-1] * [Pa] = [m]. 

        ! Note that L_w contains information about rho_uppermantle. 

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

    subroutine calc_ge_filter_2D(filt,dx,dy) 
        ! Calculate 2D Green Function

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
!        real(wp), intent(IN)  :: GE(:,:)
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy   

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp) :: ge_test_0
        real(wp) :: ge_test_1

        real(wp), parameter, dimension(42) :: rn_vals = [ 0.0,    0.011,  0.111,  1.112,  2.224, &
            3.336, 4.448,  6.672, 8.896,  11.12,  17.79,  22.24,  27.80,  33.36,  44.48,  55.60, &
            66.72,  88.96,  111.2,  133.4,  177.9,  222.4,  278.0,  333.6, 444.8,  556.0,  667.2, &
            778.4,  889.6,  1001.0, 1112.0, 1334.0, 1779.0, 2224.0, 2780.0, 3336.0, 4448.0, 5560.0, &
            6672.0, 7784.0, 8896.0, 10008.0]* 1.e3        !    # converted to meters
       
        real(wp), parameter, dimension(42) :: ge_vals = [ -33.6488, -33.64, -33.56, -32.75, -31.86, &
            -30.98, -30.12, -28.44, -26.87, -25.41,-21.80, -20.02, -18.36, -17.18, -15.71, -14.91, &
            -14.41, -13.69, -13.01,-12.31, -10.95, -9.757, -8.519, -7.533, -6.131, -5.237, -4.660, &
            -4.272,-3.999, -3.798, -3.640, -3.392, -2.999, -2.619, -2.103, -1.530, -0.292, 0.848, &
            1.676,  2.083,  2.057,  1.643]


        ! Get size of filter array and half-width
        n  = size(filt,1) 
        n2 = (n-1)/2 
       
        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_ge_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_ge_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop  
        end if 


        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center
        ! and impose correct Green function value. 

        filt = 0.
        
        do j = -n2, n2 
        do i = -n2, n2

            x  = i*dx 
            y  = j*dy 
            r  = sqrt(x**2+y**2)  

            ! Get actual index of array
            i1 = i+1+n2 
            j1 = j+1+n2 

            ! Get correct GE value for this point

            
            filt(i1,j1) = get_ge_value(r,rn_vals,ge_vals) * 1.e-12 *(dx*dy) /(9.81*max(r,dx)) ! mmr-recheck 9.81

         end do
        end do

        return 
      end subroutine calc_GE_filter_2D



      function get_ge_value(r,rn_vals,ge_vals) result(ge)

        implicit none

        real(wp), intent(IN) :: r           ! [m] Radius from point load 
        real(wp), intent(IN) :: rn_vals(:)  ! [m] Tabulated  radius values
        real(wp), intent(IN) :: ge_vals(:) ! [-] Tabulated Green function values
        real(wp) :: ge

        ! Local variables 
        integer :: k, n 
        real(wp) :: rn_now 
        real(wp) :: wt 

        n = size(rn_vals,1) 

        ! Get radius from point load
        rn_now = r   

        if (rn_now .gt. rn_vals(n)) then

            ge = 0. !ge_vals(n)

        else 

            do k = 1, n-1
                if (rn_now .ge. rn_vals(k) .and. rn_now .lt. rn_vals(k+1)) exit
            end do 

            ! Linear interpolation to get current ge value
            ge = ge_vals(k) &
                + (rn_now-rn_vals(k))/(rn_vals(k+1)-rn_vals(k))*(ge_vals(k+1)-ge_vals(k))   

        end if 

        return

      end function get_GE_value

    subroutine calc_GN_filter_2D(filt,m_earth,r_earth,dx,dy) 
        ! Calculate 2D Green Function

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
        real(wp), intent(IN)  :: m_earth
        real(wp), intent(IN)  :: r_earth
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy   

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp) :: ge_test_0
        real(wp) :: ge_test_1


        ! Get size of filter array and half-width
        n  = size(filt,1) 
        n2 = (n-1)/2 
        
        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_ge_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_ge_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop  
        end if 


        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center
        ! and impose correct Green function value. 

        filt = 0.
        
        do j = -n2, n2 
        do i = -n2, n2

            x  = i*dx 
            y  = j*dy
            
            r  = sqrt(x**2+y**2)  

            ! Get actual index of array
            i1 = i+1+n2 
            j1 = j+1+n2 

            ! Get correct GE value for this point (given by colatitude, theta)
            
            filt(i1,j1) = calc_gn_value(max(dy,r),r_earth,m_earth)* (dx*dy)

         end do
        end do

        return 

      end subroutine calc_GN_filter_2D


    function calc_gn_value(r,r_earth,m_earth) result(gn)
        !Info
        implicit none

        real(wp), intent(IN) :: r
        real(wp), intent(IN)  :: r_earth 
        real(wp), intent(IN)  :: m_earth 
        real(wp)              :: gn

        gn = r_earth / (2. * m_earth * sin (r/(2.*r_earth)) )

        return
    end function calc_gn_value
      
    
    ! === isos physics routines ======================================

    
    subroutine calc_litho_regional(w1,q1,z_bed,H_ice,z_sl,GG,rho_ice,rho_seawater,rho_uppermantle,rho_litho,g)
        ! Calculate the load on the lithosphere as
        ! distributed on an elastic plate. 

        implicit none

        real(wp), intent(INOUT) :: w1(:,:)      ! [m] Lithospheric displacement
        real(wp), intent(INOUT) :: q1(:,:)      ! [Pa] Lithospheric load
        real(wp), intent(IN)    :: z_bed(:,:)   ! [m] Bed elevation
        real(wp), intent(IN)    :: H_ice(:,:)   ! [m] Ice thickness 
        real(wp), intent(IN)    :: z_sl(:,:)    ! [m] Sea level 
        real(wp), intent(IN)    :: GG(:,:)      ! Regional filter function
        real(wp), intent(IN)    :: rho_ice
        real(wp), intent(IN)    :: rho_seawater
        real(wp), intent(IN)    :: rho_uppermantle
        real(wp), intent(IN)    :: rho_litho
        real(wp), intent(IN)    :: g 

        ! Calculate local lithospheric load and displacement first
        call calc_litho_local(w1,q1,z_bed,H_ice,z_sl,rho_ice,rho_seawater,rho_uppermantle,g)
        
        ! Convolve the estimated point load with the regional
        ! filter to obtain the distributed load w1. 

        ! recheck for stability
        call convolve_load_elastic_plate_fft(w1,q1,GG) ! recheck convol
        ! call convolve_load_elastic_plate(w1,q1,GG)

        return

    end subroutine calc_litho_regional

end module isostasy
