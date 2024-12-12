program test_isostasy

    use ncio
    
    use isostasy_defs, only : sp, dp, wp
    use fastisostasy
    use isostasy_benchmarks
    use isos_utils
    use barysealevel
    ! use ice
    
    implicit none
    
    character(len=512) :: parfldr
    character(len=512) :: outfldr
    character(len=512) :: path_par
    character(len=512) :: file_out
    character(len=512) :: file_out_extended

    character(len=56)  :: experiment
    character(len=56)  :: mantle
    character(len=56)  :: lithosphere

    real(wp) :: time, time_bp, time_init, time_end 
    real(wp) :: dtt, dt_out
    integer  :: n, nt
    integer  :: ncx, ncy, nct

    real(wp) :: r0, h0, eta

    integer  :: i, j, nx, ny, nz, slice_time
    real(wp) :: time_now
    real(wp) :: xmin, xmax, dx
    real(wp) :: ymin, ymax, dy
    real(wp) :: xcntr, ycntr
    real(wp), allocatable   :: xc(:)
    real(wp), allocatable   :: yc(:)

    real(wp), allocatable :: z_bed(:, :) 
    real(wp), allocatable :: H_ice(:, :), T_ice(:, : ,:), z_bed_ice(:, : ,:)
    real(wp), allocatable :: time_ice(:), xc_ice(:), yc_ice(:)
    real(wp), allocatable :: z_ss(:, :) 
    
    real(wp), allocatable :: mask(:, :)
    real(wp), allocatable :: z_bed_bench(:, :)

    character(len=256)  :: fldr_path, filename

    type(isos_class)    :: isos1
    type(bsl_class)     :: bsl
    ! type(ice_class)     :: ice

    ! === Define experiment to be run ====

    experiment = "test4"

    ! Tests are defined in Swierczek-Jereczek et al. (2024), GMD.
    ! Additional: "test5" = Lucía's Greenland ice-sheet load (since 15 ka)
    
    write(*,*) "experiment = ", trim(experiment)

        select case(trim(experiment))

            case("test0")
                time_init = 0.
                time_end  = 2000
                dtt       = 10.
                dt_out    = 100.
                dx = 50.e3
                dy = dx
                xmin = -3000.e3
                xmax = abs(xmin)
                ymin = xmin
                ymax = abs(ymin)
            
            case("test1")
                time_init = 0.
                time_end  = 50.e3
                dtt       = 1.0
                dt_out    = 1.e3
                dx        = 50.e3
                dy = dx
                xmin = -3000.e3
                xmax = abs(xmin)
                ymin = xmin
                ymax = abs(ymin)

            case("test2")
                time_init = 0.
                time_end  = 50.e3
                dtt       = 10.0
                dt_out    = 1.e3
                dx        = 50.e3
                dy = dx
                xmin = -3000.e3
                xmax = abs(xmin)
                ymin = xmin
                ymax = abs(ymin)

            case("test3a","test3b","test3c","test3d")
                time_init = 0.
                time_end  = 50.e3
                dtt       = 1.0
                dt_out    = 1.e3
                dx        = 50.e3
                dy = dx
                xmin = -3000.e3
                xmax = abs(xmin)
                ymin = xmin
                ymax = abs(ymin)

            case("test4")
                time_init =  -122.5e3
                time_end  =    0.25e3
                dtt       = 1.0
                dt_out    = 1.e3
                dx = 32.e3
                dy = dx
                xmin = -3040.e3
                xmax = abs(xmin)
                ymin = xmin
                ymax = abs(ymin)

            case("test5")
                time_init = 0.
                time_end  = 15.e3
                dtt       = 1.0
                dt_out    = 1.e3
                dx = 16.e3
                dy = dx
                xmin = -840.e3
                ymin = -1440.e3
                xmax = abs(xmin)
                ymax = abs(ymin)

            case("DEFAULT")
                write(*,*) 'Default values needed, stopping'
                stop

        end select

    
    ! === Define runtime information =========
    ! executable is defined in libisostasy/bin/test_isostasy.x 
    ! output directory should be predefined: output/test-isostasy

    parfldr = "par/"
    outfldr = "output/"
    path_par = trim(parfldr)//"/"//"test_isostasy_"//trim(experiment)//".nml"
    file_out = trim(outfldr)//"/"//"bedtest_"//trim(experiment)// ".nc"
    file_out_extended = trim(outfldr)//"/"//"bedtest_"//trim(experiment)// "_extended.nc"
    write(*,*) "outfldr: ",  trim(outfldr)
    write(*,*) "path_par: ", trim(path_par)
    write(*,*) "file_out: ", trim(file_out)
    write(*,*) "file_out_extended: ", trim(file_out_extended)

    write(*,*) "Initialising viscosity and rigidity fields..."
    mantle = "uniform"
    write(*,*) "viscosity field method = ", trim(mantle)
    lithosphere = "uniform"
    write(*,*) "rigidity method = ", trim(lithosphere)

    write(*,*) "time_init = ", time_init
    write(*,*) "time_end  = ", time_end
    write(*,*) "dtt       = ", dtt
    write(*,*) "dt_out    = ", dt_out

    write(*,*) "Defining grid..."
    if (mod((xmax-xmin), dx) .ne. 0.) then
        print*,'you must have an integer number of points in x-domain'
        stop
    endif
    if (mod((ymax-ymin), dy) .ne. 0.) then
        print*,'you must have an integer number of points in y-domain'
        stop
    endif
    
    nx = int( (xmax-xmin) / dx ) + 1
    ny = int( (ymax-ymin) / dy ) + 1

    allocate(xc(nx))
    allocate(yc(ny))

    do i = 1, nx 
        xc(i) = xmin + (i-1)*dx 
    end do

    do j = 1, ny
       yc(j) = ymin + (j-1)*dy
    end do

    write(*,*) "Grid info: "
    write(*,*) "dx = ", dx
    write(*,*) "dy = ", dy
    write(*,*) "nx, ny = ", nx, ny 
    write(*,*) "range(xc): ", minval(xc), maxval(xc)
    write(*,*) "range(yc): ", minval(yc), maxval(yc)
    
    ! === Define topography fields =========
    write(*,*) "Initialising topographic fields..."

    allocate(z_bed(nx, ny))
    allocate(H_ice(nx, ny))
    allocate(z_ss(nx, ny))
    allocate(z_bed_bench(nx, ny))
    allocate(mask(nx, ny))
    
    ! These inits are potentially overwritten depending on the case. See `select` below.
    z_bed       = 0.0
    H_ice       = 0.0
    z_ss         = 0.0
    z_bed_bench = z_bed

    ! Initialize bedrock model (allocate fields)
    call isos_init(isos1, path_par, "isostasy", nx, ny, dx, dy)
    call bsl_init(bsl, path_par)

    ! Define ice thickness field based on experiment being run...
    select case(trim(experiment))

        case("variable_tau") ! Constant ice thickness everywhere, with a spatially variable tau
            H_ice = 1000.0

            ! Define a mask with three different regions, which will
            ! correspond to different values of tau
            mask(1:int(nx/3.0),:) = 0.0 
            mask(int(nx/3.0)+1:2*int(nx/3.0),:) = 1.0 
            mask(2*int(nx/3.0)+1:nx,:) = 2.0 
            call isos_set_smoothed_field(isos1%domain%tau, [1.e2_wp,1.e3_wp,3.e3_wp], &
                [0.0_wp,1.0_wp,2.0_wp], mask, dx, 150.e3_wp)

        case("point_load")  ! Define ice thickness only in central grid point 
            H_ice = 0.0 
            H_ice(int((nx-1)/2),int((ny-1)/2)) = 1000.0 

        case("test0", "test1","test3a","test3b","test3c","test3d") ! ice disk of R=1000 km and H=1000 m

            r0  = 1000.0e3 ! [m]
            h0  = 1000.0   ! [m]
            eta = 1.e+21   ! [Pa s]
        
            H_ice = 0.
            xcntr = (xmax+xmin)/2.0
            ycntr = (ymax+ymin)/2.0

            do j = 1, ny
               do i = 1, nx
                  if ( (xc(i)-xcntr)**2 + (yc(j)-ycntr)**2  .le. (r0)**2 ) H_ice(i,j) = h0
               end do
            end do

            allocate(T_ice(nx, ny, 2))
            allocate(time_ice(2))
            time_ice(1) = 0.0
            time_ice(2) = 1e-9
            T_ice(:, :, 1) = 0.0_wp
            T_ice(:, :, 2) = H_ice
            z_bed(:, :) = 1e6_wp
            
        case("test2")   ! Spada et al. (2011)
         
            r0 = 6.378e6*10.*3.1416/180. ! * 0.1 ! recheck

            h0  = 1000.0   ! [m] 
            eta = 1.e+21   ! [Pa s]
        
            H_ice = 0.
            xcntr = (xmax + xmin)/2.0
            ycntr = (ymax + ymin)/2.0

            do j = 1, ny
            do i = 1, nx
                if ( (xc(i)-xcntr)**2 + (yc(j)-ycntr)**2  .le. (r0)**2 ) H_ice(i,j) = h0
            end do
            end do

            allocate(T_ice(nx, ny, 2))
            allocate(time_ice(2))
            time_ice(1) = 0.0
            time_ice(2) = 1e-9
            T_ice(:, :, 1) = 0.0_wp
            T_ice(:, :, 2) = H_ice
            z_bed(:, :) = 1e6_wp

        case("test4")  ! ICE6G_D
        ! Comment on “An Assessment of the ICE-6G_C (VM5a) Glacial Isostatic Adjustment Model” by Purcell et al.
        ! W. Richard Peltier, Donald F. Argus, Rosemarie Drummond
        ! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JB013844

            eta = 1.e+21   ! [Pa s]
            H_ice = 0.

            ! Read in H_ice
            filename = "isostasy_data/ice_history/ICE6G_D/ANT-32KM_ICE-6G_D.nc"
            nct = nc_size(filename,"time")
            ncx = nc_size(filename,"xc")
            ncy = nc_size(filename,"yc")

            if ((time_end .lt. nt) .or. (ncx .ne. nx) .or. (ncy .ne. ny)) then
                print*, 'Some dimensions do not correspond to those of the nc file.'
                stop
            endif

            allocate(T_ice(ncx, ncy, nct))
            allocate(time_ice(nct))
            call nc_read(filename, "time", time_ice, start=[1], count=[nct])
            time_ice = -1.e3 * time_ice
            call nc_read(filename, "IceT", T_ice,start=[1, 1, 1], count=[ncx, ncy, nct])
            T_ice = max(T_ice, 0.)
            H_ice = T_ice(:, :, 1)


            ! Read z_bed
            filename = "isostasy_data/topography/ANT-32KM_Latychev.nc"
            ncx = nc_size(filename,"xc")
            ncy = nc_size(filename,"yc")

            if ((ncx .ne. nx) .or. (ncy .ne. ny)) then
                print*, 'Some dimensions do not correspond to those of the nc file.'
                stop
            endif

            allocate(z_bed_ice(ncx, ncy, 1))
            call nc_read(filename,"b", z_bed_ice, start=[1, 1, 1], count=[ncx, ncy, 1])
            
            z_bed = z_bed_ice(:, :, 1)
            z_ss = 0.0_wp

        case("test5")

           ! Lucías Greenland run
            r0 = 6.378e6*10.*3.1416/180.
            h0  = 1000.0   ! [m] 
            eta = 1.e+21   ! [Pa s]
        
            write(*, *) "Reading ice .nc..."
            filename = "isostasy_data/ice_history/sims/greenland/LGM_equilibrium_15kyr.nc"
            nct = nc_size(filename,"time")
            ncx = nc_size(filename,"xc")
            ncy = nc_size(filename,"yc")

            allocate(z_bed_ice(ncx, ncy, nct))
            allocate(T_ice(ncx, ncy, nct))
            allocate(xc_ice(ncx))
            allocate(yc_ice(ncy))
            allocate(time_ice(nct))

            write(*, *) "Reading fields..."
            call nc_read(filename, "z_bed", z_bed_ice, start=[1, 1, 1], count=[ncx, ncy, nct])
            call nc_read(filename, "H_ice", T_ice, start=[1, 1, 1], count=[ncx, ncy, nct])
            call nc_read(filename, "xc", xc_ice, start=[1], count=[ncx])
            call nc_read(filename, "time", time_ice, start=[1], count=[nct])
            call nc_read(filename, "yc", yc_ice, start=[1], count=[ncy])

            H_ice = T_ice(:, :, 1)
            z_bed = z_bed_ice(:, :, 1)
            z_ss = 0.0_wp

        case DEFAULT
            write(*,*) "Error: experiment name not recognized."
            write(*,*) "experiment = ", trim(experiment)
            stop

    end select
    
    time = time_init

    ! Inititalize and write state
    call isos_init_state(isos1, z_bed, T_ice(:, :, 1), time, bsl)
    ! write(*,*) time, isos1%par%time_prognostics, isos1%par%time_diagnostics
    ! stop

    call isos_write_init(isos1, xc, yc, file_out, time_init)
    call isos_write_init_extended(isos1, file_out_extended, time_init)

    ! Determine total number of iterations to run
    nt = ceiling((time_end-time_init)/dtt) + 1

    ! Advance isostasy model
    do n = 1, nt 

        ! Update bedrock
        time = time_init + (n-1)*dtt
        call interp_2d(time_ice, T_ice, time, H_ice)
        call isos_update(isos1, H_ice, time, bsl)
        call bsl_update(bsl, time)

        ! write(*,*) "time = ", time
        ! write(*,*) "extrema H_ice: ", minval(isos1%now%Hice), maxval(isos1%now%Hice)
        ! write(*,*) "count maskgrounded: ", count(isos1%now%maskgrounded)
        ! stop

        if (mod(time-time_init, dt_out) .eq. 0.0) then  ! Write output for this timestep

            ! Calculate benchmark solutions when available and write to file
            select case(trim(experiment))

                case("test1")   ! Calculate analytical solution to elva_disk

                    call isosbench_elva_disk(z_bed_bench, r0, h0, eta, isos1%domain%dx, &
                        isos1%domain%D_lith(1,1), isos1%par%rho_ice, isos1%par%rho_uppermantle, &
                        isos1%par%g,time)

                    call isos_write_step(isos1, bsl, file_out, time, H_ice, z_bed_bench)

                case DEFAULT
                    z_bed_bench = 0.0
                    call isos_write_step(isos1, bsl, file_out, time, H_ice)
                    call isos_write_step_extended(isos1, file_out_extended, time)

            end select

        end if

        if (mod(n, 100) .eq. 1) then
            write(*,*) "time = ", time
        endif

    end do

    contains

    subroutine isos_write_init(isos, xc, yc, filename, time_init)

        implicit none

        type(isos_class), intent(IN) :: isos
        real(wp),         intent(IN) :: xc(:)
        real(wp),         intent(IN) :: yc(:)
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time_init
        
        ! Local variables
        integer :: nf
        
        ! Create the empty netcdf file
        call nc_create(filename)
        
        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename, "xc", x=xc*1e-3, units="km")
        call nc_write_dim(filename, "yc", x=yc*1e-3, units="km")
        call nc_write_dim(filename, "time", x=time_init, dx=1.0_wp, nx=1, &
            units="year", unlimited=.TRUE.)

        ! Write constant fields
        call nc_write(filename, "He_lith", isos%out%He_lith, units="km", &
            long_name="Lithosphere thickness", dim1="xc", dim2="yc" ,start=[1, 1])

        ! call nc_write(filename,"GE",isos%out%GE, units="", &
        !     long_name="Elastic Green function", dim1="xc", dim2="yc", start=[1, 1])

        ! if (isos%par%interactive_sealevel) then
        !     call nc_write(filename, "GN", isos%out%GN, units="", &
        !     long_name="SSH Green function", dim1="xc", dim2="yc", start=[1, 1])
        ! end if
        
        if (isos%par%method .eq. 2) then
            ! call nc_write(filename, "kei", isos%out%kei, units="", &
            !     long_name="Kelvin function filter", dim1="xc", dim2="yc", start=[1, 1])

            ! call nc_write(filename,"GV",isos%out%GV, units="", &
            !     long_name="Viscous Green function", dim1="xc", dim2="yc", start=[1, 1])

            call nc_write(filename, "tau", isos%out%tau, units="yr", &
                long_name="Asthenosphere relaxation timescale", &
                dim1="xc", dim2="yc", start=[1, 1])
        end if

        if (isos%par%method .eq. 3) then
            call nc_write(filename, "log10_eta_eff", log10(isos%out%eta_eff), units="Pa s", &
                long_name="Effective upper-mantle viscosity", &
                dim1="xc", dim2="yc", start=[1, 1])

            ! call nc_write(filename, "kappa", isos%out%kappa, units="", &
            !     long_name="Pseudodifferential operator in Fourier space", &
            !     dim1="xc", dim2="yc", start=[1, 1])
        end if

        return

    end subroutine isos_write_init

    subroutine isos_write_init_extended(isos, filename, time_init)

        implicit none

        type(isos_class), intent(IN) :: isos
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time_init
        
        ! Local variables
        integer :: nf
        
        ! Create the empty netcdf file
        call nc_create(filename)
        
        ! Add grid axis variables to netcdf file
        call nc_write_dim(filename, "xc", x=isos%domain%xc*1e-3, units="km")
        call nc_write_dim(filename, "yc", x=isos%domain%yc*1e-3, units="km")
        call nc_write_dim(filename, "time", x=time_init, dx=1.0_wp, nx=1, &
            units="year", unlimited=.TRUE.)

        ! Write constant fields
        call nc_write(filename, "He_lith", isos%domain%He_lith, units="km", &
            long_name="Lithosphere thickness", dim1="xc", dim2="yc" ,start=[1, 1])

        call nc_write(filename,"GE",isos%domain%GE, units="", &
            long_name="Elastic Green function", dim1="xc", dim2="yc", start=[1, 1])

        if (isos%par%interactive_sealevel) then
            call nc_write(filename, "GN", isos%domain%GN, units="", &
            long_name="SSH Green function", dim1="xc", dim2="yc", start=[1, 1])
        end if
        
        if (isos%par%method .eq. 2) then
            call nc_write(filename, "kei", isos%domain%kei, units="", &
                long_name="Kelvin function filter", dim1="xc", dim2="yc", start=[1, 1])

            call nc_write(filename,"GV",isos%domain%GV, units="", &
                long_name="Viscous Green function", dim1="xc", dim2="yc", start=[1, 1])

            call nc_write(filename, "tau", isos%domain%tau, units="yr", &
                long_name="Asthenosphere relaxation timescale", &
                dim1="xc", dim2="yc", start=[1, 1])
        end if

        if (isos%par%method .eq. 3) then
            call nc_write(filename, "log10_eta_eff", log10(isos%domain%eta_eff), units="Pa s", &
                long_name="Effective upper-mantle viscosity", &
                dim1="xc", dim2="yc", start=[1, 1])

            ! call nc_write(filename, "kappa", isos%out%kappa, units="", &
            !     long_name="Pseudodifferential operator in Fourier space", &
            !     dim1="xc", dim2="yc", start=[1, 1])
        end if

        call nc_write(filename, "maskactive", isos%domain%maskactive, units="1", &
            long_name="Active mask", dim1="xc", dim2="yc", start=[1, 1])
        return

    end subroutine isos_write_init_extended

    ! Write results to file
    subroutine isos_write_step(isos, bsl, filename, time, H_ice, z_bed_bench)

        implicit none 
        
        type(isos_class), intent(IN) :: isos
        type(bsl_class),  intent(IN) :: bsl
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time
        real(wp),         intent(IN) :: H_ice(:, :)
        real(wp),         intent(IN), optional :: z_bed_bench(:, :)

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename, "time", ncid)
        call nc_read(filename, "time", time_prev, start=[n], count=[1], ncid=ncid)
        if (abs(time-time_prev) .gt. 1e-5) n = n+1

        ! Update the time step
        call nc_write(filename,"time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        call nc_write(filename, "bsl", isos%now%bsl, units="m", long_name="Barystatic sea level", &
            dim1="time", start=[n], count=[1], ncid=ncid)
        call nc_write(filename, "A_ocean", bsl%A_ocean_now, units="m^2", long_name="Ocean area", &
            dim1="time", start=[n], count=[1], ncid=ncid)

        ! Write variables
        call nc_write(filename, "H_ice", isos%out%Hice, units="m", long_name="Ice thickness", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "z_ss", isos%out%z_ss, units="m", long_name="Sea-surface height", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename,"z_bed", isos%out%z_bed, units="m", &
            long_name="Bedrock elevation", dim1="xc", dim2="yc", dim3="time", &
            start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "dwdt", isos%out%dwdt, units="m/yr", &
            long_name="Bedrock elevation change", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "w_viscous", isos%out%w, units="m", &
            long_name="Displacement (viscous)", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "w_elastic", isos%out%we, units="m", &
            long_name="Displacement (elastic)", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "z_ss_perturbation", isos%out%dz_ss, units="m", &
            long_name="Geoid displacement", dim1="xc", dim2="yc", dim3="time", &
            start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "column_anomaly", isos%out%canom_full, units="N m^-2", &
            long_name = "Anomaly in column pressure", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "ocean_mask", isos%out%maskocean, units="1", &
            long_name = "Ocean mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "grounded_mask", isos%out%maskgrounded, units="1", &
            long_name = "Grounded mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "continent_mask", isos%out%maskcontinent, units="1", &
            long_name = "Continent mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        !# TODO: We could remove the analytical solutions alltogether and make the comparison
        ! as post-processing.
        if (present(z_bed_bench)) then
            ! Compare with benchmark solution
            call nc_write(filename, "z_bed_bench", z_bed_bench,units="m", &
                long_name="Benchmark bedrock elevation", &
                dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)
            call nc_write(filename, "err_z_bed", isos%out%w - z_bed_bench,units="m", &
                long_name="Error in bedrock elevation", &
                dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)
        end if

        call nc_close(ncid)     ! Close the netcdf file

        return 
    end subroutine isos_write_step


    ! Write results to file
    subroutine isos_write_step_extended(isos, filename, time)

        implicit none 
        
        type(isos_class), intent(IN) :: isos
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)

        ! Determine current writing time step
        n = nc_size(filename, "time", ncid)
        call nc_read(filename, "time", time_prev, start=[n], count=[1], ncid=ncid)
        if (abs(time-time_prev) .gt. 1e-5) n = n+1

        ! Update the time step
        call nc_write(filename,"time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! Write variables
        call nc_write(filename, "H_ice", isos%now%Hice, units="m", long_name="Ice thickness", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "z_ss", isos%now%z_ss, units="m", long_name="Sea-surface height", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename,"z_bed", isos%now%z_bed, units="m", &
            long_name="Bedrock elevation", dim1="xc", dim2="yc", dim3="time", &
            start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "dwdt", isos%now%dwdt, units="m/yr", &
            long_name="Bedrock elevation change", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "w_viscous", isos%now%w, units="m", &
            long_name="Displacement (viscous)", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "w_elastic", isos%now%we, units="m", &
            long_name="Displacement (elastic)", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "z_ss_perturbation", isos%now%dz_ss, units="m", &
            long_name="Geoid displacement", dim1="xc", dim2="yc", dim3="time", &
            start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "column_anomaly", isos%now%canom_full, units="N m^-2", &
            long_name = "Anomaly in column pressure", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "ocean_mask", isos%now%maskocean, units="1", &
            long_name = "Ocean mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "grounded_mask", isos%now%maskgrounded, units="1", &
            long_name = "Grounded mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "continent_mask", isos%now%maskcontinent, units="1", &
            long_name = "Continent mask", &
            dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)


        call nc_close(ncid)     ! Close the netcdf file

        return 
    end subroutine isos_write_step_extended

end program test_isostasy