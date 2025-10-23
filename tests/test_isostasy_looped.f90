program test_isostasy_looped

    use ncio
    use nml
    
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
    character(len=512) :: file_out_bsl
    character(len=512) :: file_isos_restart
    character(len=512) :: file_bsl_restart

    character(len=56)  :: experiment
    character(len=56), allocatable  :: experiments(:)
    character(len=56)  :: mantle
    character(len=56)  :: lithosphere

    real(wp) :: time, time_bp, time_init, time_end 
    real(wp) :: dtt, dt_out
    integer  :: n, nt
    integer  :: ncx, ncy, nct

    real(wp) :: r0, h0, eta

    integer  :: i, j, nx, ny, nz, slice_time, i_xp, n_xp
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

    n_xp = 2
    allocate(experiments(n_xp))
    ! experiments(1) = "test1"
    experiments(1) = "test2a"
    experiments(2) = "test2b"
    ! experiments(3) = "test3a"
    ! experiments(4) = "test3b"
    ! experiments(5) = "test3c"
    ! experiments(6) = "test3d"
    ! experiments(7) = "test4"
    ! experiments(8) = "test5"

    i_xp = 0
    do i_xp = 1, n_xp

    experiment = experiments(i_xp)

    if (allocated(xc))          deallocate(xc)
    if (allocated(yc))          deallocate(yc)
    if (allocated(z_bed))       deallocate(z_bed)
    if (allocated(H_ice))       deallocate(H_ice)
    if (allocated(T_ice))       deallocate(T_ice)
    if (allocated(time_ice))    deallocate(time_ice)
    if (allocated(z_ss))        deallocate(z_ss)
    if (allocated(mask))        deallocate(mask)
    if (allocated(z_bed_bench)) deallocate(z_bed_bench)
    if (allocated(xc_ice))      deallocate(xc_ice)
    if (allocated(yc_ice))      deallocate(yc_ice)

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

            case("test2a", "test2b")
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
    ! executable is defined in libisostasy/bin/test_isostasy_looped.x 
    ! output directory should be predefined: output/test-isostasy

    parfldr = "par/"
    outfldr = "output/"
    path_par = trim(parfldr)//"/"//"test_isostasy_"//trim(experiment)//".nml"
    file_out = trim(outfldr)//"/"//"bedtest_"//trim(experiment)// ".nc"
    file_out_extended = trim(outfldr)//"/"//"bedtest_"//trim(experiment)// "_extended.nc"
    file_out_bsl = trim(outfldr)//"/"//"bsl_"//trim(experiment)// ".nc"
    file_isos_restart = trim(outfldr)//"/"//"isos_restart_"//trim(experiment)// ".nc"
    file_bsl_restart = trim(outfldr)//"/"//"bsl_restart_"//trim(experiment)// ".nc"
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

    ! Define ice thickness field based on experiment being run...
    select case(trim(experiment))

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
            
        case("test2a", "test2b")   ! Spada et al. (2011)
         
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

    ! Initialize fastisostasy
    call bsl_init(bsl, path_par, time)
    call isos_init(isos1, path_par, "isostasy", nx, ny, dx, dy)
    call isos_init_state(isos1, z_bed, T_ice(:, :, 1), time, bsl)
    
    call isos_write_init(isos1, xc, yc, file_out, time_init)
    call isos_write_init_extended(isos1, file_out_extended, time_init)
    call bsl_write_init(bsl, file_out_bsl, time_init)

    ! Determine total number of iterations to run
    nt = ceiling((time_end-time_init)/dtt) + 1

    ! Advance isostasy model
    do n = 1, nt 

        ! Update bedrock
        time = time_init + (n-1)*dtt
        call interp_2d(time_ice, T_ice, time, H_ice)

        call bsl_update(bsl, time)
        call isos_update(isos1, H_ice, time, bsl)

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

                    call isos_write_step(isos1, file_out, time, H_ice, z_bed_bench)

                case DEFAULT
                    z_bed_bench = 0.0
                    call isos_write_step(isos1, file_out, time, H_ice)
                    call isos_write_step_extended(isos1, file_out_extended, time)
                    call bsl_write_step(bsl, file_out_bsl, time)
            end select

        end if

        if (mod(n, 100) .eq. 1) then
            write(*,*) "time = ", time, "dbsl = ", isos1%now%deltaV_bsl / bsl%A_ocean_now
        endif

    end do

    call isos_restart_write(isos1, file_isos_restart, time)
    call bsl_restart_write(bsl, file_bsl_restart, time)
    
    end do  ! loop over experiments

end program test_isostasy_looped