program test_isostasy

    use ncio
    
    use isostasy_defs, only : sp, dp, wp
    use isostasy
    use isostasy_benchmarks
    use isos_utils
    ! use ice
    
    implicit none
    
    character(len=512) :: parfldr
    character(len=512) :: outfldr
    character(len=512) :: path_par
    character(len=512) :: file_out 

    character(len=56)  :: experiment
    character(len=56)  :: visc_method
    character(len=56)  :: rigidity_method

    real(wp) :: time, time_bp 
    real(wp) :: time_init 
    real(wp) :: time_end 
    real(wp) :: dtt
    real(wp) :: dt_out 
    integer  :: n, nt
    integer  :: ncx, ncy, nct

    real(wp) :: r0, h0, eta 

    integer  :: i, j, nx, ny, nz, slice_time
    real(wp) :: time_now
    real(wp) :: xmin, xmax, dx
    real(wp) :: ymin, ymax, dy
    real(wp) :: xcntr, ycntr
    real(wp), allocatable :: xc(:)
    real(wp), allocatable :: yc(:)

    real(wp), allocatable :: z_bed_ref(:, :) 
    real(wp), allocatable :: H_ice_ref(:, :) 
    real(wp), allocatable :: z_sl_ref(:, :) 
    
    real(wp), allocatable :: z_bed(:, :) 
    real(wp), allocatable :: H_ice(:, :), T_ice(:, : ,:), z_bed_ice(:, : ,:)
    real(wp), allocatable :: time_ice(:), xc_ice(:), yc_ice(:)
    real(wp), allocatable :: z_sl(:, :) 
    
    real(wp), allocatable :: mask(:, :)
    real(wp), allocatable :: z_bed_bench(:, :)

    character(len=256)  :: fldr_path, filename

    type(isos_class)    :: isos1
    ! type(ice_class)     :: ice

    ! mmr recheck
    
    ! ! === Define runtime information =========
    ! ! executable is defined in libisostasy/bin/test_isostasy.x 
    ! ! output directory should be predefined: output/test-isostasy

    ! outfldr = "output/test-isostasy"
    ! path_par = trim(outfldr)//"/"//"test_isostasy.nml" 
    ! file_out = trim(outfldr)//"/"//"bedtest.nc"

    ! write(*,*) "outfldr: ",  trim(outfldr)
    ! write(*,*) "path_par: ", trim(path_par)
    ! write(*,*) "file_out: ", trim(file_out)

    
    ! === Define experiment to be run ====

    !experiment = "constant_thickness"
    !experiment = "variable_tau"
    !experiment = "point_load"
    !  experiment = "test0" ! 1000 m radius ice disk of 1000 m heights, solved with ELRA (case=2)
    
! tests in Swierczek-Jereczek et al. (2023), GMD.
    
    !   experiment = "test1" ! 1000 m radius ice disk of 1000 m heights, solved with ELVA Benchmark: analytical (Bueler et al. 2007); dtt = 1; dtout = 1000. ; time_end = 50.e3

    !TODO Fix analytic solution of test1
    experiment = "test2"   ! Benchmark: Spada et al. (2011) disc

    !  experiment = "test3a"  ! Gaussian reduction of lithospheric thickness at centre
    !  experiment = "test3b"  ! Gaussian increase of lithospheric thickness at centre
    !  experiment = "test3c"  ! Gaussian reduction of viscosity at centre
    !  experiment = "test3d"  ! Gaussian increase of viscosity at centre
    !  experiment = "test4"  ! ais ice6g_d 
    !  experiment = "test5"  ! Lucía's Greenland ice-sheet load (since 15 ka)
    
     
    write(*,*) "experiment = ", trim(experiment)

        select case(trim(experiment))

        case("test0")
            
                time_init = 0. 
                time_end  = 2000 
                dtt       = 10. 
                dt_out    = 100. 
                dx = 50.e3 
                xmin = -3000.e3 
                ymin = xmin !* 2
            
        case("test1")

            time_init = 0. 
            time_end  = 50.e3 
            dtt       = 1.0 
            dt_out    = 1.e3
            dx        = 50.e3 
            xmin      = -3000.e3 
            ymin      = xmin

        case("test2")

            time_init = 0. 
            time_end  = 50.e3 
            dtt       = 10.0 
            dt_out    = 1.e3
            dx        = 50.e3 
            xmin      = -3000.e3 
            ymin      = xmin

        case("test3a","test3b","test3c","test3d")
            
            time_init = 0. 
            time_end  = 50.e3 
            dtt       = 1.0 
            dt_out    = 1.e3
            dx        = 50.e3 
            xmin      = -3000.e3 
            ymin      = xmin

        case("test4")
            
            time_init =  -122.5e3 
            time_end  =     2.5e3
            dtt       = 1.0  
            dt_out    = 1.e3 
            dx = 32.e3 
            xmin = -3040.e3 
            ymin = xmin
    
        case("test5")
            
            time_init = 0. 
            time_end  = 15.e3 
            dtt       = 1.0   
            dt_out    = 1.e3 
            dx = 16.e3 
            xmin = -840.e3 
            ymin = -1440.e3

        case("DEFAULT")

            write(*,*) 'Default values needed, stopping'
            stop

        end select

    
    ! === Define runtime information =========
    ! executable is defined in libisostasy/bin/test_isostasy.x 
    ! output directory should be predefined: output/test-isostasy

    parfldr = "par/"
    outfldr = "output/test-isostasy"
    path_par = trim(parfldr)//"/"//"test_isostasy_"//trim(experiment)//".nml" 
    file_out = trim(outfldr)//"/"//"bedtest_"//trim(experiment)// ".nc"

    write(*,*) "outfldr: ",  trim(outfldr)
    write(*,*) "path_par: ", trim(path_par)
    write(*,*) "file_out: ", trim(file_out)
    
 ! === Define viscosity field to be used ====

    visc_method = "uniform"
    
    write(*,*) "viscosity field method = ", trim(visc_method)

! === Define rigidity field to be used ====

    rigidity_method = "uniform"
    
    write(*,*) "rigidity method = ", trim(rigidity_method)

    
    ! === Define simulation time ========

    write(*,*) "time_init = ", time_init 
    write(*,*) "time_end  = ", time_end 
    write(*,*) "dtt       = ", dtt 
    write(*,*) "dt_out    = ", dt_out

    ! === Define grid information ============
    ! TODO: make this more general
    xmax = abs(xmin)
    if (mod((xmax-xmin),dx).ne.0.) then
        print*,'you must have an integer number of points in x-domain'
        stop
    endif
    
    ! TODO: dy optionally different from dx
    nx      = int( (xmax-xmin) / dx ) + 1
    dy      = dx
    ymax    = abs(ymin)
    ny      = int( (ymax-ymin) / dy ) + 1
    
    if (mod((ymax-ymin),dy).ne.0.) then
        print*,'you must have an integer number of points in y-domain'
        stop
    endif

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

    allocate(z_bed_ref(nx,ny))
    allocate(H_ice_ref(nx,ny))
    allocate(z_sl_ref(nx,ny))
    allocate(z_bed(nx,ny))
    allocate(H_ice(nx,ny))
    allocate(z_sl(nx,ny))
    
    allocate(z_bed_bench(nx,ny))

    allocate(mask(nx,ny))


    z_bed_ref   = 0.0 
    H_ice_ref   = 0.0 
    z_sl_ref    = -1e3 !mmr test3a
    
    z_bed       = 0.0 
    H_ice       = 0.0  
    z_sl        = -1e3 !mmr test3a
    
    z_bed_bench = z_bed 

    write(*,*) "Initial fields defined."

    ! Initialize bedrock model (allocate fields)
    call isos_init(isos1, path_par, "isostasy", nx, ny, dx, dy)

    ! Define ice thickness field based on experiment being run...
    select case(trim(experiment))

        case("constant_thickness")
            ! Set ice thickness to a constant value everywhere

            H_ice = 1000.0

        case("variable_tau")
            ! Set ice thickness to a constant value everywhere,
           ! with a spatially variable field of tau

            H_ice = 1000.0

            ! Define a mask with three different regions, which will
            ! correspond to different values of tau
            mask(1:int(nx/3.0),:) = 0.0 
            mask(int(nx/3.0)+1:2*int(nx/3.0),:) = 1.0 
            mask(2*int(nx/3.0)+1:nx,:) = 2.0 

            ! Define tau field using the mask
!mmr           call isos_set_field(isos1%domain%tau,[1e2,1e3,3e3],[0.0_wp,1.0_wp,2.0_wp],mask,dx,sigma=150e3)

            call isos_set_smoothed_field(isos1%domain%tau, [1.e2_wp,1.e3_wp,3.e3_wp], &
                [0.0_wp,1.0_wp,2.0_wp], mask, dx, 150.e3_wp)

        case("point_load")
            ! Define ice thickness only in one grid point 

            H_ice = 0.0 
            H_ice(int((nx-1)/2),int((ny-1)/2)) = 1000.0 

         case("test0", "test1","test3a","test3b","test3c","test3d")

            ! Bueler et al. (2007): ice disk in a circle of radius 1000 km and thickness 1000 m

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

            allocate(T_ice(nx,ny,1))
            allocate(time_ice(1))

            T_ice(:, :,1) = H_ice(:, :)
            
      case("test2")

         ! Spada et al. (2011)
         
            r0 = 6.378e6*10.*3.1416/180. ! * 0.1 ! recheck

            h0  = 1000.0   ! [m] 
            eta = 1.e+21   ! [Pa s]
        
            H_ice = 0.
            xcntr = (xmax+xmin)/2.0
            ycntr = (ymax+ymin)/2.

            do j = 1, ny
               do i = 1, nx
                  if ( (xc(i)-xcntr)**2 + (yc(j)-ycntr)**2  .le. (r0)**2 ) H_ice(i,j) = h0
                end do
             end do

             allocate(T_ice(nx,ny,1))
             allocate(time_ice(1))
            
            T_ice(:, :,1) = H_ice(:, :)
            
         case("test4")

! ICE6G_D
! Comment on “An Assessment of the ICE-6G_C (VM5a) Glacial Isostatic Adjustment Model” by Purcell et al.
! W. Richard Peltier, Donald F. Argus, Rosemarie Drummond
! https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2016JB013844

            
            R0 = 6.378e6*10.*3.1416/180. 

            h0  = 1000.0   ! [m] 
            eta = 1.e+21   ! [Pa s]
        
            H_ice = 0.

            ! Read in H_ice

            filename = "/Users/montoya/work/ice_data/Antarctica/ANT-32KM/ANT-32KM_ICE-6G_D.nc"
            
            Nct = nc_size(filename,"time")
            ncx = nc_size(filename,"xc")
            ncy = nc_size(filename,"yc")

            if (time_end.lt.nt) then

               print*,'Need to increase time_end to read full data length'
               stop
               
            endif

            if (ncx.ne.nx) then
                           
               print*,'ncx not equal to nx'
               stop
               
            endif

            if (ncy.ne.ny) then
                           
               print*,'ncx not equal to nx'
               stop
               
            endif


            allocate(T_ice(ncx,ncy,nct))
            allocate(time_ice(nct))
            
            call nc_read(filename,"time",time_ice,start=[1],count=[nct])
            time_ice = -1.e3 * time_ice
            
            call nc_read(filename,"IceT",T_ice,start=[1,1,1],count=[ncx,ncy,nct])
            T_ice = max(T_ice, 0.)

            ! Initialize H_ice
                    
           H_ice = T_ice(:, :,1)
           H_ice_ref = T_ice(:, :,1)
           
           z_bed_ref = 0. ! recheck - Jan  - what do you put here? ie what are the initial conditions? LGM?


        case("test5")

           ! Lucías Greenland run
           
            r0 = 6.378e6*10.*3.1416/180. 

            h0  = 1000.0   ! [m] 
            eta = 1.e+21   ! [Pa s]
        
            ! Read in H_ice

            filename = "/Users/montoya/work/isostasy/cloned_isostasy/isostasy/input/gris_lgg/LGM_equilibrium_15kyr.nc"
            
            Nct = nc_size(filename,"time")
            ncx = nc_size(filename,"xc")
            ncy = nc_size(filename,"yc")


            if (time_end.lt.nt) then
               print*,'Need to increase time_end to read full data length'
               stop               
            endif

            if (ncx.ne.nx) then                        
              print*,'ncx not equal to nx'
               stop
            endif

            if (ncy.ne.ny) then                          
               print*,'ncx not equal to nx'
               stop
            endif
            
            allocate(z_bed_ice(ncx,ncy,nct))
            allocate(T_ice(ncx,ncy,nct))
                        
            allocate(xc_ice(ncx))
            allocate(yc_ice(ncy))
            allocate(time_ice(nct))

            call nc_read(filename,"z_bed",z_bed_ice,start=[1,1,1],count=[ncx,ncy,nct])
            call nc_read(filename,"H_ice",T_ice,start=[1,1,1],count=[ncx,ncy,nct])
            call nc_read(filename,"xc",xc_ice,start=[1],count=[ncx])
            call nc_read(filename,"time",time_ice,start=[1],count=[nct])
            call nc_read(filename,"yc",yc_ice,start=[1],count=[ncy])
            

            H_ice = T_ice(:, :,1)
            z_bed = z_bed_ice(:, :,1)
            
            H_ice_ref = T_ice(:, :,1)
            z_bed_ref = z_bed_ice(:, :,1)

        case DEFAULT

            write(*,*) "Error: experiment name not recognized."
            write(*,*) "experiment = ", trim(experiment)
            stop 

         end select
         
    ! Inititalize state
    call isos_init_state(isos1,z_bed,H_ice,z_sl,z_bed_ref,H_ice_ref,z_sl_ref,time=time_init) 

    
    ! Initialize writing output
    call isos_write_init(isos1, xc, yc, file_out, time_init)

    ! Determine total number of iterations to run
    nt = ceiling((time_end-time_init)/dtt) + 1 

    
    ! Advance isostasy model
    do n = 1, nt 

        time = time_init + (n-1)*dtt
        call interp_linear(time_ice, T_ice, time, H_ice)

        ! TODO: check if following lines can be deleted. It looks like an obvious "yes", but
        ! maybe I am missing something?)
        ! H_ice = H_ice
        ! z_bed = z_bed

        ! Update bedrock
        call isos_update(isos1, H_ice, time, z_sl)

        
        if (mod(time-time_init,dt_out) .eq. 0.0) then
!           if (mod(time,dt_out) .eq. 0.0) then
           
            ! Write output for this timestep

            ! Calculate benchmark solutions when available and write to file
            select case(trim(experiment))

            case("test1")

               ! Calculate analytical solution to elva_disk

               ! mmr: comment this to spare time; enable for test1 only
                call isosbench_elva_disk(z_bed_bench, r0, h0, eta, isos1%domain%dx, &
                    isos1%domain%D_lith(1,1), isos1%par%rho_ice, isos1%par%rho_uppermantle, &
                    isos1%par%g,time)

                call isos_write_step(isos1, file_out, time, H_ice, z_sl, z_bed_bench)

                case DEFAULT
                    z_bed_bench = 0.0
                    call isos_write_step(isos1, file_out, time, H_ice, z_sl)

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
        call nc_write_dim(filename,"xc",x=xc*1e-3,units="km")
        call nc_write_dim(filename,"yc",x=yc*1e-3,units="km")
        call nc_write_dim(filename,"time",x=time_init,dx=1.0_wp,nx=1,units="year",unlimited=.TRUE.)

        ! Write dimensions for regional filter too
        call nc_write_dim(filename, "xf", x=0, dx=1, nx=size(isos%domain%GV,1), units="pt")
        call nc_write_dim(filename, "yf", x=0, dx=1, nx=size(isos%domain%GV,2), units="pt")

        ! Write constant fields
        call nc_write(filename, "eta_eff", isos%domain%eta_eff, units="Pa s", &
            long_name="Asthenosphere effective viscosity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1])

        call nc_write(filename, "He_lith", isos%domain%He_lith, units="km", &
            long_name="Lithosphere thickness", &
            dim1="xc",dim2="yc",dim3="time",start=[1,1])

        call nc_write(filename, "D_lith", isos%domain%D_lith, units="N m", &
            long_name="Lithosphere rigidity", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1])

        ! call nc_write(filename, "tau", isos%domain%tau, units="yr", &
        !     long_name="Asthenosphere relaxation timescale", &
        !     dim1="xc",dim2="yc",start=[1,1]) 

        call nc_write(filename, "kei", isos%domain%kei, units="", &
            long_name="Kelvin function filter", dim1="xf",dim2="yf", start=[1,1])

        call nc_write(filename,"GV",isos%domain%GV, units="", &
            long_name="Regional elastic plate filter", dim1="xf", dim2="yf", start=[1,1])

        ! ! TODO recheck convol
        ! call nc_write(filename, "GN", isos%domain%GN, units="", &
        !         long_name="Geoid's elastic plate filter", dim1="xc", dim2="yc", start=[1,1])

        return

    end subroutine isos_write_init

    ! Write results to file
    subroutine isos_write_step(isos,filename,time,H_ice,z_sl,z_bed_bench)

        implicit none 
        
        type(isos_class), intent(IN) :: isos
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time
        real(wp),         intent(IN) :: H_ice(:, :) 
        real(wp),         intent(IN) :: z_sl(:, :) 
        real(wp),         intent(IN), optional :: z_bed_bench(:, :)

        ! Local variables
        integer  :: ncid, n
        real(wp) :: time_prev 

        ! Open the file for writing
        call nc_open(filename,ncid,writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_size(filename, "time", ncid)
        call nc_read(filename, "time", time_prev, start=[n], count=[1], ncid=ncid) 
        if (abs(time-time_prev) .gt. 1e-5) n = n+1 

        ! Update the time step
        call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
        
        ! Write variables
        call nc_write(filename, "H_ice", H_ice, units="m", long_name="Ice thickness", &
              dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        call nc_write(filename, "z_sl", z_sl, units="m", long_name="Sea level elevation", &
              dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        call nc_write(filename,"z_bed", isos%now%z_bed, units="m", &
            long_name="Bedrock elevation", dim1="xc", dim2="yc", dim3="time", &
            start=[1,1,n],ncid=ncid)

        call nc_write(filename, "dzbdt", isos%now%dzbdt, units="m/yr", &
            long_name="Bedrock elevation change", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        call nc_write(filename, "w_viscous", isos%now%w, units="m", &
            long_name="Displacement (viscous)", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        call nc_write(filename, "w_elastic", isos%now%we, units="m", &
            long_name="Displacement (elastic)", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        call nc_write(filename, "ssh_perturbation", isos%now%ssh_perturb, units="m", &
            long_name="Geoid displacement", dim1="xc", dim2="yc", dim3="time", &
            start=[1,1,n], ncid=ncid)

        call nc_write(filename, "column_anomaly", isos%now%canom_full, units="N m^-2", &
            long_name = "Anomaly in column pressure", &
            dim1="xc", dim2="yc", dim3="time", start=[1,1,n], ncid=ncid)

        if (present(z_bed_bench)) then 
            ! Compare with benchmark solution 
            call nc_write(filename,"z_bed_bench",z_bed_bench,units="m",long_name="Benchmark bedrock elevation", &        
                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)                                                  
            call nc_write(filename,"err_z_bed",isos%now%z_bed - z_bed_bench,units="m",long_name="Error in bedrock elevation", & 
                dim1="xc",dim2="yc",dim3="time",start=[1,1,n],ncid=ncid)  
        end if 

     
        ! Close the netcdf file
        call nc_close(ncid)

        return 

    end subroutine isos_write_step


    subroutine interp_linear(x,y,xout,yout)
        ! Simple linear interpolation of a point

        implicit none

        real(wp), dimension(:), intent(IN) :: x 
        real(wp), dimension(:, : ,:), intent(IN) :: y
        real(wp), intent(IN) :: xout
        real(wp), dimension(:, :),intent(OUT) :: yout
        integer  :: i, j, n, nout
        real(wp) :: alph

        n    = size(x)

        if (xout .lt. x(1)) then
            yout = y(:, :,1)
        else if (xout .gt. x(n)) then
            yout = y(:, :,n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do

            if (j .eq. 1) then
                yout = y(:, :,1)
            else if (j .eq. n+1) then
                yout = y(:, :,n)
            else
               alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(:, :,j-1) + alph*(y(:, :,j) - y(:, :,j-1))
!               print*,'hola', xout, alph, yout(90,90), y(90,90,j), y(90,90,j-1)
             end if
        end if

        return
    
      end subroutine interp_linear
    
end program test_isostasy

