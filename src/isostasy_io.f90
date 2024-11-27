module isostasy_io

    use isostasy_defs, only : wp, isos_class, isos_domain_class
    use nml
    use ncio

    implicit none

    private
    
    public :: isos_grid_write
    public :: isos_restart_write
    public :: isos_restart_read
    
    contains

    subroutine isos_grid_write(domain, fnm, create)

        implicit none 
        type(isos_domain_class), intent(IN) :: domain
        character(len=*),  intent(IN) :: fnm
        logical,           intent(IN) :: create

        ! Local variables 
        character(len=16) :: xnm 
        character(len=16) :: ynm 
        character(len=16) :: grid_mapping_name

        xnm = "xc"
        ynm = "yc" 
        grid_mapping_name = "crs"

        ! Create the netcdf file if desired
        if (create) then 
            ! Create the empty netcdf file
            call nc_create(fnm)

            call nc_write_dim(fnm,xnm,x=domain%xc*1e-3,units="km")
            call nc_write_dim(fnm,ynm,x=domain%yc*1e-3,units="km")

        end if

        call nc_write(fnm, "x2D", domain%x*1e-3, dim1=xnm, dim2=ynm, &
            units="km", grid_mapping=grid_mapping_name)
        call nc_write(fnm, "y2D", domain%y*1e-3, dim1=xnm, dim2=ynm, &
            units="km", grid_mapping=grid_mapping_name)
        return

    end subroutine isos_grid_write


    subroutine isos_restart_write(isos, filename, time, init)

        implicit none 

        type(isos_class),  intent(IN) :: isos
        character(len=*),  intent(IN) :: filename
        real(wp),          intent(IN) :: time
        logical, optional, intent(IN) :: init

        ! Local variables
        integer  :: ncid, n, nx, ny
        logical  :: initialize_file
        real(wp) :: time_prev
        character(len=2) :: xnm, ynm 
        real(wp), allocatable :: bsl(:, :)
        real(wp), allocatable :: bsl_ref(:, :)

        initialize_file = .TRUE. 
        if (present(init)) initialize_file = init

        xnm = "xc"
        ynm = "yc"

        nx    = isos%domain%nx
        ny    = isos%domain%ny

        allocate(bsl(nx, ny))
        allocate(bsl_ref(nx, ny))
        bsl(:, :) = isos%now%bsl
        bsl_ref(:, :) = 0.0 ! isos%ref%bsl
        write(*,*) 'Ref BSL:', isos%ref%bsl
        write(*,*) 'Ref BSL matrix:', minval(bsl_ref), maxval(bsl_ref)

        if (initialize_file) then
            ! Initialize file by writing grid info
            call isos_grid_write(isos%domain, filename, create=.TRUE.)
            call nc_write_dim(filename,"time",     x=time,dx=1.0_wp, nx=1, units="years")
            call nc_write_dim(filename,"ref_now",  x=1.0_wp, dx=1.0_wp, nx=2, units="1")
            call nc_write_attr(filename,xnm,"standard_name","projection_x_coordinate")
            call nc_write_attr(filename,ynm,"standard_name","projection_y_coordinate")
        end if
        
        ! == Begin writing data ==============================================
        
        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)
        
        if (initialize_file) then 
            ! Current time index to write will be the first and only one 
            n = 1 
        else 
            ! Determine current writing time step 
            n = nc_size(filename,"time",ncid)
            call nc_read(filename,"time",time_prev,start=[n],count=[1],ncid=ncid) 
            if (abs(time-time_prev).gt.1e-5) n = n+1 

            ! Update the time step
            call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
        end if


        call nc_write(filename, "z_bed", isos%now%z_bed, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename,"dz_ss", isos%now%dz_ss, units="m", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid,start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename,"w", isos%now%w, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "we", isos%now%we, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "bsl", bsl, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])

        call nc_write(filename, "dz_ss_ref", isos%ref%dz_ss, units="m", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "H_ice_ref", isos%ref%Hice, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "z_bed_ref", isos%ref%z_bed, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "w_ref", isos%ref%w, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "we_ref", isos%ref%we, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "bsl_ref", bsl_ref, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])

        call nc_write(filename, "log10_eta_eff", log10(isos%domain%eta_eff), units="Pa s", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "He_lith", isos%domain%He_lith, units="km", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "mask_active", isos%domain%maskactive, units="1", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        
        call nc_close(ncid)

    end subroutine isos_restart_write

    subroutine isos_restart_read(isos,filename,time)

        implicit none

        type(isos_class),  intent(INOUT) :: isos
        character(len=*),  intent(IN)    :: filename
        real(wp),          intent(IN)    :: time

        call nc_read(filename, "z_bed_ref", isos%ref%z_bed, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema z_bed_ref: ", minval(isos%ref%z_bed), maxval(isos%ref%z_bed)

        call nc_read(filename, "H_ice_ref", isos%ref%Hice, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema H_ice_ref: ", minval(isos%ref%Hice), maxval(isos%ref%Hice)

        call nc_read(filename, "dz_ss_ref", isos%ref%dz_ss, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema dz_ss_ref: ", minval(isos%ref%dz_ss), maxval(isos%ref%dz_ss)

        call nc_read(filename, "w_ref", isos%ref%w, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema w_ref: ", minval(isos%ref%w), maxval(isos%ref%w)

        call nc_read(filename, "we_ref", isos%ref%we, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema we_ref: ", minval(isos%ref%we), maxval(isos%ref%we)

        call nc_read(filename, "bsl_ref", isos%ref%w_equilibrium, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        isos%ref%bsl = sum(isos%ref%w_equilibrium) / (isos%domain%nx * isos%domain%ny)
        ! write(*,*) "bsl_ref: ", isos%ref%bsl

        ! Read current state
        call nc_read(filename, "z_bed", isos%now%z_bed, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema z_bed: ", minval(isos%now%z_bed), maxval(isos%now%z_bed)

        call nc_read(filename, "dz_ss", isos%now%dz_ss, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema dz_ss: ", minval(isos%now%dz_ss), maxval(isos%now%dz_ss)

        call nc_read(filename, "bsl", isos%now%bsl, start=[1], &
            count=[1])
        ! write(*,*) "bsl: ", isos%now%bsl

        call nc_read(filename, "w", isos%now%w, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema w: ", minval(isos%now%w), maxval(isos%now%w)

        call nc_read(filename, "we", isos%now%we, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema we: ", minval(isos%now%we), maxval(isos%now%we)

        write(*,*) "isos_restart_read:: read in restart file: ", trim(filename)

        return

    end subroutine isos_restart_read

end module isostasy_io