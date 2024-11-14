module isostasy_io

    use isostasy_defs, only : wp, isos_class, isos_domain_class
    use nml
    use ncio

    implicit none

    private
    
    public :: isos_grid_write
    public :: isos_restart_write
    
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


end module isostasy_io