module isostasy_io

    use isostasy_defs, only : wp, isos_class, isos_domain_class
    use nml
    use ncio

    implicit none

    private
    
    public :: isos_grid_write
    public :: isos_restart_write
    public :: isos_restart_read
    public :: isos_write_init
    public :: isos_write_init_extended
    public :: isos_write_step
    public :: isos_write_step_extended
    
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
        character(len=2) :: xnm, ynm 

        initialize_file = .TRUE. 
        if (present(init)) initialize_file = init

        xnm = "xc"
        ynm = "yc"

        nx    = isos%domain%nx
        ny    = isos%domain%ny

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
            n = nc_time_index(filename,"time",time,ncid)
            
            ! Update the time step
            call nc_write(filename,"time",time,dim1="time",start=[n],count=[1],ncid=ncid)
        end if


        call nc_write(filename, "z_bed", isos%now%z_bed, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename,"z_ss", isos%now%z_ss, units="m", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid,start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename,"dz_ss", isos%now%dz_ss, units="m", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid,start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename,"w", isos%now%w, units="m", dim1="xc", dim2="yc", &
            dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "we", isos%now%we, units="m", dim1="xc", dim2="yc", &
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

        call nc_write(filename, "log10_eta_eff", log10(isos%domain%eta_eff), units="Pa s", &
            dim1="xc", dim2="yc", dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
        call nc_write(filename, "He_lith", isos%domain%He_lith, units="km", dim1="xc", &
            dim2="yc", dim3="time", ncid=ncid, start=[1,1,n], count=[nx,ny,1])
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

        ! Read current state
        call nc_read(filename, "z_bed", isos%now%z_bed, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema z_bed: ", minval(isos%now%z_bed), maxval(isos%now%z_bed)

        call nc_read(filename, "z_ss", isos%now%z_ss, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema z_ss: ", minval(isos%now%z_ss), maxval(isos%now%z_ss)

        call nc_read(filename, "dz_ss", isos%now%dz_ss, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema dz_ss: ", minval(isos%now%dz_ss), maxval(isos%now%dz_ss)

        call nc_read(filename, "w", isos%now%w, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema w: ", minval(isos%now%w), maxval(isos%now%w)

        call nc_read(filename, "we", isos%now%we, start=[1,1,1], &
            count=[isos%domain%nx, isos%domain%ny, 1])
        ! write(*,*) "Extrema we: ", minval(isos%now%we), maxval(isos%now%we)

        write(*,*) "isos_restart_read:: read in restart file: ", trim(filename)

        return

    end subroutine isos_restart_read


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

        call nc_write(filename, "D_lith", isos%domain%D_lith, units="N m", &
            long_name="Lithosphere flexural rigidity", dim1="xc", dim2="yc", start=[1, 1])

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

            call nc_write(filename, "kappa", isos%domain%kappa, units="", &
                long_name="Pseudodifferential operator in Fourier space", &
                dim1="xc", dim2="yc", start=[1, 1])
        end if

        call nc_write(filename, "maskactive", isos%domain%maskactive, units="1", &
            long_name="Active mask", dim1="xc", dim2="yc", start=[1, 1])
        return

    end subroutine isos_write_init_extended

    ! Write results to file
    subroutine isos_write_step(isos, filename, time, H_ice, z_bed_bench)

        implicit none 
        
        type(isos_class), intent(IN) :: isos
        character(len=*), intent(IN) :: filename
        real(wp),         intent(IN) :: time
        real(wp),         intent(IN) :: H_ice(:, :)
        real(wp),         intent(IN), optional :: z_bed_bench(:, :)

        ! Local variables
        integer  :: ncid, n 

        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time", time, dim1="time", start=[n], count=[1], ncid=ncid)

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

        ! Open the file for writing
        call nc_open(filename, ncid, writable=.TRUE.)

        ! Determine current writing time step 
        n = nc_time_index(filename,"time",time,ncid)

        ! Update the time step
        call nc_write(filename,"time", time, dim1="time", start=[n], count=[1], ncid=ncid)
        
        ! Write variables
        call nc_write(filename, "H_ice", isos%now%Hice, units="m", long_name="Ice thickness", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "z_ss", isos%now%z_ss, units="m", long_name="Sea-surface height", &
              dim1="xc", dim2="yc", dim3="time", start=[1, 1, n], ncid=ncid)

        call nc_write(filename, "dz_ss", isos%now%dz_ss, units="m", &
            long_name="Sea-surface height change", dim1="xc", dim2="yc", dim3="time", &
            start=[1, 1, n], ncid=ncid)

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

end module isostasy_io