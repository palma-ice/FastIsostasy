module sealevel

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    use isostasy_defs, only : wp, isos_class, isos_state_class, isos_domain_class, isos_param_class
    use isos_utils
    use ncio 

    implicit none

    public :: calc_z_ss
    public :: calc_rsl
    public :: calc_bsl_constant_Aocean
    public :: calc_bsl_variable_Aocean
    public :: calc_columnanoms_load
    public :: calc_columnanoms_solidearth
    public :: calc_masks
    public :: calc_sl_contributions
    public :: calc_H_above_bsl
    public :: load_bsl_external_ts
    public :: calc_bsl_external_ts

    contains

    subroutine calc_z_ss(z_ss, bsl, z_ss_ref, dz_ss)
        implicit none
        real(wp), intent(INOUT) :: z_ss(:, :)
        real(wp), intent(IN)    :: bsl
        real(wp), intent(IN)    :: z_ss_ref(:, :)
        real(wp), intent(IN)    :: dz_ss(:, :)

        z_ss = bsl + z_ss_ref + dz_ss
        return
    end subroutine calc_z_ss

    ! Calculate the relative sea level.
    subroutine calc_rsl(state)
        implicit none
        type(isos_state_class), intent(INOUT) :: state

        state%rsl = state%z_ss - state%z_bed
        return
    end subroutine calc_rsl

    ! Calculate the barystatic sea level using a constant A_ocean.
    subroutine calc_bsl_constant_Aocean(isos)
        implicit none
        type(isos_class), intent(INOUT)     :: isos

        call calc_bsl(isos%now%bsl, isos%ref%bsl, isos%now%V_af, isos%now%V_den, &
            isos%now%V_pov, isos%par%A_ocean_pd)
        return
    end subroutine calc_bsl_constant_Aocean

    subroutine calc_bsl_variable_Aocean(isos)
        ! Warning: This subroutine is work in progress and should not be used yet.
        ! We need to figure out how to compute the bsl incrementally when A_ocean is
        ! piecewise constant.
        implicit none
        type(isos_class), intent(INOUT)     :: isos

        call interp_0d(isos%domain%bsl_vec, isos%domain%A_ocean_vec, &
            isos%now%bsl, isos%now%A_ocean)

        call calc_bsl(isos%now%bsl, isos%ref%bsl, isos%now%V_af, isos%now%V_den, &
            isos%now%V_pov, isos%now%A_ocean)
        return
    end subroutine calc_bsl_variable_Aocean

    subroutine calc_bsl(bsl, bsl_ref, V_af, V_den, V_pov, A_ocean)
        implicit none
        real(wp), intent(OUT)   :: bsl
        real(wp), intent(IN)    :: bsl_ref
        real(wp), intent(IN)    :: V_af
        real(wp), intent(IN)    :: V_den
        real(wp), intent(IN)    :: V_pov
        real(wp), intent(IN)    :: A_ocean

        bsl = bsl_ref - (V_af + V_den + V_pov) / A_ocean
        return
    end subroutine calc_bsl

    ! Calculate the barystatic sea level using a piecewise constant A_ocean.
    ! This is work in progress and should not be used yet.
    subroutine calc_bsl_pwconstant_Aocean(isos)
        implicit none
        type(isos_class), intent(INOUT)     :: isos

        isos%now%bsl = isos%now%bsl / isos%now%A_ocean
        call interp_0d(isos%domain%bsl_vec, isos%domain%A_ocean_vec, &
            isos%now%bsl, isos%now%A_ocean)
    end subroutine calc_bsl_pwconstant_Aocean

    subroutine calc_columnanoms_load(isos)
        implicit none
        type(isos_class), intent(INOUT)     :: isos 

        isos%now%canom_load(:, :) = 0
        call add_columnanom(isos%par%rho_ice, isos%now%Hice, isos%ref%Hice, isos%now%canom_load)

        ! if (isos%par%interactive_sealevel) then
        !     call add_columnanom(isos%par%rho_seawater, isos%now%rsl, isos%ref%rsl, &
        !         isos%now%canom_load, isos%now%maskocean)
        ! end if

        call maskfield(isos%now%canom_load, isos%now%canom_load, isos%domain%maskactive, &
            isos%domain%nx, isos%domain%ny)
    end subroutine calc_columnanoms_load

    ! Calculate the column anomalies of the solid Earth.
    subroutine calc_columnanoms_solidearth(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_full = isos%now%canom_load
        call add_columnanom(isos%par%rho_litho, isos%now%we, isos%ref%we, isos%now%canom_full)
        call add_columnanom(isos%par%rho_uppermantle, isos%now%w, isos%ref%w, isos%now%canom_full)
        call maskfield(isos%now%canom_full, isos%now%canom_full, isos%domain%maskactive, &
            isos%domain%nx, isos%domain%ny)
        return
    end subroutine calc_columnanoms_solidearth

    !
    subroutine add_columnanom(rho, H_now, H_ref, canom, mask)
        implicit none

        real(wp), intent(IN)    :: rho
        real(wp), intent(IN)    :: H_now(:, :)
        real(wp), intent(IN)    :: H_ref(:, :)
        real(wp), intent(INOUT) :: canom(:, :)
        logical, intent(INOUT), optional  :: mask(:, :)

        integer                 :: nx, ny
        real(wp), allocatable   :: canom_helper(:, :)

        if (present(mask)) then
            nx = size(canom, 1)
            ny = size(canom, 2)
            allocate(canom_helper(nx, ny))
            canom_helper = 0.0_wp

            call maskfield(canom_helper, rho * (H_now - H_ref), mask, nx, ny)
            canom = canom + canom_helper
        else
            canom = canom + rho * (H_now - H_ref)
        end if
        return
    end subroutine add_columnanom

    subroutine calc_mass_anom(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 

        call maskfield(isos%now%mass_anom, isos%domain%A * isos%now%canom_full, &
            isos%domain%maskactive, isos%domain%nx, isos%domain%ny)
        return
    end subroutine calc_mass_anom

    subroutine calc_H_above_bsl(state, par)
        implicit none
        type(isos_state_class), intent(INOUT)   :: state
        type(isos_param_class), intent(IN)      :: par

        state%H_above_bsl = state%Hice - state%bsl * (par%rho_seawater / par%rho_ice)
        return
    end subroutine calc_H_above_bsl

    subroutine calc_masks(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        ! write(*,*) 'Updating continent mask'
        isos%now%maskcontinent = isos%now%z_bed > 0

        ! write(*,*) 'Updating grounded mask'
        isos%now%Haf = isos%now%Hice - isos%now%rsl * (isos%par%rho_seawater / isos%par%rho_ice)
        isos%now%maskgrounded = (isos%now%Haf > 0) .and. (isos%now%Hice > 0)

        ! write(*,*) 'Updating ocean mask'
        call calc_maskocean(isos)

        return
    end subroutine calc_masks

    subroutine calc_maskocean(isos)
        implicit none
        
        type(isos_class), intent(INOUT) :: isos
        integer                         :: i, j

        do i = 1, isos%domain%nx
        do j = 1, isos%domain%ny
            if (isos%now%maskcontinent(i, j) .or. isos%now%maskgrounded(i, j)) then
                isos%now%maskocean(i, j) = .false.
            else
                isos%now%maskocean(i, j) = .true.
            endif
        end do
        end do

    end subroutine calc_maskocean

    subroutine calc_sl_contributions(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 

        isos%now%V_af = sum( (isos%now%H_above_bsl - isos%ref%H_above_bsl ) * isos%domain%A ) * &
            (isos%par%rho_ice / isos%par%rho_seawater)
        isos%now%V_den = sum((isos%now%Hice - isos%ref%Hice) * isos%par%Vden_factor * &
            isos%domain%A)
        isos%now%V_pov = 0.0 ! I think this is not needed when V_af is based on H_above_bsl
        ! write(*,*) isos%now%V_af, isos%now%V_den, isos%now%V_pov
        return
    end subroutine calc_sl_contributions


    ! === Routines for handling external bsl contributions (i.e. from ice sheets not treated within the domain) ===

    subroutine calc_bsl_external_ts(bsl_ext,ts_bsl,ts_time,time)
        ! Determine the current bsl_external value for the current time from a time seris of values.
        ! This is done by linear interpolation, where the time series is already available and
        ! passed as arguments [ts_bsl,ts_time].
        ! For paleo sea-level time series, time is usually given in units of years before present. 

        implicit none
        
        real(wp), intent(OUT) :: bsl_ext
        real(wp), intent(IN)  :: ts_bsl(:)
        real(wp), intent(IN)  :: ts_time(:)
        real(wp), intent(IN)  :: time 
        
        ! Interpolate series object to current time
        bsl_ext = interp_linear(ts_time,ts_bsl,xout=time)

        return

    end subroutine calc_bsl_external_ts

    subroutine load_bsl_external_ts(bsl,time,filename,nc_varname,conv_time_units)

        implicit none

        real(wp), intent(OUT), allocatable :: bsl(:)
        real(wp), intent(OUT), allocatable :: time(:)
        character(len=*), intent(IN) :: filename
        character(len=*), intent(IN), optional :: nc_varname
        real(wp), intent(IN), optional :: conv_time_units 

        ! Local variables
        integer :: n
        logical :: use_nc 

        n = len_trim(filename)
        if (filename(n-1:n) .eq. "nc") then
            use_nc = .TRUE.
        else
            use_nc = .FALSE.
        end if

        if (use_nc) then 

            if (.not. present(nc_varname)) then
                write(error_unit,*) "load_bsl_external_ts:: Error: No variable name provided to be loaded from a NetCDF input file. &
                &Provide the correct nc_varname as an argument to load_bsl_external_ts()."
                write(error_unit,*) "filename: ", trim(filename)
                stop
            end if

            ! Read the time series from netcdf file 
            call read_series_nc(bsl,time,filename,nc_varname)

        else

            ! Read the time series from ascii file
            call read_series(bsl,time,filename)

        end if 
        
        ! Convert time units
        if (present(conv_time_units)) then
            ! Convert time units (e.g. kiloyears to years)
            time = time*conv_time_units 
        end if 

        return

    end subroutine load_bsl_external_ts

    function interp_linear(x,y,xout) result(yout)
        ! Simple linear interpolation of a point

        implicit none 
 
        real(wp), dimension(:), intent(IN) :: x, y
        real(wp), intent(IN) :: xout
        real(wp) :: yout 
        integer  :: i, j, n, nout 
        real(wp) :: alph

        n    = size(x) 

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n 
                if (x(j) .ge. xout) exit 
            end do

            if (j .eq. 1) then 
                yout = y(1) 
            else if (j .eq. n+1) then 
                yout = y(n)
            else 
                alph = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alph*(y(j) - y(j-1))
            end if 
        end if 

        return 

    end function interp_linear
    
    subroutine read_series(var,time,filename)
        ! This subroutine will read a time series of
        ! two columns [time,var] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        real(wp), intent(INOUT), allocatable :: var(:)
        real(wp), intent(INOUT), allocatable :: time(:)
        character(len=*), intent(IN) :: filename 

        ! Local variables
        integer, parameter :: nmax = 10000
        integer :: i, stat, nt, io_unit
        character(len=256) :: str, str1 
        real(wp) :: x(nmax), y(nmax) 

        ! Open file for reading 
        open(newunit=io_unit,file=filename,status="old")

        ! Read the header in the first line: 
        read(io_unit,*,IOSTAT=stat) str

        do i = 1, nmax 
            read(io_unit,'(a100)',IOSTAT=stat) str 

            ! Exit loop if the end-of-file is reached 
            if(IS_IOSTAT_END(stat)) exit 

            str1 = adjustl(trim(str))
!            str1=str
            if ( len(trim(str1)) .gt. 0 ) then 
                if ( .not. (str1(1:1) == "!" .or. &
                            str1(1:1) == "#") ) then 
                    read(str1,*) x(i), y(i) 
                end if
            end if  
        end do 


        ! Close the file
        close(io_unit) 

        if (i .eq. nmax) then 
            write(*,*) "read_series:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data
        nt = i-1 
        
        if (allocated(time))  deallocate(time)
        if (allocated(var))   deallocate(var)

        allocate(time(nt))
        allocate(var(nt))
        
        ! Initialize variables to zero
        time  = 0.0
        var   = 0.0

        ! Store data
        time = x(1:nt) 
        var  = y(1:nt) 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(time), maxval(time)
        write(*,*) "    range var : ",minval(var),  maxval(var)
        
        return 

    end subroutine read_series

    subroutine read_series_nc(var,time,filename,varname)
        ! This subroutine will read a time series of
        ! sea level from a netcdf file.

        implicit none 

        real(wp), intent(INOUT), allocatable :: var(:)
        real(wp), intent(INOUT), allocatable :: time(:)
        character(len=*)  :: filename 
        character(len=*)  :: varname 

        integer :: nt 

        ! Allocate the time series object and store output data
        nt = nc_size(filename,"time")

        if (allocated(time))  deallocate(time)
        if (allocated(var))   deallocate(var)

        allocate(time(nt))
        allocate(var(nt))
        
        ! Initialize variables to zero
        time  = 0.0
        var   = 0.0

        ! Load data from file
        call nc_read(filename,"time",time)
        call nc_read(filename,varname,var)

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(time), maxval(time)
        write(*,*) "    range var : ",minval(var),  maxval(var)
        
        return 

    end subroutine read_series_nc

end module sealevel
