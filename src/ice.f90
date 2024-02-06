module ice

    use nml
    use ncio

    use, intrinsic :: iso_fortran_env, only: error_unit
    
    implicit none 

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = dp !mmr sp 

    real(wp), parameter :: sec_year  = 365.0*24.0*60.0*60.0   ! [s/a]
    real(wp), parameter :: pi        = 3.14159265359

    type series_type
        character(len=512) :: filename
        real(wp), allocatable :: time(:), var(:), sigma(:)
    end type

    type ice_class 
        type(series_type) :: series
        integer :: method
        real(wp) :: H_ice_const
        real(wp) :: time, H_ice, sigma
    end type

    private
    public :: ice_class
    public :: ice_init
    public :: ice_update

contains 

    subroutine ice_init(ice,filename)

        implicit none 

        type(ice_class), intent(OUT) :: ice 
        character(len=*),     intent(IN)  :: filename 

        ! Local variables 
        character(len=56) :: varname 
        logical           :: use_nc 
        integer           :: n 

        ! call nml_read(filename,"ice","method",ice%method,init=.TRUE.)
        
        ! select case(ice%method)

        !     case(0)
        !         ! Load constant sea level value from parameter file 

        !         call nml_read(filename,"ice","H_ice_const",ice%H_ice_const,init=.TRUE.)
                

        !     case(1) 
                ! ! Load transient sea-level time series 

                ! ! Determine filename from which to load sea level time series 
                ! call nml_read(filename,"ice","ice_path",ice%series%filename,init=.TRUE.)
        
                ! use_nc = .FALSE. 
                ! n = len_trim(ice%series%filename)
                ! if (ice%series%filename(n-1:n) .eq. "nc") use_nc = .TRUE. 

                ! if (use_nc) then 

                !     ! Get the variable name of interest
                !     call nml_read(filename,"ice","ice_name",varname,init=.TRUE.)
                    
                    ! Read the time series from netcdf file 
                !     call read_series_nc(ice%series,ice%series%filename,varname)
               call read_series_nc(ice%series,ice%series%filename,"T_ice")

                ! else

                !     ! Read the time series from ascii file
                !     call read_series(ice%series,ice%series%filename)

                ! end if 
        
                ! Also set H_ice_const to zero for safety 
                ice%H_ice_const = 0.0_wp 

        !     case DEFAULT 

        !         write(error_unit,*) ""
        !         write(error_unit,*) "ice_init:: Error: method not recognized."
        !         write(error_unit,*) "method = ", ice%method
        !         stop 

        ! end select
        
        return 

    end subroutine ice_init

    subroutine ice_update(ice,year_bp,H_ice_now)

        implicit none 

        type(ice_class), intent(INOUT) :: ice 
        real(wp),             intent(IN)    :: year_bp 
        real(wp), optional,   intent(IN)    :: H_ice_now 

        select case(ice%method)

            case(0)
                ! Assign sea-level constant 
                ice%time  = year_bp 
                ice%H_ice  = ice%H_ice_const
                ice%sigma = 0.0 

            case(1) 
                ice%time  = year_bp 
                ice%H_ice  = series_interp(ice%series,year_bp)
                ice%sigma = 0.0 
            
            case DEFAULT 

                write(error_unit,*) ""
                write(error_unit,*) "ice_update:: Error: method not recognized."
                write(error_unit,*) "method = ", ice%method
                stop 

        end select

        if (present(H_ice_now)) then 
            ! Overwrite any options and use argument value 

            ice%H_ice = H_ice_now 

        end if 

        return 

    end subroutine ice_update

    subroutine read_series(series,filename)
        ! This subroutine will read a time series of
        ! two columns [time,var] from an ascii file.
        ! Header should be commented by "#" or "!"
        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 

        integer, parameter :: f = 190
        integer, parameter :: nmax = 10000

        integer :: i, stat, n 
        character(len=256) :: str, str1 
        real(wp) :: x(nmax), y(nmax) 

        ! Open file for reading 
        open(f,file=filename,status="old")

        ! Read the header in the first line: 
        read(f,*,IOSTAT=stat) str

        do i = 1, nmax 
            read(f,'(a100)',IOSTAT=stat) str 

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
        close(f) 

        if (i .eq. nmax) then 
            write(*,*) "read_series:: warning: "// &
                       "Maximum length of time series reached, ", nmax
            write(*,*) "Time series in the file may be longer: ", trim(filename)
        end if 

        ! Allocate the time series object and store output data
        n = i-1 
        call series_allocate(series,n)

        series%time = x(1:n) 
        series%var  = y(1:n) 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var),  maxval(series%var)
        
        return 

    end subroutine read_series

    subroutine read_series_nc(series,filename,varname)
        ! This subroutine will read a time series of
        ! sea level from a netcdf file.

        implicit none 

        type(series_type) :: series 
        character(len=*)  :: filename 
        character(len=*)  :: varname 

        integer :: n 

        ! Allocate the time series object and store output data
        n = nc_size(filename,"time")

        call series_allocate(series,n)

        call nc_read(filename,"time",series%time)
        call nc_read(filename,varname,series%var)

        ! Convert from kiloyears to years
        ! ajr: this should not be hard coded!! 
        series%time = series%time*1e3 

        write(*,*) "read_series:: Time series read from file: "//trim(filename)
        write(*,*) "    range time: ",minval(series%time), maxval(series%time)
        write(*,*) "    range var : ",minval(series%var),  maxval(series%var)
        
        return 

    end subroutine read_series_nc

    function series_interp(series,year_bp) result(var)
        ! Wrapper for simple `interp_2d_over_time` function
        ! for series_types. 
        implicit none 

        type(series_type) :: series 
        real(wp) :: year_bp 
        real(wp) :: var 
        integer  :: nt, i 

        ! Interpolate series object
        var = interp_2d_over_time(series%time,series%var,xout=year_bp)

        return 

    end function series_interp

    function interp_2d_over_time(x,y,xout) result(yout)
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

    end function interp_2d_over_time
    
    subroutine series_allocate(series,nt)

        implicit none 

        type(series_type) :: series 
        integer :: nt 

        if (allocated(series%time))  deallocate(series%time)
        if (allocated(series%var))   deallocate(series%var)
        if (allocated(series%sigma)) deallocate(series%sigma)

        allocate(series%time(nt))
        allocate(series%var(nt))
        allocate(series%sigma(nt))
        
        ! Initialize variables to zero
        series%time  = 0.0
        series%var   = 0.0
        series%sigma = 0.0
        
        return 

    end subroutine series_allocate

end module ice
