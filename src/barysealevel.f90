module barysealevel

    use nml
    use ncio
    use, intrinsic :: iso_fortran_env, only: error_unit
    use isostasy_defs, only : sp, dp, wp, pi
    use isos_utils, only : interp_0d

    implicit none

    type series_type
        character(len=512) :: filename
        real(wp), allocatable :: time(:), var(:), sigma(:)
    end type

    type bsl_class
        character(len=56) :: method       ! "const", "file", "fastiso"
        real(wp) :: bsl_init    ! [m] Initial bsl
        real(wp) :: time        ! [yr] Current time
        real(wp) :: bsl_now     ! [m] Current bsl
        real(wp) :: A_ocean_now ! [m^2] Current ocean surface area
        real(wp) :: sigma_now
        real(wp) :: deltaV_bsl_now       ! [m^3] Volume contribution to BSL
        type(series_type) :: series ! BSL time series from file
        real(wp) :: A_ocean_pd      ! [m^2] PD A_ocean as in Goelzer et al. (2020)
        logical  :: constant_ocean_surface  ! True if ocean surface area is constant
        real(wp), allocatable :: A_ocean_vec(:)  ! [m^2] A_ocean vector for interpolation over bsl
        real(wp), allocatable :: bsl_vec(:)      ! [m] bsl vector for interpolation of A_ocean
    end type

    public :: bsl_class
    public :: bsl_init
    public :: bsl_update

contains

    subroutine bsl_init(bsl, filename)

        implicit none

        type(bsl_class),    intent(INOUT) :: bsl
        character(len=*),   intent(IN)  :: filename

        ! Local variables
        character(len=56) :: varname
        character(len=256):: A_ocean_path
        logical           :: use_nc
        integer           :: n, n_bsl_vec
        real(wp)          :: A_ocean_pd_biased

        write(*,*) "bsl_init:: reading method..."
        call nml_read(filename, "barysealevel", "method", bsl%method)
        
        select case(trim(bsl%method))

        case("const")
            write(*,*) "bsl_init:: using constant BSL."
            call nml_read(filename, "barysealevel", "bsl_init", bsl%bsl_init, init=.TRUE.)

        case("file")
            write(*,*) "bsl_init:: using BSL time series from file."

            ! Determine filename from which to load sea level time series
            call nml_read(filename, "barysealevel", "sl_path", bsl%series%filename, init=.TRUE.)
    
            use_nc = .FALSE. 
            n = len_trim(bsl%series%filename)
            if (bsl%series%filename(n-1:n) .eq. "nc") use_nc = .TRUE. 

            if (use_nc) then 

                ! Get the variable name of interest
                call nml_read(filename, "barysealevel", "sl_name", varname, init=.TRUE.)
                
                ! Read the time series from netcdf file 
                call read_series_nc(bsl%series, bsl%series%filename, varname)

            else

                ! Read the time series from ascii file
                call read_series(bsl%series, bsl%series%filename)

            end if 
    
            ! Also set bsl_init to zero for safety 
            bsl%bsl_init = 0.0_wp

        case("fastiso")

            write(*,*) "bsl_init:: using BSL from FastIsostasy domains."

            ! Load the ocean surface area from parameter file
            call nml_read(filename, "barysealevel", "bsl_init", bsl%bsl_init)
            call nml_read(filename, "barysealevel", "A_ocean_pd", bsl%A_ocean_pd)
            call nml_read(filename, "barysealevel", "A_ocean_path", A_ocean_path)
            if ((A_ocean_path .eq. "None") .or. &
                (A_ocean_path .eq. "none") .or. &
                (A_ocean_path .eq. "no")) then
                bsl%constant_ocean_surface = .TRUE.
                bsl%A_ocean_now = bsl%A_ocean_pd
            else
                bsl%constant_ocean_surface = .FALSE.
                n_bsl_vec = nc_size(A_ocean_path, "z")
                allocate(bsl%bsl_vec(n_bsl_vec))
                allocate(bsl%A_ocean_vec(n_bsl_vec))
                call nc_read(A_ocean_path, "z", bsl%bsl_vec)
                call nc_read(A_ocean_path, "A", bsl%A_ocean_vec)

                call interp_0d(bsl%bsl_vec, bsl%A_ocean_vec, bsl%bsl_now, A_ocean_pd_biased)
                bsl%A_ocean_vec = bsl%A_ocean_vec * bsl%A_ocean_pd / A_ocean_pd_biased
                call interp_0d(bsl%bsl_vec, bsl%A_ocean_vec, bsl%bsl_now, bsl%A_ocean_now)
            end if

            write(*,*) "bsl_init:: You are using the cumulated BSL contributions of FastIsostasy domains."

        case DEFAULT 

            write(error_unit,*) ""
            write(error_unit,*) "bsl_init:: Error: method not recognized."
            write(error_unit,*) "method = ", bsl%method
            stop 

        end select

        return 

    end subroutine bsl_init

    subroutine bsl_update(bsl, year_bp)

        implicit none

        type(bsl_class), intent(INOUT)  :: bsl
        real(wp), intent(IN)            :: year_bp

        select case(trim(bsl%method))
        case("const")
            ! Assign sea-level constant
            bsl%time  = year_bp
            bsl%bsl_now  = bsl%bsl_init
            bsl%sigma_now = 0.0

        case("file")
            bsl%time  = year_bp
            bsl%bsl_now  = series_interp(bsl%series, year_bp)
            bsl%sigma_now = 0.0

        case("fastiso")
            bsl%time  = year_bp
            bsl%sigma_now = 0.0

            if (bsl%constant_ocean_surface) then
                bsl%A_ocean_now = bsl%A_ocean_pd
            else
                call interp_0d(bsl%bsl_vec, bsl%A_ocean_vec, bsl%bsl_now, &
                    bsl%A_ocean_now)
            end if
            
        case DEFAULT

                write(error_unit,*) ""
                write(error_unit,*) "bsl_update:: Error: method not recognized."
                write(error_unit,*) "method = ", bsl%method
                stop 

        end select

        return
    end subroutine bsl_update

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
        ! Wrapper for simple `interp_linear` function
        ! for series_types. 
        implicit none 

        type(series_type) :: series 
        real(wp) :: year_bp 
        real(wp) :: var 
        integer  :: nt, i 

        ! Interpolate series object
        var = interp_linear(series%time,series%var,xout=year_bp)

        return 

    end function series_interp

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

end module barysealevel
