module precision

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp)
    integer,  parameter :: wp = sp 
    
    ! Tolerance settings
    real(wp), parameter :: TOL           = real(1e-5,wp)
    real(wp), parameter :: TOL_UNDERFLOW = real(1e-15,wp)

contains 

end module precision
