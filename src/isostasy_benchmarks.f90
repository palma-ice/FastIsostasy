module isostasy_benchmarks

    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit

    implicit none

    ! Internal constants
    integer,  parameter :: dp  = kind(1.d0)
    integer,  parameter :: sp  = kind(1.0)

    ! Choose the precision of the library (sp,dp) - should match isostasy definitions
    integer,  parameter :: wp = sp 

    real(wp), parameter :: pi = 3.14159265359



contains

    ! subroutine analytical_elva_disk_init()

    !     implicit none
   
    !     call calc_analytical_asthenosphere_viscous_disk_params(kappa_mod,dist2c,r, &
    !                                                     lr,kappa_min,kappa_max,dk,nx,ny,dx)                                    
    
    !     r0     = 1000.0e3 ! [m] recheck - include into routine?
    !     h0     = 1000.0   ! [m] 
    !     eta    = 1.e+21   ! [Pa s]
        
    !     ! Initialize parameters for the analytical integrand
    !     write(*,*) "    range(kappa_mod): ", minval(kappa_mod), maxval(kappa_mod)    
    !     call initialize_analytical_integrand(ana,r0,h0,D_lith_const,eta)


    !     return

    ! end subroutine analytical_elva_disk_init




end module isostasy_benchmarks