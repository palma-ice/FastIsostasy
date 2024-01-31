module solver_xlra

    use isostasy_defs, only : wp, pi
    use convolutions
    
    implicit none

    private
    
    public :: calc_asthenosphere_relax
    public :: calc_llra_equilibrium
    public :: calc_elra_equilibrium

    contains

    ! Calculate rate of change of vertical bedrock height from a relaxing asthenosphere.
    ! Allows laterally variable relaxation times (LV-LLRA and LV-ELRA)
    subroutine calc_asthenosphere_relax(dzbdt, w, w_equilibrium, tau)
        implicit none

        real(wp), intent(OUT) :: dzbdt(:, :)
        real(wp), intent(IN)  :: w(:, :)
        real(wp), intent(IN)  :: w_equilibrium(:, :)
        real(wp), intent(IN)  :: tau(:, :)

        dzbdt = (w_equilibrium-w) / tau
        return
    end subroutine calc_asthenosphere_relax

    ! Calculate the equilibrium displacement of LLRA
    subroutine calc_llra_equilibrium(w_equilibrium, canom_load, rho_uppermantle)
        implicit none

        real(wp), intent(INOUT) :: w_equilibrium(:, :)  ! [m] Lithospheric displacement
        real(wp), intent(INOUT) :: canom_load(:, :)     ! [Pa] Lithospheric load
        real(wp), intent(IN)    :: rho_uppermantle

        w_equilibrium = canom_load / rho_uppermantle
        
        return
    end subroutine calc_llra_equilibrium

    ! Calculate the equilibrium displacement of ELRA
    subroutine calc_elra_equilibrium(w_equilibrium, canom_full, GV, g)
        implicit none

        real(wp), intent(INOUT) :: w_equilibrium(:, :)  ! [m] Lithospheric displacement
        real(wp), intent(INOUT) :: canom_full(:, :)     ! [Pa] Lithospheric load
        real(wp), intent(IN)    :: GV(:, :)             ! Regional filter function
        real(wp), intent(IN)    :: g 

        ! Convolve the estimated point load with the regional
        ! filter to obtain the distributed load w. 

        ! recheck for stability
        ! call convolve_load_elastic_plate_fft(w,q1,GV) ! recheck convol
        ! TODO: use new convolution method here
        
        return

    end subroutine calc_elra_equilibrium

end module solver_xlra