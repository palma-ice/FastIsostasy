module lv_xlra

    use isostasy_defs, only : wp, pi
    use convolutions
    
    implicit none

    private
    
    public :: calc_asthenosphere_relax
    public :: calc_llra_equilibrium

    contains

    ! Calculate rate of change of vertical bedrock height from a relaxing asthenosphere.
    ! Allows laterally variable relaxation times (LV-LLRA and LV-ELRA)
    subroutine calc_asthenosphere_relax(dwdt, w, w_equilibrium, tau)
        implicit none

        real(wp), intent(OUT) :: dwdt(:, :)
        real(wp), intent(IN)  :: w(:, :)
        real(wp), intent(IN)  :: w_equilibrium(:, :)
        real(wp), intent(IN)  :: tau(:, :)

        dwdt = (w_equilibrium-w) / tau
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

end module lv_xlra