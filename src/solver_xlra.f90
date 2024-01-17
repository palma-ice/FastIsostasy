module solver_xlra

    use isostasy_defs, only : wp, pi

    implicit none

    private
    
    public :: calc_litho_local
    
    contains


    elemental subroutine calc_asthenosphere_relax(dzbdt,z_bed,z_bed_ref,w_b,tau)
        ! Calculate rate of change of vertical bedrock height
        ! from a relaxing asthenosphere.

        implicit none

        real(wp), intent(OUT) :: dzbdt 
        real(wp), intent(IN)  :: z_bed 
        real(wp), intent(IN)  :: z_bed_ref
        real(wp), intent(IN)  :: w_b        ! w_b = w1-w0
        real(wp), intent(IN)  :: tau

        dzbdt = -((z_bed-z_bed_ref) + w_b) / tau
        
        return

    end subroutine calc_asthenosphere_relax


    
    elemental subroutine calc_litho_local(w,q,z_bed,H_ice,z_sl,rho_ice,rho_seawater,rho_uppermantle,g)
        ! Calculate the local lithospheric loading from ice or ocean weight 
        ! in units of [Pa] and local equilibrium displacement w [m].

        implicit none 

        real(wp), intent(OUT) :: w
        real(wp), intent(OUT) :: q
        real(wp), intent(IN)  :: z_bed
        real(wp), intent(IN)  :: H_ice
        real(wp), intent(IN)  :: z_sl 
        real(wp), intent(IN)  :: rho_ice
        real(wp), intent(IN)  :: rho_seawater
        real(wp), intent(IN)  :: rho_uppermantle
        real(wp), intent(IN)  :: g 

        if (rho_ice*H_ice.ge.rho_seawater*(z_sl-z_bed)) then   
           
            ! Ice or land

            q = rho_ice*g*H_ice

        else
            ! Ocean

           q = rho_seawater*g*(z_sl-z_bed) 

        end if

        ! Scale to get local displacement given the load q
        w = q / (rho_uppermantle*g) 

        return 

    end subroutine calc_litho_local

end module solver_xlra