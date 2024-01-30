module solver_xlra

    use isostasy_defs, only : wp, pi
    use convolutions
    
    implicit none

    private
    
    public :: calc_asthenosphere_relax
    public :: calc_litho_local
    public :: calc_litho_regional
    
    contains

    ! TODO: this can be easily extended to LV-ELRA
    subroutine calc_asthenosphere_relax(dzbdt, z_bed, z_bed_ref, w_b, tau)
        ! Calculate rate of change of vertical bedrock height
        ! from a relaxing asthenosphere.
        implicit none

        real(wp), intent(OUT) :: dzbdt(:, :)
        real(wp), intent(IN)  :: z_bed(:, :)
        real(wp), intent(IN)  :: z_bed_ref(:, :)
        real(wp), intent(IN)  :: w_b(:, :)        ! w_b = w1-w0
        real(wp), intent(IN)  :: tau(:, :)

        dzbdt = -((z_bed-z_bed_ref) + w_b) / tau
        return
    end subroutine calc_asthenosphere_relax

    subroutine calc_litho_local(w, q, z_bed, H_ice, z_sl, rho_ice, rho_seawater, &
        rho_uppermantle, g)
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


    subroutine calc_litho_regional(w, q1, z_bed, H_ice, z_sl, GG, rho_ice, &
        rho_seawater, rho_uppermantle, rho_litho, g)
        ! Calculate the load on the lithosphere as
        ! distributed on an elastic plate. 

        implicit none

        real(wp), intent(INOUT) :: w(:,:)      ! [m] Lithospheric displacement
        real(wp), intent(INOUT) :: q1(:,:)      ! [Pa] Lithospheric load
        real(wp), intent(IN)    :: z_bed(:,:)   ! [m] Bed elevation
        real(wp), intent(IN)    :: H_ice(:,:)   ! [m] Ice thickness 
        real(wp), intent(IN)    :: z_sl(:,:)    ! [m] Sea level 
        real(wp), intent(IN)    :: GG(:,:)      ! Regional filter function
        real(wp), intent(IN)    :: rho_ice
        real(wp), intent(IN)    :: rho_seawater
        real(wp), intent(IN)    :: rho_uppermantle
        real(wp), intent(IN)    :: rho_litho
        real(wp), intent(IN)    :: g 

        ! Calculate local lithospheric load and displacement first
        ! TODO: replace with new load
        ! call calc_litho_local(w, q1, z_bed, H_ice, z_sl, rho_ice, rho_seawater, &
        !     rho_uppermantle, g)
        
        ! Convolve the estimated point load with the regional
        ! filter to obtain the distributed load w. 

        ! recheck for stability
        ! call convolve_load_elastic_plate_fft(w,q1,GG) ! recheck convol
        call convolve_load_elastic_plate(w, q1, GG)

        return

    end subroutine calc_litho_regional

end module solver_xlra