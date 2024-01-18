module sea_level
    use isostasy_defs, only : wp, isos_param_class, isos_state_class, isos_domain_class, isos_class

    ! implicit none

    private
    public :: calc_sealevel

    contains

    ! calc_sealevel()
    ! Update the sea level based on the new ice thickness field
    subroutine calc_sealevel(Hice, isos)
        implicit none
        real(wp), intent(IN)    :: Hice(:,:)
        type(isos_class), intent(INOUT) :: isos

        call calc_columnanoms_load(Hice, isos)  ! Part 1

        if (((isos%now%t - isos%ref%t) / isos%par%dt_diagnostics) .ge. isos%now%count_updates) then
            call calc_seasurfaceheight(isos)    ! Part 2
            call calc_masks(isos)               ! Part 3
            call calc_sl_contribution(isos)     ! Part 4
            isos%now%count_updates = isos%now%count_updates + 1
        endif

        return
    end subroutine calc_sealevel

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! Part 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_columnanoms_load(Hice, isos)
        implicit none
        real(wp), intent(IN)            :: Hice(:,:)
        type(isos_class), intent(INOUT)   :: isos 

        isos%now%Hice = Hice
        isos%now%Hsw = (isos%now%ssh - isos%now%z_bed) * isos%now%maskocean
        isos%now%canom_load(:, :) = 0
        call add_columnanom(isos%par%rho_ice, isos%now%Hice, isos%ref%Hice, isos%now%canom_load)
        call add_columnanom(isos%par%rho_seawater, isos%now%Hsw, isos%ref%Hsw, isos%now%canom_load)
    end subroutine calc_columnanoms_load

    !
    subroutine add_columnanom(rho, H_now, H_ref, canom)
        implicit none

        real(wp), intent(IN)    :: rho
        real(wp), intent(IN)    :: H_now(:,:)
        real(wp), intent(IN)    :: H_ref(:,:)
        real(wp), intent(INOUT)   :: canom(:,:)

        canom = canom + rho * (H_now - H_ref)
        return
    end subroutine add_columnanom

    !
    subroutine calc_columnanom_solidearth(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 
        isos%now%canom_full = isos%now%canom_load
        call add_columnanom(isos%par%rho_litho, isos%now%ue, isos%ref%ue, isos%now%canom_full)
        call add_columnanom(isos%par%rho_uppermantle, isos%now%u, isos%ref%u, isos%now%canom_full)
    end subroutine calc_columnanom_solidearth

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! Part 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_seasurfaceheight(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 
        call calc_ssh_perturbation(isos)
        isos%now%ssh = isos%ref%ssh + isos%now%ssh_perturb + isos%now%bsl
    end subroutine calc_seasurfaceheight


    subroutine calc_ssh_perturbation(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        call calc_mass_anom(isos)
        !call calc_fft_convo()
    end subroutine calc_ssh_perturbation

    !
    subroutine calc_mass_anom(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 
        isos%now%mass_anom = isos%domain%A * ( isos%now%canom_full &
            - isos%par%rho_seawater * isos%now%bsl * isos%now%maskocean * isos%ref%maskactive )
    end subroutine calc_mass_anom

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! Part 3 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_masks(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos
        call calc_maskgrounded(isos)
        call calc_maskocean(isos)
    end subroutine calc_masks

    subroutine calc_maskocean(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos
        isos%now%maskocean = ((isos%now%ssh - isos%now%z_bed) > 0) * (1 - isos%now%maskgrounded)
    end subroutine calc_maskocean

    !
    subroutine calc_maskgrounded(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos
        call calc_height_above_floatation(isos)
        isos%now%maskgrounded = isos%now%Haf > 0
    end subroutine calc_maskgrounded

    !
    subroutine calc_height_above_floatation(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%Haf = isos%now%Hice + min(isos%now%z_bed - isos%now%ssh, 0) *
            (isos%par%rho_sw / isos%par%rho_ice)
    end subroutine calc_height_above_floatation


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!! Part 4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine calc_sl_contribution(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 
        real(wp) :: V_af
        real(wp) :: density_factor
        real(wp) :: V_den
        ! real(wp) :: V_pov

        V_af = sum( (isos%now%Haf - isos%ref%Haf ) * isos%domain%A )
        density_factor = isos%par%rho_ice / isos%par%rho_water
            - isos%par%rho_ice / isos%par%rho_sw
        V_den = sum( (isos%now%Hice - isos%ref%Hice) * density_factor * isos%domain%A )
        isos%now%bsl = (V_af + V_den) / isos%par%A_ocean_pd
    end subroutine calc_sl_contribution

end module sea_level