module sealevel

    use isostasy_defs, only : wp, isos_class, isos_state_class, isos_domain_class, isos_param_class
    use isos_utils

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

end module sealevel