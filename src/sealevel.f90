module sealevel

    use isostasy_defs, only : wp, isos_class, isos_state_class, isos_domain_class, isos_param_class
    use isos_utils

    implicit none

    public :: calc_z_ss
    public :: calc_columnanoms_ice
    public :: calc_columnanoms_seawater
    public :: calc_columnanoms_lithosphere
    public :: calc_columnanoms_mantle
    public :: calc_columnanoms_load
    public :: calc_columnanoms_full
    public :: calc_mass_anom
    public :: calc_Haf
    public :: calc_masks
    public :: calc_maskcontinent
    public :: cal_maskgrounded
    public :: calc_maskocean
    public :: calc_sl_contributions

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

    subroutine calc_columnanoms_load(isos)
        implicit none
        type(isos_class), intent(INOUT)     :: isos 

        isos%now%canom_load = isos%now%canom_ice

        if (isos%par%interactive_sealevel) then
            isos%now%canom_load = isos%now%canom_load + isos%now%canom_seawater
        end if

        call maskfield(isos%now%canom_load, isos%now%canom_load, isos%domain%maskactive, &
            isos%domain%nx, isos%domain%ny)
    end subroutine calc_columnanoms_load

    ! Calculate the column anomalies of the solid Earth.
    subroutine calc_columnanoms_full(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_full = isos%now%canom_load
        isos%now%canom_full = isos%now%canom_full + isos%now%canom_litho + isos%now%canom_mantle
        call maskfield(isos%now%canom_full, isos%now%canom_full, isos%domain%maskactive, &
            isos%domain%nx, isos%domain%ny)
        return
    end subroutine calc_columnanoms_full

    subroutine calc_columnanoms_ice(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_ice = 0.0_wp
        call add_columnanom(isos%par%rho_ice, isos%now%Haf, isos%ref%Haf, isos%now%canom_ice)
    end subroutine calc_columnanoms_ice

    subroutine calc_columnanoms_seawater(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_seawater = 0.0_wp
        call add_columnanom(isos%par%rho_water, isos%now%z_ss, isos%ref%z_ss, &
            isos%now%canom_seawater, isos%now%maskocean)
    end subroutine calc_columnanoms_seawater

    subroutine calc_columnanoms_lithosphere(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_litho = 0.0_wp
        call add_columnanom(isos%par%rho_litho, isos%now%we, isos%ref%we, isos%now%canom_litho)
    end subroutine calc_columnanoms_lithosphere

    subroutine calc_columnanoms_mantle(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%canom_mantle = 0.0_wp
        call add_columnanom(isos%par%rho_uppermantle, isos%now%w, isos%ref%w, &
            isos%now%canom_mantle)
    end subroutine calc_columnanoms_mantle

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

    subroutine calc_Haf(state, par)
        implicit none
        type(isos_state_class), intent(INOUT) :: state
        type(isos_param_class), intent(IN)    :: par

        ! 1. Compute ice column thickness that corresponds to the current sea surface height
        ! (must be larger than 0)
        state%Haf = (state%z_ss - state%z_bed) * &
            (par%rho_seawater / par%rho_ice)
        where(state%Haf < 0.0_wp) state%Haf = 0.0_wp
        ! 2. Compute the difference with the actual ice column thickness
        state%Haf = state%Hice - state%Haf
        where(state%Haf < 0.0_wp) state%Haf = 0.0_wp
    end subroutine calc_Haf

    subroutine calc_masks(state)
        implicit none
        type(isos_state_class), intent(INOUT)   :: state

        call calc_maskcontinent(state)
        call cal_maskgrounded(state)
        call calc_maskocean(state)
        return
    end subroutine calc_masks

    subroutine calc_maskcontinent(state)
        implicit none
        type(isos_state_class), intent(INOUT)   :: state
        state%maskcontinent = state%z_bed > state%z_ss
    end subroutine calc_maskcontinent

    subroutine cal_maskgrounded(state)
        implicit none
        type(isos_state_class), intent(INOUT)   :: state
        state%maskgrounded = (state%Hice > 0.0_wp) .and. (state%Haf > 0.0_wp)
    end subroutine cal_maskgrounded

    subroutine calc_maskocean(state)
        implicit none
        
        type(isos_state_class), intent(INOUT)   :: state
        integer                         :: i, j, nx, ny
        nx = size(state%maskocean, 1)
        ny = size(state%maskocean, 2)

        do i = 1, nx
        do j = 1, ny
            if (state%maskcontinent(i, j) .or. state%maskgrounded(i, j)) then
                state%maskocean(i, j) = .false.
            else
                state%maskocean(i, j) = .true.
            endif
        end do
        end do

    end subroutine calc_maskocean

    subroutine calc_sl_contributions(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos
        real(wp)                          :: V_old

        V_old = isos%now%V_af + isos%now%V_den + isos%now%V_pov
        isos%now%V_af = sum( (isos%now%Haf - isos%ref%Haf ) * &
            isos%domain%A ) * (isos%par%rho_ice / isos%par%rho_seawater)
        isos%now%V_den = sum((isos%now%Hice - isos%ref%Hice) * isos%par%Vden_factor * &
            isos%domain%A)
        isos%now%V_pov = 0.0    ! I think this is tricky because sum(displacement) = 0 over
                                ! the whole Earth and V_pov depends on the exact pattern (displacement over land or sea?)

        isos%now%deltaV_bsl = (isos%now%V_af + isos%now%V_den + isos%now%V_pov) - V_old

        return
    end subroutine calc_sl_contributions

end module sealevel