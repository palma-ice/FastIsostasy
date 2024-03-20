module sealevel

    use isostasy_defs, only : wp, isos_class
    use isos_utils

    implicit none

    public :: calc_bsl_pwconstant_Aocean
    public :: calc_columnanoms_load
    public :: calc_columnanoms_solidearth
    public :: calc_masks
    public :: calc_sl_contribution

    contains

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

        if (isos%par%interactive_sealevel) then
            call add_columnanom(isos%par%rho_seawater, isos%now%rsl, isos%ref%rsl, &
                isos%now%canom_load, isos%now%maskocean)
        end if

        call maskfield(isos%now%canom_load, isos%now%canom_load, isos%domain%maskactive, &
            isos%domain%nx, isos%domain%ny)
    end subroutine calc_columnanoms_load

    !
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

    subroutine calc_masks(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos

        isos%now%maskcontinent = isos%now%rsl < 0

        ! maskgrounded
        isos%now%Haf = isos%now%Hice - isos%now%rsl * (isos%par%rho_seawater / isos%par%rho_ice)
        isos%now%maskgrounded = isos%now%Haf > 0

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

    subroutine calc_sl_contribution(isos)
        implicit none
        type(isos_class), intent(INOUT)   :: isos 
        real(wp) :: V_af
        real(wp) :: V_den
        ! real(wp) :: V_pov

        V_af = sum( (isos%now%Haf - isos%ref%Haf ) * isos%domain%A )
        V_den = sum((isos%now%Hice - isos%ref%Hice) * isos%par%Vden_factor * isos%domain%A)
        isos%now%bsl = (V_af + V_den) / isos%par%A_ocean_pd
        return
    end subroutine calc_sl_contribution

end module sealevel