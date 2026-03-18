module integrators

    use, intrinsic :: iso_fortran_env, only: error_unit
    use isostasy_defs, only : sp, dp, wp, pi, isos_class, ode_class

    implicit none

    public :: step_euler, step_rk4, step_bs32, step_tsit54

contains

    subroutine step_euler(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos

        interface
            function ode_rhs(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs
        end interface

        do while (ode%t < tf)
            ! handle last partial step cleanly
            ode%dt = min(ode%dt, tf - ode%t)

            ode%k1 = ode_rhs(isos%now%w, ode%t, isos)
            ode%x = ode%x + ode%dt * ode%k1
            isos%now%w = ode%x
            ode%t = ode%t + ode%dt
        end do

    end subroutine step_euler

    subroutine step_rk4(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos
        real(wp) :: t, h

        interface
            function ode_rhs(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs
        end interface

        do while (ode%t < tf)
            ! handle last partial step cleanly
            h = min(ode%dt, tf - ode%t)

            ode%k1 = ode_rhs(ode%x, ode%t, isos)
            ode%y1 = ode%x + 0.5*h*ode%k1
            ode%k2 = ode_rhs(ode%y1, ode%t + 0.5*h, isos)
            ode%y2 = ode%x + 0.5*h*ode%k2
            ode%k3 = ode_rhs(ode%y2, ode%t + 0.5*h, isos)
            ode%y3 = ode%x + h*ode%k3
            ode%k4 = ode_rhs(ode%y3, ode%t + h, isos)
            ode%x = ode%x + (h/6)*(ode%k1 + 2*ode%k2 + 2*ode%k3 + ode%k4)
            isos%now%w = ode%x
            ode%t = ode%t + h
        end do

    end subroutine step_rk4

    subroutine step_bs32(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos
        real(wp) :: error, tol

        interface
            function ode_rhs(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs
        end interface

        do while (ode%t < tf)
            ! handle last partial step cleanly
            ode%dt = min(ode%dt, tf - ode%t)
            ode%dt = min(ode%dt, isos%par%dt_max)

            ode%k1 = ode_rhs(ode%x, ode%t, isos)
            ode%y1 = ode%x + 0.5*ode%dt*ode%k1
            ode%k2 = ode_rhs(ode%y1, ode%t + 0.5*ode%dt, isos)
            ode%y2 = ode%x + 0.75*ode%dt*ode%k2
            ode%k3 = ode_rhs(ode%y2, ode%t + 0.75*ode%dt, isos)
            ode%y3 = ode%x + ode%dt/9 * (2*ode%k1 + 3*ode%k2 + 4*ode%k3)
            ode%k4 = ode_rhs(ode%y3, ode%t + ode%dt, isos)

            ode%x1 = ode%x + ode%dt/9 * (2*ode%k1 + 3*ode%k2 + 4*ode%k3)
            ode%x2 = ode%x + ode%dt * (7/24*ode%k1 + 1/4*ode%k2 + 1/3*ode%k3 + 3*ode%k4)

            ! Estimate error and adjust time step
            tol = isos%par%atol + isos%par%rtol * max(maxval(abs(ode%x1)), maxval(abs(ode%x2)))
            error = sqrt(sum(((ode%x2 - ode%x1) / tol)**2) / size(ode%x))
            if (error < 1) then !  .and. err_rel < isos%par%rtol
                ! Accept step
                ode%x = ode%x2
                isos%now%w = ode%x
                ode%t = ode%t + ode%dt
                ode%dt = min(ode%dt * 1/error, isos%par%dt_max)
                write(*,*) "FastIso: t =", ode%t, " dt =", ode%dt, " error =", error
            else
                ! Reject step and reduce time step
                ode%dt = ode%dt * 0.5
            end if
        end do
    end subroutine step_bs32

    subroutine step_tsit54(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos

        interface
            function ode_rhs(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs
        end interface

        real(wp) :: error, tol
        ! Coefficients for the Tsitouras 5(4) method
        real(wp), parameter :: a(6, 6) = reshape([ &
            0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            0.2_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            0.075_wp, 0.225_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            0.9777_wp, -3.7333_wp, 3.5555_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
            2.9525_wp, -11.5957_wp, 9.8228_wp, -0.2908_wp, 0.0_wp, 0.0_wp, &
            2.8462_wp, -10.7575_wp, 8.9064_wp, 0.2784_wp, -0.2735_wp, 0.0_wp], &
            [6, 6])
        real(wp), parameter :: b(6) = [0.0845_wp, 0.0_wp, 0.0_wp, 0.5179_wp, 0.1276_wp, 0.2700_wp]
        real(wp), parameter :: bh(6) = [0.0845_wp, 0.0_wp, 0.0_wp, 0.5179_wp, 0.1276_wp, 0.2700_wp - 1/5.0_wp]
        real(wp), parameter :: c(6) = [0.0_wp, 0.2_wp, 0.3_wp, 0.8_wp, 8.0_wp/9.0_wp, 1.0_wp]

        do while (ode%t < tf)
            ! handle last partial step cleanly
            ode%dt = min(ode%dt, tf - ode%t)

            ! Compute the stages
            ode%k1 = ode_rhs(ode%x, ode%t, isos)
            ode%y1 = ode%x + a(2, 1)*ode%dt*ode%k1
            ode%k2 = ode_rhs(ode%y1, ode%t + c(2)*ode%dt, isos)
            ode%y2 = ode%x + a(3, 1)*ode%dt*ode%k1 + a(3, 2)*ode%dt*ode%k2
            ode%k3 = ode_rhs(ode%y2, ode%t + c(3)*ode%dt, isos)
            ode%y3 = ode%x + a(4, 1)*ode%dt*ode%k1 + a(4, 2)*ode%dt*ode%k2 + a(4, 3)*ode%dt*ode%k3
            ode%k4 = ode_rhs(ode%y3, ode%t + c(4)*ode%dt, isos)
            ode%y4 = ode%x + a(5, 1)*ode%dt*ode%k1 + a(5, 2)*ode%dt*ode%k2 + a(5, 3)*ode%dt*ode%k3 + a(5, 4)*ode%dt*ode%k4
            ode%k5 = ode_rhs(ode%y4, ode%t + c(5)*ode%dt, isos)
            ode%y5 = ode%x + a(6, 1)*ode%dt*ode%k1 + a(6, 2)*ode%dt*ode%k2 + a(6, 3)*ode%dt*ode%k3 + a(6, 4)*ode%dt*ode%k4 + a(6, 5)*ode%dt*ode%k5
            ode%k6 = ode_rhs(ode%y5, ode%t + c(6)*ode%dt, isos)

            ! Compute the new solution and error estimate
            ode%x1 = ode%x + ode%dt * bh(1)*ode%k1 + ode%dt * bh(2)*ode%k2 + ode%dt * bh(3)*ode%k3 + &
                ode%dt * bh(4)*ode%k4 + ode%dt * bh(5)*ode%k5 + ode%dt * bh(6)*ode%k6
            ode%x2 = ode%x + ode%dt * b(1)*ode%k1 + ode%dt * b(2)*ode%k2 + ode%dt * b(3)*ode%k3 + &
                ode%dt * b(4)*ode%k4 + ode%dt * b(5)*ode%k5 + ode%dt * b(6)*ode%k6

            ! Estimate error and adjust time step
            tol = isos%par%atol + isos%par%rtol * max(maxval(abs(ode%x1)), maxval(abs(ode%x2)))
            error = sqrt(sum(((ode%x2 - ode%x1) / tol)**2) / size(ode%x))
            if (error < 1) then
                ! Accept step
                ode%x = ode%x2
                isos%now%w = ode%x
                ode%t = ode%t + ode%dt
                ode%dt = min(ode%dt * 1/error, isos%par%dt_max)
                write(*,*) "FastIso: t =", ode%t, " dt =", ode%dt, " error =", error
            else
                ! Reject step and reduce time step
                ode%dt = ode%dt * 0.5
            end if
        end do

        
    end subroutine step_tsit54

end module integrators