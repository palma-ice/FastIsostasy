module integrators

    use, intrinsic :: iso_fortran_env, only: error_unit
    use isostasy_defs, only : sp, dp, wp, pi, isos_class, ode_class

    implicit none

    public :: step_euler, step_rk4, step_bs32 !, step_tsit54

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
            isos%now%w = isos%now%w + ode%dt * ode%k1
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
            ode%t = ode%t + h
        end do

    end subroutine step_rk4

    subroutine step_bs32(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos
        real(wp) :: err_abs, err_rel

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
            err_abs = maxval(abs(ode%x2 - ode%x1))
            err_rel = err_abs / max(maxval(abs(ode%x2)), 1e-10)
            if (err_abs < isos%par%atol) then !  .and. err_rel < isos%par%rtol
                ! Accept step
                ode%x = ode%x2
                ode%t = ode%t + ode%dt
                ode%dt = min(ode%dt * (isos%par%atol/err_abs)**0.25, isos%par%dt_max)
                ! write(*,*) "Accepted: t = ", ode%t, " dt = ", ode%dt, " err_abs = ", err_abs
            else
                ! Reject step and reduce time step
                ode%dt = ode%dt * 0.5
            end if
        end do
    end subroutine step_bs32

    ! subroutine step_tsit54(ode_rhs, tf, ode, isos)
    !     implicit none
    !     real(wp), intent(in) :: tf
    !     type(ode_class), intent(inout) :: ode
    !     type(isos_class), intent(inout) :: isos
    !     real(wp) :: err_abs, err_rel

    !     real(wp) :: c2, c3, c4, c5, c6
    !     real(wp) :: a21, a31, a32, a41, a42, a43, a51, a52, a53, a54, a61, a62, a63, a64, a65
    !     real(wp) :: b1, b3, b4, b5
    !     real(wp) :: b1s, b3s, b4s, b5s, b6s

    !     c2 = 0.161
    !     c3 = 0.327
    !     c4 = 0.9
    !     c5 = 0.98
    !     c6 = 1.0

    !     a21 = 0.161
    !     a31 = 0.008468277408959993
    !     a32 = 0.3185667612541235
    !     a41 = 0.9
    !     a42 = -0.7262596873842339
    !     a43 = 0.8263596873842339
    !     a51 = 0.98
    !     a52 = 0.0
    !     a53 = 0.0
    !     a54 = 0.0
    !     a61 = 1.0
    !     a62 = 0.0
    !     a63 = 0.0
    !     a64 = 0.0
    !     a65 = 0.0

    !     b1 = 0.09646076681806523
    !     b3 = 0.001295522046139671
    !     b4 = 0.2840332026666311
    !     b5 = 0.6787682725532619

    !     b1s = 0.09171915262261698
    !     b3s = 0.0019686359636073892
    !     b4s = 0.2897153057105852
    !     b5s = 0.6734157376458712
    !     b6s = 0.002206214953488671

    !     do while (ode%t < tf)
    !         ! handle last partial step cleanly
    !         ode%dt = min(ode%dt, tf - ode%t)

    !         ode%k1 = ode_rhs(isos%now%w, ode%t, isos)
    !         ode%k2 = ode_rhs(isos%now%w + a21*ode%dt*ode%k1, ode%t + c2*ode%dt, isos)
    !         ode%k3 = ode_rhs(isos%now%w + a31*ode%dt*ode%k1 + a32*ode%dt*ode%k2, ode%t + c3*ode%dt, isos)
    !         ode%k4 = ode_rhs(isos%now%w + a41*ode%dt*ode%k1 + a42*ode%dt*ode%k2 + a43*ode%dt*ode%k3, ode%t + c4*ode%dt, isos)
    !         ode%k5 = ode_rhs(isos%now%w + a51*ode%dt*ode%k1 + a52*ode%dt*ode%k2 + a53*ode%dt*ode%k3 + a54*ode%dt*ode%k4, ode%t + c5*ode%dt, isos)
    !         ode%k6 = ode_rhs(isos%now%w + a61*ode%dt*ode%k1 + a62*ode%dt*ode%k2 + a63*ode%dt*ode%k3 + a64*ode%dt*ode%k4 + a65*ode%dt*ode%k5, ode%t + c6*ode%dt, isos)
    !         ode%x1 = isos%now%w + ode%dt * (b1*ode%k1 + b3*ode%k3 + b4*ode%k4 + b5*ode%k5)
    !         ode%x2 = isos%now%w + ode%dt * (b1s*ode%k1 + b3s*ode%k3 + b4s*ode%k4 + b5s*ode%k5 + b6s*ode%k6)

    !         isos%now%w = ode%x2
    !         ode%t = ode%t + ode%dt
    !     end do
    ! end subroutine step_tsit54

end module integrators