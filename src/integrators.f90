module integrators

    use, intrinsic :: iso_fortran_env, only: error_unit
    use isostasy_defs, only : sp, dp, wp, pi, isos_class, ode_class

    implicit none

    public :: step_euler, step_rk4, step_bs32, step_tsit54, step_rkc
    private :: spectral_radius_estimate, rkc_choose_stages, rkc_stability_boundary, &
        rkc_scalar_map, rkc_coeffs, chebyshev_table

contains

    subroutine step_euler(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

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

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

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

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

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
                ode%dt = min(ode%dt * (1/error) ** (1/isos%par%q), isos%par%dt_max)
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

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

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
                ode%dt = min(ode%dt * (1/error) ** (1/isos%par%q), isos%par%dt_max)
                write(*,*) "FastIso: t =", ode%t, " dt =", ode%dt, " error =", error
            else
                ! Reject step and reduce time step
                ode%dt = ode%dt * 0.5
            end if
        end do

        
    end subroutine step_tsit54

    ! =========================================================================
    ! step_rkc: stabilised explicit Runge-Kutta-Chebyshev method (RKC2, damped,
    ! s stages), ported from FastIsostasy.jl's RKCIntegrator (Sommeijer-Shampine-
    ! Verwer 1997). Unlike the tableau methods above, the stage count `s` is
    ! chosen every step from `dt` and an estimated spectral radius of the RHS,
    ! and the update is a three-term Chebyshev recurrence rather than a fixed
    ! stage matrix. Its real-axis stability interval grows with s^2, so it can
    ! take much larger steps than RK4/BS3/Tsit5 on stiff problems at the cost of
    ! more (cheap) stages per step, while remaining only 2nd order accurate.
    ! =========================================================================

    subroutine step_rkc(ode_rhs, tf, ode, isos)
        implicit none
        real(wp), intent(in) :: tf
        type(ode_class), intent(inout) :: ode
        type(isos_class), intent(inout) :: isos

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

        ! Sommeijer-Shampine-Verwer (1997) damped RKC2 defaults, matching
        ! FastIsostasy.jl's RKCIntegrator() default constructor.
        real(wp), parameter :: damping = 2.0_wp / 13.0_wp
        real(wp), parameter :: safety  = 1.2_wp
        integer,  parameter :: smax    = 200

        real(wp) :: lambda_max, error, tol
        integer  :: s, j
        real(wp), allocatable :: mu(:), nu(:), mutilde(:), gammatilde(:), c(:)

        ! Seed the spectral-radius estimate once per call (cf. init_integrator
        ! in FastIsostasy.jl); re-estimating every accepted step is unnecessary
        ! since the mantle relaxation operator's stiffness varies slowly.
        lambda_max = spectral_radius_estimate(ode_rhs, ode%x, ode%t, isos)

        do while (ode%t < tf)
            ! handle last partial step cleanly
            ode%dt = min(ode%dt, tf - ode%t)
            ode%dt = min(ode%dt, isos%par%dt_max)

            s = rkc_choose_stages(ode%dt, lambda_max, damping, smax, safety)
            if (allocated(mu)) deallocate(mu, nu, mutilde, gammatilde, c)
            allocate(mu(s), nu(s), mutilde(s), gammatilde(s), c(s))
            call rkc_coeffs(s, damping, mu, nu, mutilde, gammatilde, c)

            ! Rolling Chebyshev recurrence, aliased onto the existing ode work
            ! arrays (no new grid-sized allocations): k1 = F(Y_0) (constant
            ! through the step), k2 = F(Y_{j-1}) (refreshed every stage),
            ! y1/y2 = rolling Y_{j-1}/Y_{j-2}, x1 = stage scratch, x2 = error est.
            ode%k1 = ode_rhs(ode%x, ode%t, isos)
            ode%y1 = ode%x + mutilde(1) * ode%dt * ode%k1    ! Y_1
            ode%y2 = ode%x                                   ! Y_0, rolled as "Y_{-1}"

            do j = 2, s
                ode%k2 = ode_rhs(ode%y1, ode%t + c(j-1) * ode%dt, isos)
                ode%x1 = mu(j) * ode%y1 + nu(j) * ode%y2 + &
                    (1.0_wp - mu(j) - nu(j)) * ode%x + &
                    mutilde(j) * ode%dt * ode%k2 + gammatilde(j) * ode%dt * ode%k1
                ode%y2 = ode%y1
                ode%y1 = ode%x1
            end do
            ! Y_s now lives in ode%y1

            ! SSV embedded error estimate (Sommeijer-Shampine-Verwer 1997, §4):
            ! Est = 0.8*(Y_0 - Y_s) + 0.4*dt*(F(Y_0) + F(Y_s)). Stays bounded as
            ! z -> -inf (unlike a naive Euler-predictor difference), so it does
            ! not defeat the method's extended stability by over-rejecting.
            ode%k2 = ode_rhs(ode%y1, ode%t + ode%dt, isos)
            ode%x2 = 0.8_wp * (ode%x - ode%y1) + 0.4_wp * ode%dt * (ode%k1 + ode%k2)

            tol = isos%par%atol + isos%par%rtol * max(maxval(abs(ode%x)), maxval(abs(ode%y1)))
            error = sqrt(sum((ode%x2 / tol)**2) / size(ode%x))

            if (error < 1) then
                ! Accept step
                ode%x = ode%y1
                isos%now%w = ode%x
                ode%t = ode%t + ode%dt
                ode%dt = min(ode%dt * (1/error) ** (1/isos%par%q), isos%par%dt_max)
                write(*,*) "FastIso: t =", ode%t, " dt =", ode%dt, " error =", error, " s =", s
            else
                ! Reject step and reduce time step
                ode%dt = ode%dt * 0.5
            end if
        end do

        if (allocated(mu)) deallocate(mu, nu, mutilde, gammatilde, c)

    end subroutine step_rkc

    ! Nonlinear power iteration (finite-difference Jacobian-vector products, no
    ! Jacobian ever formed) estimating the dominant |eigenvalue| of a generic
    ! RHS `ode_rhs(x, t, isos)`. Same technique used by RKC/ROCK-type
    ! integrators (Sommeijer, Shampine & Verwer 1997) to size their stage count
    ! without assembling a Jacobian. Ported from FastIsostasy.jl's
    ! `spectral_radius_estimate` (src/stability.jl).
    function spectral_radius_estimate(ode_rhs, x0, t, isos, maxiter, tol) result(lambda_max)
        implicit none
        real(wp), intent(in) :: x0(:, :)
        real(wp), intent(in) :: t
        type(isos_class), intent(inout) :: isos
        integer,  intent(in), optional :: maxiter
        real(wp), intent(in), optional :: tol
        real(wp) :: lambda_max

        abstract interface
            function ode_rhs_(x, t, isos) result(dxdt)
                use isostasy_defs, only : wp, isos_class
                implicit none
                real(wp), intent(in) :: t
                real(wp), intent(in) :: x(:, :)
                type(isos_class), intent(inout) :: isos
                real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
            end function ode_rhs_
        end interface
        procedure(ode_rhs_) :: ode_rhs

        integer  :: maxiter_, iter, i, j, k, nx, ny
        real(wp) :: tol_, unorm, vnorm, h, sigma, sigma_new, eps_wp
        real(wp), allocatable :: F0(:, :), v(:, :), z(:, :), Fz(:, :)
        logical :: converged

        maxiter_ = 100
        if (present(maxiter)) maxiter_ = maxiter
        tol_ = 1.0e-2_wp
        if (present(tol)) tol_ = tol

        nx = size(x0, 1)
        ny = size(x0, 2)
        allocate(F0(nx, ny), v(nx, ny), z(nx, ny), Fz(nx, ny))

        F0 = ode_rhs(x0, t, isos)
        v = F0

        eps_wp = epsilon(1.0_wp)
        unorm = sqrt(sum(x0**2))
        vnorm = sqrt(sum(v**2))

        ! If the RHS at x0 is ~0 (e.g. a cold start at equilibrium), seed the
        ! probe direction with an alternating +/-1 pattern instead.
        if (vnorm < sqrt(eps_wp) * max(1.0_wp, unorm)) then
            k = 0
            do j = 1, ny
                do i = 1, nx
                    k = k + 1
                    if (mod(k, 2) == 1) then
                        v(i, j) = 1.0_wp
                    else
                        v(i, j) = -1.0_wp
                    end if
                end do
            end do
            vnorm = sqrt(sum(v**2))
        end if

        h = sqrt(eps_wp) * max(1.0_wp, unorm)
        sigma = 0.0_wp

        do iter = 1, maxiter_
            z = x0 + (h / vnorm) * v
            Fz = ode_rhs(z, t, isos)
            Fz = Fz - F0
            sigma_new = sqrt(sum(Fz**2)) / h
            converged = (sigma_new > 0.0_wp) .and. (abs(sigma_new - sigma) <= tol_ * sigma_new)
            sigma = sigma_new
            if (converged) exit
            v = Fz
            vnorm = sqrt(sum(v**2))
            if (vnorm < sqrt(eps_wp) * max(1.0_wp, sigma)) exit
        end do

        lambda_max = sigma
        deallocate(F0, v, z, Fz)

    end function spectral_radius_estimate

    ! Smallest stage count (clamped to [2, smax]) whose exact stability
    ! boundary covers `safety * dt * lambda_max`, seeded by the closed-form
    ! asymptotic beta(s) ~= 0.65*s^2 and refined against the exact boundary.
    function rkc_choose_stages(dt, lambda_max, damping, smax, safety) result(s)
        implicit none
        real(wp), intent(in) :: dt, lambda_max, damping, safety
        integer,  intent(in) :: smax
        integer :: s
        real(wp) :: z

        z = safety * dt * lambda_max
        if (z <= 0.0_wp) then
            s = 2
            return
        end if
        s = max(2, min(smax, ceiling(sqrt(z / 0.65_wp))))
        do while (s < smax)
            if (rkc_stability_boundary(s, damping) >= z) exit
            s = s + 1
        end do
    end function rkc_choose_stages

    ! Exact real-axis stability boundary of the s-stage damped RKC recurrence:
    ! the largest beta > 0 such that |Y_s/Y_0| <= 1 for z in [-beta, 0],
    ! evaluated by directly iterating the actual recurrence (rkc_scalar_map),
    ! not a closed-form shortcut. Found by bisection.
    function rkc_stability_boundary(s, damping) result(beta)
        implicit none
        integer,  intent(in) :: s
        real(wp), intent(in) :: damping
        real(wp) :: beta

        real(wp), allocatable :: mu(:), nu(:), mutilde(:), gammatilde(:), c(:)
        real(wp) :: lo, hi, mid, btol, eps_wp

        allocate(mu(s), nu(s), mutilde(s), gammatilde(s), c(s))
        call rkc_coeffs(s, damping, mu, nu, mutilde, gammatilde, c)

        eps_wp = epsilon(1.0_wp)
        btol = sqrt(eps_wp)
        lo = 0.0_wp
        hi = real(4 * s * s, wp)
        do while (stable(hi))
            hi = hi * 2.0_wp
        end do
        do while (hi - lo > btol * max(1.0_wp, hi))
            mid = 0.5_wp * (lo + hi)
            if (stable(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        beta = lo
        deallocate(mu, nu, mutilde, gammatilde, c)

    contains

        logical function stable(b)
            real(wp), intent(in) :: b
            stable = abs(rkc_scalar_map(mu, nu, mutilde, gammatilde, s, -b)) <= 1.0_wp + sqrt(eps_wp)
        end function stable

    end function rkc_stability_boundary

    ! Scalar-test-equation evaluation of the actual s-stage recurrence (not a
    ! closed-form shortcut): returns Y_s/Y_0 for u' = z/dt * u, dt = 1.
    function rkc_scalar_map(mu, nu, mutilde, gammatilde, s, z) result(y1)
        implicit none
        real(wp), intent(in) :: mu(:), nu(:), mutilde(:), gammatilde(:)
        integer,  intent(in) :: s
        real(wp), intent(in) :: z
        real(wp) :: y1

        real(wp) :: F0, y0, y2, Fjm1, ynext
        integer :: j

        F0 = z
        y0 = 1.0_wp
        y1 = y0 + mutilde(1) * F0
        y2 = y0
        do j = 2, s
            Fjm1 = z * y1
            ynext = mu(j) * y1 + nu(j) * y2 + (1.0_wp - mu(j) - nu(j)) * y0 + &
                mutilde(j) * Fjm1 + gammatilde(j) * F0
            y2 = y1
            y1 = ynext
        end do
    end function rkc_scalar_map

    ! Per-stage recurrence coefficients (mu, nu, mutilde, gammatilde, c) for the
    ! s-stage damped, second-order RKC method (Sommeijer-Shampine-Verwer 1997).
    ! Ported term-for-term from FastIsostasy.jl's rkc_coeffs (src/integrators.jl),
    ! itself cross-checked there against SUNDIALS' LSRKStep implementation.
    subroutine rkc_coeffs(s, damping, mu, nu, mutilde, gammatilde, c)
        implicit none
        integer,  intent(in)  :: s
        real(wp), intent(in)  :: damping
        real(wp), intent(out) :: mu(:), nu(:), mutilde(:), gammatilde(:), c(:)

        real(wp), allocatable :: Tv(:), Dv(:), Pv(:)
        real(wp) :: w0, w1, temp1, temp2, arg
        real(wp) :: b1, bj, bjm1, bjm2, a_jm1, cjm2
        integer  :: j

        w0 = 1.0_wp + damping / real(s, wp)**2
        call chebyshev_table(s, w0, Tv, Dv, Pv)   ! Tv(j) == T_j(w0), etc.

        temp1 = w0**2 - 1.0_wp
        temp2 = sqrt(temp1)
        arg   = real(s, wp) * log(w0 + temp2)
        w1 = sinh(arg) * temp1 / (cosh(arg) * real(s, wp) * temp2 - w0 * sinh(arg))

        mu = 0.0_wp
        nu = 0.0_wp
        mutilde = 0.0_wp
        gammatilde = 0.0_wp
        c = 0.0_wp

        b1 = 1.0_wp / (2.0_wp * w0)**2   ! b_0 = b_1, closed form (regularises T_1''/T_1'^2)
        mutilde(1) = b1 * w1
        c(1) = mutilde(1)

        bjm2 = b1
        bjm1 = b1
        do j = 2, s
            bj = Pv(j) / Dv(j)**2
            mu(j) = 2.0_wp * bj * w0 / bjm1
            nu(j) = -bj / bjm2
            mutilde(j) = 2.0_wp * bj * w1 / bjm1
            a_jm1 = 1.0_wp - bjm1 * Tv(j-1)      ! Tv(j-1) == T_{j-1}(w0)
            gammatilde(j) = -a_jm1 * mutilde(j)

            if (j == 2) then
                cjm2 = 0.0_wp
            else
                cjm2 = c(j-2)
            end if
            c(j) = mu(j) * c(j-1) + nu(j) * cjm2 + mutilde(j) + gammatilde(j)

            bjm2 = bjm1
            bjm1 = bj
        end do

        deallocate(Tv, Dv, Pv)

    end subroutine rkc_coeffs

    ! T_j(x), T_j'(x), T_j''(x) for every degree j = 0..s at a single point x,
    ! via the three-term Chebyshev recurrence and its first two derivatives.
    subroutine chebyshev_table(s, x, Tv, Dv, Pv)
        implicit none
        integer,  intent(in) :: s
        real(wp), intent(in) :: x
        real(wp), allocatable, intent(out) :: Tv(:), Dv(:), Pv(:)
        integer :: j

        allocate(Tv(0:s), Dv(0:s), Pv(0:s))
        Tv(0) = 1.0_wp; Dv(0) = 0.0_wp; Pv(0) = 0.0_wp
        Tv(1) = x;      Dv(1) = 1.0_wp; Pv(1) = 0.0_wp
        do j = 1, s - 1
            Tv(j+1) = 2.0_wp * x * Tv(j) - Tv(j-1)
            Dv(j+1) = 2.0_wp * Tv(j) + 2.0_wp * x * Dv(j) - Dv(j-1)
            Pv(j+1) = 4.0_wp * Dv(j) + 2.0_wp * x * Pv(j) - Pv(j-1)
        end do
    end subroutine chebyshev_table

end module integrators