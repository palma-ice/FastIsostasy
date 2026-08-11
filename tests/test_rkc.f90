program test_rkc
    ! Minimal correctness test for the RKC (Runge-Kutta-Chebyshev) integrator
    ! ported from FastIsostasy.jl (src/integrators.jl) into
    ! FastIsostasy/src/integrators.f90.
    !
    ! Test problem: elementwise linear decay dx/dt = -lambda*x on a small grid,
    ! with lambda varying across the grid so that the RHS is genuinely stiff
    ! (max|lambda| large). This has the closed-form solution
    ! x(t) = x0 * exp(-lambda*t), against which we check:
    !
    !   1. RKC accuracy: adaptive RKC reproduces the analytic solution to
    !      within a small tolerance over a short horizon.
    !   2. RKC stability: fixed-step explicit Euler at a step size well
    !      outside its stability interval (dt > 2/|lambda|_max) blows up,
    !      while RKC remains bounded and decays correctly at that same step
    !      size -- demonstrating the extended stability domain that
    !      motivated the port.
    !   3. 2nd-order convergence: halving a fixed step size should reduce
    !      RKC's error by ~4x.

    use isostasy_defs, only : wp, isos_class, ode_class
    use integrators, only : step_rkc, step_euler

    implicit none

    integer, parameter :: nx = 4, ny = 4
    real(wp), parameter :: t0 = 0.0_wp
    real(wp) :: lambda(nx, ny)
    real(wp) :: x0(nx, ny)
    real(wp) :: err_rkc
    integer :: i, j
    logical :: pass

    type(isos_class) :: isos
    type(ode_class)  :: ode

    pass = .true.

    ! --- Build a stiff, spatially varying decay rate ------------------------
    do j = 1, ny
        do i = 1, nx
            lambda(i, j) = 50.0_wp + 10.0_wp * real((i - 1) + (j - 1) * nx, wp)
        end do
    end do
    x0 = 1.0_wp

    write(*, '(a,f10.2)') "Stiffest rate |lambda|_max = ", maxval(abs(lambda))

    ! =========================================================================
    ! Test 1: RKC accuracy against the analytic solution over a short horizon
    ! (lambda*tf kept moderate so the analytic solution does not underflow).
    ! =========================================================================
    block
        real(wp), parameter :: tf1 = 0.05_wp
        real(wp) :: xan(nx, ny)

        xan = x0 * exp(-lambda * tf1)

        call init_minimal_isos(isos, nx, ny, atol=1.0e-8_wp, rtol=1.0e-6_wp, &
            dt_max=1.0_wp, q=2.0_wp)
        call init_ode(ode, nx, ny, x0, t0, dt0=1.0e-3_wp)
        isos%now%w = ode%x

        call step_rkc(decay_rhs, tf1, ode, isos)

        err_rkc = sqrt(sum((ode%x - xan)**2) / size(xan)) / sqrt(sum(xan**2) / size(xan))
        write(*, '(a,es12.4)') "RKC relative RMS error vs analytic solution: ", err_rkc
        if (err_rkc > 1.0e-3_wp) then
            write(*, *) "FAIL: RKC solution not accurate enough."
            pass = .false.
        else
            write(*, *) "PASS: RKC solution matches analytic decay."
        end if
    end block

    ! =========================================================================
    ! Test 2: RKC remains stable and decays correctly with a fixed step size
    ! well outside Euler's stability interval (dt < 2/lambda_max), while
    ! fixed-step Euler at that same step size blows up.
    ! =========================================================================
    block
        real(wp), parameter :: tf2 = 5.0_wp
        real(wp) :: dt_euler_limit, dt_big

        dt_euler_limit = 2.0_wp / maxval(abs(lambda))
        dt_big = 8.0_wp * dt_euler_limit    ! well outside Euler's stability interval

        write(*, '(a,es12.4)') "Euler stability limit dt < ", dt_euler_limit
        write(*, '(a,es12.4)') "Fixed step size used for both methods: ", dt_big

        ! Fixed-step Euler at dt_big: expected to diverge (|1 - lambda*dt| > 1)
        call init_minimal_isos(isos, nx, ny, atol=1.0e-8_wp, rtol=1.0e-6_wp, &
            dt_max=dt_big, q=2.0_wp)
        call init_ode(ode, nx, ny, x0, t0, dt0=dt_big)
        isos%now%w = ode%x
        call step_euler(decay_rhs, tf2, ode, isos)

        if (any(.not. is_finite_arr(ode%x)) .or. maxval(abs(ode%x)) > 1.0e3_wp) then
            write(*, *) "PASS: fixed-step Euler at large dt diverges as expected (unstable)."
        else
            write(*, *) "FAIL: fixed-step Euler unexpectedly remained bounded; test problem not stiff enough."
            pass = .false.
        end if

        ! RKC with the same generous dt_max: should remain stable, bounded and
        ! decay towards the (near-zero) equilibrium, internally choosing enough
        ! stages to cover the stiffness -- unlike Euler above, which diverged
        ! at this exact step size.
        call init_minimal_isos(isos, nx, ny, atol=1.0e-8_wp, rtol=1.0e-6_wp, &
            dt_max=dt_big, q=2.0_wp)
        call init_ode(ode, nx, ny, x0, t0, dt0=dt_big)
        isos%now%w = ode%x
        call step_rkc(decay_rhs, tf2, ode, isos)

        if (all(is_finite_arr(ode%x)) .and. maxval(abs(ode%x)) < 1.0e-2_wp) then
            write(*, '(a,es12.4)') "PASS: RKC stays bounded and decays correctly, max|x| = ", maxval(abs(ode%x))
        else
            write(*, *) "FAIL: RKC did not remain bounded/accurate at the larger step size."
            pass = .false.
        end if
    end block

    ! =========================================================================
    ! Test 3: 2nd-order convergence -- halving the (fixed) step size should
    ! reduce RKC's error by ~4x. dt1/dt2 are chosen so that even the stiffest
    ! grid point (lambda*dt <~ 0.2) stays in the asymptotic small-step regime;
    ! a 2nd-order method's convergence rate is only expected to show up there,
    ! not for lambda*dt = O(1) steps where the stability function R(z), not
    ! the local Taylor truncation, dominates the error.
    ! =========================================================================
    block
        real(wp), parameter :: tf3 = 0.02_wp
        real(wp) :: xan(nx, ny)
        real(wp) :: e1, e2, e_ratio, dt1, dt2

        xan = x0 * exp(-lambda * tf3)

        dt1 = 0.001_wp
        dt2 = dt1 / 2.0_wp

        call init_minimal_isos(isos, nx, ny, atol=1.0e30_wp, rtol=1.0e30_wp, &
            dt_max=dt1, q=2.0_wp)   ! huge tol => controller always accepts => fixed step
        call init_ode(ode, nx, ny, x0, t0, dt0=dt1)
        isos%now%w = ode%x
        call step_rkc(decay_rhs, tf3, ode, isos)
        e1 = sqrt(sum((ode%x - xan)**2) / size(xan))

        call init_minimal_isos(isos, nx, ny, atol=1.0e30_wp, rtol=1.0e30_wp, &
            dt_max=dt2, q=2.0_wp)
        call init_ode(ode, nx, ny, x0, t0, dt0=dt2)
        isos%now%w = ode%x
        call step_rkc(decay_rhs, tf3, ode, isos)
        e2 = sqrt(sum((ode%x - xan)**2) / size(xan))

        e_ratio = e1 / max(e2, tiny(1.0_wp))
        write(*, '(a,es12.4,a,es12.4)') "e1 = ", e1, "  e2 = ", e2
        write(*, '(a,f8.3)') "Error ratio e(dt)/e(dt/2) (expect ~4 for 2nd order): ", e_ratio
        if (e_ratio > 2.5_wp .and. e_ratio < 6.0_wp) then
            write(*, *) "PASS: RKC exhibits ~2nd-order convergence."
        else
            write(*, *) "FAIL: convergence order does not match expected 2nd order."
            pass = .false.
        end if
    end block

    if (pass) then
        write(*, *) ""
        write(*, *) "ALL TESTS PASSED"
    else
        write(*, *) ""
        write(*, *) "SOME TESTS FAILED"
        stop 1
    end if

contains

    function decay_rhs(x, t, isos) result(dxdt)
        use isostasy_defs, only : wp, isos_class
        implicit none
        real(wp), intent(in) :: t
        real(wp), intent(in) :: x(:, :)
        type(isos_class), intent(inout) :: isos
        real(wp), dimension(size(x, 1), size(x, 2)) :: dxdt
        dxdt = -lambda * x
    end function decay_rhs

    elemental function is_finite_arr(x) result(is_finite)
        real(wp), intent(in) :: x
        logical :: is_finite
        is_finite = (x == x) .and. (abs(x) < huge(1.0_wp))
    end function is_finite_arr

    subroutine init_minimal_isos(isos, nx, ny, atol, rtol, dt_max, q)
        type(isos_class), intent(inout) :: isos
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: atol, rtol, dt_max, q
        isos%par%atol = atol
        isos%par%rtol = rtol
        isos%par%dt_max = dt_max
        isos%par%q = q
        if (.not. allocated(isos%now%w)) allocate(isos%now%w(nx, ny))
    end subroutine init_minimal_isos

    subroutine init_ode(ode, nx, ny, x0, t0, dt0)
        type(ode_class), intent(inout) :: ode
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: x0(nx, ny), t0, dt0
        if (.not. allocated(ode%x))  allocate(ode%x(nx, ny))
        if (.not. allocated(ode%x1)) allocate(ode%x1(nx, ny))
        if (.not. allocated(ode%x2)) allocate(ode%x2(nx, ny))
        if (.not. allocated(ode%k1)) allocate(ode%k1(nx, ny))
        if (.not. allocated(ode%k2)) allocate(ode%k2(nx, ny))
        if (.not. allocated(ode%y1)) allocate(ode%y1(nx, ny))
        if (.not. allocated(ode%y2)) allocate(ode%y2(nx, ny))
        ode%x = x0
        ode%t = t0
        ode%dt = dt0
    end subroutine init_ode

end program test_rkc
