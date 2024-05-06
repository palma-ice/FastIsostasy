
module lv_elva

    use, intrinsic :: iso_c_binding
    use isostasy_defs, only : sp, dp, wp, pi, isos_domain_class, isos_param_class
    use finite_differences
    use isos_utils

    implicit none
    include 'fftw3.f03'

    private
    
    public :: calc_lvelva
    public :: calc_layerboundaries
    public :: layered_viscosity
    public :: calc_effective_viscosity
    public :: calc_fft_backward_r2r
    public :: calc_fft_forward_r2r
    public :: calc_fft_backward_c2r
    public :: calc_fft_forward_r2c
    public :: calc_kappa
    public :: calc_beta
    public :: convenient_calc_kappa

    contains

    ! Calculate vertical displacement rate (viscous part) on rectangular domain.
    subroutine calc_lvelva(dwdt, w, canom_full, maskactive, g, nu, D_lith, eta,&
        kappa, nx, ny, dx_matrix, dy_matrix, sec_per_year, forward_plan, backward_plan)

        implicit none
        
        real(wp), intent(INOUT) :: dwdt(:, :)
        real(wp), intent(IN)    :: w(:, :)
        real(wp), intent(IN)    :: canom_full(:, :)
        logical,  intent(IN)    :: maskactive(:, :)
        real(wp), intent(IN)    :: g
        real(wp), intent(IN)    :: nu
        real(wp), intent(IN)    :: D_lith(:, :) 
        real(wp), intent(IN)    :: eta(:, :)   ! [Pa s] Viscosity, eta=1e21 by default. 
        real(wp), intent(IN)    :: kappa(:, :)
        integer, intent(IN)     :: nx, ny
        real(wp), intent(IN)    :: sec_per_year
        real(wp), intent(IN)    :: dx_matrix(:, :)
        real(wp), intent(IN)    :: dy_matrix(:, :)
        type(c_ptr), intent(IN) :: forward_plan
        type(c_ptr), intent(IN) :: backward_plan

        real(wp), allocatable :: p(:, :)
        real(dp), allocatable :: f(:, :)
        real(dp), allocatable :: f_hat(:, :)
        real(dp), allocatable :: dwdt_hat(:, :)

        real(dp), allocatable :: dwdt_dp(:,:) 

        real(wp), allocatable :: w_x(:, :)
        real(wp), allocatable :: w_xy(:, :)
        real(wp), allocatable :: w_xx(:, :)
        real(wp), allocatable :: w_yy(:, :)

        real(wp), allocatable :: Mxy(:, :)
        real(wp), allocatable :: Mxy_x(:, :)
        real(wp), allocatable :: Mxy_xy(:, :)

        real(wp), allocatable :: Mxx(:, :) 
        real(wp), allocatable :: Mx_xx(:, :)
        real(wp), allocatable :: Myy(:, :)
        real(wp), allocatable :: My_yy(:, :)

        allocate(p(nx, ny))
        allocate(f(nx, ny))
        allocate(f_hat(nx, ny))
        allocate(dwdt_hat(nx, ny))

        allocate(dwdt_dp(nx, ny))

        allocate(w_x(nx, ny))
        allocate(w_xy(nx, ny))
        allocate(w_xx(nx, ny))
        allocate(w_yy(nx, ny))

        allocate(Mxx(nx, ny))
        allocate(Myy(nx, ny))

        allocate(Mx_xx(nx, ny))
        allocate(My_yy(nx, ny))
        allocate(Mxy_x(nx, ny))
        allocate(Mxy_xy(nx, ny))

        ! Finite differences
        call calc_derivative_x(w_x, w, dx_matrix, nx, ny)
        call calc_derivative_y(w_xy, w_x, dy_matrix, nx, ny)
        call calc_derivative_xx(w_xx, w, dx_matrix, nx, ny)
        call calc_derivative_yy(w_yy, w, dy_matrix, nx, ny)

        ! Ventsel and Krauthammer (2001): Thin Plates and Shells.
        ! Theory, Analysis, and Applications. Eq (2.13, 2.23)
        Mxx = -D_lith*(w_xx + nu*w_yy)
        Myy = -D_lith*(w_yy + nu*w_xx)
        Mxy = -D_lith * (1.0-nu) * w_xy

        ! Finite differences
        call calc_derivative_x(Mxy_x, Mxy, dx_matrix, nx, ny)
        call calc_derivative_y(Mxy_xy, Mxy_x, dy_matrix, nx, ny)
        call calc_derivative_xx(Mx_xx, Mxx, dx_matrix, nx, ny)
        call calc_derivative_yy(My_yy, Myy, dy_matrix, nx, ny)

        call maskfield(p, -g * canom_full, maskactive, nx, ny)
        f = (p + Mx_xx + 2.0_wp * Mxy_xy + My_yy) / (2.0_wp * eta)
        call calc_fft_forward_r2r(forward_plan, f, f_hat)
        dwdt_hat = f_hat / kappa

        ! write(*,*) sum(dwdt_hat), sum(f_hat)
        call calc_fft_backward_r2r(backward_plan, dwdt_hat, dwdt_dp)
        call apply_zerobc_at_corners_dp(dwdt_dp, nx, ny)
        dwdt = dwdt_dp 
        
        ! Rate of viscous asthenosphere uplift per unit time (seconds)
        dwdt = dwdt * sec_per_year  !  [m/s] x [s/a] = [m/a]

        return
    end subroutine calc_lvelva
    
    subroutine calc_layerboundaries(layer_boundaries, He_lith, layer_boundaries_vec)
        implicit none
        real(wp), intent(INOUT) :: layer_boundaries(:, :, :)
        real(wp), intent(IN)    :: He_lith(:, :)
        real(wp), intent(IN)    :: layer_boundaries_vec(:)
        integer                 :: k, nl

        nl = size(layer_boundaries_vec)
        layer_boundaries = 0.0_wp
        layer_boundaries(:, :, 1) = He_lith

        do k = 2, nl
            layer_boundaries(:, :, k) = layer_boundaries_vec(k)
        end do
        return
    end subroutine calc_layerboundaries

    subroutine layered_viscosity(eta, eta_vec)
        implicit none
        real(wp), intent(INOUT) :: eta(:, :, :)
        real(wp), intent(IN)    :: eta_vec(:)
        integer                 :: k, nl

        nl = size(eta_vec)
        do k = 1, nl
            eta(:, :, k) = eta_vec(k)
        end do

        return
    end subroutine layered_viscosity

    subroutine calc_effective_viscosity(eta_eff, eta, dx, dy, layer_boundaries)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:, :)
        real(wp), intent(IN)     :: eta(:, :, :)
        real(wp), intent(IN)     :: dx, dy
        real(wp), intent(IN)     :: layer_boundaries(:, :, :)

        real(wp) :: Lx, Ly, L
        real(wp) :: kappa
        real(wp), allocatable ::  eta_c(:, :)
        real(wp), allocatable ::  eta_ratio(:, :)
        real(wp), allocatable ::  inv_ratio(:, :)
        real(wp), allocatable ::  R(:, :)
        real(wp), allocatable ::  dz(:, :)
        real(wp), allocatable ::  c(:, :)
        real(wp), allocatable ::  s(:, :)
       
        integer  :: i, j, k, nx, ny, nl

        nx = size(eta_eff, 1)
        ny = size(eta_eff, 2)
        nl = size(eta, 3)
        if (nl .ne. size(layer_boundaries, 3)) then
            write(*,*) "Number of levels in eta and in boundaries do not coincide."
            stop
        end if

        allocate(eta_c(nx, ny))
        allocate(eta_ratio(nx, ny))
        allocate(inv_ratio(nx, ny))
        allocate(R(nx, ny))
        allocate(dz(nx, ny))
        allocate(c(nx, ny))
        allocate(s(nx, ny))

        if (nl .eq. 1) then
            eta_eff = eta(:, :, 1)

        else if (nl .ge. 1) then
            Lx = dx * (nx-1)
            Ly = dy * (ny-1)
            L = (Lx + Ly) / 2.0
            kappa = 2*pi/L

            ! Start with n-th layer: viscous half space
            eta_eff(:, :) = eta(:, :, nl)
            
            do k = nl, 2
                write(*,*) k
                dz = layer_boundaries(:, :, k) - layer_boundaries(:, :, k-1)
                eta_c = eta(:, :, k-1)
                eta_ratio = eta_c / eta_eff
                inv_ratio = 1. / eta_ratio

                c = cosh(dz*kappa)
                s = sinh(dz*kappa)

                R = (2.0 * eta_ratio * c * s + (1-eta_ratio**2) * (dz*kappa)**2 + &
                    (eta_ratio*s)**2 + c**2 ) / ((eta_ratio + inv_ratio)* c * s + &
                    (eta_ratio - inv_ratio)*dz*kappa + s**2 + c**2)
                eta_eff = R * eta_eff
            end do

        else
            print*,'Number of levels is wrong: nl = ', nl
            stop
        endif

        return
    end subroutine calc_effective_viscosity

    subroutine calc_beta(beta, kappa, D_lith, rho_uppermantle, g) 
        ! Calculate analytical solution as in Bueler et al 2007 (eq 11)
                    
        real(wp), intent(OUT)  :: beta(:, :) 
        real(wp), intent(IN)   :: kappa(:, :)
        real(wp), intent(IN)   :: D_lith(:, :)
        real(wp), intent(IN)   :: rho_uppermantle
        real(wp), intent(IN)   :: g 
        
        ! Local variables
        integer  :: i, j, nx, ny 
        integer  :: ip, iq, ic, jc 

        beta = 0.0 
        nx = size(beta, 1)
        ny = size(beta, 2)
        ic = (nx-1)/2 + 1
        jc = (ny-1)/2 + 1

        do i = 1, nx
            if (i.le.ic) then 
                ip = i-1
            else
                ip = nx-i+1
            end if
            do j = 1, ny
                if (j.le.jc) then
                    iq = j-1  
                else
                    iq = ny-j+1
                end if
                beta(i,j)   = rho_uppermantle*g + D_lith(i,j)*kappa(i,j)**4
            end do
        end do

        return
    end subroutine calc_beta

    subroutine convenient_calc_kappa(domain)
        implicit none
        type(isos_domain_class), intent(INOUT)  :: domain

        call calc_kappa(domain%kappa, domain%nx, domain%ny, domain%dx, domain%dy)
        return
    end subroutine convenient_calc_kappa

    subroutine calc_kappa(kappa, nx, ny, dx, dy)
        implicit none
        integer, intent(IN)     :: nx, ny
        real(wp), intent(IN)    :: dx, dy
        real(wp), intent(OUT)   :: kappa(:, :)

        integer :: i, j, ic, jc
        real(wp):: mu_x, mu_y
        real(wp):: p, q

        mu_x = 2._wp * pi / ((nx-1) * dx)    ! 2 * pi / Lx
        mu_y = 2._wp * pi / ((ny-1) * dy)    ! 2 * pi / Ly
        ic = (nx-1)/2 + 1
        jc = (ny-1)/2 + 1

        kappa = 0.0
        do i = 1, nx
            
            if (i .le. ic) then 
                p = mu_x * (i-1)
            else
                p = mu_x * (nx-i+1)
            end if

            do j = 1, ny

                if (j .le. jc) then
                    q = mu_y * (j-1)
                else
                    q = mu_y * (ny-j+1)
                end if

                kappa(i, j)  = (p*p + q*q)**0.5

            end do
        end do
        kappa(1,1) = (kappa(1,2) + kappa(2,1)) / 2.0

        return
    end subroutine calc_kappa

end module lv_elva