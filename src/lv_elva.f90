
module lv_elva

    use, intrinsic :: iso_c_binding
    use isostasy_defs, only : wp, pi, isos_domain_class, isos_param_class
    use finite_differences
    use isos_utils

    implicit none
    include 'fftw3.f03'

    private
    
    public :: calc_lvelva
    public :: calc_effective_viscosity_3layer_channel
    public :: calc_effective_viscosity_3d
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
        real(wp), allocatable :: f(:, :)
        real(wp), allocatable :: f_hat(:, :)
        real(wp), allocatable :: dwdt_hat(:, :)

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
        call calc_fft_backward_r2r(backward_plan, dwdt_hat, dwdt)
        call apply_zerobc_at_corners(dwdt, nx, ny)

        ! Rate of viscous asthenosphere uplift per unit time (seconds)
        dwdt = dwdt * sec_per_year  !  [m/s] x [s/a] = [m/a]

        return
    end subroutine calc_lvelva

    ! TODO: this subroutine should be removed
    subroutine calc_effective_viscosity_3layer_channel(eta_eff, visc_c, thck_c, He_lith, &
        n_lev, dx, dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:, :)
        real(wp), intent(IN)     :: visc_c
        real(wp), intent(IN)     :: thck_c
        real(wp), intent(IN)     :: He_lith !(:, :) 
        integer,  intent(IN)     :: n_lev
        real(wp), intent(IN)     :: dx, dy
        
        real(wp) :: Lx, Ly, L
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin
        real(wp), allocatable :: xc(:), yc(:)

        real(wp), allocatable ::  R(:, :)
        real(wp), allocatable ::  eta(:, :,:)
        real(wp), allocatable ::  eta_ratio(:, :)
        real(wp), allocatable ::  eta_c(:, :)
        real(wp), allocatable ::  dz(:, :,:)
        real(wp), allocatable ::  dz_c(:, :)
        real(wp), allocatable ::  eta_ratiom1(:, :)
        real(wp), allocatable ::  c(:, :)
        real(wp), allocatable ::  s(:, :)
        real(wp), allocatable :: kappa(:, :)
       
        integer  :: i, j, k, nx, ny

        nx = size(eta_eff,1)
        ny = size(eta_eff,2) 
        
        allocate(xc(nx))
        allocate(yc(ny))
        allocate(R(nx, ny))
        allocate(eta_ratio(nx, ny))
        allocate(eta_ratiom1(nx, ny))
        allocate(c(nx, ny))
        allocate(s(nx, ny))
        allocate(eta(nx,ny,n_lev))
        allocate(eta_c(nx, ny))
        allocate(dz_c(nx, ny))
        allocate(dz(nx,ny,n_lev))
        allocate(kappa(nx, ny))

        if (n_lev.lt.2) then

           print*,'n_lev should be at least 2'
           stop
           
        else if (n_lev.gt.3) then
           print*,'Option n_lev > 3 not enabled for viscosity yet'
           stop

        else if (n_lev.eq.3) then

            do i = 1, nx
            xc(i) = dx*(i-1)
            end do
            xmin = xc(1)
            xmax = xc(nx)

            do j = 1, ny
            yc(j) = dy*(j-1)
            enddo
            ymin = yc(1)
            ymax = yc(ny)

            
            xcntr = (xmax+xmin)/2.0
            ycntr = (ymax+ymin)/2.0

            Lx = xmax - xmin
            Ly = ymax - ymin
            L = (Lx + Ly) / 2.0

            kappa = 2*pi/L 
                    
            eta(:, :,1)  = 0.
            eta(:, :,2)  = visc_c      
            eta(:, :,3)  = eta_eff

            dz(:, :,1)  = He_lith*1.e3
            dz(:, :,2)  = thck_c*1.e3     ![m]
            dz(:, :,3)  = 2000.*1.e3      ![m]
        
            ! Start with n-th layer: viscous half space
            
            eta_eff(:, :) = eta(:, :,n_lev)
            
            do k = 1, n_lev-1

            eta_c = eta(:, :,n_lev-1)
            dz_c = dz(:, :,n_lev-k+1)

            eta_ratio = eta_c/eta_eff
            eta_ratiom1 = 1./eta_ratio

            c = cosh(dz_c*kappa)
            s = sinh(dz_c*kappa)

            R = (2.0 * eta_ratio * c * s + (1-eta_ratio**2) * (dz_c*kappa)**2 + &
                (eta_ratio*s)**2 + c**2 ) / ((eta_ratio + eta_ratiom1)* c * s + &
                (eta_ratio - eta_ratiom1)*dz_c*kappa + s**2 + c**2)
            
            eta_eff = R*eta_c
            
            end do

        else

            print*,'n_lev = ', n_lev
        
            stop

        endif
        ! move this out for symmetry     eta_eff    = eta_eff  * (1.5/(1. + nu))         ! [Pa s]
     
        return
    end subroutine calc_effective_viscosity_3layer_channel

    subroutine calc_effective_viscosity_3d(eta_eff, eta, dx, dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:, :)
        real(wp), intent(IN)     :: eta(:, :,:)
        real(wp), intent(IN)     :: dx, dy
        
        real(wp) :: Lx, Ly, L
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin
        real(wp), allocatable :: xc(:), yc(:)

        real(wp), allocatable ::  R(:, :)
        real(wp), allocatable ::  eta_ratio(:, :)
        real(wp), allocatable ::  eta_c(:, :)
        real(wp), allocatable ::  dz(:, :,:)
        real(wp), allocatable ::  dz_c(:, :)
        real(wp), allocatable ::  eta_ratiom1(:, :)
        real(wp), allocatable ::  c(:, :)
        real(wp), allocatable ::  s(:, :)
        real(wp), allocatable :: kappa(:, :)
       
        integer  :: i, j, k, nx, ny, n_lev

        nx = size(eta_eff,1)
        ny = size(eta_eff,2)
        n_lev = size(eta,3)
        
        allocate(xc(nx))
        allocate(yc(ny))
        allocate(R(nx, ny))
        allocate(eta_ratio(nx, ny))
        allocate(eta_ratiom1(nx, ny))
        allocate(c(nx, ny))
        allocate(s(nx, ny))
        allocate(eta_c(nx, ny))
        allocate(dz_c(nx, ny))
        allocate(dz(nx,ny,n_lev))
        allocate(kappa(nx, ny))

        dz = 100.0*1.e3

        if (n_lev.lt.2) then

           print*,'n_lev should be at least 2'
           stop
           
!        else if (n_lev.gt.3) then
!           print*,'Option n_lev > 3 not enabled for viscosity yet'
!           stop

        else if (n_lev.ge.3) then

            do i = 1, nx
                xc(i) = dx*(i-1)
            end do
            xmin = xc(1)
            xmax = xc(nx)

            do j = 1, ny
                yc(j) = dy*(j-1)
            enddo
            ymin = yc(1)
            ymax = yc(ny)

            
            xcntr = (xmax+xmin)/2.0
            ycntr = (ymax+ymin)/2.0

            Lx = xmax - xmin
            Ly = ymax - ymin
            L = (Lx + Ly) / 2.0

            kappa = 2*pi/L 
                    
        
            ! Start with n-th layer: viscous half space
            eta_eff(:, :) = eta(:, :,n_lev)
            
            do k = 1, n_lev-1

                eta_c = eta(:, :,n_lev-1)
                dz_c = 100.*1.e3 !dz(:, :,n_lev-k+1)
                
                eta_ratio = eta_c/eta_eff
                eta_ratiom1 = 1./eta_ratio

                c = cosh(dz_c*kappa)
                s = sinh(dz_c*kappa)

                R = (2.0 * eta_ratio * c * s + (1-eta_ratio**2) * (dz_c*kappa)**2 + &
                    (eta_ratio*s)**2 + c**2 ) / ((eta_ratio + eta_ratiom1)* c * s + &
                    (eta_ratio - eta_ratiom1)*dz_c*kappa + s**2 + c**2)
                
                eta_eff = R*eta_c
            
            end do

        else
            print*,'n_lev = ', n_lev
            stop
        endif

        return
    end subroutine calc_effective_viscosity_3d

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