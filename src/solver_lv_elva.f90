
module solver_lv_elva

    use, intrinsic :: iso_c_binding
    use isostasy_defs, only : wp, pi, isos_domain_class, isos_param_class
    use finite_differences
    use isos_utils

    implicit none
    include 'fftw3.f03'

    private
    
    public :: calc_lvelva_viscous_square
    public :: calc_effective_viscosity_3layer_channel
    public :: calc_effective_viscosity_3d
    public :: calc_fft_backward_r2r
    public :: calc_fft_forward_r2r
    public :: calc_fft_backward_c2r
    public :: calc_fft_forward_r2c
    public :: calc_kappa
    public :: calc_beta

    contains

    subroutine calc_lvelva_viscous_square(dzbdt, w, canom_full, nu, mu, D_lith, eta, &
        kappa, nsq, par, domain)
    
        ! Extend a given domain [nx,ny] so that it is square based on the largest dimension
        ! of original data. Solve calc_asthenosphere_viscous on the square arrays, then 
        ! extract solution onto original domain [nx,ny].
      
        implicit none

        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: w(:,:)
        real(wp), intent(IN)    :: canom_full(:,:)
        real(wp), intent(IN)    :: nu, mu
        real(wp), intent(IN)    :: D_lith(:,:)
        real(wp), intent(IN)    :: eta(:,:)
        real(wp), intent(INOUT) :: kappa(:,:)
        integer, intent(IN)     :: nsq
        type(isos_param_class)  :: par
        type(isos_domain_class) :: domain

        ! Local variables
        real(wp), allocatable :: sq_dzbdt(:,:)
        real(wp), allocatable :: sq_w(:,:)
        real(wp), allocatable :: sq_canom_full(:,:)
        real(wp), allocatable :: sq_D_lith(:,:)
        real(wp), allocatable :: sq_eta(:,:)    

        ! Step 0: determine size of square array and allocate variables
        allocate(sq_dzbdt(nsq,nsq))
        allocate(sq_w(nsq,nsq))
        allocate(sq_canom_full(nsq,nsq))
        allocate(sq_D_lith(nsq,nsq))
        allocate(sq_eta(nsq,nsq))

        ! Step 1: populate variables on a square grid
        sq_dzbdt = 0.0 
        call extend_array(sq_dzbdt, dzbdt, fill_with="mirror", val=0.0_wp)
        call extend_array(sq_w, w, fill_with="mirror", val=0.0_wp)
        call extend_array(sq_canom_full, canom_full, fill_with="mirror", val=0.0_wp)
        call extend_array(sq_D_lith, D_lith, fill_with="mirror", val=0.0_wp)
        call extend_array(sq_eta, eta, fill_with="mirror", val=0.0_wp)

        ! Step 2: solve
        call calc_lv_asthenosphere_viscous(sq_dzbdt, sq_w, sq_canom_full, nu, mu, sq_D_lith, &
            sq_eta, kappa, nsq, nsq, par, domain)

        ! Step 3: get solution on original grid
        call reduce_array(dzbdt, sq_dzbdt)
        call reduce_array(w, sq_w)
        
        return

    end subroutine calc_lvelva_viscous_square


    ! Calculate vertical displacement rate (viscous part) on rectangular domain.
    subroutine calc_lv_asthenosphere_viscous(dzbdt, u, canom_full, nu, mu, D_lith, eta, &
        kappa, nx, ny, par, domain)

        implicit none

        real(wp), parameter :: epsilon = 1.e-2 
        
        real(wp), intent(INOUT) :: dzbdt(:,:)
        real(wp), intent(INOUT) :: u(:,:)
        real(wp), intent(INOUT) :: canom_full(:,:)
        real(wp), intent(IN)    :: nu, mu
        real(wp), intent(IN)    :: D_lith(:,:) 
        real(wp), intent(IN)    :: eta(:,:)   ! [Pa s] Viscosity, eta=1e21 by default. 
        real(wp), intent(INOUT) :: kappa(:,:)
        integer, intent(IN)     :: nx, ny
        type(isos_param_class)  :: par
        type(isos_domain_class) :: domain

        real(wp), allocatable :: f(:,:)
        real(wp), allocatable :: f_hat(:,:)

        real(wp), allocatable :: prod(:,:)
        real(wp), allocatable :: prod_hat(:,:)

        real(wp), allocatable :: dudt(:,:)
        real(wp), allocatable :: u_x(:,:)
        real(wp), allocatable :: u_y(:,:)

        real(wp), allocatable :: u_xx(:,:)
        real(wp), allocatable :: u_xy(:,:)
        real(wp), allocatable :: u_yy(:,:)

        real(wp), allocatable :: Mx(:,:) 
        real(wp), allocatable :: My(:,:)

        real(wp), allocatable :: Mxy(:,:)
        real(wp), allocatable :: Mxy_x(:,:)
        real(wp), allocatable :: Myx_y(:,:)

        real(wp), allocatable :: Mx_xx(:,:)
        real(wp), allocatable :: My_yy(:,:)
        
        real(wp), allocatable :: Mxy_xy(:,:)

        allocate(f(nx,ny))
        allocate(f_hat(nx,ny))

        allocate(prod(nx,ny))
        allocate(prod_hat(nx,ny))

        allocate(dudt(nx,ny))

        allocate(u_x(nx,ny))
        allocate(u_y(nx,ny))

        allocate(u_xx(nx,ny))
        allocate(u_xy(nx,ny))
        allocate(u_yy(nx,ny))

        allocate(Mx(nx,ny))
        allocate(My(nx,ny))

        allocate(Mxy_x(nx,ny))
        allocate(Myx_y(nx,ny))
        
        allocate(Mx_xx(nx,ny))
        allocate(My_yy(nx,ny))
        allocate(Mxy_xy(nx,ny))

        ! Finite differences
        call calc_derivative_xx(u_xx, u, domain%dx, nx, ny)
        call calc_derivative_yy(u_yy, u, domain%dy, nx, ny)
        call calc_derivative_x(u_x, u, domain%dx, nx, ny)
        call calc_derivative_y(u_xy, u_x, domain%dy, nx, ny)

        ! Ventsel and Krauthammer (2001): Thin Plates and Shells.
        ! Theory, Analysis, and Applications. Eq (2.13, 2.23)
        Mx = -D_lith*(u_xx + nu*u_yy)
        My = -D_lith*(u_yy + nu*u_xx)
        Mxy = -D_lith*(1.0-nu)*u_xy

        ! Finite differences
        call calc_derivative_xx(Mx_xx, Mx, domain%dx, nx, ny)
        call calc_derivative_yy(My_yy, My, domain%dy, nx, ny)
        call calc_derivative_x(Mxy_x, Mxy, domain%dx, nx, ny)
        call calc_derivative_y(Mxy_xy, Mxy_x, domain%dy, nx, ny)

        f = (canom_full + Mx_xx + 2.0*Mxy_xy + My_yy) / (2. * eta)
        call calc_fft_forward_r2r(domain%forward_fftplan_r2r, f, f_hat)

        prod_hat = f_hat / kappa / mu
        call calc_fft_backward_r2r(domain%backward_fftplan_r2r, prod_hat, prod)
        dudt = prod  ! [m/s]

        dudt  = dudt - 0.25 * (dudt(1,1) + dudt(nx,ny) + dudt(1,ny) + dudt(nx,1))

        ! Rate of viscous asthenosphere uplift per unit time (seconds)
        dzbdt = dudt * par%sec_per_year  !  [m/s] x [s/a] = [m/a] 

        deallocate(f)
        deallocate(f_hat)
        deallocate(prod)
        deallocate(prod_hat)
        deallocate(dudt)
        
        deallocate(u_x)
        deallocate(u_y)
        
        deallocate(u_xx)
        deallocate(u_xy)
        deallocate(u_yy)

        deallocate(Mx_xx)
        deallocate(Mxy_xy)
        deallocate(My_yy)

        return

    end subroutine calc_lv_asthenosphere_viscous

    subroutine calc_effective_viscosity_3layer_channel(eta_eff, visc_c, thck_c, He_lith, &
        n_lev, dx, dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:,:)
        real(wp), intent(IN)     :: visc_c
        real(wp), intent(IN)     :: thck_c
        real(wp), intent(IN)     :: He_lith !(:,:) 
        integer,  intent(IN)     :: n_lev
        real(wp), intent(IN)     :: dx, dy
        
        real(wp) :: Lx, Ly, L
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin
        real(wp), allocatable :: xc(:), yc(:)

        real(wp), allocatable ::  R(:,:)
        real(wp), allocatable ::  eta(:,:,:)
        real(wp), allocatable ::  eta_ratio(:,:)
        real(wp), allocatable ::  eta_c(:,:)
        real(wp), allocatable ::  dz(:,:,:)
        real(wp), allocatable ::  dz_c(:,:)
        real(wp), allocatable ::  eta_ratiom1(:,:)
        real(wp), allocatable ::  c(:,:)
        real(wp), allocatable ::  s(:,:)
        real(wp), allocatable :: kappa(:,:)
       
        integer  :: i, j, k, nx, ny

        nx = size(eta_eff,1)
        ny = size(eta_eff,2) 
        
        allocate(xc(nx))
        allocate(yc(ny))
        allocate(R(nx,ny))
        allocate(eta_ratio(nx,ny))
        allocate(eta_ratiom1(nx,ny))
        allocate(c(nx,ny))
        allocate(s(nx,ny))
        allocate(eta(nx,ny,n_lev))
        allocate(eta_c(nx,ny))
        allocate(dz_c(nx,ny))
        allocate(dz(nx,ny,n_lev))
        allocate(kappa(nx,ny))

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
                    
            eta(:,:,1)  = 0.
            eta(:,:,2)  = visc_c      
            eta(:,:,3)  = eta_eff

            dz(:,:,1)  = He_lith*1.e3
            dz(:,:,2)  = thck_c*1.e3     ![m]
            dz(:,:,3)  = 2000.*1.e3      ![m]
        
            ! Start with n-th layer: viscous half space
            
            eta_eff(:,:) = eta(:,:,n_lev)
            
            do k = 1, n_lev-1

            eta_c = eta(:,:,n_lev-1)
            dz_c = dz(:,:,n_lev-k+1)

            eta_ratio = eta_c/eta_eff
            eta_ratiom1 = 1./eta_ratio

            c = cosh(dz_c*kappa)
            s = sinh(dz_c*kappa)

            R = (2.0 * eta_ratio * c * s + (1-eta_ratio**2) * (dz_c*kappa)**2 + (eta_ratio*s)**2 + c**2 )/&
                    ((eta_ratio + eta_ratiom1)* c * s + (eta_ratio - eta_ratiom1)*dz_c*kappa + s**2 + c**2)
            
            eta_eff = R*eta_c
            
            end do

        else

            print*,'n_lev = ', n_lev
        
            stop

        endif
        ! move this out for symmetry     eta_eff    = eta_eff  * (1.5/(1. + nu))         ! [Pa s]

        
        deallocate(xc)
        deallocate(yc)        
        deallocate(R)
        deallocate(eta_ratio)
        deallocate(eta_ratiom1)
        deallocate(c)
        deallocate(s)
        deallocate(eta)
        deallocate(eta_c)
        deallocate(dz_c)
        deallocate(dz)
        deallocate(kappa)
     
        return
        
    end subroutine calc_effective_viscosity_3layer_channel

    subroutine calc_effective_viscosity_3d(eta_eff,eta,dx,dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:,:)
        real(wp), intent(IN)     :: eta(:,:,:)
        real(wp), intent(IN)     :: dx, dy
        
        real(wp) :: Lx, Ly, L
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin
        real(wp), allocatable :: xc(:), yc(:)

        real(wp), allocatable ::  R(:,:)
        real(wp), allocatable ::  eta_ratio(:,:)
        real(wp), allocatable ::  eta_c(:,:)
        real(wp), allocatable ::  dz(:,:,:)
        real(wp), allocatable ::  dz_c(:,:)
        real(wp), allocatable ::  eta_ratiom1(:,:)
        real(wp), allocatable ::  c(:,:)
        real(wp), allocatable ::  s(:,:)
        real(wp), allocatable :: kappa(:,:)
       
        integer  :: i, j, k, nx, ny, n_lev

        nx = size(eta_eff,1)
        ny = size(eta_eff,2)

        n_lev = size(eta,3)
        
        allocate(xc(nx))
        allocate(yc(ny))        
        allocate(R(nx,ny))
        allocate(eta_ratio(nx,ny))
        allocate(eta_ratiom1(nx,ny))
        allocate(c(nx,ny))
        allocate(s(nx,ny))
        allocate(eta_c(nx,ny))
        allocate(dz_c(nx,ny))
        allocate(dz(nx,ny,n_lev))
        allocate(kappa(nx,ny))

! mmr recheck this 

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
            
            eta_eff(:,:) = eta(:,:,n_lev)
            
            do k = 1, n_lev-1

            eta_c = eta(:,:,n_lev-1)
            dz_c = 100.*1.e3 !dz(:,:,n_lev-k+1)
            
            eta_ratio = eta_c/eta_eff
            eta_ratiom1 = 1./eta_ratio

            c = cosh(dz_c*kappa)
            s = sinh(dz_c*kappa)

            R = (2.0 * eta_ratio * c * s + (1-eta_ratio**2) * (dz_c*kappa)**2 + (eta_ratio*s)**2 + c**2 )/&
                    ((eta_ratio + eta_ratiom1)* c * s + (eta_ratio - eta_ratiom1)*dz_c*kappa + s**2 + c**2)
            
            eta_eff = R*eta_c
            
            end do

        else

            print*,'n_lev = ', n_lev
        
            stop

        endif
     
        deallocate(xc)
        deallocate(yc)
        deallocate(R)
        deallocate(eta_ratio)
        deallocate(eta_ratiom1)
        deallocate(c)
        deallocate(s)
        deallocate(eta_c)
        deallocate(dz_c)
        deallocate(dz)
        deallocate(kappa)
     
        return
    end subroutine calc_effective_viscosity_3d

    subroutine calc_beta(beta, kappa, mu, D_lith, rho_uppermantle, g) 
        ! Calculate analytical solution as in Bueler et al 2007 (eq 11)
                    
        real(wp), intent(OUT)  :: beta(:,:) 
        real(wp), intent(OUT)  :: mu      
        real(wp), intent(IN)   :: kappa(:,:)
        real(wp), intent(IN)   :: D_lith(:,:)
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
                beta(i,j)   = rho_uppermantle*g + D_lith(i,j)*(mu**4)*kappa(i,j)**4
            end do
        end do

        return
    end subroutine calc_beta

    subroutine convenient_calc_kappa(domain)
        implicit none
        type(isos_domain_class), intent(INOUT)  :: domain

        call calc_kappa(domain%kappa, domain%nx, domain%ny)
        return
    end subroutine convenient_calc_kappa

    subroutine calc_kappa(kappa, nx, ny)

        integer, intent(IN)     :: nx, ny
        real(wp), intent(OUT)   :: kappa(:,:)

        integer :: i, j, ic, jc, ip, iq

        kappa = 0.0

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
              kappa(i,j)  = (ip*ip + iq*iq)**0.5
           end do
        end do
        kappa(1,1) = (kappa(1,2) + kappa(2,1)) / 2.0

        return
    end subroutine calc_kappa



    ! http://www.fftw.org/fftw3_doc/The-Discrete-Hartley-Transform.html

    ! The discrete Hartley transform (DHT) is an invertible linear transform closely
    ! related to the DFT. In the DFT, one multiplies each input by cos - i * sin
    ! (a complex exponential), whereas in the DHT each input is multiplied by simply cos +
    ! sin. Thus, the DHT transforms n real numbers to n real numbers, and has the convenient
    ! property of being its own inverse. In FFTW, a DHT (of any positive n) can be specified
    ! by an r2r kind of FFTW_DHT.

    ! Like the DFT, in FFTW the DHT is unnormalized, so computing a DHT of size n followed
    ! by another DHT of the same size will result in the original array multiplied by n.

    ! The DHT was originally proposed as a more efficient alternative to the DFT for real
    ! data, but it was subsequently shown that a specialized DFT (such as FFTW’s r2hc or r2c
    ! transforms) could be just as fast. In FFTW, the DHT is actually computed by
    ! post-processing an r2hc transform, so there is ordinarily no reason to prefer it from
    ! a performance perspective. However, we have heard rumors that the DHT might be the
    ! most appropriate transform in its own right for certain applications, and we would be
    ! very interested to hear from anyone who finds it useful.

    subroutine calc_fft_forward_r2r(plan, in, out)

        implicit none 
        type(c_ptr), intent(IN)     :: plan
        real(wp), intent(INOUT)     :: in(:,:)
        real(wp), intent(INOUT)     :: out(:,:)

        call fftw_execute_r2r(plan, in, out)

        return
      
    end subroutine calc_fft_forward_r2r

    subroutine calc_fft_backward_r2r(plan, in, out)

        implicit none 

        type(c_ptr), intent(IN)     :: plan
        real(wp), intent(INOUT)     :: in(:,:)
        real(wp), intent(INOUT)     :: out(:,:)

        integer(kind=4)             :: m, n

        m = size(in, 1)
        n = size(in, 2)

        call fftw_execute_r2r(plan, in, out)
        out = out / (m*n*1.)

        return
    end subroutine calc_fft_backward_r2r


    !!!!!!
    ! http://www.fftw.org/fftw3_doc/Real_002ddata-DFTs.html
    ! Fftw computes an unnormalized transform: computing an r2c
    ! followed by a c2r transform (or vice versa) will result in the
    ! original data multiplied by the size of the transform (the
    ! product of the logical dimensions). An r2c transform produces
    ! the same output as a FFTW_FORWARD complex DFT of the same input,
    ! and a c2r transform is correspondingly equivalent to
    ! FFTW_BACKWARD.
      
    subroutine calc_fft_forward_r2c(plan, in, out)

        implicit none 

        type(c_ptr), intent(IN)     :: plan
        real(wp), intent(INOUT)     :: in(:,:)
        complex(wp), intent(INOUT)  :: out(:,:)

        call fftw_execute_dft_r2c(plan, in, out)

        return
    end subroutine calc_fft_forward_r2c

    subroutine calc_fft_backward_c2r(plan, in, out)
        implicit none

        type(c_ptr), intent(IN)    :: plan
        complex(wp), intent(INOUT) :: in(:,:)
        real(wp), intent(INOUT)    :: out(:,:)

        integer(kind=4)            :: m,n

        m = size(in,1)
        n = size(in,2)
        call fftw_execute_dft_c2r(plan, in, out)
        out = out / (m*n*1.)

        return
    end subroutine calc_fft_backward_c2r

end module solver_lv_elva