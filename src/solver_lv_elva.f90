
module solver_lv_elva

    use isostasy_defs, only : sp, dp, wp, pi 
!    use solver_elva !, only : fill_with
    
    implicit none

    private
    public :: calc_lv_asthenosphere_viscous_square
    public :: calc_gaussian_viscosity
    public :: calc_gaussian_rigidity
    public :: calc_effective_viscosity_3layer_channel
    public :: calc_effective_viscosity_3d
    public :: extend_array 
    public :: reduce_array
    public :: calc_asthenosphere_viscous_params
    public :: calc_fft_backward_r2r
    Public :: calc_fft_forward_r2r
    public :: calc_fft_backward_c2c
    public :: calc_fft_forward_c2c

  contains

    
    subroutine calc_lv_asthenosphere_viscous_square(dzbdt,w,w_el,q,nu,D_lith,eta,rho_a,rho_l,g,dx,dt)
    
        ! Extend a given domain [nx,ny] so that it is square based on the 
        ! largest dimension of original data. 
        ! Solve calc_asthenosphere_viscous on the square arrays, then 
        ! extract solution onto original domain [nx,ny].
      
        implicit none

        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: w(:,:)
        real(wp), intent(INOUT) :: w_el(:,:)    
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: nu  
        real(wp), intent(IN)    :: D_lith(:,:)
        real(wp), intent(IN)    :: eta(:,:)   
        real(wp), intent(IN)    :: rho_a
        real(wp), intent(IN)    :: rho_l 
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer :: i, j, nx, ny
        integer :: nsq
        real(wp), allocatable :: sq_dzbdt(:,:)
        real(wp), allocatable :: sq_w(:,:)
        real(wp), allocatable :: sq_w_el(:,:)    
        real(wp), allocatable :: sq_q(:,:)
        real(wp), allocatable :: sq_D_lith(:,:)
        real(wp), allocatable :: sq_eta(:,:)    

        ! Step 0: determine size of square array and allocate variables

        nx  = size(dzbdt,1)
        ny  = size(dzbdt,2)
        nsq = max(nx,ny) 
        
        allocate(sq_dzbdt(nsq,nsq))
        allocate(sq_w(nsq,nsq))
        allocate(sq_w_el(nsq,nsq))
        allocate(sq_q(nsq,nsq))
        allocate(sq_D_lith(nsq,nsq))
        allocate(sq_eta(nsq,nsq))

        ! Step 1: populate variables on a square grid
        
        sq_dzbdt = 0.0 
        call extend_array(sq_dzbdt,dzbdt,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_w,w,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_w_el,w_el,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_q,q,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_D_lith,D_lith,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_eta,eta,fill_with="mirror",val=0.0_wp)

        ! Step 2: solve

        call calc_lv_asthenosphere_viscous(sq_dzbdt,sq_w,sq_w_el,sq_q,nu,sq_D_lith,sq_eta,rho_a,rho_l,g,dx,dt)

        ! Step 3: get solution on original grid

        call reduce_array(dzbdt,sq_dzbdt)
        call reduce_array(w,sq_w)
        call reduce_array(w_el,sq_w_el)
        
        return

    end subroutine calc_lv_asthenosphere_viscous_square


   subroutine calc_lv_asthenosphere_viscous(dzbdt,u,u_el,q,nu,D_lith,eta,rho_a,rho_l,g,dx,dt)
        ! Calculate rate of change of vertical bedrock height
        ! from a viscous half-space asthenosphere.
        ! Contains viscous component only.

        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        real(wp), parameter :: epsilon = 1.e-2 
        
        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(IN)    :: u(:,:)
        real(wp), intent(IN)    :: u_el(:,:)
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: nu
        real(wp), intent(IN)    :: D_lith(:,:) 
        real(wp), intent(IN)    :: eta(:,:)   ! [Pa s] Viscosity, eta=1e21 by default. 
        real(wp), intent(IN)    :: rho_a
        real(wp), intent(IN)    :: rho_l
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer  :: l, m, i, j
        real(wp) :: sec_per_year, dy ! mmr recheck - dy should be input
        logical  :: fft_r2r, fft_c2c, finite_difs,  finite_difs_1, finite_difs_2
        
        real(wp), allocatable :: kappa(:,:)
        real(wp), allocatable :: kappa_p(:,:)
        real(wp), allocatable :: kappa_q(:,:)
        real(wp), allocatable :: beta(:,:)
        real(wp) :: mu 

        real(wp), allocatable :: f(:,:)        
        real(wp), allocatable :: f_hat(:,:)    

        real(wp), allocatable :: prod(:,:)     
        real(wp), allocatable :: prod_hat(:,:)
        
        real(wp), allocatable :: p(:,:)

        real(wp), allocatable :: u_hat(:,:)
        
        real(wp), allocatable :: dudt(:,:)     
        real(wp), allocatable :: u_x(:,:)     
        real(wp), allocatable :: u_y(:,:)     

        real(wp), allocatable :: u_xx(:,:)
        real(wp), allocatable :: u_xy(:,:)
        real(wp), allocatable :: u_yy(:,:)
        real(wp), allocatable :: u_yx(:,:)

        real(wp), allocatable :: u_xx_hat(:,:)
        real(wp), allocatable :: u_xy_hat(:,:)
        real(wp), allocatable :: u_yy_hat(:,:)

        real(wp), allocatable :: Mx(:,:) 
        real(wp), allocatable :: My(:,:)

        real(wp), allocatable :: Mxy(:,:)
        real(wp), allocatable :: Mxy_x(:,:)
        real(wp), allocatable :: Myx_y(:,:)

        real(wp), allocatable :: Mx_hat(:,:)
        real(wp), allocatable :: My_hat(:,:)
        real(wp), allocatable :: Mxy_hat(:,:)

        real(wp), allocatable :: Mx_xx(:,:)
        real(wp), allocatable :: My_yy(:,:)
        
        real(wp), allocatable :: Mxy_xy(:,:)
        real(wp), allocatable :: Mxy_yx(:,:)

        real(wp), allocatable :: Mx_xx_hat(:,:)
        real(wp), allocatable :: Mxy_xy_hat(:,:)
        real(wp), allocatable :: My_yy_hat(:,:)


        character :: boundaries

        sec_per_year = 3600.0*24.0*365.0 ![s]

        dy = dx 

        l = size(dzbdt,1)
        m = size(dzbdt,2)


        allocate(kappa(l,m))
        allocate(kappa_p(l,m))
        allocate(kappa_q(l,m))
        allocate(beta(l,m))

        
        allocate(u_hat(l,m))
        allocate(f(l,m))
        allocate(f_hat(l,m))
        
        allocate(prod(l,m))
        allocate(prod_hat(l,m))
        
        allocate(dudt(l,m))
        
        allocate(u_x(l,m))
        allocate(u_y(l,m))
        
        allocate(u_xx(l,m))
        allocate(u_xy(l,m))
        allocate(u_yy(l,m))
        allocate(u_yx(l,m))


        allocate(u_xx_hat(l,m))
        allocate(u_xy_hat(l,m))
        allocate(u_yy_hat(l,m))


        allocate(Mx(l,m))
        allocate(My(l,m))
        
        allocate(Mxy_x(l,m))
        allocate(Myx_y(l,m))
        

        allocate(Mx_xx(l,m))
        allocate(My_yy(l,m))        
        allocate(Mxy_xy(l,m))
       
        allocate(Mx_xx_hat(l,m))
        allocate(Mxy_xy_hat(l,m))
        allocate(My_yy_hat(l,m))

        allocate(Mx_hat(l,m))
        allocate(Mxy_hat(l,m))
        allocate(My_hat(l,m))

        allocate(p(l,m))
        
       
        ! Initialize

        call calc_asthenosphere_viscous_params(kappa,kappa_p,kappa_q,beta,mu,D_lith,rho_a,g,dx)  ! recheck this belongs out of here (once)

!mmr recheck - try with finite differences too

        ! Finite differences

        finite_difs = .true.

        if (finite_difs) then

           call calc_derivative_xx(u_xx,u,dx)
           call calc_derivative_yy(u_yy,u,dy)
           
           call calc_derivative_x(u_x,u,dx)
           call calc_derivative_y(u_xy,u_x,dy)

        else
           
           call calc_fft_forward_r2r(u,u_hat)

           u_xx_hat = kappa_p**2 * mu**2 * u_hat
           u_yy_hat = kappa_q**2 * mu**2 * u_hat
           u_xy_hat = kappa * mu**2 * u_hat

           call calc_fft_backward_r2r(u_xx_hat,u_xx)
           call calc_fft_backward_r2r(u_yy_hat,u_yy)
           call calc_fft_backward_r2r(u_xy_hat,u_xy)
         
        endif
        
        !  Ventsel and Krauthammer (2001): Thin Plates and Shells.
        ! Theory, Analysis, and Applications. Eq (2.13, 2.23)

        Mx = -D_lith*(u_xx + nu*u_yy) 
        My = -D_lith*(u_yy + nu*u_xx) 
        Mxy = -D_lith*(1.0-nu)*u_xy    

         ! Finite differences

        if (finite_difs) then

           call calc_derivative_xx(Mx_xx,Mx,dx)
           call calc_derivative_yy(My_yy,My,dy)

           call calc_derivative_x(Mxy_x,Mxy,dx)
           call calc_derivative_y(Mxy_xy,Mxy_x,dy)
           
        else 

           call calc_fft_forward_r2r(Mx,Mx_hat)
           call calc_fft_forward_r2r(My,My_hat)
           call calc_fft_forward_r2r(Mxy,Mxy_hat)

           Mx_xx_hat = kappa_p**2 * mu**2 * Mx_hat
           My_yy_hat = kappa_q**2 * mu**2 * My_hat
           Mxy_xy_hat = kappa * mu**2 * Mxy_hat

           call calc_fft_backward_r2r(Mx_xx_hat,Mx_xx)
           call calc_fft_backward_r2r(My_yy_hat,My_yy)
           call calc_fft_backward_r2r(Mxy_xy_hat,Mxy_xy)
                         
     endif
  
     ! Viscous + elastic
!     p =  - q - rho_a*g*u -rho_l*g*u_el

     f = - q - rho_a*g*u + Mx_xx + 2.0*Mxy_xy + My_yy  - rho_l*g*u_el
     
     kappa(1,1) = (kappa(1,2) + kappa(2,1)) / 2.0 ! 1.         

! original case (Jans, as in main.pdf)
!        
!        call calc_fft_forward_r2r(f,f_hat)
!        prod_hat = f_hat/kappa/mu
!        call calc_fft_backward_r2r(prod_hat,prod)
!        dudt = prod/(2.0*eta)  ! [m/s]

        ! new case!
        call calc_fft_forward_r2r(f/(2.*eta),f_hat)
        prod_hat = f_hat/kappa/mu
        call calc_fft_backward_r2r(prod_hat,prod)
        dudt = prod  ! [m/s]

        dudt  = dudt - 0.25*(dudt(1,1)+dudt(l,m)+dudt(1,m)+dudt(l,1)) 
        
! Rate of viscous asthenosphere uplift per unit time (seconds)
        
        dzbdt = dudt * sec_per_year  !  [m/s] x [s/a] = [m/a] 



        
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
        deallocate(u_yx)

!        deallocate(Mx_x)
!        deallocate(My_y)

        deallocate(Mx_xx)
        deallocate(Mxy_xy)
        deallocate(My_yy)

        deallocate(p)

        
        return

    end subroutine calc_lv_asthenosphere_viscous

      
   subroutine calc_derivative_x(dudx,u,dx)
      ! Get simple horizontal derivatives wrt x
  
        implicit none

        real(wp), intent(OUT) :: dudx(:,:) 
        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: dx 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1

        nx = size(u,1)
        ny = size(u,2) 

        dudx = 0.0
            
        do j = 1, ny 
           do i = 2, nx-1

              im1 = i-1
              ip1 = i+1

              dudx(i,j) = (u(ip1,j)-u(im1,j))/(2.0*dx)

           end do
        end do

        dudx(nx,:) = (u(nx,:)-u(nx-1,:))/dx
        dudx(1,:) = (u(2,:)-u(1,:))/dx

        return

      end subroutine calc_derivative_x

         subroutine calc_derivative_y(dudy,u,dy)
      ! Get simple horizontal derivatives
  
        implicit none

        real(wp), intent(OUT) :: dudy(:,:) 
        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: dy
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        dudy = 0.0

            
        do j = 2, ny-1  
           do i = 1, nx

              jm1 = j-1
              jp1 = j+1

              dudy(i,j) = (u(i,jp1)-u(i,jm1))/(2.0*dy)

           end do
        end do

        dudy(:,ny) = (u(:,ny)-u(:,ny-1))/dy
        dudy(:,1) = (u(:,2)-u(:,1))/dy
                
        return

      end subroutine calc_derivative_y

      subroutine calc_derivative_xx(uxx,u,dx)
      ! Get simple horizontal derivatives
  
        implicit none

        real(wp), intent(OUT) :: uxx(:,:) 
        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: dx 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1

        nx = size(u,1)
        ny = size(u,2) 

        uxx = 0.0
            
        do j = 1, ny  
           do i = 2, nx-1

              im1 = i-1
              ip1 = i+1

              uxx(i,j) = (u(ip1,j)- 2*u(i,j) + u(im1,j))/(dx*dx)

           end do
        end do

        uxx(1,:) = (u(3,:)-2*u(2,:)+u(1,:))/(dx*dx)
        uxx(nx,:) = (u(nx,:) - 2*u(nx-1,:) + u(nx-2,:))/(dx*dx)
                        
        return

      end subroutine calc_derivative_xx

      subroutine calc_derivative_yy(uyy,u,dy)
      ! Get simple horizontal derivatives
  
        implicit none

        real(wp), intent(OUT) :: uyy(:,:) 
        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: dy 
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: jm1, jp1 

        nx = size(u,1)
        ny = size(u,2) 

        uyy = 0.0
            
        do j = 2, ny-1  
           do i = 1, nx


              jm1 = j-1
              jp1 = j+1

              uyy(i,j) = (u(i,jp1) - 2*u(i,j) + u(i,jm1))/(dy*dy)

           end do
        end do

        uyy(:,1) = (u(:,3)-2*u(:,2)+u(:,1))/(dy*dy)
        uyy(:,ny) = (u(:,ny) - 2*u(:,ny-1) + u(:,ny-2))/(dy*dy)
                        
        return

      end subroutine calc_derivative_yy


        subroutine calc_horizontal_derivatives_2D(dudx,dudy,u,dx,dy)
      ! Get simple horizontal derivatives
  
        implicit none

        real(wp), intent(OUT) :: dudx(:,:) 
        real(wp), intent(OUT) :: dudy(:,:) 
        real(wp), intent(IN)  :: u(:,:) 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        
        ! Local variables 
        integer :: i, j, nx, ny 
        integer :: im1, ip1, jm1, jp1 

        nx = size(dudx,1)
        ny = size(dudx,2) 

        dudx = 0.0
        dudy = 0.0

            
        do j = 2, ny-1  
           do i = 2, nx-1

              im1 = i-1
              ip1 = i+1
              jm1 = j-1
              jp1 = j+1

              dudx(i,j) = (u(ip1,j)-u(im1,j))/(2.0*dx)
              dudy(i,j) = (u(i,jp1)-u(i,jm1))/(2.0*dy)


           end do
        end do

        dudx(nx,:) = (u(nx,:)-u(nx-1,:))/dx
        dudx(1,:) = (u(2,:)-u(1,:))/dx

        dudy(:,ny) = (u(:,ny)-u(:,ny-1))/dy
        dudy(:,1) = (u(:,2)-u(:,1))/dy
                
        return

      end subroutine calc_horizontal_derivatives_2D


      subroutine calc_gaussian_viscosity(eta,eta_0,sign,dx,dy)

        real(wp), intent(OUT) :: eta(:,:)
        real(wp), intent(IN) :: eta_0, sign, dx, dy
        real(wp) :: Lx, Ly, L, det_eta_sigma, f
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin

        real(wp), allocatable :: xc(:), yc(:)

        real(wp) :: eta_sigma(2,2)
        real(wp) :: eta_sigma_m1(2,2)

        integer  :: i, j, nx, ny

        nx = size(eta,1)
        ny = size(eta,2)

        allocate(xc(nx))
        allocate(yc(ny))

        
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


        Lx = (xmax - xmin)/2.
        Ly = (ymax - ymin)/2.


        eta_sigma_m1 = reshape ([ (4.0/Lx)**2,  0.0_wp, 0.0_wp,  (4.0/Ly)**2], shape = shape(eta_sigma_m1))

        do i = 1, nx
           do j = 1, ny
              f = exp(-0.5*dot_product([xc(i),yc(j)]-[xcntr,ycntr], matmul(eta_sigma_m1, [xc(i),yc(j)]-[xcntr,ycntr]) ))              
              eta(i,j) = eta_0 * 10**(sign * f)
           enddo
        enddo

        eta = exp(log10(1.e21 / eta)) * eta
        
        return
        
      end subroutine calc_gaussian_viscosity

      subroutine calc_gaussian_rigidity(He_lith,He_lith_0, He_lith_1,sign,dx,dy)

        real(wp), intent(IN) :: He_lith_0, He_lith_1,sign, dx, dy
        real(wp), intent(OUT) :: He_lith(:,:)
        real(wp) :: Lx, Ly, L, det_He_lith_sigma
        real(wp) :: xcntr, ycntr, xmax, xmin, ymax, ymin

        real(wp), allocatable :: xc(:), yc(:)

        real(wp) :: He_lith_sigma(2,2)
        real(wp) :: He_lith_sigma_m1(2,2)

        integer  :: i, j, nx, ny

        nx = size(He_lith,1)
        ny = size(He_lith,2)

        allocate(xc(nx))
        allocate(yc(ny))

        
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

        Lx = (xmax - xmin)/2.
        Ly = (ymax - ymin)/2.

        He_lith_sigma_m1 = reshape ([ (4.0/Lx)**2,  0.0_wp, 0.0_wp,  (4.0/Ly)**2], shape = shape(he_lith_sigma_m1))

        do i = 1, nx
           do j = 1, ny
              He_lith(i,j) = He_lith_0 +   &
                   sign*He_lith_1 * exp(-0.5*dot_product([xc(i),yc(j)]-[xcntr,ycntr], matmul(He_lith_sigma_m1, [xc(i),yc(j)]-[xcntr,ycntr]) ))
           enddo
        enddo
        
        return
        
      end subroutine calc_gaussian_rigidity

        
      subroutine calc_effective_viscosity_3layer_channel(eta_eff,visc_c,thck_c,He_lith,n_lev,nu,dx,dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:,:)
        real(wp), intent(IN)     :: visc_c
        real(wp), intent(IN)     :: thck_c
        real(wp), intent(IN)     :: He_lith !(:,:) 
        integer,  intent(IN)     :: n_lev               
        real(wp), intent(IN)     :: nu, dx, dy
        
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
     
     eta_eff    = eta_eff  * (1.5/(1. + nu))         ! [Pa s]

        
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

      subroutine calc_effective_viscosity_3d(eta_eff,eta,nu,dx,dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:,:)
        real(wp), intent(IN)     :: eta(:,:,:)
        real(wp), intent(IN)     :: nu, dx, dy
        
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
     
     eta_eff    = eta_eff  * (1.5/(1. + nu))         ! [Pa s]

        
        deallocate(xc)
        deallocate(yc)        
        deallocate(R)
        deallocate(eta_ratio)
        deallocate(eta_ratiom1)
        deallocate(c)
        deallocate(s)
!        deallocate(eta)
        deallocate(eta_c)
        deallocate(dz_c)
        deallocate(dz)
        deallocate(kappa)
     
        return
        
      end subroutine calc_effective_viscosity_3d



    subroutine calc_asthenosphere_viscous_params(kappa,kappa_p,kappa_q,beta,mu,D_lith,rho_a,g,dx) 
        ! Calculate parameters needed for elastic lithosphere viscous asthenosphere (ELVA)                                        
        ! solution as in Bueler et al 2007 (eq 11)
                    
      real(wp), intent(OUT)  :: kappa(:,:)
      real(wp), intent(OUT)  :: kappa_p(:,:)
      real(wp), intent(OUT)  :: kappa_q(:,:)
        real(wp), intent(OUT)  :: beta(:,:) 
        real(wp), intent(OUT)  :: mu      
        real(wp), intent(IN)   :: D_lith(:,:)
        real(wp), intent(IN)   :: rho_a 
        real(wp), intent(IN)   :: g 
        real(wp), intent(IN)   :: dx
        
        ! Local variables
        integer  :: i, j, nx, ny 
        integer  :: ip, iq, ic, jc 
        real(wp) :: xd, yd

        nx = size(kappa,1)
        ny = size(kappa,2)

        ! Calculate mu

        mu = 2.*pi/((nx-1)*dx)

        ! Calculate kappa and beta

        kappa = 0.0
        kappa_p = 0.0
        kappa_q = 0.0
        beta  = 0.0 

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
                kappa_p(i,j) = ip 
                kappa_q(i,j) = iq 
                beta(i,j)   = rho_a*g + D_lith(i,j)*(mu**4)*kappa(i,j)**4
            end do
        end do

        return
      
    end subroutine calc_asthenosphere_viscous_params

      
      subroutine calc_kappa(field,kappa)

        ! Calculate kappa and beta

        real(wp), intent(IN)  :: field(:,:)
        real(wp), allocatable, intent(OUT) :: kappa(:,:)

        integer :: i, j, ic, jc, ip, iq, nx, ny

        nx = size(field,1)
        ny = size(field,2)

        allocate(kappa(nx,ny))
        
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


        end subroutine calc_kappa

    subroutine make_fft_plans(in,plan_fwd,plan_bck) 
 
        use, intrinsic :: iso_c_binding
        implicit none 
        include 'fftw3.f03' 

        real(wp), intent(IN)       :: in(:,:)
        type(c_ptr), intent(OUT)   :: plan_fwd, plan_bck
        complex (dp), allocatable  :: in_aux(:,:)
        complex (dp), allocatable  :: out_aux(:,:)
        integer (kind = 4)         :: l, m


        l    = size(in,1)
        m    = size(in,2)

        if(l.ne.m) then
            write(*,*) "Dimensions do not match, stopping now"
            stop
        end if

        in_aux = in

! r2r      
        plan_fwd = fftw_plan_dft_2d(l,m,in_aux,out_aux,-1,FFTW_ESTIMATE)
        plan_bck = fftw_plan_dft_2d(l,m,out_aux,in_aux,+1,FFTW_ESTIMATE)

        return

    end subroutine make_fft_plans

    subroutine calc_fft_forward_r2r(in,out)

        use, intrinsic :: iso_c_binding
        implicit none 
        include 'fftw3.f03'  

        real(wp), intent(IN)       :: in(:,:)
        real(wp), intent(OUT)      :: out(:,:)

        real(wp), allocatable      :: rec(:,:)
        real(dp), allocatable      :: in_aux(:,:)
        real(dp), allocatable      :: out_aux(:,:) 
        real(dp), allocatable      :: rec_aux(:,:)
        real(dp)                   :: dx, cc
        type(c_ptr)                :: plan
        integer(kind=4)            :: m,n
        logical                    :: print_check


        m = size(in,1)
        n = size(in,2)

        print_check = .false.

        allocate(in_aux(m,n))
        allocate(out_aux(m,n))
        allocate(rec(m,n))
        allocate(rec_aux(m,n))


! http://www.fftw.org/fftw3_doc/The-Discrete-Hartley-Transform.html
      
!      The discrete Hartley transform (DHT) is an invertible linear
!      transform closely related to the DFT. In the DFT, one
!      multiplies each input by cos - i * sin (a complex exponential),
!      whereas in the DHT each input is multiplied by simply cos +
!      sin. Thus, the DHT transforms n real numbers to n real numbers,
!      and has the convenient property of being its own inverse. In
!      FFTW, a DHT (of any positive n) can be specified by an r2r kind
!      of FFTW_DHT.

!      Like the DFT, in FFTW the DHT is unnormalized, so computing a
!      DHT of size n followed by another DHT of the same size will
!      result in the original array multiplied by n.

!      The DHT was originally proposed as a more efficient alternative
!      to the DFT for real data, but it was subsequently shown that a
!      specialized DFT (such as FFTW’s r2hc or r2c transforms) could
!      be just as fast. In FFTW, the DHT is actually computed by
!      post-processing an r2hc transform, so there is ordinarily no
!      reason to prefer it from a performance perspective.5 However,
!      we have heard rumors that the DHT might be the most appropriate
!      transform in its own right for certain applications, and we
!      would be very interested to hear from anyone who finds it
!      useful.


        in_aux = in 
        plan = fftw_plan_r2r_2d(m,n,in_aux,out_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)

        call fftw_execute_r2r(plan, in_aux, out_aux)
        out = real(out_aux/sqrt(m*n*1.))

        if (print_check) then
            call r4mat_print_some ( m, n, in, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the original data:' )
            plan =  fftw_plan_r2r_2d(m,n,out_aux,rec_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)
            call fftw_execute_r2r(plan, out_aux, rec_aux)
            rec_aux = rec_aux/(m*n)
            rec = real(rec_aux,wp) 
            call fftw_destroy_plan(plan)
            call r4mat_print_some (m, n, rec, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the recovered data:' )
        end if

        deallocate(in_aux)
        deallocate(out_aux)
        deallocate(rec)
        deallocate(rec_aux)
      
        return
      
    end subroutine calc_fft_forward_r2r

    subroutine calc_fft_backward_r2r(in,out)

        use, intrinsic :: iso_c_binding
        implicit none 
        include 'fftw3.f03'  

        real(wp), intent(IN)       :: in(:,:)
        real(wp), intent(OUT)      :: out(:,:)

        real(wp), allocatable      :: rec(:,:)
        real(dp), allocatable      :: in_aux(:,:)
        real(dp), allocatable      :: out_aux(:,:) 
        real(dp), allocatable      :: rec_aux(:,:)
        real(dp)                   :: dx, cc
        type(c_ptr)                :: plan
        integer(kind=4)            :: m,n
        logical                    :: print_check

        m = size(in,1)
        n = size(in,2)


        allocate(in_aux(m,n))
        allocate(out_aux(m,n))
        allocate(rec(m,n))
        allocate(rec_aux(m,n))


        in_aux = in
        plan = fftw_plan_r2r_2d(m,n,in_aux,out_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)

        call fftw_execute_r2r(plan, in_aux, out_aux)
        out = real(out_aux/sqrt(m*n*1.)) 

        if (print_check) then
            call r4mat_print_some ( m, n, in, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the original data:' )
            plan =  fftw_plan_r2r_2d(m,n,out_aux,rec_aux,FFTW_DHT,FFTW_DHT,FFTW_ESTIMATE)
            call fftw_execute_r2r(plan, out_aux, rec_aux)
            rec_aux = rec_aux/(m*n)
            rec = real(rec_aux,wp) 
            call fftw_destroy_plan(plan)
            call r4mat_print_some ( m, n, rec, n/2-2, m/2-2, n/2+2, m/2+2, '  Part of the recovered data:' )
        end if

        deallocate(in_aux)
        deallocate(out_aux)
        deallocate(rec)
        deallocate(rec_aux)

        return
    
    end subroutine calc_fft_backward_r2r

    
    ! http://www.fftw.org/fftw3_doc/Real_002ddata-DFTs.html
    ! FFTW computes an unnormalized transform: computing an r2c
    ! followed by a c2r transform (or vice versa) will result in the
    ! original data multiplied by the size of the transform (the
    ! product of the logical dimensions). An r2c transform produces
    ! the same output as a FFTW_FORWARD complex DFT of the same input,
    ! and a c2r transform is correspondingly equivalent to
    ! FFTW_BACKWARD.
    
 
    subroutine calc_fft_forward_c2c(in,out)

        use, intrinsic :: iso_c_binding
        implicit none 
        include 'fftw3.f03'  

        real(wp), intent(IN)       :: in(:,:)
        real(wp), intent(OUT)      :: out(:,:)

        real(wp), allocatable      :: rec(:,:)
        complex(dp), allocatable   :: in_aux(:,:)   
        complex(dp), allocatable   :: out_aux(:,:)     
        real(dp), allocatable      :: rec_aux(:,:)
        real(dp)                   :: dx, cc
        type(c_ptr)                :: plan
        integer(kind=4)            :: m, n

        m = size(in,1)
        n = size(in,2)

      
        allocate(in_aux(m,n))
        allocate(out_aux(m,n))
        allocate(rec(m,n))
        allocate(rec_aux(m,n))


        in_aux = in
        plan = fftw_plan_dft_2d(m,n,in_aux,out_aux,-1,FFTW_ESTIMATE)

        call fftw_execute_dft(plan, in_aux, out_aux)                 
        call fftw_destroy_plan(plan)
        out = real(out_aux/sqrt(m*n*1.)) 


        deallocate(in_aux)
        deallocate(out_aux)
        deallocate(rec)
        deallocate(rec_aux)

        return
      
    end subroutine calc_fft_forward_c2c

    subroutine calc_fft_backward_c2c(in,out)

        use, intrinsic :: iso_c_binding
        implicit none 
        include 'fftw3.f03'  

        real(wp), intent(IN)       :: in(:,:)
        real(wp), intent(OUT)      :: out(:,:)

        real(wp), allocatable      :: rec(:,:)
        complex(dp), allocatable   :: in_aux(:,:)   
        complex(dp), allocatable   :: out_aux(:,:)     
        real(dp), allocatable      :: rec_aux(:,:)
        real(dp)                   :: dx, cc
        type(c_ptr)                :: plan
        integer(kind=4)            :: m, n

        m = size(in,1)
        n = size(in,2)


        allocate(in_aux(m,n))
        allocate(out_aux(m,n))
        allocate(rec(m,n))
        allocate(rec_aux(m,n))


        in_aux = in
        plan = fftw_plan_dft_2d(m,n,in_aux,out_aux,1,FFTW_ESTIMATE)

        call fftw_execute_dft(plan, in_aux, out_aux)                 
        call fftw_destroy_plan(plan)
        out = real(out_aux/sqrt(m*n*1.)) 


        deallocate(in_aux)
        deallocate(out_aux)
        deallocate(rec)
        deallocate(rec_aux)

        return

    end subroutine calc_fft_backward_c2c



    ! ===== ARRAY SIZING FUNCTIONS ==============================

    subroutine extend_array(ve,v,fill_with,val)
        ! Given an array ve with dimensions >= those of v, 
        ! populate ve with values of v, centered in the larger array.
        ! Fill extra cells with value of user's choice. 

        implicit none

        real(wp), intent(INOUT) :: ve(:,:) 
        real(wp), intent(IN)    :: v(:,:)
        character(len=*), intent(IN) :: fill_with
        real(wp), intent(IN), optional :: val

        ! Local variables
        integer  :: i, j, nx, ny, nx1, ny1
        integer  :: i0, i1, j0, j1
        real(wp) :: fill_value 

        nx = size(v,1)
        ny = size(v,2) 

        nx1 = size(ve,1)
        ny1 = size(ve,2)

        ! Determine fill value to be used

        select case(trim(fill_with))

            case("zero","zeros")

                fill_value = 0.0

            case("mean")

                fill_value = sum(v) / real(nx*ny,wp)

            case("val")

                if (.not. present(val)) then
                    write(*,*) "extend_array:: Error: for fill_with='val', the optional &
                    &argument val must be provided and is missing right now."
                    stop
                end if

                fill_value = val

            case("mirror")
                ! Set fill_value to zero, it will not be used 

                fill_value = 0.0 

            case DEFAULT
                write(*,*) "extend_array:: Error: choice of 'fill_with' not recognized."
                write(*,*) "fill_with = ", trim(fill_with)
                stop
        
        end select

        ! Initialize extended array to correct fill value

        ve = fill_value 

        ! Fill in the actual values centered within the extended array

        if (nx .eq. nx1) then
            i0 = 1
        else
            i0 = floor( (nx1-nx)/2.0 )
        end if 
        i1 = i0+nx-1

        if (ny .eq. ny1) then
            j0 = 1
        else
            j0 = floor( (ny1-ny)/2.0 )
        end if 
        j1 = j0+ny-1

        ve(i0:i1,j0:j1) = v


        if (trim(fill_with) .eq. "mirror") then
            ! Populate extended array region with mirrored points 

            ! Left
            if (i0 .gt. 1) then 
                ve(1:i0-1,j0:j1) = v(i0-1:1:-1,:)
            end if
            ! Right
            if (i1 .lt. nx1) then 
                ve(i1+1:nx1,j0:j1) = v(nx:nx-(nx1-i1)+1:-1,:)
            end if
            
            ! Bottom
            if (j0 .gt. 1) then 
                ve(i0:i1,1:j0-1) = v(:,j0-1:1:-1)
            end if
            ! Top
            if (j1 .lt. ny1) then 
                ve(i0:i1,j1+1:ny1) = v(:,ny:ny-(ny1-j1)+1:-1)
            end if
            
            ! To do - populate the corners too. For now ignore, and keep 
            ! extended array values equal to fill_value=0.0. 

        end if 

        return
    
    end subroutine extend_array

    subroutine reduce_array(v,ve)
        ! Given an array of dimensions >= those of v
        ! extract centered values of ve of interest,
        ! discarding remaining values around the borders.

        implicit none

        real(wp), intent(INOUT) :: v(:,:)
        real(wp), intent(IN)    :: ve(:,:) 
        
        ! Local variables
        integer  :: i, j, nx, ny, nx1, ny1
        integer  :: i0, i1, j0, j1
        real(wp) :: fill_value 

        nx = size(v,1)
        ny = size(v,2) 

        nx1 = size(ve,1)
        ny1 = size(ve,2)

        ! Fill in the actual values from those centered within the extended array

        if (nx .eq. nx1) then
            i0 = 1
        else
            i0 = floor( (nx1-nx)/2.0 )
        end if 
        i1 = i0+nx-1

        if (ny .eq. ny1) then
            j0 = 1
        else
            j0 = floor( (ny1-ny)/2.0 )
        end if 
        j1 = j0+ny-1

        v = ve(i0:i1,j0:j1)

        return
    
    end subroutine reduce_array
    
 end module solver_lv_elva


    
