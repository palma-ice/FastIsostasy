!mmr2

module solver_lv_elva

    use isostasy_defs, only : sp, dp, wp, pi 
    use solver_elva !, only : fill_with
    
    implicit none

    private
    public :: calc_lv_asthenosphere_viscous_square
!    public :: calc_lv_asthenosphere_viscous
    public :: calc_gaussian_viscosity
  contains


    
    subroutine calc_lv_asthenosphere_viscous_square(dzbdt,w,q,nu,D_lith,eta,rho_a,g,dx,dt)
    
        ! Extend a given domain [nx,ny] so that it is square based on the 
        ! largest dimension of original data. 
        ! Solve calc_asthenosphere_viscous on the square arrays, then 
        ! extract solution onto original domain [nx,ny].
      
        implicit none

        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: w(:,:)    
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: nu  
        real(wp), intent(IN)    :: D_lith
        real(wp), intent(IN)    :: eta(:,:) !mmr2                  ! [Pa s] Viscosity, eta=1e21 by default.
        real(wp), intent(IN)    :: rho_a 
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer :: i, j, nx, ny
        integer :: nsq
        real(wp), allocatable :: sq_dzbdt(:,:)
        real(wp), allocatable :: sq_w(:,:)    
        real(wp), allocatable :: sq_q(:,:)    

        ! Step 0: determine size of square array and allocate variables

        nx  = size(dzbdt,1)
        ny  = size(dzbdt,2)
        nsq = max(nx,ny)
        
        allocate(sq_dzbdt(nsq,nsq))
        allocate(sq_w(nsq,nsq))
        allocate(sq_q(nsq,nsq))

        ! Step 1: populate variables on a square grid
        
        sq_dzbdt = 0.0 
        call extend_array(sq_w,w,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_q,q,fill_with="mirror",val=0.0_wp)

        ! Step 2: solve

! recheck note: we could give the case and choose here        call calc_asthenosphere_viscous(sq_dzbdt,sq_w,sq_q,D_lith,eta,rho_a,g,dx,dt)
        call calc_lv_asthenosphere_viscous(sq_dzbdt,sq_w,sq_q,nu,D_lith,eta,rho_a,g,dx,dt)

        ! Step 3: get solution on original grid

        call reduce_array(dzbdt,sq_dzbdt)
        call reduce_array(w,sq_w)
        
        return

    end subroutine calc_lv_asthenosphere_viscous_square


   subroutine calc_lv_asthenosphere_viscous(dzbdt,u,q,nu,D_lith,eta,rho_a,g,dx,dt)
        ! Calculate rate of change of vertical bedrock height
        ! from a viscous half-space asthenosphere.
        ! Contains viscous component only.

        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        real(wp), parameter :: epsilon = 1.e-2 ! hereiam 1.e-2 ! hereiam
        
        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: u(:,:) 
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: nu
        real(wp), intent(IN)    :: D_lith
        real(wp), intent(IN)    :: eta(:,:)   ! [Pa s] Viscosity, eta=1e21 by default. ! recheck - laterally variable
        real(wp), intent(IN)    :: rho_a 
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer  :: l, m, i, j
        real(wp) :: sec_per_year, dy ! mmr2 recheck dy should be input
        logical  :: fft_r2r, fft_c2c        
        
        real(wp), allocatable :: kappa(:,:)
        real(wp), allocatable :: kappa_p(:,:)
        real(wp), allocatable :: kappa_q(:,:)
        real(wp), allocatable :: beta(:,:)
        real(wp) :: mu 

        real(wp), allocatable :: f(:,:)        
        real(wp), allocatable :: f_hat(:,:)    
        real(wp), allocatable :: prod(:,:)     
        real(wp), allocatable :: prod_hat(:,:) 
        real(wp), allocatable :: mask(:,:)     

        real(wp), allocatable :: u0(:,:)
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
        
        real(wp), allocatable :: mxx(:,:)
        real(wp), allocatable :: mxy(:,:)
        real(wp), allocatable :: myy(:,:)

        real(wp), allocatable :: mxx_hat(:,:)
        real(wp), allocatable :: mxy_hat(:,:)
        real(wp), allocatable :: myy_hat(:,:)


        real(wp), allocatable :: mm(:,:)
        real(wp), allocatable :: mm_x(:,:)
        real(wp), allocatable :: mm_y(:,:)
        real(wp), allocatable :: mm_xx(:,:)
        real(wp), allocatable :: mm_xy(:,:)
        real(wp), allocatable :: mm_yx(:,:)
        real(wp), allocatable :: mm_yy(:,:)


        real(wp), allocatable :: mm_xx_hat(:,:)
        real(wp), allocatable :: mm_xy_hat(:,:)
        real(wp), allocatable :: mm_yy_hat(:,:)


        character :: boundaries

        sec_per_year = 3600.0*24.0*365.0 ![s]

        dy = dx !mmr2 recheck dy should be input

        l = size(dzbdt,1)
        m = size(dzbdt,2)


        allocate(kappa(l,m))
        allocate(kappa_p(l,m))
        allocate(kappa_q(l,m))
        allocate(beta(l,m))

        
        allocate(u0(l,m))
        allocate(u_hat(l,m))
        allocate(f(l,m))
        allocate(f_hat(l,m))
        allocate(prod(l,m))
        allocate(prod_hat(l,m))
        allocate(mask(l,m))
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


        allocate(mm(l,m))
        allocate(mm_x(l,m))
        allocate(mm_y(l,m))


        allocate(mm_xx(l,m))
        allocate(mm_xy(l,m))
        allocate(mm_yy(l,m))


        allocate(mxx_hat(l,m))
        allocate(mxy_hat(l,m))
        allocate(myy_hat(l,m))

        allocate(mm_xx_hat(l,m))
        allocate(mm_xy_hat(l,m))
        allocate(mm_yy_hat(l,m))
        
       
        ! Initialize

        u0 = u

 !       print*,'hola u0', u0(101,101)

        
        call calc_asthenosphere_viscous_params(kappa,kappa_p,kappa_q,beta,mu,D_lith,rho_a,g,dx)  ! recheck this belongs out of here (once)

!!!!!!!!!!!!!!!
!        print*,'Tell Jan I think we can make it much simpler because we have a prognostic equation'
!hereiam        stop
!!!!!!!!!!!!!!!
        
        
! Calculate first derivatives wrt x and y
        
!        call calc_horizontal_derivatives_2D(u_x,u_y,u,dx,dy)    ! recheck this calculation

! Calculate second derivatives wrt x and y. Note u_xy should be ideally identical to u_yx
        
!        call calc_horizontal_derivatives_2D(u_xx,u_xy,u_x,dx,dy)
!        call calc_horizontal_derivatives_2D(u_yx,u_yy,u_y,dx,dy)
        

        call calc_fft_forward_r2r(u,u_hat)

        u_xx_hat = kappa_p**2 * mu**2 * u_hat
        u_yy_hat = kappa_q**2 * mu**2 * u_hat
        u_xy_hat = kappa * mu**2 * u_hat
        
        call calc_fft_backward_r2r(u_xx_hat,u_xx)
        call calc_fft_backward_r2r(u_yy_hat,u_yy)
        call calc_fft_backward_r2r(u_xy_hat,u_xy)
         

        
!  Ventsel and Krauthammer (2001): Thin Plates and Shells. Theory, Analysis, and Applications

        
        mxx = -D_lith*(u_xx + nu*u_yy) 
        myy = -D_lith*(u_yy + nu*u_xx) 
        mxy = -D_lith*(1.0-nu)*u_xy    
 
        
        ! mm = mxx
        ! call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        ! call calc_horizontal_derivatives_2D(mm_xx,mm_xy,mm_x,dx,dy)

        ! mm_xy  = 0.
        ! mm = mxy
        ! call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        ! call calc_horizontal_derivatives_2D(mm_xx,mm_xy,mm_x,dx,dy)

        ! mm = myy
        ! call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        ! call calc_horizontal_derivatives_2D(mm_yx,mm_yy,mm_y,dx,dy)

        call calc_fft_forward_r2r(mxx,mxx_hat)
        call calc_fft_forward_r2r(myy,myy_hat)
        call calc_fft_forward_r2r(mxy,mxy_hat)

        mm_xx_hat = kappa_p**2 * mu**2 * mxx_hat
        mm_yy_hat = kappa_q**2 * mu**2 * myy_hat
        mm_xy_hat = kappa * mu**2 * mxy_hat
        
        call calc_fft_backward_r2r(mm_xx_hat,mm_xx)
        call calc_fft_backward_r2r(mm_yy_hat,mm_yy)
        call calc_fft_backward_r2r(mm_xy_hat,mm_xy)
                         
        
        f = -q - rho_a*g*u + mm_xx + 2.0*mm_xy + mm_yy ! recheck mmr2 note sign criteria has changed wrt ELVA
        

        call calc_fft_forward_r2r(f,f_hat)

        kappa(1,1) = 1.         
        prod_hat = f_hat/kappa/mu

        call calc_fft_backward_r2r(prod_hat,prod)

       
        dudt = prod/(2.0*eta)  ! [m/s]
        
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

        deallocate(mm)
        deallocate(mm_x)
        deallocate(mm_y)

        deallocate(mm_xx)
        deallocate(mm_xy)
        deallocate(mm_yy)

        
        return

    end subroutine calc_lv_asthenosphere_viscous


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


      ! subroutine calc_effective_viscosity(eta_eff,eta_ratio,nlevels)

      !   real(wp), intent(OUT)    :: eta_eff(:,:)
      !   real(wp), intent(IN)     :: eta_ratio(:,:)
      !   integer(wp), intent(IN)  :: nlevels

      !   ! Local variables 
      !   integer :: i, j, nx, ny, l

      !   nx = size(eta_ratio,1)
      !   ny = size(eta_ratio,2) 

      !   eta_eff(:,:,nlevels) = eta(:,:,
      !   do l = 1, nlevels

      !      eta_ratio =  eta(:,:,l)/ eta(:,:,l+1)/
      !      calc_effective_viscosity_scaling(r,eta_ratio,tc,kappa,level)

      !      eta_eff(k) = eta_eff
      !   enddo
        
      ! end subroutine calc_effective_viscosity


      subroutine calc_gaussian_viscosity(eta,dx,dy)

        real(wp), intent(IN) :: dx, dy
        real(wp), intent(INOUT) :: eta(:,:)
        real(wp) :: Lx, Ly, L, det_eta_sigma
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


        Lx = xmax - xmin
        Ly = ymax - ymin
        L = (Lx + Ly) / 2.0

        eta_sigma_m1 = reshape ([ 0.5*(4./L)**2,  0.0_wp &
             , 0.0_wp,  0.5*(4./L)**2], shape = shape(eta_sigma_m1))
        
        det_eta_sigma = (L/4.)**2

        do i = 1, nx
           do j = 1, ny
              eta(i,j) =  1.e+21 * exp(-0.5*dot_product([xc(i),yc(j)]-[xcntr,ycntr], matmul(eta_sigma_m1, [xc(i),yc(j)]-[xcntr,ycntr]) )) ! &/ (2.0*pi*det_eta_sigma**0.5)
!              if (xc(i).eq.xcntr.and.yc(j).eq.ycntr) then
!                 print*,'hola', eta(i,j)/1.e21
!                 stop
!               endif
           enddo
        enddo

        
        return
        
      end subroutine calc_gaussian_viscosity
        
      subroutine calc_effective_viscosity_scaling(r,eta_ratio,tc,kappa,level)

        implicit none

        real(wp), intent(OUT)    :: r(:,:) 
        real(wp), intent(IN)     :: eta_ratio(:,:)
        real(wp), intent(IN)     :: tc(:,:) 
        real(wp), intent(IN)     :: kappa 
        integer(wp), intent(IN)  :: level
       
        real(wp), allocatable ::  eta_ratiom1(:,:)
        real(wp), allocatable ::  c(:,:)
        real(wp), allocatable ::  s(:,:)

        ! Local variables 
        integer :: i, j, nx, ny 

        nx = size(eta_ratio,1)
        ny = size(eta_ratio,2) 

        allocate(eta_ratiom1(nx,ny))
        allocate(c(nx,ny))
        allocate(s(nx,ny))
        
        c = cosh(tc*kappa)
        s = sinh(tc*kappa)
        eta_ratiom1 = 1./eta_ratio
        
        r = (2.0 * eta_ratio * s + (1-eta_ratio**2) * (tc*kappa)**2 + eta_ratio**2 * s**2 + c**2 )/&
             (eta_ratio + eta_ratiom1)* c * s + (eta_ratio + eta_ratiom1)*(tc*kappa) + s**2 + c**2

        return
        
      end subroutine calc_effective_viscosity_scaling
      
 end module solver_lv_elva

 !mmr2


    
