
module solver_lv_elva

    use isostasy_defs, only : sp, dp, wp, pi 
    use solver_elva !, only : fill_with
    
    implicit none

    private
    public :: calc_lv_asthenosphere_viscous_square
    public :: calc_gaussian_viscosity
    public :: calc_gaussian_rigidity
    public :: calc_effective_viscosity
  contains


    
    subroutine calc_lv_asthenosphere_viscous_square(dzbdt,w,w_el,q,nu,D_lith,eta,rho_a,g,dx,dt)
    
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

        ! Step 0: determine size of square array and allocate variables

        nx  = size(dzbdt,1)
        ny  = size(dzbdt,2)
        nsq = max(nx,ny)
        
        allocate(sq_dzbdt(nsq,nsq))
        allocate(sq_w(nsq,nsq))
        allocate(sq_w_el(nsq,nsq))
        allocate(sq_q(nsq,nsq))

        ! Step 1: populate variables on a square grid
        
        sq_dzbdt = 0.0 
        call extend_array(sq_w,w,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_w_el,w_el,fill_with="mirror",val=0.0_wp)
        call extend_array(sq_q,q,fill_with="mirror",val=0.0_wp)

        ! Step 2: solve

        call calc_lv_asthenosphere_viscous(sq_dzbdt,sq_w,sq_w_el,sq_q,nu,D_lith,eta,rho_a,g,dx,dt)

        ! Step 3: get solution on original grid

        call reduce_array(dzbdt,sq_dzbdt)
        call reduce_array(w,sq_w)
        call reduce_array(w_el,sq_w_el)
        
        return

    end subroutine calc_lv_asthenosphere_viscous_square


   subroutine calc_lv_asthenosphere_viscous(dzbdt,u,u_el,q,nu,D_lith,eta,rho_a,g,dx,dt)
        ! Calculate rate of change of vertical bedrock height
        ! from a viscous half-space asthenosphere.
        ! Contains viscous component only.

        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        real(wp), parameter :: epsilon = 1.e-2, rho_l = 2600.0 ! recheck - input in namelist
        
        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: u(:,:)
        real(wp), intent(IN)    :: u_el(:,:)
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: nu
        real(wp), intent(IN)    :: D_lith(:,:) 
        real(wp), intent(IN)    :: eta(:,:)   ! [Pa s] Viscosity, eta=1e21 by default. 
        real(wp), intent(IN)    :: rho_a 
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer  :: l, m, i, j
        real(wp) :: sec_per_year, dy ! mmr recheck - dy should be input
        logical  :: fft_r2r, fft_c2c, finite_difs        
        
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


        real(wp), allocatable :: p(:,:)


        character :: boundaries

        sec_per_year = 3600.0*24.0*365.0 ![s]

        dy = dx !mmr - recheck dy should be input

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

        allocate(p(l,m))
        
       
        ! Initialize

        u0 = u

        call calc_asthenosphere_viscous_params(kappa,kappa_p,kappa_q,beta,mu,D_lith,rho_a,g,dx)  ! recheck this belongs out of here (once)

!mmr recheck - try with finite differences too
        
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
 
        ! Finite differences

        finite_difs = .false.

        if (finite_difs) then

            
        mm = mxx
        call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        call calc_horizontal_derivatives_2D(mm_xx,mm_xy,mm_x,dx,dy)

        mm_xy = 0.
        mm = mxy
        call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        call calc_horizontal_derivatives_2D(mm_xx,mm_xy,mm_x,dx,dy)

        mm = myy
        call calc_horizontal_derivatives_2D(mm_x,mm_y,mm,dx,dy)
        call calc_horizontal_derivatives_2D(mm_yx,mm_yy,mm_y,dx,dy)

        else 

           ! mmr Calculate second derivatives with FFT
           
        call calc_fft_forward_r2r(mxx,mxx_hat)
        call calc_fft_forward_r2r(myy,myy_hat)
        call calc_fft_forward_r2r(mxy,mxy_hat)

        mm_xx_hat = kappa_p**2 * mu**2 * mxx_hat
        mm_yy_hat = kappa_q**2 * mu**2 * myy_hat
        mm_xy_hat = kappa * mu**2 * mxy_hat
        
        call calc_fft_backward_r2r(mm_xx_hat,mm_xx)
        call calc_fft_backward_r2r(mm_yy_hat,mm_yy)
        call calc_fft_backward_r2r(mm_xy_hat,mm_xy)
                         
     endif
  
     ! Viscous + elastic

     p =  - q - rho_a*g*u -rho_l*g*u_el

! mmr spare
     ! viscous only
! test1
!     f = - q - rho_a*g*u + mm_xx + 2.0*mm_xy + mm_yy
! viscous + elastic     
     ! test2, test3
     f = - q - rho_a*g*u -rho_l*g*u_el + mm_xx + 2.0*mm_xy + mm_yy

        
!     print*,'hola pressure', p(61,61)/9.81
!     stop
     
        kappa(1,1) = (kappa(1,2) + kappa(2,1)) / 2. ! 1.         

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

        deallocate(mm)
        deallocate(mm_x)
        deallocate(mm_y)

        deallocate(mm_xx)
        deallocate(mm_xy)
        deallocate(mm_yy)

        deallocate(p)

        
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

        
      subroutine calc_effective_viscosity(eta_eff,visc_c,thck_c,He_lith,n_lev,nu,dx,dy)

        implicit none

        real(wp), intent(INOUT)  :: eta_eff(:,:)
        real(wp), intent(IN)     :: visc_c
        real(wp), intent(IN)     :: thck_c
        real(wp), intent(IN)     :: He_lith(:,:) 
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
        
      end subroutine calc_effective_viscosity







      
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
      
 end module solver_lv_elva


    