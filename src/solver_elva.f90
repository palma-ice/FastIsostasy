module solver_elva

    use isostasy_defs, only : sp, dp, wp, pi 

    implicit none

    private
    public :: calc_asthenosphere_viscous_square
    public :: calc_asthenosphere_viscous

contains

    subroutine calc_asthenosphere_viscous_square(dzbdt,w,q,D_lith,eta,rho_a,g,dx,dt)
    
        ! Extend a given domain [nx,ny] so that it is square based on the 
        ! largest dimension of original data. 
        ! Solve calc_asthenosphere_viscous on the square arrays, then 
        ! extract solution onto original domain [nx,ny]. 

        implicit none

        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: w(:,:)    
        real(wp), intent(IN)    :: q(:,:)    
        real(wp), intent(IN)    :: D_lith
        real(wp), intent(IN)    :: eta                  ! [Pa s] Viscosity, eta=1e21 by default.
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
        call extend_array(sq_w,w,fill_with="val",val=0.0_wp)
        call extend_array(sq_q,q,fill_with="val",val=0.0_wp)

        ! Step 2: solve

        call calc_asthenosphere_viscous(sq_dzbdt,sq_w,sq_q,D_lith,eta,rho_a,g,dx,dt)

        ! Step 3: get solution on original grid

        call reduce_array(dzbdt,sq_dzbdt)
        call reduce_array(w,sq_w)
        
        return

    end subroutine calc_asthenosphere_viscous_square

    subroutine calc_asthenosphere_viscous(dzbdt,w,q,D_lith,eta,rho_a,g,dx,dt)
        ! Calculate rate of change of vertical bedrock height
        ! from a viscous half-space asthenosphere.
        ! Contains viscous component only.

        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'

        real(wp), intent(OUT)   :: dzbdt(:,:)
        real(wp), intent(INOUT) :: w(:,:)
        real(wp), intent(IN)    :: q(:,:)
        real(wp), intent(IN)    :: D_lith
        real(wp), intent(IN)    :: eta                  ! [Pa s] Viscosity, eta=1e21 by default.
        real(wp), intent(IN)    :: rho_a 
        real(wp), intent(IN)    :: g 
        real(wp), intent(IN)    :: dx
        real(wp), intent(IN)    :: dt
        
        ! Local variables
        integer  :: l, m, i, j
        real(wp) :: dt_sec
        logical  :: fft_r2r, fft_c2c

        real(wp), allocatable :: kappa(:,:)
        real(wp), allocatable :: beta(:,:)
        real(wp) :: mu 

        real(wp), allocatable :: w0(:,:)
        real(wp), allocatable :: q_hat(:,:)
        real(wp), allocatable :: w_hat(:,:)

        real(c_double), allocatable :: data_in(:,:)
        real(c_double), allocatable :: data_inbis(:,:)
        complex(c_double_complex), allocatable :: data_out(:,:)

        complex(dp), allocatable :: q_hat_c(:,:)
        complex(dp), allocatable :: w_hat_c(:,:)

        real(wp), allocatable    :: w_hat_c_re(:,:)
        real(wp), allocatable    :: w_hat_c_im(:,:)      

        
        integer, parameter :: nd = 2
        
        dt_sec = dt * 3600*24*365 ! [s] 

        l = size(dzbdt,1)
        m = size(dzbdt,2)

        allocate(kappa(l,m))
        allocate(beta(l,m))

        allocate(w0(l,m))
        allocate(q_hat(l,m))
        allocate(w_hat(l,m))
        allocate(q_hat_c(l,m/2+1)) 
        allocate(w_hat_c(l,m/2+1))

        allocate(w_hat_c_re(l,m/2+1))
        allocate(w_hat_c_im(l,m/2+1))

        ! Step 0: initialize helper variables

        call calc_asthenosphere_viscous_params(kappa,beta,mu,D_lith,rho_a,g,dx) 

        !  Initialize 
        
        w0 = w

        q_hat    = 0.0
        w_hat    = 0.0
        w_hat_c  = 0.0
        dzbdt    = 0.0

        fft_r2r = .true.
        fft_c2c = .false.
      
        if (fft_r2r) then

            ! fft of load
            ! 

            call calc_fft_forward_r2r(q,q_hat)

            ! fft of displacement

            call calc_fft_forward_r2r(w,w_hat)

            ! displacement fft at timestep n+1    

            if (dt.gt.0)  w_hat = ( ( 2.*eta*mu*kappa - (dt_sec/2.)*beta)*w_hat + dt_sec*q_hat)/(2*eta*mu*kappa + (dt_sec/2.)*beta)

            ! Inverse fft to obtain displacement at n+1

            call calc_fft_backward_r2r(w_hat,w)

        else if (fft_c2c) then

            ! fft of load

            call calc_fft_forward_c2c(q,q_hat)

            ! fft of displacement

            call calc_fft_forward_c2c(w,w_hat)

            ! displacement fft at timestep n+1    

            if (dt.gt.0)  w_hat = ( ( 2.*eta*mu*kappa - (dt_sec/2.)*beta)*w_hat + dt_sec*q_hat)/(2*eta*mu*kappa + (dt_sec/2.)*beta)

            ! Inverse fft to obtain displacement at n+1

            call calc_fft_backward_c2c(w_hat,w)

        else  

            write(*,*) "Error, you have to choose an option to calculate the FFTs."
            stop
      
        end if
         
        ! Impose boundary conditions 

        w  = w - 0.25*(w(1,1)+w(l,m)+w(1,m)+w(l,1)) 

        ! Rate of viscous asthenosphere uplift

        if (dt.gt.0.)  dzbdt = -(w-w0)/dt

        return

    end subroutine calc_asthenosphere_viscous

    subroutine calc_asthenosphere_viscous_params(kappa,beta,mu,D_lith,rho_a,g,dx) 
        ! Calculate parameters needed for elastic lithosphere viscous asthenosphere (ELVA)                                        
        ! solution as in Bueler et al 2007 (eq 11)
                    
        real(wp), intent(OUT)  :: kappa(:,:)
        real(wp), intent(OUT)  :: beta(:,:) 
        real(wp), intent(OUT)  :: mu      
        real(wp), intent(IN)   :: D_lith
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
                beta(i,j)   = rho_a*g + D_lith*(mu**4)*kappa(i,j)**4
            end do
        end do

        return
      
    end subroutine calc_asthenosphere_viscous_params

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
!mmr      call fftw_destroy_plan(plan)
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
!mmr      call fftw_destroy_plan(plan)
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

end module solver_elva