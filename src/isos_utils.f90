module isos_utils

    use, intrinsic :: iso_c_binding
    use isostasy_defs, only : wp, pi, isos_domain_class, isos_param_class, &
        isos_state_class, isos_output_class

    implicit none
    include 'fftw3.f03'

    private
    
    public :: midindex

    public :: calc_density_correction_factor
    public :: calc_flexural_lengthscale
    public :: calc_homogeneous_rigidity
    public :: calc_heterogeneous_rigidity

    public :: apply_zerobc_at_corners

    public :: calc_fft_backward_r2r
    public :: calc_fft_forward_r2r
    public :: calc_fft_backward_c2r
    public :: calc_fft_forward_r2c

    public :: copy_state
    public :: cropdomain2output
    public :: cropstate2output
    public :: calc_cropindices
    public :: extendice2isostasy

    public :: interp_0d_over_time
    public :: interp_2d_over_time

    public :: maskfield
    public :: isos_set_field
    public :: isos_set_smoothed_field
    public :: smooth_gauss_2D
    public :: gauss_values
    public :: calc_gaussian_rigidity
    public :: calc_gaussian_viscosity
    
    contains

    ! ===== COMPUTE SECONDARY PHYSICAL CONSTANTS FROM PRIMARY ONES ======

    ! See Goelzer et al. (2020)
    subroutine calc_density_correction_factor(par)
        implicit none
        type(isos_param_class), intent(INOUT)   :: par

        par%Vden_factor = par%rho_ice / par%rho_water - par%rho_ice / par%rho_seawater
        return
    end subroutine calc_density_correction_factor

    ! Calculate the flexural length scale (Coulon et al, 2021, Eq. in text after Eq. 3)
    ! Note: will be on the order of 100km
    subroutine calc_flexural_lengthscale(L_w, D_lith_const, rho_uppermantle, g)
        implicit none
        real(wp), intent(INOUT) ::L_w
        real(wp), intent(IN)    :: D_lith_const, rho_uppermantle, g

        L_w = (D_lith_const / (rho_uppermantle*g))**0.25
        return
    end subroutine calc_flexural_lengthscale


    ! Calculate flexural rigidity based on effective elastic thickness of the lithosphere
    ! (He_lith), Young's modulus and Poisson's ratio. See Coulon et al. (2021) Eq. 5 & 6.
    subroutine calc_homogeneous_rigidity(D_lith, E, He_lith, nu)
        implicit none
        real(wp), intent(OUT) :: D_lith
        real(wp), intent(IN)  :: E
        real(wp), intent(IN)  :: He_lith
        real(wp), intent(IN)  :: nu

        D_lith = (E*1e9) * (He_lith*1e3)**3 / (12.0 * (1.0-nu**2))
        return
    end subroutine calc_homogeneous_rigidity

    subroutine calc_heterogeneous_rigidity(D_lith, E, He_lith, nu)
        implicit none
        real(wp), intent(OUT)   :: D_lith(:, :)
        real(wp), intent(IN)    :: E
        real(wp), intent(IN)    :: He_lith(:, :)
        real(wp), intent(IN)    :: nu

        D_lith = (E*1e9) * (He_lith*1e3)**3 / (12.0*(1.0-nu**2))

        return
    end subroutine calc_heterogeneous_rigidity

    ! ===== BC FUNCTIONS ==============================

    subroutine apply_zerobc_at_corners(x, nx, ny)
        real(wp), intent(INOUT) :: x(:, :)
        integer, intent(IN)     :: nx, ny

        x = x - 0.25 * (x(1,1) + x(nx,ny) + x(1,ny) + x(nx,1))
        return
    end subroutine apply_zerobc_at_corners


    ! ===== FFT Functions ==============================

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
        real(wp), intent(INOUT)     :: in(:, :)
        real(wp), intent(INOUT)     :: out(:, :)

        call fftw_execute_r2r(plan, in, out)

        return
      
    end subroutine calc_fft_forward_r2r

    subroutine calc_fft_backward_r2r(plan, in, out)

        implicit none 

        type(c_ptr), intent(IN)     :: plan
        real(wp), intent(INOUT)     :: in(:, :)
        real(wp), intent(INOUT)     :: out(:, :)

        integer(kind=4)            :: nx, ny
        nx = size(in,1)
        ny = size(in,2)

        call fftw_execute_r2r(plan, in, out)
        out = out / (nx * ny * 1.)

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
        real(wp), intent(INOUT)     :: in(:, :)
        complex(wp), intent(INOUT)  :: out(:, :)

        call fftw_execute_dft_r2c(plan, in, out)

        return
    end subroutine calc_fft_forward_r2c

    subroutine calc_fft_backward_c2r(plan, in, out)
        implicit none

        type(c_ptr), intent(IN)    :: plan
        complex(wp), intent(INOUT) :: in(:, :)
        real(wp), intent(INOUT)    :: out(:, :)

        integer(kind=4)            :: nx, ny
        nx = size(in,1)
        ny = size(in,2)

        call fftw_execute_dft_c2r(plan, in, out)
        out = out / (nx * ny * 1._wp)

        return
    end subroutine calc_fft_backward_c2r

    ! ===== MISC ==============================

    subroutine copy_state(ref, now)
        implicit none
        type(isos_state_class), intent(INOUT)   :: ref
        type(isos_state_class), intent(IN)      :: now

        ref%z_bed           = now%z_bed
        ref%dwdt           = now%dwdt
        ref%q               = now%q
        ref%w               = now%w
        ref%w_equilibrium   = now%w_equilibrium
        ref%we              = now%we
        ref%cplx_out_aux    = now%cplx_out_aux

        ref%Haf             = now%Haf
        ref%Hice            = now%Hice
        ref%Hseawater       = now%Hseawater

        ref%ssh             = now%ssh
        ref%ssh_perturb     = now%ssh_perturb
        ref%canom_load      = now%canom_load
        ref%canom_full      = now%canom_full
        ref%mass_anom       = now%mass_anom

        ref%maskocean       = now%maskocean
        ref%maskgrounded    = now%maskgrounded
        ref%maskcontinent   = now%maskcontinent
        return
    end subroutine copy_state

    ! ===== ARRAY SIZING FUNCTIONS ==============================

    subroutine cropdomain2output(output, domain)
        implicit none
        type(isos_output_class), intent(INOUT)  :: output
        type(isos_domain_class), intent(IN)     :: domain
        integer                                 :: i1, i2, j1, j2

        i1 = domain%icrop1
        i2 = domain%icrop2
        j1 = domain%jcrop1
        j2 = domain%jcrop2

        output%He_lith = domain%He_lith(i1:i2, j1:j2)
        output%D_lith = domain%D_lith(i1:i2, j1:j2)
        output%eta_eff = domain%eta_eff(i1:i2, j1:j2)
        output%tau = domain%tau(i1:i2, j1:j2)
        output%kappa = domain%kappa(i1:i2, j1:j2)

        output%kei = domain%kei(i1:i2, j1:j2)
        output%GE = domain%GE(i1:i2, j1:j2)
        output%GV = domain%GV(i1:i2, j1:j2)
        output%GN = domain%GN(i1:i2, j1:j2)

        return
    end subroutine cropdomain2output

    subroutine cropstate2output(output, now, i1, i2, j1, j2)
        implicit none
        type(isos_output_class), intent(INOUT)  :: output
        type(isos_state_class), intent(IN)      :: now
        integer, intent(IN)                     :: i1, i2, j1, j2

        output%Hice = now%Hice(i1:i2, j1:j2)
        output%ssh = now%ssh(i1:i2, j1:j2)
        output%z_bed = now%z_bed(i1:i2, j1:j2)
        output%dwdt = now%dwdt(i1:i2, j1:j2)
        output%w = now%w(i1:i2, j1:j2)
        output%we = now%we(i1:i2, j1:j2)
        output%ssh_perturb = now%ssh_perturb(i1:i2, j1:j2)
        output%canom_full = now%canom_full(i1:i2, j1:j2)

        return
    end subroutine cropstate2output

    subroutine calc_cropindices(icrop1, icrop2, jcrop1, jcrop2, nx, ny)
        implicit none
        integer, intent(INOUT)  :: icrop1, icrop2, jcrop1, jcrop2
        integer, intent(IN)     :: nx, ny
        integer                 :: pad

        if (nx .eq. ny) then
            icrop1 = 1
            icrop2 = nx
            jcrop1 = 1
            jcrop2 = ny
        else if (nx < ny) then
            if ( mod(ny - nx, 2) .eq. 0) then
                pad = (ny - nx) / 2
            else
                pad = (ny - nx + 1) / 2
            end if

            icrop1 = pad
            icrop2 = ny-pad
            jcrop1 = 1
            jcrop2 = ny
        else
            if ( mod(nx - ny, 2) .eq. 0) then
                pad = (nx - ny) / 2
            else
                pad = (nx - ny + 1) / 2
            end if

            jcrop1 = pad
            icrop1 = 1
            icrop2 = nx
            jcrop2 = nx-pad
        end if

        return
    end subroutine calc_cropindices

    subroutine extendice2isostasy(now, z_bed, H_ice, ssh, i1, i2, j1, j2)
        implicit none
        type(isos_state_class), intent(INOUT)   :: now
        real(wp), intent(IN)                    :: z_bed(:, :)
        real(wp), intent(IN)                    :: H_ice(:, :)
        real(wp), intent(IN)                    :: ssh(:, :)
        integer, intent(IN)                     :: i1, i2, j1, j2

        now%z_bed = 0.0
        now%z_bed(i1:i2, j1:j2) = z_bed
        now%Hice = 0.0
        now%Hice(i1:i2, j1:j2) = H_ice
        now%ssh = 0.0
        now%ssh(i1:i2, j1:j2) = ssh
        return
    end subroutine extendice2isostasy


    ! ===== INTERPOLATION FUNCTIONS ==============================

    ! Linear interpolation of scalar value over time
    subroutine interp_0d_over_time(x, y, xout, yout)

        implicit none

        real(wp), dimension(:), intent(IN) :: x 
        real(wp), dimension(:), intent(IN) :: y
        real(wp), intent(IN) :: xout
        real(wp), intent(OUT) :: yout
        integer  :: j, n
        real(wp) :: alpha

        n    = size(x)

        if (xout .lt. x(1)) then
            yout = y(1)
        else if (xout .gt. x(n)) then
            yout = y(n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do

            if (j .eq. 1) then
                yout = y(1)
            else if (j .eq. n+1) then
                yout = y(n)
            else
                alpha = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(j-1) + alpha*(y(j) - y(j-1))
             end if
        end if

        return
    end subroutine interp_0d_over_time

    ! Simple linear interpolation of 2D field over time
    subroutine interp_2d_over_time(x, y, xout, yout)

        implicit none

        real(wp), dimension(:), intent(IN) :: x 
        real(wp), dimension(:, : ,:), intent(IN) :: y
        real(wp), intent(IN) :: xout
        real(wp), dimension(:, :),intent(OUT) :: yout
        integer  :: j, n
        real(wp) :: alpha

        n    = size(x)

        if (xout .lt. x(1)) then
            yout = y(:, :,1)
        else if (xout .gt. x(n)) then
            yout = y(:, :,n)
        else
            do j = 1, n
                if (x(j) .ge. xout) exit
            end do

            if (j .eq. 1) then
                yout = y(:, :,1)
            else if (j .eq. n+1) then
                yout = y(:, :,n)
            else
                alpha = (xout - x(j-1)) / (x(j) - x(j-1))
                yout = y(:, :,j-1) + alpha*(y(:, :,j) - y(:, :,j-1))
             end if
        end if

        return
    end subroutine interp_2d_over_time

    ! ===== MASKING FUNCTIONS ==============================

    ! Mask the input field
    subroutine maskfield(out, in, mask, nx, ny)
        implicit none 

        real(wp), intent(INOUT) :: out(:, :) 
        real(wp), intent(IN)    :: in(:, :)
        logical, intent(IN)     :: mask(:, :)
        integer, intent(IN)     :: nx, ny

        integer     :: i, j

        do i = 1, nx
            do j = 1, ny
                if (mask(i,j)) then
                    out(i,j) = in(i,j)
                else
                    out(i,j) = 0.0_wp
                endif
            enddo
        enddo

        return
    end subroutine maskfield

    subroutine isos_set_field(var, var_values, mask_values, mask)
        ! Impose multiple var values according to the corresponding
        ! locations provided in a mask. 
        ! Additionally impose Gaussian smoothing via the
        ! smoothing radius sigma.

        implicit none 

        real(wp), intent(OUT) :: var(:, :) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:, :) 

        ! Local variables
        integer :: j, n 

        ! Determine how many values should be assigned
        n = size(var_values,1)

        ! Initially set var=0 everywhere
        var = 0.0_wp 

        ! Loop over unique var values and assign them
        ! to correct regions of domain.
        do j = 1, n 

            where(mask .eq. mask_values(j)) var = var_values(j)

        end do
        
        return
    end subroutine isos_set_field

    subroutine isos_set_smoothed_field(var, var_values, mask_values, mask, dx, sigma)
        ! Impose multiple var values according to the corresponding
        ! locations provided in a mask. 
        ! Additionally impose Gaussian smoothing via the
        ! smoothing radius sigma.

        implicit none 

        real(wp), intent(OUT) :: var(:, :) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:, :) 
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: sigma

        ! Safety check 
        if (sigma .le. dx) then 
            write(*,*) "isos_set_field:: Error: sigma must be larger than dx."
            write(*,*) "dx    = ", dx 
            write(*,*) "sigma = ", sigma 
            stop 
        end if 

        call isos_set_field(var, var_values, mask_values, mask)
        
        ! Apply Gaussian smoothing as desired
        call smooth_gauss_2D(var, dx,sigma)
        
        return
    end subroutine isos_set_smoothed_field

    subroutine smooth_gauss_2D(var, dx, sigma, mask_apply, mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:, :)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: sigma  
        logical,    intent(IN), optional :: mask_apply(:, :) 
        logical,    intent(IN), optional :: mask_use(:, :) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2, k 
        real(wp), allocatable :: filter0(:, :), filter(:, :) 
        real(wp), allocatable :: var_old(:, :) 
        logical,  allocatable :: mask_apply_local(:, :) 
        logical,  allocatable :: mask_use_local(:, :) 

        real(wp), allocatable :: var_ext(:, :), var_ref_ext(:, :) 

        nx    = size(var,1)
        ny    = size(var,2)

        ! Determine half-width of filter as 2-sigma
        n2 = ceiling(2.0_wp*sigma / dx)

        ! Get total number of points for filter window in each direction
        n = 2*n2+1
        
        allocate(var_old(nx,ny))
        allocate(mask_apply_local(nx,ny))
        allocate(mask_use_local(nx,ny))
        allocate(filter0(n,n))
        allocate(filter(n,n))

        allocate(var_ext(-n2:nx+n2,-n2:ny+n2))
        allocate(var_ref_ext(-n2:nx+n2,-n2:ny+n2))
        
        ! Check whether mask_apply is available 
        if (present(mask_apply)) then 
            ! use mask_use to define neighborhood points
            
            mask_apply_local = mask_apply 

        else
            ! Assume that everywhere should be smoothed

            mask_apply_local = .TRUE.
        
        end if

        ! Check whether mask_use is available 
        if (present(mask_use)) then 
            ! use mask_use to define neighborhood points
            
            mask_use_local = mask_use 

        else
            ! Assume that mask_apply also gives the points to use for smoothing 

            mask_use_local = mask_apply_local
        
        end if

        ! Calculate default 2D Gaussian smoothing kernel
        filter0 = gauss_values(dx,dx,sigma=sigma,n=n)

        var_old = var 

        var_ref_ext = -9999.0 

        var_ref_ext(1:nx,1:ny) = var 
        do i = 0, -n2, -1
            k = -i+1
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 
        do i = nx+1, nx+n2 
            k = nx + ((nx+1)-i)
            var_ref_ext(i,:) = var_ref_ext(k,:) 
        end do 

        do j = 0, -n2, -1
            k = -j+1
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 
        do j = ny+1, ny+n2 
            k = ny + ((ny+1)-j)
            var_ref_ext(:,j) = var_ref_ext(:,k) 
        end do 

        if (count(var_ref_ext .eq. -9999.0) .gt. 0) then 
            write(*,*) "Missing points!"
            stop 
        end if 

        do j = 1, ny
        do i = 1, nx


            !if (mask_apply_local(i,j)) then 
                ! Apply smoothing to this point 

                ! Limit filter input to neighbors of interest
                filter = filter0 
                !where(.not. mask_use_local(i-n2:i+n2,j-n2:j+n2) ) filter = 0.0

                ! If neighbors are available, normalize and perform smoothing  
                if (sum(filter) .gt. 0.0) then 
                    filter = filter/sum(filter)
                    var_ext(i,j) = sum(var_ref_ext(i-n2:i+n2,j-n2:j+n2)*filter) 
                end if  

            !end if 

        end do 
        end do 

        ! Get variable on normal grid 
        var = var_ext(1:nx,1:ny)

        return
    end subroutine smooth_gauss_2D

    function gauss_values(dx,dy,sigma,n) result(filt)
        ! Calculate 2D Gaussian smoothing kernel
        ! https://en.wikipedia.org/wiki/Gaussian_blur

        implicit none 

        real(wp), intent(IN) :: dx 
        real(wp), intent(IN) :: dy 
        real(wp), intent(IN) :: sigma 
        integer,  intent(IN) :: n 
        real(wp) :: filt(n,n) 

        ! Local variables 
        real(wp) :: x, y  
        integer  :: n2, i, j, i1, j1  

        if (mod(n,2) .ne. 1) then 
            write(*,*) "gauss_values:: error: n can only be odd."
            write(*,*) "n = ", n 
        end if 

        n2 = (n-1)/2 

        do j = -n2, n2 
        do i = -n2, n2 
            x = i*dx 
            y = j*dy 

            i1 = i+1+n2 
            j1 = j+1+n2 
            filt(i1,j1) = 1.0/(2.0*pi*sigma**2)*exp(-(x**2+y**2)/(2*sigma**2))

        end do 
        end do 
        
        ! Normalize to ensure sum to 1
        filt = filt / sum(filt)

        return
    end function gauss_values


    subroutine calc_gaussian_viscosity(eta,eta_0,sign,dx,dy)

        real(wp), intent(OUT) :: eta(:, :)
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
        real(wp), intent(OUT) :: He_lith(:, :)
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
                He_lith(i,j) = He_lith_0 + sign*He_lith_1 * &
                    exp(-0.5*dot_product([xc(i),yc(j)]-[xcntr,ycntr], &
                    matmul(He_lith_sigma_m1, [xc(i),yc(j)]-[xcntr,ycntr]) ))
            enddo
        enddo
        
        return
    end subroutine calc_gaussian_rigidity

    function midindex(n) result(n2)
        implicit none
        integer, intent(IN) :: n
        integer             :: n2

        if (mod(n,2) .eq. 0) then
            n2 = n/2
        else
            n2 = (n-1) / 2
        end if

    end function midindex

end module isos_utils