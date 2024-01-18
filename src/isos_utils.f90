module isos_utils

    use isostasy_defs, only : wp, pi

    implicit none

    private
    
    public :: extend_array
    public :: reduce_array
    public :: isos_set_field
    public :: isos_set_smoothed_field
    public :: smooth_gauss_2D
    public :: gauss_values
    public :: calc_gaussian_rigidity
    public :: calc_gaussian_viscosity
    
    contains



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

    subroutine isos_set_field(var, var_values, mask_values, mask)
        ! Impose multiple var values according to the corresponding
        ! locations provided in a mask. 
        ! Additionally impose Gaussian smoothing via the
        ! smoothing radius sigma.

        implicit none 

        real(wp), intent(OUT) :: var(:,:) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:,:) 

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

        real(wp), intent(OUT) :: var(:,:) 
        real(wp), intent(IN)  :: var_values(:)
        real(wp), intent(IN)  :: mask_values(:)
        real(wp), intent(IN)  :: mask(:,:) 
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
        call smooth_gauss_2D(var,dx=dx,sigma=sigma)
        
        return

    end subroutine isos_set_smoothed_field

    subroutine smooth_gauss_2D(var,dx,sigma,mask_apply,mask_use)
        ! Smooth out a field to avoid noise 
        ! mask_apply designates where smoothing should be applied 
        ! mask_use   designates which points can be considered in the smoothing filter 

        implicit none

        real(wp),   intent(INOUT) :: var(:,:)      ! [nx,ny] 2D variable
        real(wp),   intent(IN)    :: dx 
        real(wp),   intent(IN)    :: sigma  
        logical,    intent(IN), optional :: mask_apply(:,:) 
        logical,    intent(IN), optional :: mask_use(:,:) 

        ! Local variables
        integer  :: i, j, nx, ny, n, n2, k 
        real(wp), allocatable :: filter0(:,:), filter(:,:) 
        real(wp), allocatable :: var_old(:,:) 
        logical,  allocatable :: mask_apply_local(:,:) 
        logical,  allocatable :: mask_use_local(:,:) 

        real(wp), allocatable :: var_ext(:,:), var_ref_ext(:,:) 

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

end module isos_utils