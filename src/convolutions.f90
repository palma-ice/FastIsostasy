
module convolutions

    use isostasy_defs, only : wp, pi, isos_domain_class
    use isos_utils

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    private
    
    public :: convenient_calc_convolution_indices
    public :: calc_convolution_indices
    public :: convolve_load_elastic_plate
    public :: convenient_samesize_fftconvolution
    public :: samesize_fftconvolution

    contains

    subroutine convenient_calc_convolution_indices(domain)
        implicit none
        type(isos_domain_class), intent(INOUT)  :: domain

        call calc_convolution_indices(domain%i1, domain%i2, domain%j1, domain%j2, &
            domain%convo_offset, domain%ny, domain%nx)
        return
    end subroutine convenient_calc_convolution_indices

    subroutine calc_convolution_indices(i1, i2, j1, j2, convo_offset, nx, ny)
        implicit none
        integer, intent(INOUT)  :: i1, i2, j1, j2, convo_offset
        integer, intent(IN)     :: nx, ny

        if ( mod(nx,2).eq.0 ) then
            i1 = nx/2
        else
            i1 = int(nx/2) + 1
        endif
        i2 = 2*nx - 1 - int(nx/2)
 
        if ( mod(ny,2).eq.0 ) then
            j1 = ny/2
        else
            j1 = int(ny/2) + 1
        endif
        j2 = 2*ny - 1 - int(ny/2)

        if ( mod(ny - nx, 2).eq.0 ) then
            convo_offset = (ny - nx) / 2
        else
            convo_offset = (ny - nx - 1) / 2
        endif

        return
    end subroutine calc_convolution_indices

    subroutine convolve_load_elastic_plate(w1,q1,GG)
        ! Spread the load q1 [Pa] from each point in the grid
        ! via the regional Green's function scaling GG [m N-1]

        implicit none

        real(wp), intent(OUT) :: w1(:,:)        ! [m] Lithospheric displacement
        real(wp), intent(IN)  :: q1(:,:)        ! [Pa] Lithospheric loading
        real(wp), intent(IN)  :: GG(:,:)        ! Regional scaling filter

        ! Local variables
        !integer :: ip, jp, lpx, lpy
        real(wp), allocatable :: q1_ext(:,:)
        real(wp), allocatable :: w_reg(:,:)

        integer :: i, j, nx ,ny, nr

        nx = size(w1,1)
        ny = size(w1,2)

        ! Size of regional neighborhood 
        nr = (size(GG,1)-1)/2 

        ! Populate load on extended grid
        allocate(q1_ext(1-nr:nx+nr,1-nr:ny+nr))

        ! First fill in main grid points with current point load
        q1_ext(1:nx,1:ny) = q1 

        ! Populate the extended grid points
        do i = 1, nx
            q1_ext(i,1-nr:0)=q1_ext(i,1)
            q1_ext(i,ny+1:ny+nr)=q1_ext(i,ny)
        end do
        do j = 1, ny
            q1_ext(1-nr:0,j)=q1_ext(1,j)
            q1_ext(NX+1:NX+nr,j)=q1_ext(nx,j)
        end do
        
        ! Populate the extended grid corner points     
        q1_ext(1-nr:0,1-nr:0)         = q1_ext(1,1)
        q1_ext(1-nr:0,ny+1:ny+nr)     = q1_ext(1,ny)
        q1_ext(nx+1:nx+nr,1-nr:0)     = q1_ext(nx,1)
        q1_ext(nx+1:nx+nr,ny+1:ny+nr) = q1_ext(nx,ny)

        ! ----- allocation de w_reg  et de croix -----------

        allocate(w_reg(-nr:nr,-nr:nr))

        w_reg = 0.
        
        do j = 1, ny
            do i = 1, nx
                ! Apply the neighborhood scaling to the deflection in the neighborhood
                w_reg = GG * q1_ext(i-nr:i+nr,j-nr:j+nr)

                ! Sum to get total deflection at current point due to all neighbors
                w1(i,j) = sum(w_reg)
            end do
        end do

        return

    end subroutine convolve_load_elastic_plate

    subroutine convenient_samesize_fftconvolution(out, in1, in2, domain)

        real(wp), intent(OUT)   :: out(:,:)
        real(wp), intent(IN)    :: in1(:,:)
        real(wp), intent(IN)    :: in2(:,:)
        type(isos_domain_class), intent(In)  :: domain

        call samesize_fftconvolution(out, in1, in2, domain%i1, domain%i2, &
            domain%j1, domain%j2, domain%nx, domain%ny, &
            domain%forward_fftplan_r2r, domain%backward_fftplan_r2r)
        return
    end subroutine convenient_samesize_fftconvolution
    
    ! Compute "out" as the (fft-based) convolution of "in1" and "in2"
    subroutine samesize_fftconvolution(out, in1, in2, i1, i2, j1, j2, nx, ny, &
        forward_plan, backward_plan)

        implicit none

        real(wp), intent(OUT)   :: out(:,:)
        real(wp), intent(IN)    :: in1(:,:)
        real(wp), intent(IN)    :: in2(:,:)
        integer, intent(IN)     :: i1, i2, j1, j2, nx, ny

        ! Local variables
        real(wp), allocatable :: out_ext(:,:)
        real(wp), allocatable :: in1_ext(:,:)
        real(wp), allocatable :: in2_ext(:,:)
        complex(wp), allocatable :: out_ext_hat(:,:)
        complex(wp), allocatable :: in1_ext_hat(:,:)
        complex(wp), allocatable :: in2_ext_hat(:,:)
        
        type(c_ptr), intent(IN)     :: forward_plan
        type(c_ptr), intent(IN)     :: backward_plan

        real(wp), allocatable :: w_reg(:,:)

        integer :: i, j

        ! Populate load on extended grid
        allocate(    out_ext(1:2*nx-1,1:2*ny-1))
        allocate(    in1_ext(1:2*nx-1,1:2*ny-1))
        allocate(    in2_ext(1:2*nx-1,1:2*ny-1))
        allocate(out_ext_hat(1:2*nx-1,1:2*ny-1))
        allocate(in1_ext_hat(1:2*nx-1,1:2*ny-1))
        allocate(in2_ext_hat(1:2*nx-1,1:2*ny-1))

        ! Pad with zeros
        in1_ext = 0.
        in2_ext = 0.
        in1_ext(1:nx,1:ny) = in1
        in2_ext(1:nx,1:ny) = in2

        ! TODO: we can precompute this
        call calc_fft_forward_r2c(forward_plan, in1_ext, in1_ext_hat)
        call calc_fft_forward_r2c(forward_plan, in2_ext, in2_ext_hat)

        ! Multiply both
        out_ext_hat =  in1_ext_hat * in2_ext_hat

        ! Invert product
        call calc_fft_backward_c2r(backward_plan, out_ext_hat, out_ext)

        ! TODO: check if normalisation is correctly accounted for in the new version.
        ! normalisation is needed because of the FFTW definition
        out(1:nx,1:ny) = out_ext(i1:i2, j1:j2) ! * ((2*nx-1) * (2*ny-1)) ** 0.5  ! sqrt

        return
    end subroutine samesize_fftconvolution

end module convolutions