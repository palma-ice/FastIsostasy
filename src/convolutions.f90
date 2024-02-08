
module convolutions

    use isostasy_defs, only : wp, pi, isos_domain_class
    use isos_utils

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    private
    
    public :: convenient_calc_convolution_indices
    public :: calc_convolution_indices
    public :: precomputed_fftconvolution
    public :: precompute_kernel

    contains

    subroutine convenient_calc_convolution_indices(domain)
        implicit none
        type(isos_domain_class), intent(INOUT)  :: domain

        call calc_convolution_indices(domain%i1, domain%i2, domain%j1, domain%j2, &
            domain%offset, domain%ny, domain%nx)
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

    ! Compute `out` as the (fft-based) convolution of a precomputed `kernel` and
    ! input field `in`.
    subroutine precomputed_fftconvolution(out, kernel, in, i1, i2, j1, j2, &
        offset, nx, ny, forward_plan, backward_plan)

        implicit none

        real(wp),    intent(OUT)    :: out(:, :)
        complex(wp), intent(IN)     :: kernel(:, :)
        real(wp), intent(IN)        :: in(:, :)
        integer, intent(IN)         :: i1, i2, j1, j2, nx, ny, offset
        type(c_ptr), intent(IN)     :: forward_plan
        type(c_ptr), intent(IN)     :: backward_plan

        ! Local variables
        real(wp),    allocatable    :: helper_real(:, :)
        complex(wp), allocatable    :: helper_cplx(:, :)

        ! Populate load on extended grid
        allocate(helper_real(2*nx-1, 2*ny-1))
        allocate(helper_cplx(2*nx-1, 2*ny-1))

        ! Zero-padded FFT
        helper_real = 0.
        helper_real(1:nx, 1:ny) = in
        call calc_fft_forward_r2c(forward_plan, helper_real, helper_cplx)

        ! Compute and invert product
        helper_cplx =  kernel * helper_cplx
        call calc_fft_backward_c2r(backward_plan, helper_cplx, helper_real)

        call apply_zerobc_at_corners(helper_real, 2*nx-1, 2*ny-1)
        out(1:nx, 1:ny) = helper_real(i1+offset:i2+offset, j1-offset:j2-offset)
        return
    end subroutine precomputed_fftconvolution

    subroutine precompute_kernel(plan, kernel, fftkernel, nx, ny)
        implicit none

        type(c_ptr), intent(IN)     :: plan
        real(wp),    intent(IN)     :: kernel(:, :)
        complex(wp), intent(INOUT)  :: fftkernel(:, :)
        real(wp), allocatable       :: extended_kernel(:, :)
        integer, intent(INOUT)      :: nx, ny

        allocate(extended_kernel(2*nx-1, 2*ny-1))
        extended_kernel = 0.0_wp
        extended_kernel(1:nx, 1:ny) = kernel
        call calc_fft_forward_r2c(plan, extended_kernel, fftkernel)

        return
    end subroutine precompute_kernel

end module convolutions