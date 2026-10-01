
module convolutions

    use isostasy_defs, only : sp, dp, wp, pi, isos_domain_class
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
            domain%offset, domain%nx, domain%ny)
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
    ! input field `in`, using the FFTW work arrays of `domain`.
    subroutine precomputed_fftconvolution(out, kernel, in, domain)

        implicit none

        real(wp),    intent(OUT)            :: out(:, :)
        complex(dp), intent(IN)             :: kernel(:, :)
        real(wp),    intent(IN)             :: in(:, :)
        type(isos_domain_class), intent(IN) :: domain

        associate(nx => domain%nx, ny => domain%ny, offset => domain%offset, &
            i1 => domain%i1, i2 => domain%i2, j1 => domain%j1, j2 => domain%j2)

        ! Zero-padded FFT
        domain%fft_r = 0.0_dp
        domain%fft_r(1:nx, 1:ny) = in
        call calc_fft_forward_r2c(domain%forward_dftplan_r2c, domain%fft_r, domain%fft_c)

        ! Compute and invert product
        domain%fft_c = kernel * domain%fft_c
        call calc_fft_backward_c2r(domain%backward_dftplan_c2r, domain%fft_c, domain%fft_r)

        call apply_zerobc_at_corners_dp(domain%fft_r, 2*nx-1, 2*ny-1)
        out(1:nx, 1:ny) = domain%fft_r(i1+offset:i2+offset, j1-offset:j2-offset)

        end associate
        return
    end subroutine precomputed_fftconvolution

    ! FFT of the zero-padded `kernel`, as needed by precomputed_fftconvolution.
    function precompute_kernel(kernel, domain) result(fftkernel)
        implicit none

        real(wp), intent(IN)                :: kernel(:, :)
        type(isos_domain_class), intent(IN) :: domain
        complex(dp)                         :: fftkernel(size(domain%fft_c, 1), &
                                                         size(domain%fft_c, 2))

        domain%fft_r = 0.0_dp
        domain%fft_r(1:domain%nx, 1:domain%ny) = kernel
        call calc_fft_forward_r2c(domain%forward_dftplan_r2c, domain%fft_r, domain%fft_c)
        fftkernel = domain%fft_c

        return
    end function precompute_kernel

end module convolutions