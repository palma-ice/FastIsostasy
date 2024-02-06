
module convolutions

    use isostasy_defs, only : wp, pi, isos_domain_class
    use isos_utils

    use, intrinsic :: iso_c_binding
    implicit none
    include 'fftw3.f03'

    private
    
    public :: convenient_calc_convolution_indices
    public :: calc_convolution_indices
    public :: convenient_samesize_fftconvolution
    public :: samesize_fftconvolution
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

    subroutine convenient_samesize_fftconvolution(out, in1, in2, domain)

        real(wp), intent(OUT)   :: out(:, :)
        real(wp), intent(IN)    :: in1(:, :)
        real(wp), intent(IN)    :: in2(:, :)
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

        real(wp), intent(OUT)   :: out(:, :)
        real(wp), intent(IN)    :: in1(:, :)
        real(wp), intent(IN)    :: in2(:, :)
        integer, intent(IN)     :: i1, i2, j1, j2, nx, ny

        ! Local variables
        real(wp), allocatable       :: out_ext(:, :)
        real(wp), allocatable       :: in1_ext(:, :)
        real(wp), allocatable       :: in2_ext(:, :)
        complex(wp), allocatable    :: out_ext_hat(:, :)
        complex(wp), allocatable    :: in1_ext_hat(:, :)
        complex(wp), allocatable    :: in2_ext_hat(:, :)
        
        type(c_ptr), intent(IN)     :: forward_plan
        type(c_ptr), intent(IN)     :: backward_plan

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
        in1_ext(1:nx, 1:ny) = in1
        in2_ext(1:nx, 1:ny) = in2

        !# TODO: we can precompute this
        call calc_fft_forward_r2c(forward_plan, in1_ext, in1_ext_hat)
        call calc_fft_forward_r2c(forward_plan, in2_ext, in2_ext_hat)

        ! Multiply both
        out_ext_hat =  in1_ext_hat * in2_ext_hat

        ! Invert product
        call calc_fft_backward_c2r(backward_plan, out_ext_hat, out_ext)
        out(1:nx,1:ny) = out_ext(i1:i2, j1:j2)
        ! (i1+offset:i2+offset, j1-offset:j2-offset)

        return
    end subroutine samesize_fftconvolution

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
        ! write(*,*) "helper real 1: ", sum(helper_real),  "helper cplx 1: ", sum(helper_cplx)

        ! Multiply both
        helper_cplx =  kernel * helper_cplx

        ! Invert product
        call calc_fft_backward_c2r(backward_plan, helper_cplx, helper_real)
        ! write(*,*) "helper real 2: ", sum(helper_real),  "helper cplx 2: ", sum(helper_cplx)

        out(1:nx,1:ny) = helper_real(i1+offset:i2+offset, j1-offset:j2-offset)
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