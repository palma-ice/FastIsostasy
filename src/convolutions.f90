
module convolutions

    use isostasy_defs, only : wp, pi

    implicit none

    private
    
    public :: convolve_load_elastic_plate
    public :: convolve_load_elastic_plate_fft
    public :: calc_geoid_displacement

    contains

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

    subroutine convolve_load_elastic_plate_fft(w1, q1, GG, i1, i2, j1, j2, nx, ny)
        ! Spread the load q1 [Pa] from each point in the grid
        ! via the regional Green's function scaling GG [m N-1]

        implicit none

        real(wp), intent(OUT)   :: w1(:,:)        ! [m] Lithospheric displacement
        real(wp), intent(IN)    :: q1(:,:)        ! [Pa] Lithospheric loading
        real(wp), intent(IN)    :: GG(:,:)        ! Regional scaling filter
        integer, intent(IN)     :: i1, i2, j1, j2, nx, ny

        ! Local variables
        real(wp), allocatable :: q1_ext(:,:)
        real(wp), allocatable :: GG_ext(:,:)
        real(wp), allocatable :: w_ext(:,:)
        complex(wp), allocatable :: q1_ext_hat(:,:)
        complex(wp), allocatable :: GG_ext_hat(:,:)
        complex(wp), allocatable :: w_ext_hat(:,:)
        
        real(wp), allocatable :: w_reg(:,:)

        integer :: i, j, nr

        ! Size of regional neighborhood 
        nr = size(GG, 1) !(size(GG,1)-1)/2

        if ((nr.ne.nx).or.(nr.ne.ny)) then
           print*,'nr, nx, ny =', nr, nx, ny
           stop
        endif


        ! Populate load on extended grid
        allocate(    q1_ext(1:2*nx-1,1:2*ny-1))
        allocate(    GG_ext(1:2*nx-1,1:2*ny-1))
        allocate(GG_ext_hat(1:2*nx-1,1:2*ny-1))
        allocate(q1_ext_hat(1:2*nx-1,1:2*ny-1))
        allocate(     w_ext(1:2*nx-1,1:2*ny-1))
        allocate( w_ext_hat(1:2*nx-1,1:2*ny-1))

        ! Pad with zeros
        q1_ext = 0.
        GG_ext = 0.
        
        ! First fill in main grid points with current point load
        q1_ext(1:nx,1:ny) = q1 

        ! First fill in main grid points with current point load
        GG_ext(1:nx,1:ny) = GG

        call calc_fft_forward_r2c(GG_ext,GG_ext_hat)

        ! ! Fourier transfor q1_ext

        call calc_fft_forward_r2c(q1_ext,q1_ext_hat)

        ! Multiply both

        w_ext_hat =  q1_ext_hat * GG_ext_hat

        ! Invert product

        call calc_fft_backward_c2r(w_ext_hat,w_ext)

        ! normalisation is needed because of the FFTW definition
        w1(1:nx,1:ny) = w_ext(i1:i2, j1:j2) * ((2*nx-1) * (2*ny-1)) ** 0.5  ! sqrt

        return
    end subroutine convolve_load_elastic_plate_fft

    subroutine calc_geoid_displacement(ssh_perturb, q, gn, i1, i2, j1, j2, nx, ny)

        implicit none
      
        real(wp), intent(INOUT)  :: ssh_perturb(:,:)
        real(wp), intent(IN)     :: GN(:,:)
        real(wp), intent(IN)     :: q(:,:)
        integer, intent(IN)     :: i1, i2, j1, j2, nx, ny

    ! recheck convol     
    !      call convolve_load_elastic_plate(ssh_perturb,q,GN)
        call convolve_load_elastic_plate_fft(ssh_perturb, q, GN, i1, i2, j1, j2, nx, ny)
        
        ssh_perturb  = ssh_perturb - 0.25 * (ssh_perturb(1,1) + ssh_perturb(nx,ny) & 
            + ssh_perturb(1,ny) + ssh_perturb(nx,1))

    end subroutine calc_geoid_displacement

end module convolutions