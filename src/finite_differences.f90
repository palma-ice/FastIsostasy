module finite_differences

    use isostasy_defs, only : wp, pi 

    implicit none

    private
    
    public :: calc_derivative_x
    public :: calc_derivative_y
    public :: calc_derivative_xx
    public :: calc_derivative_yy
    
    contains

    ! Get simple horizontal derivatives wrt x
    subroutine calc_derivative_x(dudx, u, dx, nx, ny)
        implicit none
        real(wp), intent(OUT)   :: dudx(:, :) 
        real(wp), intent(IN)    :: u(:, :) 
        real(wp), intent(IN)    :: dx(:, :)
        integer,  intent(IN)    :: nx, ny
        
        ! Local variables 
        integer :: i, j, im1, ip1

        dudx = 0.0

        do j = 1, ny 
           do i = 2, nx-1
              im1 = i-1
              ip1 = i+1
              dudx(i, j) = (u(ip1,j)-u(im1,j))/(2*dx(i, j))
           end do
        end do

        dudx(nx,:) = (u(nx,:)-u(nx-1,:))/dx(nx,:)
        dudx(1,:) = (u(2,:)-u(1,:))/dx(1,:)

        return
    end subroutine calc_derivative_x

    subroutine calc_derivative_y(dudy,u,dy,nx,ny)  
        implicit none
        real(wp), intent(OUT)   :: dudy(:, :)
        real(wp), intent(IN)    :: u(:, :)
        real(wp), intent(IN)    :: dy(:, :)
        integer,  intent(IN)    :: nx, ny

        ! Local variables 
        integer :: i, j, jm1, jp1

        dudy = 0.0

        do j = 2, ny-1  
           do i = 1, nx
              jm1 = j-1
              jp1 = j+1
              dudy(i,j) = (u(i,jp1)-u(i,jm1))/(2*dy(i, j))
           end do
        end do

        dudy(:,ny) = (u(:,ny)-u(:,ny-1))/dy(:,ny)
        dudy(:,1) = (u(:,2)-u(:,1))/dy(:,1)

        return
    end subroutine calc_derivative_y

    ! Get simple horizontal derivatives
    subroutine calc_derivative_xx(uxx,u,dx,nx,ny)
        implicit none
        real(wp), intent(OUT)   :: uxx(:, :)
        real(wp), intent(IN)    :: u(:, :)
        real(wp), intent(IN)    :: dx(:, :)
        integer,  intent(IN)    :: nx, ny

        ! Local variables 
        integer :: i, j, im1, ip1

        uxx = 0.0

        do j = 1, ny  
           do i = 2, nx-1
              im1 = i-1
              ip1 = i+1
              uxx(i,j) = (u(ip1,j)- 2*u(i,j) + u(im1,j))/(dx(i, j)*dx(i, j))
           end do
        end do

        uxx(1,:) = (u(3,:)-2*u(2,:)+u(1,:))/(dx(1,:)*dx(1,:))
        uxx(nx,:) = (u(nx,:) - 2*u(nx-1,:) + u(nx-2,:))/(dx(1,:)*dx(1,:))

        return
    end subroutine calc_derivative_xx

    subroutine calc_derivative_yy(uyy,u,dy,nx,ny)
        implicit none
        real(wp), intent(OUT)   :: uyy(:, :)
        real(wp), intent(IN)    :: u(:, :)
        real(wp), intent(IN)    :: dy(:, :)
        integer,  intent(IN)    :: nx, ny

        ! Local variables 
        integer :: i, j, jm1, jp1

        uyy = 0.0
            
        do j = 2, ny-1  
           do i = 1, nx
              jm1 = j-1
              jp1 = j+1
              uyy(i,j) = (u(i,jp1) - 2*u(i,j) + u(i,jm1))/(dy(i, j)*dy(i, j))
           end do
        end do

        uyy(:,1) = (u(:,3)-2*u(:,2)+u(:,1))/(dy(:,1)*dy(:,1))
        uyy(:,ny) = (u(:,ny) - 2*u(:,ny-1) + u(:,ny-2))/(dy(:,ny)*dy(:,ny))

        return
    end subroutine calc_derivative_yy

end module finite_differences