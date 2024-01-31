module green_functions

    use isostasy_defs, only : wp, pi

    implicit none

    private
    
    public :: calc_greens_function_scaling
    public :: calc_ge_filter_2D
    public :: get_ge_value
    public :: calc_gn_value
    public :: calc_GN_filter_2D
    
    contains


    ! The Green's function (Eq. 3 of Coulon et al, 2021) gives displacement G in [m]
    ! as a function of the distance r from the point load P_b [Pa]. 

    ! Here G0 is calculated, which is G without including the point load. G0 has units
    ! of [m N-1]. This can then be multiplied with the actual magnitude of the
    ! point load to obtain G.
    ! G = G0 * P_b = [m N-1] * [Pa] = [m]. 

    ! Note that L_w contains information about rho_uppermantle. 
    subroutine calc_greens_function_scaling(G0, kei2D, L_w, D_lith, dx, dy)


        implicit none

        real(wp), intent(OUT) :: G0(:,:) 
        real(wp), intent(IN)  :: kei2D(:,:) 
        real(wp), intent(IN)  :: L_w 
        real(wp), intent(IN)  :: D_lith 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        
        G0 = -L_w**2 / (2.0*pi*D_lith) * kei2D * (dx*dy)

        return

    end subroutine calc_greens_function_scaling

    subroutine calc_GE_filter_2D(filt, dx, dy)
        ! Calculate 2D Green Function

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp) :: ge_test_0
        real(wp) :: ge_test_1

        real(wp), parameter, dimension(42) :: rn_vals = [ 0.0,    0.011,  0.111,  1.112,  2.224, &
            3.336, 4.448,  6.672, 8.896,  11.12,  17.79,  22.24,  27.80,  33.36,  44.48,  55.60, &
            66.72,  88.96,  111.2,  133.4,  177.9,  222.4,  278.0,  333.6, 444.8,  556.0,  667.2, &
            778.4,  889.6,  1001.0, 1112.0, 1334.0, 1779.0, 2224.0, 2780.0, 3336.0, 4448.0, 5560.0, &
            6672.0, 7784.0, 8896.0, 10008.0]* 1.e3        !    # converted to meters
       
        real(wp), parameter, dimension(42) :: ge_vals = [ -33.6488, -33.64, -33.56, -32.75, -31.86, &
            -30.98, -30.12, -28.44, -26.87, -25.41,-21.80, -20.02, -18.36, -17.18, -15.71, -14.91, &
            -14.41, -13.69, -13.01,-12.31, -10.95, -9.757, -8.519, -7.533, -6.131, -5.237, -4.660, &
            -4.272,-3.999, -3.798, -3.640, -3.392, -2.999, -2.619, -2.103, -1.530, -0.292, 0.848, &
            1.676,  2.083,  2.057,  1.643]


        ! Get size of filter array and half-width
        n  = size(filt, 1)
        write(*,*) "n = ", n
        n2 = (n-1) / 2 
       
        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_ge_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_ge_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop
        end if 


        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center
        ! and impose correct Green function value. 

        filt = 0.
        
        do j = -n2, n2
        do i = -n2, n2

            x  = i*dx
            y  = j*dy
            r  = sqrt(x**2+y**2)

            ! Get actual index of array
            i1 = i+1+n2
            j1 = j+1+n2

            ! Get correct GE value for this point
            filt(i1,j1) = get_ge_value(r, rn_vals, ge_vals) * 1.e-12 *(dx*dy) /(9.81*max(r,dx))
            ! TODO: recheck 9.81

        end do
        end do

        return 
      end subroutine calc_GE_filter_2D

      function get_ge_value(r,rn_vals,ge_vals) result(ge)

        implicit none

        real(wp), intent(IN) :: r           ! [m] Radius from point load 
        real(wp), intent(IN) :: rn_vals(:)  ! [m] Tabulated  radius values
        real(wp), intent(IN) :: ge_vals(:) ! [-] Tabulated Green function values
        real(wp) :: ge

        ! Local variables 
        integer :: k, n 
        real(wp) :: rn_now 
        real(wp) :: wt 

        n = size(rn_vals,1) 

        ! Get radius from point load
        rn_now = r   

        if (rn_now .gt. rn_vals(n)) then

            ge = 0. !ge_vals(n)

        else 

            do k = 1, n-1
                if (rn_now .ge. rn_vals(k) .and. rn_now .lt. rn_vals(k+1)) exit
            end do 

            ! Linear interpolation to get current ge value
            ge = ge_vals(k) &
                + (rn_now-rn_vals(k))/(rn_vals(k+1)-rn_vals(k))*(ge_vals(k+1)-ge_vals(k))   

        end if 

        return

      end function get_GE_value

    ! Compute Green function (Coulon et al. 2021) to determine the perturbation
    ! of the sea-surface height
    subroutine calc_GN_filter_2D(filt,m_earth,r_earth,dx,dy) 
        ! Calculate 2D Green Function

        implicit none 

        real(wp), intent(OUT) :: filt(:,:) 
        real(wp), intent(IN)  :: m_earth
        real(wp), intent(IN)  :: r_earth
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy   

        ! Local variables 
        integer  :: i, j, i1, j1, n, n2
        real(wp) :: x, y, r

        real(wp) :: ge_test_0
        real(wp) :: ge_test_1


        ! Get size of filter array and half-width
        n  = size(filt,1) 
        n2 = (n-1)/2 
        
        ! Safety check
        if (size(filt,1) .ne. size(filt,2)) then 
            write(*,*) "calc_ge_filt:: error: array 'filt' must be square [n,n]."
            write(*,*) "size(filt): ", size(filt,1), size(filt,2)
            stop
        end if 

        ! Safety check
        if (mod(n,2) .ne. 1) then 
            write(*,*) "calc_ge_filt:: error: n can only be odd."
            write(*,*) "n = ", n
            stop  
        end if 


        ! Loop over filter array in two dimensions,
        ! calculate the distance from the center
        ! and impose correct Green function value. 

        filt = 0.
        
        do j = -n2, n2 
        do i = -n2, n2

            x  = i*dx 
            y  = j*dy
            
            r  = sqrt(x**2+y**2)  

            ! Get actual index of array
            i1 = i+1+n2 
            j1 = j+1+n2 

            ! Get correct GE value for this point (given by colatitude, theta)
            
            filt(i1,j1) = calc_gn_value(max(dy,r),r_earth,m_earth)* (dx*dy)

         end do
        end do

        return 

      end subroutine calc_GN_filter_2D


    function calc_gn_value(r,r_earth,m_earth) result(gn)
        !Info
        implicit none

        real(wp), intent(IN) :: r
        real(wp), intent(IN)  :: r_earth 
        real(wp), intent(IN)  :: m_earth 
        real(wp)              :: gn

        gn = r_earth / (2. * m_earth * sin (r/(2.*r_earth)) )

        return
    end function calc_gn_value

end module green_functions