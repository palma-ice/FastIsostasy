module green_functions

    use isostasy_defs, only : wp, pi
    use isos_utils

    implicit none

    private
    
    public :: calc_viscous_green
    public :: calc_elastic_green
    public :: calc_z_ss_green
    public :: get_ge_value
    public :: calc_gn_value
    
    contains


    ! The Green's function (Eq. 3 of Coulon et al, 2021) gives displacement G in [m]
    ! as a function of the distance r from the point load P_b [Pa].

    ! Here GV is calculated, which is G without including the point load. GV has units
    ! of [m N-1]. This can then be multiplied with the actual magnitude of the
    ! point load to obtain G.
    ! G = GV * P_b = [m N-1] * [Pa] = [m]. 

    ! Note that L_w contains information about rho_uppermantle. 
    subroutine calc_viscous_green(GV, kei2D, L_w, D_lith, dx, dy)
        implicit none

        real(wp), intent(OUT) :: GV(:, :) 
        real(wp), intent(IN)  :: kei2D(:, :) 
        real(wp), intent(IN)  :: L_w 
        real(wp), intent(IN)  :: D_lith 
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy
        
        GV = -L_w**2 / (2.0*pi*D_lith) * kei2D * (dx*dy)

        return
    end subroutine calc_viscous_green

    subroutine calc_elastic_green(filt, dx, dy)
        ! Calculate 2D Green Function

        implicit none 

        real(wp), intent(OUT) :: filt(:, :) 
        real(wp), intent(IN)  :: dx
        real(wp), intent(IN)  :: dy

        ! Local variables 
        integer  :: i, j, i1, j1, nx, ny, mx, my
        real(wp) :: x, y, r

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
        nx  = size(filt, 1)
        ny  = size(filt, 2)
        mx = midindex(nx)
        my = midindex(ny)

        ! Fill matrix of Green function.
        filt = 0.
        do j = 1, ny
        do i = 1, nx

            x  = (i-1-mx)*dx
            y  = (j-1-mx)*dy
            r  = sqrt(x**2+y**2)

            ! Get correct GE value for this point
            filt(i, j) = get_ge_value(r, rn_vals, ge_vals) * 1.e-12 * (dx*dy) / &
                (9.81*max(r, dx))
        end do
        end do

        return 
    end subroutine calc_elastic_green

    ! Compute Green function (Coulon et al. 2021) to determine the z_ss perturbation.
    subroutine calc_z_ss_green(filt, m_earth, r_earth, dx, dy) 
        implicit none 

        real(wp), intent(OUT) :: filt(:, :) 
        real(wp), intent(IN)  :: m_earth
        real(wp), intent(IN)  :: r_earth
        real(wp), intent(IN)  :: dx 
        real(wp), intent(IN)  :: dy   

        ! Local variables 
        integer  :: i, j, i1, j1, nx, ny, mx, my
        real(wp) :: x, y, r

        ! Get size of filter array and half-width
        nx  = size(filt, 1)
        ny  = size(filt, 2)
        mx = midindex(nx)
        my = midindex(ny)

        ! Fill matrix of Green function.
        filt = 0.
        do j = 1, ny
        do i = 1, nx

            x  = (i-1-mx)*dx
            y  = (j-1-mx)*dy
            r  = sqrt(x**2+y**2)

            ! Get correct GE value for this point (given by colatitude, theta)
            filt(i, j) = calc_gn_value(max(dy,r), r_earth, m_earth)
        end do
        end do

        return
    end subroutine calc_z_ss_green


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


    function calc_gn_value(r, r_earth, m_earth) result(gn)
        implicit none
        real(wp), intent(IN) :: r
        real(wp), intent(IN)  :: r_earth 
        real(wp), intent(IN)  :: m_earth 
        real(wp)              :: gn

        gn = r_earth / (2. * m_earth * sin (r/(2.*r_earth)) )

        return
    end function calc_gn_value

end module green_functions