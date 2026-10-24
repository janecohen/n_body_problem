module initialize

    implicit none

contains

    subroutine get_position_components(position, dim, R)

        ! ---------------------------------------------------------------
        ! get_position_components(): generates random coordinates for a particle 
        ! within a circle or sphere of radius R
        ! ---------------------------------------------------------------

        real(kind=8), intent(in) :: R
        integer, intent(in) :: dim
        real(kind=8), intent(out) :: position(dim)
        real(kind=8) :: mag

        do 

            call random_number(position)
            position = (2.0d0*position - 1.0d0)*R ! values in [-R,R]

            mag = SQRT(SUM(position**2))

            if (mag <= R) exit
                
        end do

    end subroutine get_position_components

    subroutine get_velocity_components(velocity, position, dim, omega)

        ! ---------------------------------------------------------------
        ! get_velocity_components(): generates velocity components for a particle
        ! ---------------------------------------------------------------

        integer, intent(in) :: dim
        real(kind=8), intent(in) :: position(dim), omega
        real(kind=8), intent(out) :: velocity(dim)

        velocity = 0.0d0

        velocity(1) = -omega*position(2) 
        velocity(2) = omega*position(1)

    end subroutine get_velocity_components

    subroutine barycentric_velocity(v_bar, vel, N, dim)

        ! ---------------------------------------------------------------
        ! barycentric_velocity(): calculates the barycentric velocity of a system
        ! ---------------------------------------------------------------
        
        integer, intent(in) :: N, dim
        real(kind=8), intent(in) :: vel(N,dim)
        real(kind=8), intent(out) :: v_bar(dim)

        v_bar = SUM(vel, dim=1) / N

    end subroutine barycentric_velocity

end module initialize