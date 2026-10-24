module forces
    
    implicit none

contains

    subroutine force_vector(pos_i, pos_j, F, N, dim, eps)

        ! ---------------------------------------------------------------
        ! force_vector(): calculates the force between two particles
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim
        real(kind=8), intent(in) :: pos_i(dim), pos_j(dim), eps
        real(kind=8), intent(out) :: F(dim)
        real(kind=8) :: dist

        dist = SQRT(SUM((pos_j - pos_i)**2))

        F = (pos_j - pos_i) / (N**2 * (dist**2 + eps**2)**1.5d0)

    end subroutine force_vector

    subroutine acceleration_vector(Acc, pos, N, dim, i, mass, eps)

        ! ---------------------------------------------------------------
        ! acceleration_vector(): calculates the acceleration of particle i due to all 
        ! other particles
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim, i
        real(kind=8), intent(in) :: pos(N,dim), mass, eps
        real(kind=8), intent(out) :: Acc(dim)
        integer :: j
        real(kind=8) :: F(dim), sum_forces(dim)

        sum_forces = 0.0d0

        do j=1, N ! iterate through each particle acting on particle i

            if (i == j) cycle

            call force_vector(pos(i,:), pos(j,:), F, N, dim, eps)

            sum_forces = sum_forces + F

        end do  

        Acc = sum_forces/mass
     
    end subroutine acceleration_vector

    subroutine potential_energy(Epot, pos, N, dim, eps)

        ! ---------------------------------------------------------------
        ! potential_energy(): calculates the total potential energy of the system
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim
        real(kind=8), intent(in) :: pos(N,dim), eps
        real(kind=8), intent(out) :: Epot
        integer :: i, j
        real(kind=8) :: Epot_i

        Epot = 0.0d0

        do i=1,N
            Epot_i = 0.0d0

            do j=1,N
                if (i == j) cycle
                Epot_i = Epot_i + 1.0d0 / (N**2 * SQRT(SUM((pos(j,:) - pos(i,:))**2) +eps**2))
            end do

            Epot = Epot + Epot_i
        end do

        Epot = -Epot/2.0d0 ! account for double effect of looping i,j

    end subroutine potential_energy

    subroutine kinetic_energy(Ekin, vel, N, dim)

        ! ---------------------------------------------------------------
        ! kinetic_energy(): calculates the total kinetic energy of the system
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim
        real(kind=8), intent(in) :: vel(N,dim)
        real(kind=8), intent(out) :: Ekin
        integer :: i

        Ekin = 0.0d0

        do i=1,N
            Ekin = Ekin + SUM(vel(i,:)**2) / N
        end do

        Ekin = Ekin / 2.0d0

    end subroutine kinetic_energy

    subroutine angular_momentum(L, pos, vel, N, dim)

        ! ---------------------------------------------------------------
        ! angular_momentum(): calculates the total angular momentum of the system
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim
        real(kind=8), intent(in) :: pos(N,dim), vel(N,dim)
        real(kind=8), intent(out) :: L(3)
        real(kind=8) :: prod(3)
        integer :: i

        L = 0.0d0

        do i=1,N
            call cross_product(pos(i,:),vel(i,:), prod, dim)
            L = L + prod / N
        end do

    end subroutine angular_momentum

    subroutine cross_product(a, b, c, dim)

        ! ---------------------------------------------------------------
        ! cross_product(): calculates the cross product of two vectors
        ! ---------------------------------------------------------------

        implicit none
        real(kind=8), intent(in) :: a(dim), b(dim)
        integer, intent(in) :: dim
        real(kind=8), intent(out) :: c(3)
        real(kind=8) :: a3(3), b3(3)

        a3 = 0.0d0
        b3 = 0.0d0

        a3(:dim) = a
        b3(:dim) = b

        c(1) = a3(2)*b3(3) - a3(3)*b3(2)
        c(2) = a3(3)*b3(1) - a3(1)*b3(3)
        c(3) = a3(1)*b3(2) - a3(2)*b3(1)

    end subroutine cross_product

end module forces