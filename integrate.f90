module integrate

    use initialize
    use forces

    implicit none

contains

    subroutine leapfrog(pos, vel, dt, N, dim, iterations, mass, eps, save_filename, quantity_filename)

        ! ---------------------------------------------------------------
        ! leapfrog(): integrates the equations of motion using the leapfrog method
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim, iterations
        real(kind=8), intent(inout) :: pos(N,dim), vel(N,dim)
        real(kind=8), intent(in) :: dt, mass, eps
        character(len=*), intent(in) :: save_filename, quantity_filename
        integer :: t, i
        real(kind=8) :: Acc(N,dim), Epot, Ekin, L(3), v_bar(dim)

        open(unit=10, file=save_filename, status="replace")
        open(unit=11, file=quantity_filename, status="replace")

        ! Save initial conditions
        do i= 1, N
            write(10, *) pos(i,:), vel(i,:)
        end do

        ! calculate initial conserved quantities
        call potential_energy(Epot, pos, N, dim, eps)
        call kinetic_energy(Ekin, vel, N, dim)
        call angular_momentum(L, pos, vel, N, dim)
        call barycentric_velocity(v_bar, vel, N, dim)

        ! save quantities at t = 0
        write(11, *) Epot, Ekin, L, v_bar

        ! initial acceleration
        do i=1,N
            call acceleration_vector(Acc(i,:), pos, N, dim, i, mass, eps)
        end do

        do t=1,iterations ! loop through all time iterations

            ! half kick
            vel = vel + dt/2.0d0 * Acc         

            ! drift
            pos = pos + dt * vel

            ! half kick
            do i=1,N
                call acceleration_vector(Acc(i,:), pos, N, dim, i, mass, eps)
            end do

            vel = vel + dt/2.0d0 * Acc


            ! sava data to file 
            do i= 1, N
                write(10, *) pos(i,:), vel(i,:)
            end do

            ! calculate energy 
            call potential_energy(Epot, pos, N, dim, eps)
            call kinetic_energy(Ekin, vel, N, dim)
            call angular_momentum(L, pos, vel, N, dim)
            call barycentric_velocity(v_bar, vel, N, dim)

            ! save data to file
            write(11, *) Epot, Ekin, L, v_bar

        end do

        ! close text files
        close(10)
        close(11)

    end subroutine leapfrog

    subroutine rk4(pos, vel, dt, N, dim, iterations, mass, eps, save_filename, quantity_filename)

        ! ---------------------------------------------------------------
        ! rk4(): integrates the equations of motion using the 4th order Runge-Kutta method
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim, iterations
        real(kind=8), intent(inout) :: pos(N,dim), vel(N,dim)
        real(kind=8), intent(in) :: dt, mass, eps
        real(kind=8) :: k1(N,dim), k2(N,dim), k3(N,dim), k4(N,dim)
        real(kind=8) :: l1(N,dim), l2(N,dim), l3(N,dim), l4(N,dim), pos_temp(N,dim)
        character(len=*), intent(in) :: save_filename, quantity_filename
        integer :: t, i
        real(kind=8) :: Epot, Ekin, L(3), v_bar(dim)

        open(unit=10, file=save_filename, status="replace")
        open(unit=11, file=quantity_filename, status="replace")

        ! save initial conditions
        do i= 1, N
            write(10, *) pos(i,:), vel(i,:)
        end do

        ! calculate initial conserved quantities
        call potential_energy(Epot, pos, N, dim, eps)
        call kinetic_energy(Ekin, vel, N, dim)
        call angular_momentum(L, pos, vel, N, dim)
        call barycentric_velocity(v_bar, vel, N, dim)

        ! save quantities at t = 0
        write(11, *) Epot, Ekin, L, v_bar

        do t=1,iterations ! loop through all time iterations

            k1 = vel
            
            do i=1,N
                call acceleration_vector(l1(i,:), pos, N, dim, i, mass, eps)
            end do

            k2 = vel + dt/2.0d0 * l1

            pos_temp = pos + dt/2.0d0 * k1

            do i=1,N
                call acceleration_vector(l2(i,:), pos_temp, N, dim, i, mass, eps)
            end do

            k3 = vel + dt/2.0d0 * l2

            pos_temp = pos + dt/2.0d0 * k2

            do i=1,N
                call acceleration_vector(l3(i,:), pos_temp, N, dim, i, mass, eps)
            end do

            k4 = vel + dt * l3

            pos_temp = pos + dt * k3

            do i=1,N
                call acceleration_vector(l4(i,:), pos_temp, N, dim, i, mass, eps)
            end do 

            ! update position and velocity
            pos = pos + dt * 1.0d0/6.0d0 * (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)
            vel = vel + dt * 1.0d0/6.0d0 * (l1 + 2.0d0*l2 + 2.0d0*l3 + l4)

            ! sava data to file 
            do i= 1, N
                write(10, *) pos(i,:), vel(i,:)
            end do

            ! calculate conserved quantities
            call potential_energy(Epot, pos, N, dim, eps)
            call kinetic_energy(Ekin, vel, N, dim)
            call angular_momentum(L, pos, vel, N, dim)
            call barycentric_velocity(v_bar, vel, N, dim)

            ! save data to file
            write(11, *) Epot, Ekin, L, v_bar

        end do

        close(10)
        close(11)
        
    end subroutine rk4


    subroutine rk8(pos, vel, dt, N, dim, iterations, mass, eps, save_filename, quantity_filename)

        ! ---------------------------------------------------------------
        ! rk8(): integrates the equations of motion using the 8th order Runge-Kutta method
        ! ---------------------------------------------------------------

        integer, intent(in) :: N, dim, iterations
        real(kind=8), intent(inout) :: pos(N,dim), vel(N,dim)
        real(kind=8), intent(in) :: dt, mass, eps
        character(len=*), intent(in) :: save_filename, quantity_filename
        integer, parameter :: stages = 12
        integer :: t, i, j, s, ios
        real(kind=8) :: A(stages,stages), B(stages), k(stages,N,dim), p(stages,N,dim)
        real(kind=8) :: pos_temp(N,dim), vel_temp(N,dim)
        real(kind=8) :: Epot, Ekin, L(3), v_bar(dim)

        ! NOTE on variables: using p instead of l for RK8 coefficient as L is already used for angular momentum

        ! read RK8 coefficients
        open(unit=12, file="dop853_coefficients.txt", status="old", action="read", iostat=ios)

        ! check file opened with no error

        if (ios /= 0) then
            print *, 'cannot open dop853_coefficients.txt'
            stop
        end if

        ! populate A and B coefficient lists
        do i=1,stages
            read(12,*) A(i,:)
        end do

        read(12,*) B
        close(12)

        open(unit=10, file=save_filename, status="replace")
        open(unit=11, file=quantity_filename, status="replace")

        ! save initial conditions
        do i=1,N
            write(10,*) pos(i,:), vel(i,:)
        end do

        ! calculate initial conserved quantities
        call potential_energy(Epot, pos, N, dim, eps)
        call kinetic_energy(Ekin, vel, N, dim)
        call angular_momentum(L, pos, vel, N, dim)
        call barycentric_velocity(v_bar, vel, N, dim)

        write(11,*) Epot, Ekin, L, v_bar

        ! save quantities at t = 0
        do t=1,iterations

            ! evaluate all RK8 stages
            do s=1,stages

                pos_temp = pos
                vel_temp = vel

                do j=1,s-1
                    pos_temp = pos_temp + dt*A(s,j)*k(j,:,:)
                    vel_temp = vel_temp + dt*A(s,j)*p(j,:,:)
                end do

                ! position derivative = velocity
                k(s,:,:) = vel_temp

                ! velocity derivative = acceleration
                do i=1,N
                    call acceleration_vector(p(s,i,:), pos_temp, N, dim, i, mass, eps)
                end do

            end do

            ! update position and velocity
            do s=1,stages
                pos = pos + dt*B(s)*k(s,:,:)
                vel = vel + dt*B(s)*p(s,:,:)
            end do

            ! sava data to file 
            do i=1,N
                write(10,*) pos(i,:), vel(i,:)
            end do

            ! calculate conserved quantities
            call potential_energy(Epot, pos, N, dim, eps)
            call kinetic_energy(Ekin, vel, N, dim)
            call angular_momentum(L, pos, vel, N, dim)
            call barycentric_velocity(v_bar, vel, N, dim)

            write(11,*) Epot, Ekin, L, v_bar

        end do

        close(10)
        close(11)

    end subroutine rk8



end module integrate