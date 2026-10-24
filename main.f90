program main

    use initialize
    use forces
    use integrate

    implicit none

    integer :: N, dim, iterations, iter, simulation_type, nseed, integration_type, integration_direction
    integer, allocatable :: seed(:)
    real(kind=8), parameter :: R = 1.0d0
    real(kind=8) :: M_solar, M_total, m, omega, dt, eps
    real(kind=8), allocatable :: pos(:,:), vel(:,:), v_bar(:)
    character(len=50) :: output_filename, quantity_filename, output_backwards_filename, quantity_backwards_filename
    character(len=50) :: simulation_name
    real(kind=8) :: start_time, end_time, elapsed_time

    ! ---------------------------------------------------------------
    ! user prompt - simulation type
    ! ---------------------------------------------------------------

    print *, "Select simulation:"
    print *, "1 : 2D N-body"
    print *, "2 : 3D N-body"
    print *, "3 : two-body circular orbit"
    read(*,*) simulation_type

    ! ---------------------------------------------------------------
    ! define simulation parameters
    ! ---------------------------------------------------------------

    select case (simulation_type)

    case (1)
        simulation_name = "2D N-body"
        N = 100
        dim = 2
        iterations = 1000
        dt = 0.001d0
        eps = 0.01d0
        omega = 1.0d0
        output_filename = "data/data_2d.txt"
        quantity_filename = "data/quantity_2d.txt"
        output_backwards_filename = "data/data_2d_backwards.txt"
        quantity_backwards_filename = "data/quantity_2d_backwards.txt"

    case (2)
        simulation_name = "3D N-body"
        N = 100
        dim = 3
        iterations = 1000
        dt = 0.001d0
        eps = 0.01d0
        omega = 1.0d0
        output_filename = "data/data_3d.txt"
        quantity_filename = "data/quantity_3d.txt"
        output_backwards_filename = "data/data_3d_backwards.txt"
        quantity_backwards_filename = "data/quantity_3d_backwards.txt"

    case (3)
        simulation_name = "Two-body circular orbit"
        N = 2
        dim = 2
        iterations = 1000
        dt = 0.001d0
        eps = 0.0d0
        omega = 1.0d0
        output_filename = "data/data_two_body.txt"
        quantity_filename = "data/quantity_two_body.txt"
        output_backwards_filename = "data/data_two_body_backwards.txt"
        quantity_backwards_filename = "data/quantity_two_body_backwards.txt"

    case default
        print *, "invalid simulation type. valid selections: 1,2,3"
        stop 1

    end select

    ! ---------------------------------------------------------------
    ! user prompt - integration type
    ! ---------------------------------------------------------------
    print *, "Select integration type:"
    print *, "1 : leapfrog"
    print *, "2 : rk4"
    print *, "3 : rk8"
    read(*,*) integration_type

    ! ---------------------------------------------------------------
    ! user prompt - integartion direction
    ! ---------------------------------------------------------------

    print *, "Select integration direction:"
    print *, "1 : forward only"
    print *, "2 : forward and backwards"
    read(*,*) integration_direction

    ! ---------------------------------------------------------------
    ! allocate arrays and set variables
    ! ---------------------------------------------------------------

    ! define array sizes
    allocate(pos(N,dim))
    allocate(vel(N,dim))
    allocate(v_bar(dim))

    ! define variable values
    M_solar = 1.989d0*10d0**30 ! one solar mass
    M_total = 10d0**11 * M_solar ! total mass of system
    m = 1.0d0/N ! mass of each body (dimensionless)

    ! ---------------------------------------------------------------
    ! initialize positions and velocities
    ! ---------------------------------------------------------------
    
    if (simulation_type==3) then

        ! two-body circular orbit

        ! initial positions
        pos(1,:) = [-0.5d0, 0.0d0]
        pos(2,:) = [ 0.5d0, 0.0d0]

        ! initial velocities
        vel(1,:) = [0.0d0, -0.5d0]
        vel(2,:) = [0.0d0,  0.5d0]

    else

        ! initialize fixed seed for debugging
        call random_seed(size=nseed)
        allocate(seed(nseed))
        seed = 12347
        call random_seed(put=seed) 

        ! populate pos: generate random coordinates for each particle
        do iter=1,N
            call get_position_components(pos(iter,:), dim, R)
        end do

        ! populate vel: use position and angular momentum for each particle
        do iter=1,N
            call get_velocity_components(vel(iter,:), pos(iter,:), dim, omega)
        end do

        ! calculate barycentric velocity and subtract from velocity
        call barycentric_velocity(v_bar, vel, N, dim)

        do iter=1,N
            vel(iter,:) = vel(iter,:) - v_bar
        end do
 
    end if

    ! ---------------------------------------------------------------
    ! print simulation information
    ! ---------------------------------------------------------------
    print *
    print *, "simulation type:  ", simulation_name
    print *, "number of bodies:", N
    print *, "dimensions:", dim
    print *, "iterations:", iterations
    print *, "time step:", dt
    print *, "total time:", (iterations)*dt
    print *, "force regularization, eps:", eps
    print *, "angular velocity:", omega
    print *, "mass:", m

    ! start timer
    call cpu_time(start_time)

    ! ---------------------------------------------------------------
    ! begin integration
    ! ---------------------------------------------------------------
    select case (integration_type)

    case (1)
        print *, "integration type:  leapfrog"
        print *, "simulation is running..."
        call leapfrog(pos, vel, dt, N, dim, iterations, m, eps, output_filename, quantity_filename)
    case (2)
        print *, "integration type:  rk4"
        print *, "simulation is running..."
        call rk4(pos, vel, dt, N, dim, iterations, m, eps, output_filename, quantity_filename)
    case (3)
        print *, "integration type:  rk8"
        print *, "simulation is running..."
        call rk8(pos, vel, dt, N, dim, iterations, m, eps, output_filename, quantity_filename)

    case default
    print *, "invalid integration type. valid selections: 1,2,3"
    stop 1

    end select

    ! ---------------------------------------------------------------
    ! begin backwards integration - if applicable
    ! ---------------------------------------------------------------

    if (integration_direction == 2) then
        select case (integration_type)

        case (1)
            call leapfrog(pos, vel, -dt, N, dim, iterations, m, eps, output_backwards_filename, quantity_backwards_filename)
        case (2)
            call rk4(pos, vel, -dt, N, dim, iterations, m, eps, output_backwards_filename, quantity_backwards_filename)
        case (3)
            call rk8(pos, vel, -dt, N, dim, iterations, m, eps, output_backwards_filename, quantity_backwards_filename)
        end select
        
    end if

    print *
    print *, "simulation ran successfully..."
    print *, "data saved to:  ", output_filename

    ! stop timer
    call cpu_time(end_time)

    ! ---------------------------------------------------------------
    ! calculate and print elapsed time
    ! ---------------------------------------------------------------

    elapsed_time = end_time - start_time

    print *, "Execution time (s): ", elapsed_time, "for N =", N, "particles, ", iterations, "iterations."

end program main