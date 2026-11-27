PROGRAM VERLET_2
  IMPLICIT NONE
  !DECLARE integers, variables and double precision
  INTEGER :: i, k, nk, n, input_unit
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND (p=13, r=300)
  !tau, time scalars (one value for the whole simulation in this case force is constant but we didnt mention it in scalar, because we are generalising it, which will be used for next problem
  REAL (KIND=wp) :: tau, time
  
  !mass, position, velocity,force - array. all particle has its own value in three dimension
  REAL (KIND=wp), ALLOCATABLE :: m(:)
  REAL (KIND=wp), ALLOCATABLE :: x(:,:), v(:,:)
  REAL (KIND=wp), ALLOCATABLE :: f(:,:), fnext(:,:)

  !open input and read simulation paramters 
  OPEN(UNIT=11, FILE="input.txt", STATUS="old", ACTION="read")

  !read numbermof time steps, time intreval and number of particles
  READ(UNIT=11, FMT=*) nk, tau
  READ(UNIT=11, FMT=*) n
  
  !ALLOCATE memory of array according to the particle)
  ALLOCATE(m(n))
  ALLOCATE(x(n,3), v(n,3), f(n,3), fnext(n,3))

  !read initial velocity positon
  DO i = 1, n
    READ(UNIT=11, FMT =*) m(i), x(i,1), x(i,2), x(i,3), v(i,1), v(i,2), v(i,3)
  END DO
  !close the file once it read
  CLOSE(11)

  !assign the force for all the particle
  f(i,1) = 0.0_wp
  f(i,2) = 0.1_wp
  f(i,3) = 0.0_wp
  
  !verlet algorithm start here
  !outer loop: reapeat this calculation for each time steps
  DO k = n, nk
    !inner loop: calculate this for each particle
    DO i = 1, n
      !calculate new position, we are using dynamic array (:,:)
      x(i,:) = x(i,:) + tau*v(i,:) + (tau**2/(2.0_wp*m(i)))*f(i,:)
  
 
      !calculate new velocities, we are using dynamic array (:,:)
      v(i,:) = v(i,:) + (tau/(2.0_wp*m(i)))*(f(i,:)+fnext(i,:))

      f(i,:) = fnext(i,:)
      
    END DO
  END DO
  !calculate total simulation time
  time = nk*tau
  
  !print the final time 
  PRINT*, "final time (s):", time
  DO i = 1, n
    x(i,2) = x(i,2) + v(i, 2) + (1/2*(f(i,2)/m(i)))* time**2
    v(i,2) = v(i,2) + (f(i,2)/m(i))*time**2
  END DO
  !loop -  print new velocity and position for all the particles
  DO i = 1, n
    PRINT*, "Particle:", i
    PRINT*, "Final Position:", x(i,:)
    PRINT*, "Final velocity:", v(i,:)
    PRINT*, "PARTICLE POSITION AT Y-DIMENSION:", x(i,2)
    PRINT*, "PARTICLE VELOCITY AT Y-DIMENSION:", v(i,2)
  END DO

  ! deallocate the memory array after its used
  DEALLOCATE(m, x, f, v, fnext)
 

  !end the program 
  
  END PROGRAM VERLET_2
