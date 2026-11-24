PROGRAM VERLET
  IMPLICIT NONE
  INTEGER :: K
  INTEGER, PARAMETER :: nsteps = 600

  REAL :: mass, tau, time
  REAL :: fx, fy, fz
  REAL :: fx_new, fy_new, fz_new
  REAL :: x, y, z
  REAL :: vx, vy,vz
  REAL :: x_new, y_new, z_new
  REAL :: vx_new, vy_new, vz_new
  ! k = iteration, tau = 0.2s, time = time taken, fx, fy, fz = initial force, x, y, z = inital position, vx, vy, vz = current velocity, x_new,.... = new position, vx_new,... = new velocity (assigning the variables)
  mass = 1.0D0
  tau = 0.2D0
  fx = 0D0
  fy = 0.1D0
  fz = 0.0D0
  x = 0.0D0
  y = 0.0D0
  z = 0.0D0
  vx = 0.0D0
  vy = 0.0D0
  vz = 0.0D0
  ! inital condition, we have to give values for those variable
  !verlet velocity algorithm
  DO K = 1, nsteps
  !calculate new postion
  x_new = x + tau*vx + (tau**2/(2.0D0*mass))*fx
  y_new = y + tau*vy + (tau**2/(2.0D0*mass))*fy
  z_new = z + tau*vz + (tau**2/(2.0D0*mass))*fz

  !calculate the force at the new time, here force is constant
  fx_new = fx
  fy_new = fy
  fz_new = fz

  !calculate the new velocity
  vx_new = vx + (tau/(2.0D0*mass))*(fx+fx_new)
  vy_new = vy + (tau/(2.0D0*mass))*(fy+fy_new)
  vz_new = vz + (tau/(2.0D0*mass))*(fz+fz_new)

  !update the old postion and velocity
  x = x_new
  y = y_new
  z = z_new

  vx = vx_new
  vy = vy_new
  vz = vz_new

  END DO

  !calculate the final simulaton time
  time = nsteps*tau

  !display the final result
  PRINT*, "Final time (s):", time

  PRINT*, "Final postion (m):"
  PRINT*, "X=", x
  PRINT*, "Y=", y
  PRINT*, "Z=", z
  PRINT*, "Final velocity (m/s):"
  PRINT*, "Vx =", vx
  PRINT*, "Vy =", vy
  PRINT*, "Vz =", vz
END PROGRAM VERLET
 
