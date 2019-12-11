
SUBROUTINE basin_1d(dx, N_x, dt, T, fileout)
  IMPLICIT NONE
  USE constants
  INTEGER, INTENT(IN) :: N_x
  REAL(DP), INTENT(OUT) :: bc_0, bc_N
  REAL(DP), INTENT(IN) :: dx, dt, upper_time
  CHARACTER(LEN = 20), INTENT(IN):: fileout

  CALL rossby_solver(dx, N_x, dt, T, fileout)

  bc_0 = 0.0   !Boundary condition
  bc_N = 0.0   !Boundary condition
  
END SUBROUTINE basin_1d

SUBROUTINE initial_condition_1d(int_psi, int_zeta)
  IMPLICIT NONE
  USE constants
  INTEGER :: i  
  REAL(DP) :: eps
  REAL(DP), INTENT(IN) :: dx, dt, upper_time, int_psi, int_zeta
  REAL(DP), INTENT(IN) :: bc_N, bc_0
  REAL(DP), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: psi_0, zeta_0

  eps = 1.0E-10
  
  IF (abs(int_psi(1) - bc_0) > eps .and. abs(int_psi(N_x)-bc_N) > eps) THEN
     print*, ini_psi(1) "," ini_psi(N_x)
     print*, "Error: Chosen initial condition does not satisfy BC!"
     print*, "Terminating program..." 
     Exit
  END IF

  ALLOCATE (psi_0(N_x), zeta_0(N_x))

  DO i = 1, N_x
     psi_0(i) = int_psi(i)
     zeta_0(i) = int_zeta(i)
  END DO
  
END SUBROUTINE initial_condition_1d

!* Function that solves the vorticity equation as a set of two coupled PDE's. Uses Forward Euler
! * time-stepping for the vorticity and centered differences for the the streamfunction. Solves the
! * 1D Poisson equation at every time step using a tridiagonal linear algebra solver (Thomas algorithm)*/

SUBROUTINE basin_euler
  IMPLICIT NONE
  USE constants
  INTEGER :: n, i
  REAL(DP) :: alpha, dxdx, t
  REAL(DP), INTENT(IN) ::  dt, dx, bc_0, bc_N, upper_time
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: diag, rhs_tridiag
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_prev, zeta_curr
  REAL(DP), DIMENSION(:), INTENT(IN) :: psi_0, zeta_0 

  !!display_model_config()

  alpha = dt/(2.0*dx)
  dxdx  = dx*dx
  t     = 0.0
  
  !! allocate previous and current states
  ALLOCATE (psi_prev(N_x),psi_curr(N_x), zeta_prev(N_x), zetz_curr(N_x), diag(N_x-2), rhs_tridiag(N_x-2))
  
  !set initial and boundary conditions
  DO i = 1, N_x
     psi_prev(i)  = psi_0(i)
     zeta_prev(i) = zeta_0(i)
  END DO

  psi_curr(1)    = bc_0
  psi_curr(N_x)  = bc_N
  zeta_curr(1)   = zeta_prev(1)
  zeta_curr(N_x) = zeta_prev(N_x)

  CALL write_state_to_file(0.0, psi_prev) !! Write initial condition to file
  
  DO n = 2, upper_time-1

     ! forward euler for vorticity
     DO i = 2, N_x-1
       zeta_curr(i) = zeta_prev(i) = alpha*(psi_prev(i+1)-psi_prev(i-1))
     END DO

     ! solve poisson equation for new streamfunction
     DO i = 2, N_x-1
        rhs_tridiag(i-1) = -dxdx*zeta_curr(i)   
     END DO

     ! interior solver by poisson
     CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_curr+1)

     DO i = 2, N_x-1
        psi_prev(i) = psi_curr(i)
        zeta_prev(i) = zeta_curr(i)
     END DO

     t = t+dt

     IF (n%20 .eq. 0) THEN
        CALL write_state_to_file(t, psi_curr)
     END IF 
  END DO
     
   DEALLOCATE(psi_prev, psi_curr, zeta_prev, zeta_curr, diag, rhs_tridiag)

END SUBROUTINE basin_euler


SUBROUTINE basin_leapfrog
  IMPLICIT NONE
  USE constants
  INTEGER :: n, i
  REAL(DP) :: alpha, gamm, dxdx, t
  REAL(DP), INTENT(IN) :: dt, dx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: diag, rhs_tridiag
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_prev, zeta_pp, zeta_curr
  REAL(DP), DIMENSION(:), INTENT(IN) :: psi_0, zeta_0

  !! display_model_config()

  alpha = dt/(2.0*dx)
  gramm = dt/dx 
  dxdx  = dx*dx
  t     = 0.0

  ALLOCATE (psi_prev(N_x), psi_curr(N_x), zeta_prev(N_x), zeta_pp(N_x), zeta_curr(N_x), diag(N_x-2), rhs_tridiag(N_x-2))
  
  DO i = 1, N_x
     psi_prev(i) = psi_0(i)
     zeta_prev(i) = zeta_0(i)
  END DO

  psi_curr(1)    = bc_0
  psi_curr(N_x)  = bc_N
  zeta_curr(1)   = zeta_pp(1)
  zeta_curr(N_x) = zeta_pp(N_x)

  CALL write_state_to_file(0.0, psi_prev);   !Write initial condition to file
  
  ! Initial Euler
  DO n = 2, N_x-1
     zeta_prev(i) = zeta_0(i) - alpha*(psi_0(i+1)-psi_0(i-1))
  END DO

  DO i = 2, N_x-1
     rhs_tridiag(i-1) = -dxdx*zeta_prev(i)
  END DO

  CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_prev+1)
  
  ! time loop for vorticity using leapfrog
  DO n = 3, upper_time-1

     ! forward leapfrog for vorticity
     DO i = 2, N_x-1
        zeta_curr(i) = zeta_pp(i) - gamm*(psi_prev(i+1) - psi_prev(i-1))
     END DO
 
     !Solve poisson equation for new streamfunction
     DO i = 2, N_x-1
        rhs_tridiag(i-1) = -dx*dx*zeta_curr(i)
     END DO
    
     CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_curr+1)

     DO i = 2, N_x-1
        psi_prev(i)  = psi_curr(i)
        zeta_pp(i)   = zeta_prev(i)
        zeta_prev(i) = zeta_curr(i)
     END DO

     t = t+dt

     IF (n%20 .eq. 0) THEN
        CALL write_state_to_file(t, psi_curr)
     END IF
  END DO

  DEALLOCATE (psi_prev, psi_curr, zeta_prev, zeta_pp, zeta_curr, diag, rhs_tridiag)

END SUBROUTINE basin_leapfrog

!Poisson equation solver
SUBROUTINE tridiagonal(b, y, N, solution)
  INTEGER :: i 
  INTEGER, INTENT(IN) :: N
  REAL(DP), DIMENSION(:), INTENT(IN) :: b, y, solution
 
  b(1) = 2.0

  ! forward substitution
  DO i = 2, N
     b(i) = (i + 2)/REAL(DP)(i + 1)
     y(i) = y(i) + (y(i-1)/b(i-1));
  ENDDO

  ! backward substitution
  solution(N) = y(N)/b(N)
  DO i = N-1, 1, -1
     solution(i) = (y(i)+solution(i+1))/b(i)
  ENDDO

END SUBROUTINE tridiagonal


SUBROUTINE poision_jacobi(g, bc_0y, bc_1y, bc_x0, bc_x1, dx, dy, N_x, N_y, max_iter, f)
  INTEGER :: i, j, iter
  INTEGER, INTENT(IN) :: N_x, N_y, max_iter
  REAL(DP) :: dxdx, dydy, dxdxdydy, dxdx_plus_dydy, eps, diff
  REAL(DP), INTENT(IN) :: g, bc_0y, bc_1y, bc_x0, yyyyyyyyyy, dx, dy, f
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: f_tmp

  dxdx = dx*dx
  dydy = dy*dy
  dxdxdydy = dxdx*dydy
  dxdx_plus_dydy = 2*(dxdx+dydy)

  DO j = 1, N_y
     f(1+j) = bc_0y(j)
     f((N_x-1)*N_y + j) = bc_1y(j)
  ENDDO

  DO i = 1, N_x
     f(i*N_y + 0) = bc_x0(i)
     f(i*N_y + (N_y-1)) = bc_x1(i)
  ENDDO

   iter = 0
   diff = 1.0E+20
   eps  = 1.0E-6
   ALLOCATE (f_tmp, N_x*N_y)

   DO WHILE (iter <= max_iter .and. abs(diff) > eps) 
      diff = 0.0
      
      DO i = 1, N_x
         DO j = 1, N_y
            f_tmp(i*N_y+j) = f(i*N_y+j)
         ENDDO
      ENDDO

      DO i = 2, N_x-1
         DO j = 2, N_y-1
            f(i*N_y + j) = (dydy*(f_tmp((i+1)*N_y + j) + f_tmp((i-1)*N_y + j)) + dxdx*(f_tmp(i*N_y + (j+1)) + f_tmp(i*N_y + (j-1))) - dxdxdydy*g(i*N_y + j))/dxdx_plus_dydy;
                diff += f(i*N_y + j) - f_tmp(i*N_y + j);
         ENDDO
      ENDDO

      iter = iter+1
   END DO
   DEALLOCATE(f_tmp)


END SUBROUTINE poisson_jacobi
