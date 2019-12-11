MODULE constants
  INTEGER,  PARAMETER, PUBLIC :: dp = KIND(1.0D0)
  INTEGER,  PARAMETER, PUBLIC :: dpc = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
  REAL(DP), PARAMETER, PUBLIC :: pi=4.D0*DATAN(1.D0)
END MODULE constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!PROGRAM rossby_wave
!  USE constants
!  IMPLICIT NONE

!  CALL simulation_1d()
!  CALL simulation_2d()
!  CALL funciton(array, n)

!END PROGRAM rossby_wave

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM simulation_1d
  USE constants
  IMPLICIT NONE
  INTEGER :: N_x, i
  CHARACTER(len = 30) :: output_02
  REAL(DP) :: L, dx, upper_time, dt, sigma, x
  REAL(DP), ALLOCATABLE, DIMENSION(:):: ini_psi, ini_zeta, ini_psi_gauss, ini_zeta_gauss
    
  L   = 1.0               !Lenght of domain
  dx  = 1.0/40.0          !Spatial step in x-direction
  N_x = INT(L/dx +1)           !Number of spatial point in x-direction
  upper_time  = 200.0     !Upper time limit
  dt  = 0.01              !Time step
  sigma = 0.1             !Gaussian initial condition

  ALLOCATE (ini_psi(N_x), ini_zeta(N_x), ini_psi_gauss(N_x), ini_zeta_gauss(N_x))
  
 ! ini_psi(1) = 0.0
 ! ini_zeta(1) = 0.0
 ! ini_psi_gauss(1) = exp(-((0-0.5)/sigma)*((0-0.5)/sigma))
 ! ini_zeta_gauss(1) = (4.0*((0 - 0.5)/(sigma*sigma))*((0 - 0.5)/(sigma*sigma)) -2.0/(sigma*sigma))*exp(-(((0-0.5)/sigma)**2))
   
  DO i = 1, N_x
     x = (i-1)*dx
     ini_psi(i) = sin(4.0*pi*x)
     ini_zeta(i) = -16.0*pi*pi*sin(4.0*pi*x)
     ini_psi_gauss(i) = exp(-((x-0.5)/sigma)*((x-0.5)/sigma))
     ini_zeta_gauss(i) = (4.0*((x - 0.5)/(sigma*sigma))*((x - 0.5)/(sigma*sigma)) -2.0/(sigma*sigma))*exp(-(((x-0.5)/sigma)**2))
     !print*, x, ini_psi(i)
  END DO

  output_02 = "sin_basin_euler.dat"
  CALL basin_1d(output_02)
  CALL initial_condition_1d(ini_psi, ini_zeta, N_x)
  !CALL basin_leapfrog(N_x)
  CALL basin_euler(N_x)

  DEALLOCATE (ini_psi, ini_zeta, ini_psi_gauss, ini_zeta_gauss)
END PROGRAM simulation_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE basin_1d(fileout)
  USE constants
  INTEGER :: N_x
  REAL(DP) :: bc_0, bc_N
  REAL(DP) :: dx, dt, upper_time
  CHARACTER(LEN = 30), INTENT(IN):: fileout

  !CALL rossby_solver(dx, N_x, dt, upper_time, fileout)
  OPEN(UNIT=6,FILE=fileout) 

  bc_0 = 0.0   !Boundary condition
  bc_N = 0.0   !Boundary condition
  
END SUBROUTINE basin_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE initial_condition_1d(ini_psi, ini_zeta, N_x)
  USE constants
  INTEGER :: i 
  INTEGER :: N_x 
  REAL(DP) :: eps
  REAL(DP), DIMENSION(N_x), INTENT(IN) :: ini_psi, ini_zeta
  REAL(DP) :: bc_N, bc_0
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_0, zeta_0

  eps = 1.0E-10
  DO i = 1, N_x
  !print*, ini_psi(i), ini_zeta(i)
  ENDDO
  IF (abs(ini_psi(1) - bc_0) .lt. eps .and. abs(ini_psi(N_x)-bc_N) .lt. eps) THEN
     print*, ini_psi(1) , ini_psi(N_x)
     print*, "Error: Chosen initial condition does not satisfy BC!"
     print*, "Terminating program..." 
     STOP
  END IF

  ALLOCATE (psi_0(N_x), zeta_0(N_x))

  DO i = 1, N_x
     psi_0(i) = ini_psi(i)
     zeta_0(i) = ini_zeta(i)
   !  print*, psi_0(i), zeta_0(i)
  END DO

  
END SUBROUTINE initial_condition_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve the vorticity equation using Forward Euler
! time-stepping for the vorticity and centered differences for the the streamfunction. 
! Solves the 1D Poisson equation using a tridiagonal linear algebra solver (Thomas algorithm)

SUBROUTINE basin_euler(N_x)
  USE constants
  INTEGER :: n, i
  INTEGER :: N_X
  REAL(DP) :: alpha, dxdx, t, upper_time
  REAL(DP) ::  dt, dx, bc_0, bc_N
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: diag, rhs_tridiag
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, zeta_prev, zeta_curr
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_curr
  REAL(DP), DIMENSION(N_x) :: psi_0, zeta_0
  !REAL(DP), DIMENSION(:) :: zeta_0 

  !!display_model_config()

  alpha = dt/(2.0*dx)
  dxdx  = dx*dx
  t     = 0.0
  
  !! allocate previous and current states
  ALLOCATE (psi_prev(N_x),psi_curr(N_x), zeta_prev(N_x), zeta_curr(N_x), diag(N_x-2), rhs_tridiag(N_x-2))
  
  !set initial and boundary conditions
  DO i = 1, N_x
     psi_prev(i)  = psi_0(i)
     zeta_prev(i) = zeta_0(i)
     !print*, psi_prev(i)
  END DO

  psi_curr(1)    = bc_0
  psi_curr(N_x)  = bc_N
  zeta_curr(1)   = zeta_prev(1)
  zeta_curr(N_x) = zeta_prev(N_x)

!  CALL write_file(t, psi_prev, N_x) !! Write initial condition to file
  
  DO n = 2, int(upper_time)
   
     ! forward euler for vorticity
     DO i = 2, N_x-1
       zeta_curr(i) = zeta_prev(i) - alpha*(psi_prev(i+1)-psi_prev(i-1))
     END DO

     ! solve poisson equation for new streamfunction
     DO i = 2, N_x-1
        rhs_tridiag(i-1) = -dxdx*zeta_curr(i)   
     END DO

     ! interior solver by poisson
     CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_curr)

     DO i = 2, N_x-1
        psi_prev(i) = psi_curr(i)
        zeta_prev(i) = zeta_curr(i)
        print*, t, psi_curr(i)
     END DO

     t = t+dt
    

    ! IF (MOD(n,20) .eq. 0) THEN
       ! print*, dt
        !CALL write_file(t, psi_curr, N_x)
    ! END IF 

  END DO
     
   DEALLOCATE(psi_prev, psi_curr, zeta_prev, zeta_curr, diag, rhs_tridiag)

END SUBROUTINE basin_euler

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE basin_leapfrog(N_x)
  USE constants
  INTEGER :: n, i
  INTEGER ::  N_x
  REAL(DP) :: alpha, gamm, dxdx, t
  REAL(DP) :: dt, dx
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: diag, rhs_tridiag
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_prev, zeta_pp, zeta_curr
  REAL(DP), DIMENSION(N_x) :: psi_0, zeta_0

  !! display_model_config()

  alpha = dt/(2.0*dx)
  gamm = dt/dx 
  dxdx  = dx*dx
  t     = 0.0

  ALLOCATE (psi_prev(N_x), psi_curr(N_x), zeta_prev(N_x), zeta_pp(N_x), zeta_curr(N_x), diag(N_x-2), rhs_tridiag(N_x-2))
  
  DO i = 1, N_x
     psi_prev(i) = psi_0(i)
     zeta_pp(i) = zeta_0(i)
  END DO

  psi_curr(1)    = bc_0
  psi_curr(N_x)  = bc_N
  zeta_curr(1)   = zeta_pp(1)
  zeta_curr(N_x) = zeta_pp(N_x)

  CALL write_file(t, psi_prev, N_x);   !Write initial condition to file
  
 
  DO n = 2, N_x-1
     zeta_prev(i) = zeta_0(i) - alpha*(psi_0(i+1)-psi_0(i-1))
  END DO

  DO i = 2, N_x-1
     rhs_tridiag(i-1) = -dxdx*zeta_prev(i)
  END DO

  CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_prev)
  
  ! time loop for vorticity using leapfrog
  DO n = 2, int(upper_time)

     ! forward leapfrog for vorticity
     DO i = 2, N_x-1
        zeta_curr(i) = zeta_pp(i) - gamm*(psi_prev(i+1) - psi_prev(i-1))
        rhs_tridiag(i-1) = -dx*dx*zeta_curr(i)
     END DO
    
     CALL tridiagonal(diag, rhs_tridiag, N_x-2, psi_curr)

     DO i = 2, N_x-1
        psi_prev(i)  = psi_curr(i)
        zeta_pp(i)   = zeta_prev(i)
        zeta_prev(i) = zeta_curr(i)
     END DO

     t = t+dt

     !IF (MOD(n,20) .eq. 0) THEN
        CALL write_file(t, psi_curr, N_x)
     !END IF
  END DO

  DEALLOCATE (psi_prev, psi_curr, zeta_prev, zeta_pp, zeta_curr, diag, rhs_tridiag)

END SUBROUTINE basin_leapfrog

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Poisson equation solver

SUBROUTINE tridiagonal(b, y, N, solution)
  USE constants
  INTEGER :: i 
  INTEGER, INTENT(IN) :: N
  REAL(DP), DIMENSION(N) :: b, y
  REAL(DP), DIMENSION(N), INTENT(OUT) :: solution
 
  b(1) = 2.0
  

  ! forward substitution
  DO i = 2, N
     b(i) = (i + 2)/REAL((i + 1),DP)
     y(i) = y(i) + (y(i-1)/b(i-1));
  ENDDO

  ! backward substitution
  solution(N) = y(N)/b(N)
  DO i = N-1, 1, -1
     solution(i) = (y(i)+solution(i+1))/b(i)
  ENDDO

END SUBROUTINE tridiagonal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_file(t, psi, n)
  USE constants
  INTEGER :: i, n
  !INTEGER, INTENT(IN) :: N_x
  REAL(DP), DIMENSION(n), INTENT(IN) :: psi
  REAL(DP), INTENT(IN) :: t

!  write (6,'(6F12.6)') t
  DO i = 1, n
     write(6,'(6F12.6)') t, psi(i)
  ENDDO


END SUBROUTINE write_file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


