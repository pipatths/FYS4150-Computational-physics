SUBROUTINE basin_2d (dx, dy, N_x, N_y, dt, T, fileout)
  REAL(DP) :: dx, dy, dt, T
  INTEGER :: N_x, N_y
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: bc_0y, bc_1y, bc_x0, bc_x1

  CALL rossby_solver(dx,dy, N_x, Ny, dt, T, fileout)
  ALLOCATE(bc_0y(N_y), bc_1y(N_y), bc_x0(N_x), bc_x1(N_x))
  
  DO i = 1, N_x
     bc_x0(i) = 0.0 !No flow at southern boundary
     bc_x1(i) = 0.0 !NO flow at Northern boundary
  ENDDO
  DO j = 1, N_y
     bc_0y(j) = 0.0 !NO flow at western boundary
     bc_1y(j) = 0.0 !NO flow at eastern boundary
  ENDDO

END SUBROUTINE basin_2d

SUBROUTINE intial_condition(ini_psi, ini_zeta)
  INTEGER :: i, j, N_y, N_x
  REAL(DP), INTENT(IN) :: bc_x0, bc_x1, bc_0y, bc_1y
  REAL(DP) :: int_psi, int_zeta, psi_0, zeta_0, eps

  eps = 1.0E-10
  DO i = 1, N_x
     IF (abs(ini_psi(i*N_y) - bc_x0(i)) > eps .and. abs(ini_psi(i*N_y+N_y-1)-bc_x1(i))>eps) THEN
        print*, ini_psi(1) " , " ini_psi(N_x-1)
        print*, "Error: Chosen initial condition does not satisfy BC!" 
        print*,  "Terminating program..." 
        exit
     ENDIF
  ENDDO

  DO j = 1, N_y)
     IF (abs(ini_psi(0+j) - bc_0y(j)) > eps .and. abs(ini_psi((N_x-1)*N_y+j)-bc_1y(j)) > eps) THEN
        print*, ini_psi(1) " , " ini_psi(N_x-1)
        print*, "Error: Chosen initial condition does not satisfy BC!" 
        print*,  "Terminating program..." 
        exit
     ENDIF
 
  ENDDO

  DO i = 1, N_x
     DO j = 1, N_y
        psi_0(i*N_y+j) = ini_psi(i*N_y+j)
        zeta_0(i*N_y+j) = ini_zeta(i*N_y+j)
     ENDDO 
  ENDDO
END SUBROUTINE intial_condition


SUBROUTINE basin_leapfrog
this->display_model_config();
  INTEGER :: N_x, N_y, i, j, n
  REAL(DP) :: alpha, gamm, t
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_pp, zeta_prev, zeta_curr
  REAL(DP), INTENT(IN) :: dt, dx

  alpha = dt/(2.0*dx)
  gamm  = dt/dx
  t = 0.0

  ALLOCATE (psi_prev(N_x*N_y), psi_curr(N_x*N_y), zeta_prev(N_x*N_y), zeta_pp(N_x*N_y), zeta_curr(N_x*N_y))

! Initial condition and boundary conditions

  DO i = 1, N_x
     DO j = 1, N_y
        psi_prev(i*N_y+j) = psi_0(i*N_y+j)
        zeta_pp(i*N_y+j)  = zeta_0(i*N_y+j)

        psi_curr(0+j) = bc_0y(j)
        psi_curr((N_x-1)*N_y+j) = bc_1y(j)
        zeta_curr(0+j) = zeta_pp(0+j)
        zeta_curr((N_x-1)*N_y+j) = zeta_pp((N_x-1)*N_y+j)
     ENDDO
     
     psi_curr(i*N_y+0) = bc_x0(i)
     psi_curr(i*N_y+N_y-1) = bc_x1(i)
     zeta_curr(i*N_y+0) = zeta_pp(i*N_y+0)
     zeta_curr(i*N_y+N_y-1) = zeta_pp(i*N_y+N_y-1) 
  ENDDO
 this->write_state_to_file(0.0, psi_prev);   // Write initial condition to file

!Initial euler 

  DO i = 2, N_x-1
     DO j = 2, N_y-1
        zeta_prev(i*N_y+j) = zeta_0(i*N_y+j) - alpha*(psi_0((i+1)*N_y+1) - psi_0((i-1)*N_y+j))
     ENDDO    
  ENDDO

  CALL poisson_jacobi(zeta_prev, bc_0y, bc_1y, bc_x0, bc_x1, dx, dy, N_x, N_y, psi_prev)

   this->write_state_to_file(0.0, psi_prev);

! loop over time using leapfrog for vorticity

  DO n = 3, t<T
     ! forward time for vorticity
     DO i = 2, N_x-1
        DO j = 2, N_y-1
           zeta_curr(i*N_y+j) = zeta_pp(i*N_y+j) - gamm*(psi_prev((i+1)*N_y+j) - psi_prev((i-1)*N_y+j))
        ENDDO
     ENDDO
  ! Solve 2D poisson and update streamfunction
     CALL poisson_jacobi(zeta_curr, bc_0y, bc_1y, bc_x0, bc_x1, dx, dy, N_x, N_y, 50, psi_curr);

     DO i = 2, N_x-1
        DO j = 2, N_y-1
           psi_prev(i*N_y+j)  = psi_curr(i*N_y+j)
           zeta_pp(i*N_y+j)   = zeta_prev(i*N_y+j) 
           zeta_prev(i*N_y+j) = zeta_curr(i*N_y+j)
        ENDDO
     ENDDO
     t = t+dt
     IF (n%50 == 0) THEN
        this->write_state_to_file(t, psi_curr);
     ENDIF

  ENDDO

  DEALLOCATE(psi_prev, psi_curr, zeta_pp, zeta_prev, zeta_curr)
 
! DEALLOCATE(bc_0y, bc_1y, bc_x0, bc_x1)
  

END SUBROUTINE basin_leapfrog
