periodic_solver_1d::periodic_solver_1d(double dx, int N_x,
                        double dt, double T, std::string fileout) : rossby_solver(dx, N_x, dt, T, fileout)

SUBROUTINE initial_condition(ini_psi, ini_zeta)
  INTEGER :: i
  REAL(DP) :: psi_0, zeta_0, eps = 1.0E-10
  REAL(DP), INTENT(IN) :: ini_psi, ini_psi

  IF ((abs(init_psi(1) - init_psi(N_x)) > eps) || (abs(init_zeta(1) - init_zeta(N_x)) > eps)) THEN
  print*, "Error: Chosen initial condition does not satisfy periodic BC!"
  print*, "Terminating program..." 
  exit(EXIT_FAILURE);
  END IF

  DO i = 1, N_x
     psi_0(i) = ini_psi(i)
     zeta_0(i) = ini_zeta(i) 
  ENDDO
  
END SUBROUTINE initial_condition

SUBROUTINE periodic_euler(ini_psi, ini_zeta)
  INTEGER :: i, n
  REAL(DP) :: alpha, dxdx
  REAL(DP), INTENT(IN) :: dt, t, T
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_prev, zetz_curr, rhs_poisson
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A
 !! display_model_cofig
  ALLOCATE (psi_prev(N_x), psi_curr(N_x), zeta_prev(N_x), zetz_curr(N_x), rhs_poisson(N_x-1))

  ALLOCATE (A(N_x-1,N_x-1)

  CALL initialize_periodic_matrix(A, N_x-1, N_x-1); 

!! Initial condition
  DO i = 1, N_x
     psi_prev(i) = psi_0(i)
     zeta_prev(i) = zeta_0(i)
  ENDDO
  this->write_state_to_file(0.0, psi_prev.memptr());   !! Write initial condition to file
  !! time loop using leapfrog for voricity
  DO n = 2, t<T
     DO i = 2, N_x-1
        zeta_curr(i) = zeta_prev(i) - alpha*(psi_prev(i+1) - psi_prev(i-1))
     ENDDO
     zeta_curr(0) = zeta_prev(1) - alpha*(psi_prev(2) - psi_prev(N_x-2))
     zeta_curr(N_x-1) = zeta_curr(1)

     DO i = 1, N_x-1
        rhs_poisson(i) = -dxdx*zeta_curr(i)
     ENDDO
     psi_curr = solve(A, rhs_poison)   !! Solve 1D Poisson eq. to get psi at current time-step
     psi_curr(N_x-1) = psi_curr(0)

     DO i = 1, N_x
        psi_prev(i) = psi_curr(i)
        zeta_prev(i) = zeta_curr(i)
     ENDDO
     t = t+dt

     IF (n%20 ==0) THEN
       this->write_state_to_file(t, psi_curr.memptr());
     END IF
  ENDDO
END SUBROUTINE initial_condition


SUBROUTINE periodic_leapfrog(ini_psi, ini_zeta)
  INTEGER :: i, n
  REAL(DP) :: alpha,gamm, dxdx, t
  REAL(DP), INTENT(IN) :: dt, t, T
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: psi_prev, psi_curr, zeta_pp,  zeta_prev, zetz_curr, rhs_poisson
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: A
 !! display_model_cofig

  alpha = dt/(2.0*dx)
  gamm  = dt/dx
  dxdx  = dx*dx
  t     = 0.0
  ALLOCATE (psi_prev(N_x), psi_curr(N_x), zeta_pp(N_x), zeta_prev(N_x), zetz_curr(N_x), rhs_poisson(N_x-1))

  ALLOCATE (A(N_x-1,N_x-1)

  CALL initialize_periodic_matrix(A, N_x-1, N_x-1); 

!! Initial condition
  DO i = 1, N_x
     psi_prev(i) = psi_0(i)
     zeta_prev(i) = zeta_0(i)
  ENDDO

  this->write_state_to_file(0.0, psi_prev.memptr());   !! Write initial condition to file


  !! time loop using leapfrog for voricity
  DO n = 2, N_x-1
     zeta_prev(i) = zeta_0(i) - alpha*(psi_0(i+1) - psi_0(i-1))
  ENDDO
  !!boundary condition
  zeta_prev(1) = this->zeta_0(1) - alpha*(this->psi_0(2) - psi_0(N_x-2));           
  zeta_prev(N_x-1) = zeta_prev(1);    

  DO  i = 1, N_x-1
      rhs_poisson(i) = -dxdx*zeta_prev(i)
  ENDDO

  psi_prev = arma::solve(A, rhs_poisson)
  pri_prev(N_x-1) = psi_prev(1)


! time loop using leapfrog for vorticity

  DO n = 3, t<T
     DO i = 1, N_x-1
        zeta_curr(i) = zeta_pp(i) - gamma*(psi_prev(i+1) - psi_prev(i-1))
     ENDDO

     !! Boundary condition
     zeta_curr(1) = zeta_pp(1) - gamma*(psi_prev(2) - psi_prev(N_x-2))
     zeta_curr(N_x) = zeta_curr[1]; 

     !! solve posson equation
     DO i = 1, N_x-1
        rhs_poisson(i) = -dxdx*zeta_curr(i)
     ENDDO

     psi_curr = arma::solve(A, rhs_poisson)
     psi_curr(N_x-1) = psi_curr(0)

     DO i = 1, N_x
        psi_prev(i) = psi_curr(i)
        zeta_pp(i)  = zeta_prev(i)
        zeta_prve(i) = zeta_curr(i)
     ENDDO
     t = t+dt

     IF (n%20 == 0) THEN
        this->write_state_to_file(t, psi_curr.memptr());
     END IF
  ENDDO
END SUBROUTINE periodic_leapfrog



SUBROUTINE initial_periodic_matrix(arma::mat&A, n_rows, n_cols)
  INTEGER :: n_rows, num_cols, i, j
  REAL(DP), ALLOCATABLE, DIMENSION(:,:), INTENT(IN) :: A
 
  DO i = 1, n_rows
     DO j = 1, n_cols
        IF (i == j) THEN
           A(i,j) = 2.0
        ELSE IF (abs(i-j) == 1 ) THEN
           A(i,j) = -1.0
        END IF
     ENDDO
  ENDDO
  A(0, n_cols-1) = -1.0
  A(n_row-1, 0)  = -1.0
END SUBROUTINE periodic_leapfrog










