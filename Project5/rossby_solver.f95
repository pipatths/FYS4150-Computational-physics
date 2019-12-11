
SUBROUTINE rossby_solver(dx, N_x, dt, T, fileout)
  INTEGER :: N_x
  REAL :: dx, dt, T
  REAL, ALLOCATABLE, DIMENSION(:) :: psi_0, zeta_0
  system_dim = "1D"
  dx = dx

  ALLOCATE (psi_0(N_x), zeta(N_x))
 
  this->outfile.open(fileout);
    this->outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);

END SUBROUTINE rossby_solver

SUBROUTINE rossby_solver_2d(dx, dy, N_x, N_y, dt, T, fileout)
  INTEGER :: N_x, N_y
  REAL :: dx, dy, dt, T
  REAL, ALLOCATABLE, DIMENSION(:) :: psi_0, zeta_0

  ALLOCATE (psi_0(N_x*N_y), zeta(N_x*N_y))
  
      this->outfile.open(fileout);
    this->outfile << std::setiosflags(std::ios::showpoint | std::ios::uppercase);
END SUBROUTINE rossby_solver_2d

SUBROUTINE intial_condition(psi_holder, zeta_holder)
  INTEGER :: i, j 
  REAL(DP) ::psi_holder, zeta_holder
  IF (system_dim.compare"1D" == 0) THEN
     DO i = 1, N_x
        psi_holder(i) = psi_0(i)
        zeta_holder(i) = zeta_0(i)
     ENDDO
  ELSE IF ( system_dim.compare("2d") ==0) THEN
     DO i = 1, N_x
        DO j = 0, N_y
           psi_holder(i*N_y+j) = psi_0(i*N_y+j)
           zeta_holder(i*N_y+j) = zeta_0(i*N_y+j)
        ENDDO
     ENDDO
  END IF
END SUBROUTINE intial_condition


SUBROUTINE write_file(t, psi)
  this->outfile << std::setw(15) << std::setprecision(8) << t;
  IF (system_dim.compare("1D") == 0) THEN
     DO i = 1, N_x
        this->outfile << std::setw(15) << std::setprecision(8) << psi[i];
     ENDDO
  ELSE IF (system_dim.compare("2D") ==0) THEN
     DO i = 1, N_x 
        DO j = 1, N_y
           this->outfile << std::setw(15) << std::setprecision(8) << psi[i*this->N_y + j];
        ENDDO
     ENDDO
  ENDIF
  this->outfile << std::endl;

END SUBROUTINE write_file



SUBROUTINE display_model
  print*, "Simulation parameters:"
  print*, "dx = " dx
  
  IF (system_dim.compare("2D") ==0) THEN
     print*, "dy = " dy
  END IF   
  print*, "N_x = " N_x

  IF (system_dim.compare("2D") == 0) THEN
     print*, "N_y = " N_y
  END IF
  print*, "dt = " dt
  print*, "T = " T
  print*, "Output in " fileout
END SUBROUTINE display_model
