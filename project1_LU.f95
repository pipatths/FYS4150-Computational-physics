MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0)
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
END MODULE constants
USE constants
IMPLICIT NONE
CONTAINS
!----------------------------------------------------------------------------  
!---------------------------------------------------------------------------- 
  SUBROUTINE derivative(n, x, computed_derivative, u_exact, h, f, i, &
             exact_derivative)
    USE constants
    INTEGER, INTENT(IN) :: n
    INTEGER  :: i
    REAL(DP), DIMENSION(n), INTENT(INOUT) :: computed_derivative, x, &
                                             u_exact, f, exact_derivative
    REAL(DP), INTENT(OUT) :: h
    x(1) = 0
    x(n) = 1     
    u_exact(1) = 0
    u_exact(n) = 0
    h = (x(n)-x(1))/(n+1)
    f(1) = 100.*EXP(-10.*x(1))
    f(n) = 100.*EXP(-10.*x(n))
    computed_derivative(1) = (100.)*EXP(-10.*x(1))
    computed_derivative(n) = (100.)*EXP(-10.*x(n))
    DO i=2, n-1
       x(i) = x(1)+((i-1)*h)
       u_exact(i) = 1.-x(i)+(x(i)*EXP(-10.))-EXP(-10.*x(i))
       f(i) = 100.*EXP(-10.*x(i))  
    ENDDO  
    DO i=2, n-1
       computed_derivative(i) = -(u_exact(i-1)+u_exact(i+1)-(2*u_exact(i)))/(h*h)
    ENDDO
    DO i = 1, n
       exact_derivative(i) = -100.*EXP(-10.*x(i))
    ENDDO

  END SUBROUTINE derivative
!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------   
  SUBROUTINE matrixa(n, i, j, Amat, g, h, f)
    INTEGER, INTENT(IN)  :: n
    INTEGER  :: i, j
    REAL(DP) :: h
    REAL(DP), DIMENSION(n) :: f
    REAL(DP), DIMENSION(n), INTENT(OUT) :: g
    REAL(DP), DIMENSION(n,n), INTENT(OUT) :: Amat 
    DO j=1, n
       DO i=1, n
          g(i) = h*h*f(i)
          IF (i==j) THEN
             Amat(i,j) = 2
          ELSE IF (i==j+1) THEN
             Amat(i,j) = -1
          ELSE IF (i==j-1) THEN
             Amat(i,j) = -1
          ELSE 
             Amat(i,j) = 0 
          END IF
       ENDDO
    ENDDO
    
  END SUBROUTINE matrixa
!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------
  SUBROUTINE lu_decomposition (n, Amat, g, u_lu)
!============================================================
! Inverse matrix
! Method: Based LU decomposition for Ax=b
! Old Dominion University
!---------------------------------------------------------- 
   INTEGER, INTENT(IN) :: n
   INTEGER :: i, j, K
   REAL(DP), DIMENSION(n,n) :: Amat 
   REAL(DP), DIMENSION(n,n) :: L, U, c
   REAL(DP), DIMENSION(n) :: d, g
   REAL(DP), DIMENSION(n), INTENT(INOUT) :: u_lu
   REAL(DP) :: coeff

   L=0.0
   U=0.0
   d=0.0

! Forward elimination
   DO k=1, n-1
      DO i=k+1,n 
         coeff=Amat(i,k)/Amat(k,k) 
         L(i,k) = coeff 
         DO j=k+1,n 
            Amat(i,j) = Amat(i,j)-coeff*Amat(k,j) 
         ENDDO
      ENDDO
   ENDDO
! Create L and U matrices
   DO i=1,n
      L(i,i) = 1.0
   ENDDO

   DO j=1,n
      DO i=1,j 
         U(i,j) = Amat(i,j)
      ENDDO
   ENDDO
! Compute inverse matrix c

      d(1) = g(1)
      DO i=2,n 
         d(i)=g(i) 
         DO j=1,i-1 
            d(i) = d(i) - L(i,j)*d(j)
         ENDDO
      ENDDO
! Back substitution
      u_lu(n)= g(n)/U(n,n)
      DO i = n-1,1,-1 
         u_lu(i) = g(i) 
         DO j=n,i+1,-1
            u_lu(i)=u_lu(i)-U(i,j)*u_lu(j) 
         ENDDO
         u_lu(i) = u_lu(i)/U(i,i)
      ENDDO
      u_lu(1) = 0.
   
  END SUBROUTINE lu_decomposition
!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------
  SUBROUTINE relative_error(i, n, u_lu_error, u_lu, u_exact)
    INTEGER  :: i
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n), INTENT(OUT) :: u_lu_error
    REAL(DP), DIMENSION(n), INTENT(IN) :: u_lu, u_exact
  
    u_lu_error(1) = 0.
    u_lu_error(n) = 0.
    DO i= 2, n-1
      u_lu_error(i) = log10(abs((u_lu(i)-u_exact(i))/u_exact(i)))
     ENDDO
  END SUBROUTINE

END MODULE functions

!----------------------------------------------------------------------------  
!---------------------------------------------------------------------------- 
PROGRAM second_derivative 
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 
  INTEGER :: n, i, j
  REAL(DP) :: h
  REAL(DP), ALLOCATABLE, DIMENSION(:)  :: computed_derivative, x, u_exact, &
                                          f, g, exact_derivative, &
                                          u_lu_error, u_lu
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Amat
  !  read in input data from screen 
  !WRITE(*,*) 'Read number of steps'
  !READ(*,*) n
  n = 10
  ! open file to write results on
  OPEN(UNIT=7,FILE='lu10.dat')

  !  allocate space in memory for the one-dimensional arrays                       
  ALLOCATE(computed_derivative(n), x(n), u_exact(n), f(n), g(n), &
           exact_derivative(n), u_lu_error(n), u_lu(n))
  ALLOCATE(Amat(n,n))             
  ! compute the second derivative
  computed_derivative = 0.0_dp; x = 0.0_dp; u_exact = 0.0_dp; &
  f = 0.0_dp; g = 0.0_dp; h = 0.0_dp; Amat = 0.0_dp; &
  exact_derivative = 0.0_dp; u_lu_error = 0.0_dp; u_lu = 0.0_dp;

  CALL derivative(n, x, computed_derivative, u_exact, h, f, i, &
                  exact_derivative)
  CALL matrixa(n, i, j, Amat, g, h, f)
  CALL lu_decomposition (n, Amat, g, u_lu)
  CALL relative_error(i, n, u_lu_error, u_lu, u_exact)

  !  Then we print the results to file  
  WRITE(7, '(t1,a,t19,a,t37,a,t55,a)') [character(20) :: 'x', &
        'u_exact', 'u-LU', 'LU error'] 
  
  DO i=1,  n
    
     WRITE(7,'(E16.10,2X,E16.10,2X,E16.10,2X,E16.10)') x(i), &
           u_exact(i), u_lu(i), u_lu_error(i)
     
  ENDDO
    
  DEALLOCATE(computed_derivative, x, u_exact, f, g, &
             exact_derivative, u_lu_error, u_lu)

END PROGRAM second_derivative


