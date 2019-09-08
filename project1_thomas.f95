MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
END MODULE constants

MODULE functions
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
       x(i) = x(1)+((i)*h)
       u_exact(i) = 1-x(i)+(x(i)*EXP(-10.))-EXP(-10*x(i))
       f(i) = 100.*EXP(-10.*x(i))  
    ENDDO  
    DO i=2, n-1
       computed_derivative(i) = -(u_exact(i-1)+u_exact(i+1)-(2*u_exact(i)))/(h*h)
    ENDDO
    DO i = 1, n
       exact_derivative(i) = -100*EXP(-10*x(i))
    ENDDO

  END SUBROUTINE derivative
!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------   
  SUBROUTINE matrixa (n, i, j, Amat, g, h, f)
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
  SUBROUTINE thomas(Amat,u_thomas, n, g, g_new_thomas)
!============================================================
! Solutions to a system of tridiagonal linear equations C*x=b
! Method: the Thomas method
! Alex G. November 2009
!-----------------------------------------------------------
! input ...
! c(n,3) - array of coefficients for matrix C
! b(n)   - vector of the right hand coefficients b
! n      - number of equations
! output ...
! x(n)   - solutions
! comments ...
! the original arrays c(n,3) and b(n) will be destroyed 
! during the calculation
!===========================================================
    IMPLICIT NONE
    INTEGER :: n, i, j, k
    REAL(DP), DIMENSION(n), INTENT(INOUT) :: g
    REAL(DP), DIMENSION(n), INTENT(OUT) :: u_thomas, g_new_thomas 
    REAL(DP), DIMENSION(n,n), INTENT(INOUT) :: Amat
    REAL(DP) :: coeff, sum_a
 
  !step 1: forward elimination
    DO k = 1, n-1
       DO i = k+1, n
          coeff=Amat(i,k)/Amat(k,k)
          DO j = k+1, n
             Amat(i,j)=Amat(i,j)-coeff*Amat(k,j)
          ENDDO
          g_new_thomas(i)=g(i)-coeff*g(k)      
       ENDDO
    ENDDO

  !step 2: back substitution
    u_thomas(n) = g_new_thomas(n)/Amat(n,n)
    DO i= n-1, 1, -1
       sum_a = 0.0
       DO j= i+1, 1, n
          sum_a = sum_a + Amat(i,j)*u_thomas(j)
       ENDDO
       u_thomas(i) = (g_new_thomas(i)-sum_a)/Amat(i,i)
    ENDDO
  END SUBROUTINE thomas
!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------
  SUBROUTINE relative_error(u_thomas_error, u_exact, i, u_thomas, n)
    INTEGER  :: i
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n), INTENT(OUT):: u_thomas_error
    REAL(DP), DIMENSION(n), INTENT(IN) :: u_thomas, u_exact

    u_thomas_error(1) = 0.
    u_thomas_error(n) = 0.
    DO i= 2, n-1
       u_thomas_error(i) = log10(abs((u_thomas(i)-u_exact(i))/u_exact(i)))
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
  REAL(DP), ALLOCATABLE, DIMENSION(:)  :: computed_derivative, x, u_exact, &
                                          f, g, exact_derivative, u_thomas, &
                                          g_new_thomas, u_thomas_error
  REAL(DP) :: h
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Amat
  !  read in input data from screen 
  !WRITE(*,*) 'Read number of steps'
  !READ(*,*) n
  n = 10
  ! open file to write results on
  OPEN(UNIT=7,FILE='thomas10.dat')

  !  allocate space in memory for the one-dimensional arrays                       
  ALLOCATE(computed_derivative(n), x(n), u_exact(n), f(n), &
           g(n), exact_derivative(n), u_thomas(n), &
           g_new_thomas(n), u_thomas_error(n))
  ALLOCATE(Amat(n,n))             
  ! compute the second derivative
  computed_derivative = 0.0_dp; x = 0.0_dp; u_exact = 0.0_dp; &
  f = 0.0_dp; g = 0.0_dp; h = 0.0_dp; u_thomas = 0.0_dp; Amat = 0.0_dp; 
  g_new_thomas  = 0.0_dp; u_thomas_error = 0.0_dp;

  CALL derivative(n, x, computed_derivative, u_exact, h, f, i, &
             exact_derivative)
  CALL matrixa(n, i, j, Amat, g, h, f)
  CALL thomas(Amat,u_thomas, n, g, g_new_thomas)
  CALL relative_error(u_thomas_error, u_exact, i, u_thomas, n)

  !  Then we print the results to file  
  WRITE(7, '(t1,a,t19,a,t37,a,t55,a)') [character(20) :: 'x', 'u_exact', &
        'u-Thomas', 'u-Thomas error'] 
  
  DO i=1,  n
     WRITE(7,'(E16.10,2X,E16.10,2X,E16.10,2X,E16.10,2X,E16.10)') x(i), u_exact(i), &
           u_thomas(i), u_thomas_error(i)
  ENDDO
    
  DEALLOCATE (computed_derivative, x, u_exact, f, g, exact_derivative, Amat, &
              u_thomas, g_new_thomas, u_thomas_error)

END PROGRAM second_derivative


