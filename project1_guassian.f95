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
       exact_derivative(i) = 100.*EXP(-10.*x(i))
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
  SUBROUTINE gaussian_elimination(a,b,c,g,u_gauss,n,b_new, g_new, h, i)
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n), INTENT(IN) :: g
    REAL(DP), DIMENSION(n) :: a, b, c, b_new, g_new
    REAL(DP), DIMENSION(n), INTENT(OUT) :: u_gauss
    REAL(DP), INTENT(IN) :: h
    INTEGER :: i, j     

    u_gauss(1) = 0.
    u_gauss(n) = 0.

    DO i = 1, n
       a(i) = -1
       b(i) = 2
       c(i) = -1
    ENDDO
! forward substituion
    g_new(1) = g(1)
    b_new(1) = b(1)
    DO i = 2, n
       b_new(i) = b(i)-(c(i-1)*a(i-1)/b(i-1))
       g_new(i) = g(i)-(g(i-1)*a(i-1)/b(i-1))
    ENDDO

! backward substitution
    u_gauss(n) = g_new(n)/b_new(n)
    DO i = n-1, 2, -1
       u_gauss(i) = g_new(i)-(u_gauss(i+1)*c(i)/b_new(i))
    ENDDO  

  END SUBROUTINE gaussian_elimination

!----------------------------------------------------------------------------  
!----------------------------------------------------------------------------
  SUBROUTINE relative_error(derivative_error, u_gauss_error, &
                            computed_derivative, exact_derivative, u_gauss, &
                            u_exact, i, n)
    INTEGER  :: i
    INTEGER, INTENT(IN) :: n
    REAL(DP), DIMENSION(n) :: derivative_error, u_gauss_error, u_thomas_error
    REAL(DP), DIMENSION(n), INTENT(IN) :: computed_derivative, exact_derivative, &
                                          u_gauss, u_exact

    derivative_error(1) = 0.0
    derivative_error(n) = 0.0                                      
    DO i= 2, n-1
       derivative_error(i) = log10(abs((computed_derivative(i)-exact_derivative(i))/exact_derivative(i)))
    ENDDO

    u_gauss_error(1) = 0.0
    u_gauss_error(n) = 0.0 
    DO i= 2, n-1
       derivative_error(i) = log10(ABS((computed_derivative(i)-exact_derivative(i))/exact_derivative(i)))
       u_gauss_error(i) = log10(ABS((u_gauss(i)-u_exact(i))/u_exact(i)))
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
                                          f, a, b, c, g, u_gauss, b_new, &
                                          g_new, exact_derivative, &
                                          derivative_error, u_gauss_error
  REAL(DP) :: h
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Amat
  !  read in input data from screen 
  !WRITE(*,*) 'Read number of steps'
  !READ(*,*) n
  n = 1000
  ! open file to write results on

  OPEN(UNIT=7,FILE='gaussian1000.dat')

  !  allocate space in memory for the one-dimensional arrays                       
  ALLOCATE(computed_derivative(n), x(n), u_exact(n), f(n), &
           a(n), b(n), c(n), g(n), u_gauss(n), b_new(n), g_new(n), &
           exact_derivative(n), derivative_error(n), u_gauss_error(n))
  ALLOCATE(Amat(n,n))             
  ! compute the second derivative
  computed_derivative = 0.0_dp; x = 0.0_dp; u_exact = 0.0_dp; &
  f = 0.0_dp; a = 0.0_dp; b = 0.0_dp; c = 0.0_dp; &
  g = 0.0_dp; u_gauss = 0.0_dp; b_new = 0.0_dp; g_new = 0.0_dp; &
  h = 0.0_dp; Amat = 0.0_dp; 

  CALL derivative(n, x, computed_derivative, u_exact, h, f, i, &
                  exact_derivative)
  CALL matrixa(n, i, j, Amat, g, h, f)
  CALL gaussian_elimination(a,b,c,g,u_gauss,n,b_new, g_new, h, i)
  CALL relative_error(derivative_error, u_gauss_error, &
                      computed_derivative, exact_derivative, u_gauss, &
                      u_exact, i, n)

  !  Then we print the results to file  
  WRITE(7, '(t1,a,t19,a,t37,a,t55,a,t73,a,t91,a,t109,a)') &
        [character(20) :: 'x', 'u-gaussian', 'u-exact', 'Gaussian error', &
        'Computed deriv', 'exact deriv', 'deriv error'] 
  
  DO i=1,  n
    
     WRITE(7,'(E16.10,2X,E16.10,2X,E16.10,2X,E16.10,2X,E16.10,2X,E16.10,2X,E16.10)')&
           x(i), u_gauss(i), u_exact(i), u_gauss_error(i), computed_derivative(i), &
           exact_derivative(i), derivative_error(i)
     
  ENDDO
    
  DEALLOCATE(computed_derivative, x, u_exact, f, a, b, c, g, u_gauss, &
             b_new, g_new, exact_derivative, Amat, derivative_error, u_gauss_error) 
END PROGRAM second_derivative


