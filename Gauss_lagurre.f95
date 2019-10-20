!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
  REAL(DP), PARAMETER, PUBLIC :: pi=4.D0*DATAN(1.D0)
END MODULE constants
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
MODULE functions
USE constants
IMPLICIT NONE
CONTAINS

!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

DOUBLE PRECISION FUNCTION int_function(r1, r2, theta1, theta2, phi1, phi2)
  REAL(DP), INTENT(IN) :: r1, r2, theta1, theta2, phi1, phi2 
  REAL(DP) :: alf, ZERO
  REAL(DP) :: exp1, exp2, deno
    alf = 2.0
    ZERO = 1E-12
    exp1 = -(2*alf-1)*(r1+r2)
    deno = (r1**2 + r2**2 - 2*r1*r2*(COS(theta1)*COS(theta2) + SIN(theta1)*SIN(theta2)*(COS(phi1-phi2)))) 

    IF (deno < ZERO) THEN
       int_function = 0.
    ELSE

       int_function = EXP(exp1)*SIN(theta1)*SIN(theta2)/SQRT(deno)
    ENDIF
END FUNCTION int_function

!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

  SUBROUTINE gauleg(x1,x2,x,w,n)
    IMPLICIT NONE
    INTEGER :: i, j, m, n
    REAL(DP) :: eps, x1, x2, x, w 
    DIMENSION :: x(n), w(n) 
    PARAMETER (eps=3.D-14)
    REAL(DP) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    DO i=1,m
       z1=0.
       z=COS(3.141592654d0*(i-.25d0)/(n+.5d0))
       DO WHILE ( ABS(z-z1) > EPS)
          p1=1.
          p2=0.
          DO j=1,n
             p3=p2
             p2=p1
             p1=((2.*j-1.)*z*p2-(j-1.)*p3)/j
          ENDDO
          pp=n*(z*p1-p2)/(z*z-1.)
          z1=z
          z=z-p1/pp
       ENDDO
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.*xl/((1.-z*z)*pp*pp)
       w(n+1-i)=w(i)
    ENDDO
  END SUBROUTINE gauleg
END MODULE functions

!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

PROGRAM integral
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: grid
  INTEGER :: h, i, j, k, l, m, c
  INTEGER :: begin_time, end_time
  REAL(DP) :: elapsed_time, exact
  REAL(DP) :: ZERO
  REAL(DP) :: alf, int_laguerre
  REAL(DP) :: theta_min, theta_max, phi_min, phi_max
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: r1, wr1
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: theta, wtheta, phi, wphi

  ALLOCATE(grid(5))
  grid(1) = 10
  grid(2) = 15
  grid(3) = 20
  grid(4) = 25
  grid(5) = 30

  ZERO = 1.0E-10
  alf   = 2.0
  theta_min = 0.
  theta_max = pi
  phi_min = 0.
  phi_max = 2.*pi

  OPEN(UNIT=5,FILE='gauss_leguerre_mpi.dat')   
  DO c = 1, 5
  ALLOCATE(r1(grid(c)), wr1(grid(c)))
  ALLOCATE(theta(grid(c)), wtheta(grid(c)), phi(grid(c)), wphi(grid(c)))
 
  CALL system_clock(begin_time)
  CALL gauss_laguerre(r1, wr1,grid(c),alf)
  CALL gauleg(theta_min,theta_max,theta, wtheta,grid(c))
  CALL gauleg(phi_min,phi_max,phi,wphi,grid(c))

 int_laguerre = 0.
  DO h = 1, grid(c)
     DO i = 1, grid(c)
        DO j = 1, grid(c)
           DO k = 1, grid(c)
              DO l = 1, grid(c)
                 DO m = 1, grid(c)
                 int_laguerre = int_laguerre+wr1(h)*wr1(i)*wtheta(j)*wtheta(k)*wphi(l)*wphi(m)* &
                                int_function(r1(h),r1(i),theta(j),theta(k),phi(l),phi(m))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  
  CALL system_clock(end_time)
 
  elapsed_time =  (end_time-begin_time)/1000.
  exact = 5*3.141592654d0*3.141592654d0/16/16
  
  WRITE(5,'(I3.2, 2X, F10.6, 2X, F10.6, 2X, F10.6, 2X, F10.6, 2x, F10.3)')&
        grid(c), int_laguerre, exact, int_laguerre-exact,&
        (int_laguerre-exact)/exact, elapsed_time

  DEALLOCATE(r1, wr1, theta, wtheta, phi, wphi) 
  ENDDO
  DEALLOCATE(grid)
END PROGRAM integral


