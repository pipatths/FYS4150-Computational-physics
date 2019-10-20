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

DOUBLE PRECISION FUNCTION int_function(x1, y1, z1, x2, y2, z2)
  REAL(DP), INTENT(IN) :: x1, y1, z1, x2, y2, z2 
  REAL(DP) :: alf, ZERO
  REAL(DP) :: exp1, exp2, deno
    alf = 2.0
    ZERO = 1E-6
    exp1 = -2*alf*sqrt(x1**2+y1**2+z1**2)
    exp2 = -2*alf*sqrt(x2**2+y2**2+z2**2)
    deno = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    IF (deno < ZERO) THEN
       int_function = 0.
    ELSE

       int_function = EXP(exp1+exp2)/deno
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


END MODULE
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
PROGRAM integral
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: grid
  INTEGER :: h, i, j, k, l, m, c
  REAL(DP) :: a, b, int_gauss, alf, ZERO
  INTEGER :: begin_time, end_time
  REAL(DP) :: elapsed_time, exact
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x, w

 ! WRITE(*,*) 'Read number of integration points'
 ! READ(*,*) n
 ! WRITE(*,*) 'Read lower limit of integration'
 ! READ(*,*) a
 ! WRITE(*,*) 'Read upper limit of integration'
 ! READ(*,*) b
  a = -2.
  b = 2.
  ALLOCATE(grid(5))
  grid(1) = 10
  grid(2) = 15
  grid(3) = 20
  grid(4) = 25
  grid(5) = 30

  OPEN(UNIT=5,FILE='gauss_legedre_mpi.dat') 
  DO c = 1, 5
  ALLOCATE(x(grid(c)), w(grid(c)))
  CALL system_clock(begin_time)
  CALL gauleg (a, b, x, w, grid(c))

 int_gauss = 0.
  DO h = 1, grid(c)
     DO i = 1, grid(c)
        DO j = 1, grid(c)
           DO k = 1, grid(c)
              DO l = 1, grid(c)
                 DO m = 1, grid(c)
                    int_gauss = int_gauss+w(h)*w(i)*w(j)*w(k)*w(l)*w(m) * int_function(x(h),x(i),x(j),x(k),x(l),x(m))
           
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
        grid(c), int_gauss, exact, int_gauss-exact, (int_gauss-exact)/exact, elapsed_time
    
  DEALLOCATE(x, w) 
  ENDDO
  DEALLOCATE(grid)
END PROGRAM integral

