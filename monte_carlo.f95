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

DOUBLE PRECISION FUNCTION int_function(x)
  REAL(DP), DIMENSION(6), INTENT(IN) :: x
  REAL(DP) :: alf, ZERO
  REAL(DP) :: exp1, exp2, deno
    alf = 2.0
    ZERO = 1E-6
    exp1 = -2*alf*sqrt(x(1)**2+x(2)**2+x(3)**2)
    exp2 = -2*alf*sqrt(x(4)**2+x(5)**2+x(6)**2)
    deno = sqrt((x(1)-x(4))**2+(x(2)-x(5))**2+(x(3)-x(6))**2)

    IF (deno < ZERO) THEN
       int_function = 0.
    ELSE

       int_function = EXP(exp1+exp2)/deno
    ENDIF   
END FUNCTION int_function
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

 REAL(DP) FUNCTION ran0(idum)
    IMPLICIT NONE
    INTEGER :: idum,ia,im,iq,ir,mask,k
    REAL(DP) :: am
    PARAMETER (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836,mask=123459876)

    idum=ieor(idum,MASK)
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    IF  (idum < 0) idum=idum+IM
    ran0=am*idum
    idum=ieor(idum,MASK)

  END FUNCTION ran0
END MODULE
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
PROGRAM integral
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: grid

  INTEGER :: i, j, d, c
  INTEGER :: idum
  REAL(DP) :: int_MC, int_MC2, temp_MC, variance
  REAL(DP) :: length, jacobi_det, stdev
  INTEGER :: begin_time, end_time
  REAL(DP) :: elapsed_time, exact

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x

  ALLOCATE(grid(5))
  grid(1) = 10000
  grid(2) = 100000
  grid(3) = 1000000
  grid(4) = 10000000
  grid(5) = 100000000

  d = 6
  idum = -1
  length = 3.
  int_MC = 0.
  int_MC2 = 0.
  jacobi_det = (2.*length)**d

  OPEN(UNIT=5,FILE='Brute_force_MC.dat')   
  DO c = 1, 5

  ALLOCATE (x(d))

  CALL system_clock(begin_time)
  DO i = 1, grid(c)
     DO j = 1, 6
        x(j) = -length+2*length*ran0(idum)
     ENDDO
     temp_MC = int_function(x)
     int_MC  = int_MC+temp_MC
     int_MC2 = int_MC2+temp_MC*temp_MC
  ENDDO
  int_MC   = jacobi_det*int_MC/grid(c)
  int_MC2  = jacobi_det*int_MC2/grid(c)
  variance = int_MC2-int_MC*int_MC
  stdev    = jacobi_det*SQRT(ABS(variance)/grid(c))
  
  CALL system_clock(end_time)

  elapsed_time =  (end_time-begin_time)/1000.
  exact = 5*3.141592654d0*3.141592654d0/16/16
  

  WRITE(5,'(I15.2, 2X, F10.6, 2X, F10.6, 2X, F10.6, 2X, F10.6, 2x, F10.3)')&
        grid(c), int_MC, exact, int_MC-exact, (int_MC-exact)/exact, elapsed_time
   
  DEALLOCATE(x) 
  ENDDO
  DEALLOCATE(grid)
END PROGRAM integral

