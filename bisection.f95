MODULE constants
  INTEGER,  PARAMETER :: dp = KIND(1.0D0)
  INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0))
  REAL(DP), PARAMETER, PUBLIC ::  truncation=1.0E-10
  REAL(DP), PARAMETER, PUBLIC :: pi=4.D0*DATAN(1.D0)
END MODULE constants

MODULE functions
USE constants
IMPLICIT NONE
CONTAINS
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
  SUBROUTINE matrixa (n, h_step, rho, non_diag, diag ,Amat, poli)
    INTEGER :: n, poli
    INTEGER :: i, j
    REAL(DP) :: h_step, non_diag, diag
    REAL(DP), DIMENSION(n) :: rho
    REAL(DP), DIMENSION(n,n), INTENT(OUT) :: Amat 
    poli = 5
    rho(1) = 0.
    rho(n) = 1.
    h_step = rho(n)/n
    diag = 2.0/(h_step*h_step)
    non_diag = -1.0/(h_step*h_step)
    
    DO i =2, n-1
       rho(i) = rho(1)+((i-1)*h_step)
    ENDDO

!!  Set the tridiagomal matrix
    Amat(1,1) = diag
    Amat(1,2) = non_diag
    DO j=1, n
       DO i=1, n
          IF (i==j) THEN
             Amat(i,j) = diag
          ELSE IF (i==j+1) THEN
             Amat(i,j) = non_diag
          ELSE IF (i==j-1) THEN
             Amat(i,j) = non_diag
          ELSE 
             Amat(i,j) = 0.0 
          END IF
       ENDDO
    ENDDO
    
  END SUBROUTINE matrixa

!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
SUBROUTINE bisection(a, b, c, f, diag, non_diag, pr, poli, Amat)
  REAL(DP) :: eps, a, b, c, f
  REAL(DP), INTENT(IN) :: non_diag, diag
  INTEGER:: poli
  INTEGER :: i
  REAL(DP), DIMENSION(poli), INTENT(INOUT) :: pr
  REAL(DP), DIMENSION(poli,poli), INTENT(INOUT) :: Amat
  eps = 1.0e-08
  a = 0
  b = 3

   DO i = 1, poli
      pr(i) = FindDet(Amat,i)
   ENDDO

15  IF (f(a)*f(b) .LT. 0) THEN
       c=(a+b)/2.0
    ELSE
       WRITE (*,*)"Try with another values of a and b"
       RETURN
    END IF
 
    IF (f(a)*f(c) .LT. 0) THEN
       b=c
    ELSE
       a=c
    END IF

    IF (ABS(b-a) .GT. eps) GO TO 15

    WRITE(*,*)"The root is",c
   
END SUBROUTINE bisection

!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
REAL(DP) FUNCTION f(x, diag, non_diag, pr, poli)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: poli 
  REAL(DP) :: x
  REAL(DP), INTENT(IN) :: diag, non_diag
  REAL(DP), DIMENSION(poli), INTENT(IN) :: pr
  f= (diag-x)*pr(poli-1)-non_diag**2*pr(poli-2)
end function
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
REAL(DP) FUNCTION FindDet(Amat, poli)
    IMPLICIT NONE
    REAL(DP), DIMENSION(poli,poli), INTENT(INOUT) :: Amat
    REAL(DP) :: m, temp
    INTEGER :: i, j, k, l
    INTEGER, INTENT(IN) :: poli
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, poli-1
        IF (Amat(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, poli
                IF (Amat(i,k) .NE. 0) THEN
                    DO j = 1, poli
                        temp = Amat(i,j)
                        Amat(i,j)= Amat(k,j)
                        Amat(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, poli
            m = Amat(j,k)/Amat(k,k)
            DO i = k+1, poli
                Amat(j,i) = Amat(j,i) - m*Amat(k,i)
            END DO
        END DO
    END DO
    
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, poli
        FindDet = FindDet * Amat(i,i)
    END DO
    
END FUNCTION FindDet
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

END MODULE
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
PROGRAM eigen_problem
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 

  INTEGER :: n, i, j, k, l, poli
!  INTEGER, POINTER :: i, j
!  INTEGER, TARGET :: k, L, 
  REAL(DP) :: h_step, non_diag, diag, eps, a, b, c, m, temp, x
  LOGICAL :: DetExists
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rho, pr
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Amat 
!  REAL(DP), ALLOCATABLE, DIMENSION(:,:) ::matrix

  !  read in input data from screen 
  WRITE(*,*) 'Read number of steps'
  READ(*,*) n
  ! open file to write results on

  !OPEN(UNIT=7,FILE='eigenvalues.dat')
  OPEN(UNIT=5,FILE='rho.dat')
  OPEN(UNIT=3,FILE='eigenvector.dat')
  OPEN(UNIT=7,FILE='amat.dat')


  !  allocate space in memory for the one-dimensional arrays                       
  ALLOCATE(Amat(n,n), rho(n), pr(poli))             
  ! compute the second derivative
  h_step = 0.0_dp; non_diag = 0.0_dp; diag = 0.0_dp; &
  eps = 0.0_dp; a = 0.0_dp; b = 0.0_dp; c = 0.0_dp; &
  m  = 0.0_dp; temp = 0.0_dp; rho  = 0.0_dp; &
  x = 0.0_dp; pr = 0.0_dp; Amat = 0.0_dp;

  CALL matrixa(n, h_step, rho, non_diag, diag ,Amat, poli)
!  CALL maxoffdiag(Amat, k, l, n, max_off)
!  CALL jacobi(n, R, eps, iter, Amat, max_off)
!  CALL analytical(n, exact, diag, non_diag, x, Amat, error)
!  CALL sort(n,x)
  CALL bisection(a, b, c, f, diag, non_diag, pr, poli, Amat)

!  print *, h_step, non_diag, diag
!  print *, iter, max_off



!  DO i=1, n
!     WRITE(5,'(E16.10, 2X, E16.10, 2X, E16.10, 2X, E16.10)') rho(i), x(i), &
!           exact(i), x(i)-exact(i)
!     WRITE(3,'(100g10.3)') (R(i,j), j=1,n)
!     WRITE(7,'(100g10.3)') (Amat(i,j), j=1,n)
!  ENDDO
    
  DEALLOCATE(Amat, rho, pr) 
END PROGRAM eigen_problem

