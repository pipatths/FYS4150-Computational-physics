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
  SUBROUTINE matrixa (n, h_step, rho, non_diag, diag ,Amat)
    INTEGER :: n
    INTEGER :: i, j
    REAL(DP) :: h_step, non_diag, diag
    REAL(DP), DIMENSION(n) :: rho
    REAL(DP), DIMENSION(n,n), INTENT(OUT) :: Amat 

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
SUBROUTINE jacobi (n, R, eps, iter, Amat, max_off)
  INTEGER :: i, j, n, k, l, itmax
  INTEGER, INTENT(OUT) :: iter
  REAL(DP), DIMENSION(n,n), INTENT(INOUT) :: Amat
  REAL(DP), DIMENSION(n,n), INTENT(OUT) :: R
  REAL(DP) :: eps, max_off, tau, t, c, s  
!  REAL(DP), INTENT(IN) :: max_off

  iter = 0
  eps = 1.0e-08
  itmax = 1000000
  DO i = 1, n
     DO j = 1, n
        R(i,j) = 0.0
        R(i,i) = 1.0
     ENDDO
  ENDDO


  DO WHILE (max_off.GT.eps .AND. iter.LT.itmax)
     iter = iter+1
10   CALL maxoffdiag (Amat, k, l, n, max_off)
     CALL rotate(Amat, R, k, l, n,tau, t, c, s, *10)
  ENDDO
  
END SUBROUTINE jacobi
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

SUBROUTINE maxoffdiag (Amat, k, l, n, max_off)
  INTEGER :: n
  INTEGER :: i, j
  INTEGER, INTENT(OUT) :: k, l
  REAL(DP), DIMENSION(n,n), INTENT(INOUT) :: Amat
  REAL(DP), INTENT(OUT):: max_off
  max_off = 0.0
  DO i = 1, n
     DO j = i+1, n
        IF ((Amat(i,j)**2).GT.max_off) THEN
           max_off = Amat(i,j)**2
           l = i
           k = j
        END IF
     ENDDO
  ENDDO

END SUBROUTINE maxoffdiag
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

SUBROUTINE rotate(Amat, R, k, l, n,tau, t, c, s, *)
  INTEGER :: k, l, n, i, j
  REAL(DP), DIMENSION(n,n), INTENT(INOUT) :: Amat
  REAL(DP), DIMENSION(n,n), INTENT(OUT) :: R
  REAL(DP) :: s, c, t, tau, a_kk, a_ll, a_ik, a_il, r_ik, r_il
   
  IF (Amat(k,l).NE.0.0) THEN
     tau = (Amat(l,l)-Amat(k,k))/(2*Amat(k,l))
     IF (tau.GT.0) THEN
        t = 1.0/(tau+SQRT(1.0+tau*tau))
     ELSE
        t = -1.0/(-tau+SQRT(1.0+tau*tau))
     END IF
     c = 1/SQRT(1+t*t)
     s = c*t
  ELSE
     c = 1.0
     s = 0.0
  END IF
  a_kk = Amat(k,k)
  a_ll = Amat(l,l)

  Amat(k,k) = c*c*a_kk-2.0*c*s*Amat(k,l)+s*s*a_ll
  Amat(l,l) = s*s*a_kk+2.0*c*s*Amat(k,l)+c*c*a_ll
  Amat(k,l) = 0.0
  Amat(l,K) = 0.0

  DO i = 1, n-1
     IF (i.NE.k .AND. i.NE.l) THEN
        a_ik = Amat(i,k)
        a_il = Amat(i,l)
        Amat(i,k) = c*a_ik-s*a_il
        Amat(k,i) = Amat(i,k)
        Amat(i,l) = c*a_il+s*a_ik
        Amat(l,i) = Amat(i,l)
     END IF
     r_ik = R(i,k)
     r_il = R(i,l)
     R(i,k) = c*r_ik-s*r_il
     R(i,l) = c*r_il-s*r_ik
  ENDDO
  RETURN 10
END SUBROUTINE rotate
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

SUBROUTINE analytical(n, exact, diag, non_diag, x, &
	              Amat, error)
  INTEGER :: i, j, n
  REAL(DP), DIMENSION(n), INTENT(OUT) :: exact, x, error
  REAL(DP), INTENT(IN) :: diag, non_diag
  REAL(DP), DIMENSION(n,n), INTENT(IN) :: Amat
 
  DO i = 1,n
     exact(i) = diag+2*non_diag*cos(i*pi/(n+1)) 
     x(i) = Amat(i,i)
     error(i) = x(i)-exact(i)
  ENDDO


END SUBROUTINE analytical
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 

SUBROUTINE sort(n,x)
  INTEGER :: i, j, n
  REAL(DP), DIMENSION(n), INTENT(INOUT) :: x
;  REAL(DP), DIMENSION(n):: eigval
  REAL(DP) :: a
 
  DO j = 2, n
     a = x(j)
     DO i = j-1, 1, -1
        IF (x(i).LE.a) go to 20
            x(i+1) = x(i)
     ENDDO
     i = 0
20   x(i+1) = a
  ENDDO   

END SUBROUTINE sort


END MODULE
!---------------------------------------------------------------------------- 
!---------------------------------------------------------------------------- 
PROGRAM eigen_problem
  USE constants
  USE functions
  IMPLICIT NONE
  ! declarations of variables 

  INTEGER :: n, i, j, k, l, iter, itmax
!  INTEGER, POINTER :: i, j
!  INTEGER, TARGET :: k, L
  REAL(DP) :: h_step, non_diag, diag, eps, max_off, &
              s, c, t, tau, a_kk, a_ll, a_ik, a_il, &
              r_ik, r_il, y

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: rho, x, exact, eig_val, error
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: Amat, R

  !  read in input data from screen 
  WRITE(*,*) 'Read number of steps'
  READ(*,*) n
  ! open file to write results on

  !OPEN(UNIT=7,FILE='eigenvalues.dat')
  OPEN(UNIT=5,FILE='rho.dat')
  OPEN(UNIT=3,FILE='eigenvector.dat')
  OPEN(UNIT=7,FILE='amat.dat')


  !  allocate space in memory for the one-dimensional arrays                       
  ALLOCATE(Amat(n,n), R(n,n), rho(n), x(n), exact(n), &
           eig_val(n), error(n))             
  ! compute the second derivative
  h_step = 0.0_dp; non_diag = 0.0_dp; diag = 0.0_dp; &
  eps = 0.0_dp; max_off = 0.0_dp; s = 0.0_dp; c = 0.0_dp; &
  t = 0.0_dp; tau = 0.0_dp; a_kk = 0.0_dp; a_ll = 0.0_dp; &
  a_ik = 0.0_dp; a_il = 0.0_dp; r_ik = 0.0_dp; r_il = 0.0_dp; &
  x = 0.0_dp; exact = 0.0_dp; error = 0.0_dp;

  CALL matrixa(n, h_step, rho, non_diag, diag, Amat)
  CALL maxoffdiag(Amat, k, l, n, max_off)
  CALL jacobi(n, R, eps, iter, Amat, max_off)
  CALL analytical(n, exact, diag, non_diag, x, Amat, error)
  CALL sort(n,x)

  print *, h_step, non_diag, diag
  print *, iter, max_off
  


  DO i=1, n
     WRITE(5,'(E16.10, 2X, E16.10, 2X, E16.10, 2X, E16.10, 2X, E16.10)') rho(i), x(i), &
           exact(i), x(i)-exact(i), (x(i)-exact(i))/exact(i)
     WRITE(3,'(100g10.3)') (R(i,j), j=1,n)
     WRITE(7,'(100g10.3)') (Amat(i,j), j=1,n)
  ENDDO
    
  DEALLOCATE(Amat, R, rho, x, exact, error) 
END PROGRAM eigen_problem

