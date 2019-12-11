SUBROUTINE tridiag_general(a, b, c, y, N, solution)
  USE constants
  
  INTEGER :: n, i
  REAL(DP) :: a, b, c, y, solution
  
  !Forward substitution
  DO i = 1, n
     b(i) = b(i) - a(i-1)*c(i-1)/b(i-1)    !Eliminating lower diagonal
     y(i) = y(i) - (a(i-1)/b(i-1))*y(i-1)  
  ENDDO
  
  !Backward substitution
  solution(n-1) = y(n-1)/b(n-1)     ! obtain final element of solution

  DO i = n-2, 1, -1
     solution(i) =(y(i) -c(i)*solution(i+1))/b(i)
  ENDDO
 
END SUBROUTINE tridiag_general


SUBROUTINE tridiag_ferrari(b, y, n, solution)
  INTEGER :: i, n
  REAL(DP) :: b, y, solution
  b(1) = 2.0

  ! forward substitution
  DO i = 2,n 
     b(i) = (i+2)/(REAL(i+1))
     y(i) = y(i) + (y(i-1)/b(i-1))
  ENDDO

  !backward substitution
  solution(n-1) = y(n-1/b(n-1)

  DO i = n-2, 0, -1
     solution(i) = (y(i) + solution(i+1))/b(i)

  ENDDO
                                                                        

END SUBROUTINE tridiag_ferrari


SUBROUTINE poision_jacobi(g, bc_0y, bc_1y, bc_x0, bc_x1, dx, dy, N_x, N_y, max_iter, f)
  INTEGER :: i, j, iter
  INTEGER, INTENT(IN) :: N_x, N_y, max_iter
  REAL(DP) :: dxdx, dydy, dxdxdydy, dxdx_plus_dydy, eps, diff
  REAL(DP), INTENT(IN) :: g, bc_0y, bc_1y, bc_x0, yyyyyyyyyy, dx, dy, f
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: f_tmp

  dxdx = dx*dx
  dydy = dy*dy
  dxdxdydy = dxdx*dydy
  dxdx_plus_dydy = 2*(dxdx+dydy)

  DO j = 1, N_y
     f(1+j) = bc_0y(j)
     f((N_x-1)*N_y + j) = bc_1y(j)
  ENDDO

  DO i = 1, N_x
     f(i*N_y + 0) = bc_x0(i)
     f(i*N_y + (N_y-1)) = bc_x1(i)
  ENDDO

   iter = 0
   diff = 1.0E+20
   eps  = 1.0E-6
   ALLOCATE (f_tmp, N_x*N_y)

   DO WHILE (iter <= max_iter .and. abs(diff) > eps) 
      diff = 0.0
      
      DO i = 1, N_x
         DO j = 1, N_y
            f_tmp(i*N_y+j) = f(i*N_y+j)
         ENDDO
      ENDDO

      DO i = 2, N_x-1
         DO j = 2, N_y-1
            f(i*N_y + j) = (dydy*(f_tmp((i+1)*N_y + j) + f_tmp((i-1)*N_y + j)) + dxdx*(f_tmp(i*N_y + (j+1)) + f_tmp(i*N_y + (j-1))) - dxdxdydy*g(i*N_y + j))/dxdx_plus_dydy;
                diff += f(i*N_y + j) - f_tmp(i*N_y + j);
         ENDDO
      ENDDO

      iter = iter+1
   END DO
   DEALLOCATE(f_tmp)


END SUBROUTINE poisson_jacobi


SUBROUTINE  poisson_jacobi_periodic (g, dx, dy, N_x, N_y, max_iter, f)
  INTEGER :: i, j, iter
  INTEGER, INTENT(IN) :: N_x, N_y, max_iter
  REAL(DP) :: dxdx, dydy, dxdxdydy, dxdx_plus_dydy, eps, diff
  REAL(DP), INTENT(IN) :: g, bc_0y, bc_1y, bc_x0, yyyyyyyyyy, dx, dy, f
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: f_tmp

  dxdx = dx*dx
  dydy = dy*dy
  dxdxdydy = dxdx*dydy
  dxdx_plus_dydy = 2*(dxdx+dydy)

  inter = 0
  diff = 1.0E+20
  eps  = 1.0E-6
  
  ALLOCATE (f_tmp(N_x*N_y))

  DO WHILE (iter<= max_it .and. abs(diff) > eps)
     diff = 0.0
    
     DO i = 1, N_x
        DO j = 1, N_y
           f_tmp(i*N_y+j) = f(i*N_y+j)
        ENDDO
     ENDDO
   
     DO i = 1, N_x
         DO j = 1, N_y

            ! at x = 0, y = 0
            IF (i == 0 .and. j == 0) THEN
               f(0*N_y + 0) = (dydy*(f_tmp(1*N_y + 0) + f_tmp((N_x-2)*N_y + 0)) + dxdx*(f_tmp(0*N_y + 1) + f_tmp(0*N_y + N_y-2)) - dxdxdydy*g(0*N_y + 0))/dxdx_plus_dydy

            ! at x = 1, y = 0
            ELSE IF (i == N_x-1 .and. j ==0) THEN
               f((N_x-1)*N_y + 0) = (dydy*(f_tmp(1*N_y + 0) + f_tmp((N_x-2)*N_y + 0)) + dxdx*(f_tmp((N_x-1)*N_y + 1) + f_tmp((N_x-1)*N_y + (N_y-2))) - dxdxdydy*g*(N_x-1)*N_y + 0))/dxdx_plus_dydy
            
            ! at x = 0, y = 1
            ELSE IF (i == 0 .and. j == N_y-1) THEN
               f(0*N_y + (N_y-1)) = (dydy*(f_tmp(1*N_y + (N_y-1)) + f_tmp((N_x-2)*N_y + (N_y-1))) + dxdx*(f_tmp(0*N_y + 1) + f_tmp(0*N_y + (N_y-2))) - dxdxdydy*g(0*N_y + (N_y-1)))/dxdx_plus_dydy
 
            ! at x = 1, y = 1
            ELSE IF (i == N_x-1 .and. j == N_y-1) THEN
               f((N_x-1)*N_y + (N_y-1)) = (dydy*(f_tmp(1*N_y + (N_y-1)) + f_tmp((N_x-2)*N_y + (N_y-1))) + dxdx*(f_tmp((N_x-1)*N_y + 1) + f_tmp((N_x-1)*N_y + N_y-2)) - dxdxdydy*g((N_x-1)*N_y + (N_y-1)))/dxdx_plus_dydy

            ! at x = 0
            ELSE IF (i == 0 .and. j .ne. 0 .and. j .ne. N_y-1) THEN
               f(0*N_y + j) = (dydy*(f_tmp(1*N_y + j) + f_tmp((N_x-2)*N_y + j)) + dxdx*(f_tmp(0*N_y + (j+1)) + f_tmp(0*N_y + (j-1))) - dxdxdydy*g(0*N_y + j))/dxdx_plus_dydy

            ! at x = 1
            ELSE IF (i == N_x-1 .and. j .ne. 0 .and. j .ne. N_y-1) THEN
               f((N_x-1)*N_y + j) = (dydy*(f_tmp(1*N_y + j) + f_tmp((N_x-2)*N_y + j)) + dxdx*(f_tmp((N_x-1)*N_y + (j+1)) + f_tmp((N_x-1)*N_y + (j-1))) - dxdxdydy*g((N_x-1)*N_y + j))/dxdx_plus_dydy

	    ! at y = 0
            ELSE IF (j == 0 && i .ne. 0 && i .ne. N_x-1) THEN
               f(i*N_y + 0) = (dydy*(f_tmp((i+1)*N_y + 0) + f_tmp((i-1)*N_y + 0)) + dxdx*(f_tmp(i*N_y + 1) + f_tmp(i*N_y + (N_y-2))) - dxdxdydy*g(i*N_y + 0))/dxdx_plus_dydy

            ! at y = 1
            ELSE IF (j == N_y-1 && i .ne. 0 && i .ne. N_x-1) THEN
               f(i*N_y + (N_y-1)) = (dydy*(f_tmp((i+1)*N_y + (N_y-1)) + f_tmp((i-1)*N_y + N_y-1)) + dxdx*(f_tmp(i*N_y + 1) + f_tmp(i*N_y + (N_y-2))) - dxdxdydy*g(i*N_y + (N_y-1)))/dxdx_plus_dydy
 
            ! interior
            ELSE
               f(i*N_y + j) = (dydy*(f_tmp[(i+1)*N_y + j) + f_tmp((i-1)*N_y + j)) + dxdx*(f_tmp(i*N_y + (j+1)) + f_tmp(i*N_y + (j-1))) - dxdxdydy*g(i*N_y + j))/dxdx_plus_dydy

            ENDIF
            
            diff = diff+f(i*N_y+j) - f_tmp(i*N_y+j)
              
         ENDDO
      ENDDO

      iter = iter+1



  ENDDO
  DEALLOCATE (f_tmp)

END SUBROUTINE poisson_jacobi_periodic
