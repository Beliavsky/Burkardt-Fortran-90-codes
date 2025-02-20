program main

!*****************************************************************************80
!
!! power_method_test() tests power_method().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER_METHOD_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test the POWER_METHOD library.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POWER_METHOD_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 uses POWER_METHOD on the Fibonacci2 matrix.
!
!  Discussion:
!
!    This matrix, despite having a single dominant eigenvalue, will generally
!    converge only very slowly under the power method.  This has to do with
!    the fact that the matrix has only 3 eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 100

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) cos_x1x2
  real ( kind = rk ) ctime
  real ( kind = rk ) ctime1
  real ( kind = rk ) ctime2
  integer i
  integer it_max
  integer it_num
  real ( kind = rk ) lambda
  real ( kind = rk ) phi
  real ( kind = rk ) sin_x1x2
  real ( kind = rk ) tol
  real ( kind = rk ) x(n)
  real ( kind = rk ) x2(n)

  call fibonacci2 ( n, a )

  call random_number ( harvest = x(1:n) )

  it_max = 500
  tol = 0.00000001D+0

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Use POWER_METHOD on the Fibonacci2 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol

  call cpu_time ( ctime1 )

  call power_method ( n, a, x, it_max, tol, lambda, it_num )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,f14.10)'  ) '  Estimated eigenvalue   = ', lambda
  write ( *, '(a,f14.10)'  ) '  Correct value          = ', phi
  write ( *, '(a,g14.6)'   ) '  ||Error||              = ', abs ( lambda - phi )
!
!  X2 is the exact eigenvector.
!  Computing it "backwards" this way avoids overflow when N is large.
!
  x2(n) = 1.0D+00
  do i = n, 2, -1
    x2(i-1) = x2(i) / phi
  end do

  x2(1:n) = x2(1:n) / sqrt ( sum ( x2(1:n)**2 ) )
!
!  The sine of the angle between X and X2 is a measure of error.
!
  cos_x1x2 = dot_product ( x(1:n), x2(1:n) )
  sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.10)' ) &
    '  Sine of angle between true and estimated vectors = ', sin_x1x2

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 uses POWER_METHOD2 on the Fibonacci2 matrix.
!
!  Discussion:
!
!    This matrix, despite having a single dominant eigenvalue, will generally
!    converge only very slowly under the power method.  This has to do with
!    the fact that the matrix has only 3 eigenvectors.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) cos_x1x2
  real ( kind = rk ) ctime
  real ( kind = rk ) ctime1
  real ( kind = rk ) ctime2
  integer i
  integer it_max
  integer it_num
  complex ( kind = ck ) lambda
  real ( kind = rk ) phi
  real ( kind = rk ) sin_x1x2
  real ( kind = rk ) tol
  complex ( kind = ck ) v(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) x2(n)

  call fibonacci2 ( n, a )

  call random_number ( harvest = x(1:n) )

  it_max = 500
  tol = 0.00000001D+0

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Use POWER_METHOD2 on the Fibonacci2 matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol

  call cpu_time ( ctime1 )

  call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )
  x = real ( v, kind = rk )
  x = x / sqrt ( sum ( x(1:n)**2 ) )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue   = ', lambda
  write ( *, '(a,f14.10)'  ) '  Correct value          = ', phi
  write ( *, '(a,g14.6)'   ) '  ||Error||              = ', abs ( lambda - phi )
!
!  X2 is the exact eigenvector.
!  Computing it "backwards" this way avoids overflow when N is large.
!
  x2(n) = 1.0D+00
  do i = n, 2, -1
    x2(i-1) = x2(i) / phi
  end do

  x2(1:n) = x2(1:n) / sqrt ( sum ( x2(1:n)**2 ) )
!
!  The sine of the angle between X and X2 is a measure of error.
!
  cos_x1x2 = dot_product ( x(1:n), x2(1:n) )
  sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) )

  write ( *, '(a)' ) ' '
  write ( *, '(a,f14.10)' ) &
    '  Sine of angle between true and estimated vectors = ', sin_x1x2

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 uses POWER_METHOD2 on the TRIS matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 July 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 100

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) beta
  real ( kind = rk ) ctime
  real ( kind = rk ) ctime1
  real ( kind = rk ) ctime2
  real ( kind = rk ) gamma
  integer i
  integer it_max
  integer it_num
  complex ( kind = ck ) lambda
  complex ( kind = ck ) lambda_max
  complex ( kind = ck ) lambda_vec(n)
  real ( kind = rk ) tol
  complex ( kind = ck ) v(n)
  real ( kind = rk ) x(n)
!
!  If ALPHA * GAMMA is negative, this matrix will have complex eigenvalues.
!
  alpha = -1.0D+00
  beta = 10.0D+00
  gamma = +8.0D+0

  call tris ( n, n, alpha, beta, gamma, a )

  call random_number ( harvest = x(1:n) )

  it_max = 4000
  tol = 0.00000001D+0

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Use POWER_METHOD2 on the TRIS (tridiagonal scalar) matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Matrix order N         = ', n
  write ( *, '(a,i8)'      ) '  Maximum iterations     = ', it_max
  write ( *, '(a,g14.6)'   ) '  Error tolerance        = ', tol
!
!  Estimate the solution.
!
  call cpu_time ( ctime1 )

  call power_method2 ( n, a, x, it_max, tol, lambda, v, it_num )

  x = real ( v, kind = rk )
  x = x / sqrt ( sum ( x(1:n)**2 ) )

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime1

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)'      ) '  Number of iterations   = ', it_num
  write ( *, '(a,g14.6)'   ) '  CPU time               = ', ctime
  write ( *, '(a,2f14.10)' ) '  Estimated eigenvalue   = ', lambda
!
!  Get the exact solution.
!
  call tris_eigenvalues ( n, alpha, beta, gamma, lambda_vec )

  lambda_max = lambda_vec(1)
  do i = 2, n
    if ( abs ( lambda_max ) < abs ( lambda_vec(i) ) ) then
      lambda_max = lambda_vec(i)
    end if
  end do

  write ( *, '(a,2f14.10)' ) '  Correct max eigenvalue = ', lambda_max
  write ( *, '(a,g14.6)'   ) '  ||Error||              = ', abs ( lambda - lambda_max )

  return
end
