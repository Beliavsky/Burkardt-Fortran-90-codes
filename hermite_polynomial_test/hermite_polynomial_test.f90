program main

!*****************************************************************************80
!
!! hermite_polynomial_test() tests hermite_polynomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  integer e
  integer p

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hermite_polynomial_test():'
  write ( *, '(a)' ) '  Fortran90 version.'
  write ( *, '(a)' ) '  Test hermite_polynomial().'

  call hermite_polynomial_test01 ( )
  call hermite_polynomial_test02 ( )
  call hermite_polynomial_test03 ( )
  call hermite_polynomial_test04 ( )
  call hermite_polynomial_test05 ( )
  call hermite_polynomial_test06 ( )
  call hermite_polynomial_test07 ( )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test08 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test08 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test09 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test09 ( p, e )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test10 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test10 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test11 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test11 ( p, e )

  p = 5
  b = 0.0D+00
  call hermite_polynomial_test12 ( p, b )

  p = 5
  b = 1.0D+00
  call hermite_polynomial_test12 ( p, b )

  p = 5
  e = 0
  call hermite_polynomial_test13 ( p, e )

  p = 5
  e = 1
  call hermite_polynomial_test13 ( p, e )

  call hermite_polynomial_test14 ( )

  call hermite_polynomial_test15 ( )

  call h_to_monomial_matrix_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hermite_polynomial_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine h_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! h_to_monomial_matrix_test() tests h_to_monomial_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 6

  real ( kind = rk8 ) h(n)
  real ( kind = rk8 ) H2M(n,n)
  integer i
  real ( kind = rk8) IDENT(n,n )
  real ( kind = rk8 ) M2H(n,n)
  real ( kind = rk8 ) m(n)

  write ( *,'(a)' ) ''
  write ( *,'(a)' ) 'h_to_monomial_matrix_test():'
  write ( *,'(a)' ) '  h_to_monomial_matrix() returns a matrix which'
  write ( *,'(a)' ) '  converts a physicist''s Hermite polynomial to monomial form.'
  write ( *,'(a)' ) ''

  call h_to_monomial_matrix ( n, H2M )
  call monomial_to_h_matrix ( n, M2H )
  IDENT = matmul ( H2M, M2H )
  call r8mat_print ( n, n, IDENT, '  Product H2M * M2H:' )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Monomial versions of first Hermite polynomials:'
  do i = 0, n - 1
    h(1:n) = 0.0 
    h(i+1) = 1.0
    m = matmul ( H2M, h )
    call r8poly_print ( i, m, '' )
  end do

  return
end
subroutine hermite_polynomial_test01 ( )

!*****************************************************************************80
!
!! hermite_polynomial_test01() tests h_polynomial_value().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk ) e
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ), allocatable :: fx2_vec(:)
  integer n
  real ( kind = rk ) x
  real ( kind = rk ) xvec(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hermite_polynomial_test01():'
  write ( *, '(a)' ) '  h_polynomial_values() stores values of'
  write ( *, '(a)' ) '  the physicist''s Hermite polynomials.'
  write ( *, '(a)' ) '  H_POLYNOMIAL_VALUE evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X           H(N,X)                    H(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call h_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    xvec(1) = x
    allocate ( fx2_vec(n+1) )
    call h_polynomial_value ( 1, n, xvec, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test02 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST02 tests HE_POLYNOMIAL_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk ) e
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ), allocatable :: fx2_vec(:)
  integer n
  real ( kind = rk ) x
  real ( kind = rk ) xvec(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hermite_polynomial_test02():'
  write ( *, '(a)' ) '  HE_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the probabilist''s Hermite polynomials.'
  write ( *, '(a)' ) '  HE_POLYNOMIAL_VALUE evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X          He(N,X)' // &
    '                   He(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call he_polynomial_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    xvec(1) = x
    allocate ( fx2_vec(1:n+1) )
    call he_polynomial_value ( 1, n, xvec, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )
    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test03 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST03 tests HF_FUNCTION_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk ) e
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ), allocatable :: fx2_vec(:)
  integer n
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST03:'
  write ( *, '(a)' ) '  HF_FUNCTION_VALUES stores values of'
  write ( *, '(a)' ) '  the Hermite function Hf(n,x).'
  write ( *, '(a)' ) '  HF_FUNCTION_VALUE evaluates the function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                        Tabulated                 Computed'
  write ( *, '(a)' ) '     N        X          Hf(N,X)' // &
    '                   Hf(N,X)                     Error'

  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hf_function_values ( n_data, n, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2_vec(n+1) )
    call hf_function_value ( 1, n, x, fx2_vec )
    fx2 = fx2_vec(n+1)
    deallocate ( fx2_vec )

    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f12.6,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, x, fx1, fx2, e

  end do

  return
end
subroutine hermite_polynomial_test04 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST04 tests H_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  real ( kind = rk ), allocatable :: hz(:,:)
  character ( len = 80 ) title
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST04:'
  write ( *, '(a)' ) '  H_POLYNOMIAL_ZEROS computes the zeros of H(n,x)'
  write ( *, '(a)' ) '  Check by calling H_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call h_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for H(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( hz(degree,0:degree) )
    call h_polynomial_value ( degree, degree, z, hz )
    write ( title, '(a,i1,a)' ) '  Evaluate H(', degree, ',z):'
    call r8vec_print ( degree, hz(1:degree,degree), title )

    deallocate ( hz )
    deallocate ( z )

  end do

  return
end
subroutine hermite_polynomial_test05 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST05 tests HE_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  real ( kind = rk ), allocatable :: hz(:,:)
  character ( len = 80 ) title
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST05:'
  write ( *, '(a)' ) '  HE_POLYNOMIAL_ZEROS computes the zeros of He(n,x)'
  write ( *, '(a)' ) '  Check by calling HE_POLYNOMIAL there.'

  do degree = 1, 5

    allocate ( z(1:degree) )
    call he_polynomial_zeros ( degree, z )
    write ( title, '(a,i1,a)' ) '  Computed zeros for He(', degree, ',z):'
    call r8vec_print ( degree, z, title )

    allocate ( hz(degree,0:degree) )
    call he_polynomial_value ( degree, degree, z, hz )
    write ( title, '(a,i1,a)' ) '  Evaluate He(', degree, ',z):'
    call r8vec_print ( degree, hz(1:degree,degree), title )

    deallocate ( hz )
    deallocate ( z )

  end do

  return
end
subroutine hermite_polynomial_test06 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST06 tests H_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  real ( kind = rk ), allocatable :: f(:)
  integer n
  real ( kind = rk ) q
  real ( kind = rk ) q_exact
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST06:'
  write ( *, '(a)' ) '  H_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with H(n,x)'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call h_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -oo < X < +00 ) X^E exp(-X^2) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   E       Q_Estimate      Q_Exact'
  write ( *, '(a)' ) ' '

  allocate ( f(1:n) )

  do e = 0, 2 * n - 1
    if ( e == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**e
    end if
    q = dot_product ( w(1:n), f(1:n) )
    call h_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine hermite_polynomial_test07 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST07 tests HE_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  real ( kind = rk ), allocatable :: f(:)
  integer n
  real ( kind = rk ) q
  real ( kind = rk ) q_exact
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST07:'
  write ( *, '(a)' ) '  HE_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with He(n,x)'

  n = 7
  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call he_quadrature_rule ( n, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral ( -oo < X < +00 ) X^E exp(-0.5*X^2) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   E       Q_Estimate      Q_Exact'
  write ( *, '(a)' ) ' '

  allocate ( f(1:n) )

  do e = 0, 2 * n - 1
    if ( e == 0 ) then
      f(1:n) = 1.0D+00
    else
      f(1:n) = x(1:n)**e
    end if
    q = dot_product ( w(1:n), f(1:n) )
    call he_integral ( e, q_exact )
    write ( *, '(2x,i2,2x,g14.6,2x,g14.6)' ) e, q, q_exact
  end do

  deallocate ( f )
  deallocate ( w )
  deallocate ( x )

  return
end
subroutine hermite_polynomial_test08 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST08 tests HN_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, real ( kind = rk ) B, the coefficient of X in the exponential factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST08'
  write ( *, '(a)' ) '  Compute a normalized physicist''s Hermite exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hn(I,X) Hn(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hn_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test09 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST09 tests HN_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, integer E, the exponent of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST09'
  write ( *, '(a)' ) '  Compute a normalized physicist''s Hermite power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hn(I,X) Hn(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized physicist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hn_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test10 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST10 tests HEN_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, real ( kind = rk ) B, the coefficient of X in the exponential factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST10'
  write ( *, '(a)' ) '  Compute a normalized probabilist''s Hermite exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hen(I,X) Hen(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hen(I,X) = normalized probabilist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hen_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test11 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST11 tests HEN_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, integer E, the exponent of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST11'
  write ( *, '(a)' ) '  Compute a normalized probabilist''s Hermite power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hen(I,X) Hen(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hn(I,X) = normalized probabilist''s Hermite polynomial of degree I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hen_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test12 ( p, b )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST12 tests HF_EXPONENTIAL_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, real ( kind = rk ) B, the coefficient of X in the exponential factor.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST12'
  write ( *, '(a)' ) '  Compute a Hermite function exponential product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) exp(B*X) Hf(I,X) Hf(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hf(I,X) = Hermite function of "degree" I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponential argument coefficient B = ', b

  allocate ( table(0:p,0:p) )

  call hf_exponential_product ( p, b, table )

  call r8mat_print ( p + 1, p + 1, table, '  Exponential product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test13 ( p, e )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST13 tests HF_POWER_PRODUCT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 February 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P, the maximum degree of the polynomial 
!    factors.
!
!    Input, integer E, the exponent of X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e
  integer p
  real ( kind = rk ), allocatable :: table(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST13'
  write ( *, '(a)' ) '  Compute a Hermite function power product table.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Tij = integral ( -oo < X < +oo ) X^E Hf(I,X) Hf(J,X) exp(-X*X) dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  where Hf(I,X) = Hermite function of "degree" I.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Maximum degree P = ', p
  write ( *, '(a,g14.6)' ) '  Exponent of X, E = ', e

  allocate ( table(0:p,0:p) )

  call hf_power_product ( p, e, table )

  call r8mat_print ( p + 1, p + 1, table, '  Power product table:' )

  deallocate ( table )

  return
end
subroutine hermite_polynomial_test14 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST14 tests H_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST14'
  write ( *, '(a)' ) &
    '  H_POLYNOMIAL_COEFFICIENTS determines the physicist''s Hermite polynomial coefficients.'

  call h_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  H(', i, ',x) ='
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( c(i,j) == 0.0D+00 ) then

      else if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x^', j
      end if
    end do
  end do
 
  return
end
subroutine hermite_polynomial_test15 ( )

!*****************************************************************************80
!
!! HERMITE_POLYNOMIAL_TEST15 tests HE_POLYNOMIAL_COEFFICIENTS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLYNOMIAL_TEST15'
  write ( *, '(a)' ) &
    '  HE_POLYNOMIAL_COEFFICIENTS determines the probabilist''s Hermite polynomial coefficients.'

  call he_polynomial_coefficients ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  He(', i, ',x) ='
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( c(i,j) == 0.0D+00 ) then

      else if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x^', j
      end if
    end do
  end do
 
  return
end
