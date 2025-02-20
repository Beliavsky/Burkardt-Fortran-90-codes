program main

!*****************************************************************************80
!
!! gegenbauer_polynomial_test() tests gegenbauer_polynomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'gegenbauer_polynomial_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test gegenbauer_polynomial().'

  call gegenbauer_alpha_check_test ( )
  call gegenbauer_ek_compute_test ( )
  call gegenbauer_integral_test ( )
  call gegenbauer_polynomial_value_test ( )
  call gegenbauer_ss_compute_test ( )
  call gegenbauer_to_monomial_matrix_test ( )

  call imtqlx_test ( )

  call monomial_to_gegenbauer_matrix_test ( )

  call r8_hyper_2f1_test ( )
  call r8_psi_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'gegenbauer_polynomial_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine gegenbauer_alpha_check_test ( )

!*****************************************************************************80
!
!! gegenbauer_alpha_check_test() compares gegenbauer_alpha_check().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  logical check
  integer n
  real ( kind = rk ) r

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_alpha_check_test():'
  write ( *, '(a)' ) '  gegenbauer_alpha_check() checks that ALPHA is legal.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       ALPHA   Check?'
  write ( *, '(a)' ) ''

  do n = 1, 10

    call random_number ( harvest = r )
    alpha = -5.0D+00 + r * 10.0D+00
    call gegenbauer_alpha_check ( alpha, check )
    write ( *, '(2x,f10.4,7x,l1)' ) alpha, check

  end do

  return
end
subroutine gegenbauer_ek_compute_test ( )

!*****************************************************************************80
!
!! gegenbauer_ek_compute_test() tests gegenbauer_ek_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  integer i
  integer n
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( : ) :: x

  alpha = 0.50D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_ek_compute_test():'
  write ( *, '(a)' ) '  gegenbauer_ek_compute() computes a Gauss-Gegenbauer rule;'
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Using parameter ALPHA = ', alpha
  write ( *, '(a,g14.6,a,g14.6)' ) '  Integration interval is [-1,+1].'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '                 W                         X'
  write ( *, '(a)' ) ''

  do n = 1, 10

    allocate ( w(n) )
    allocate ( x(n) )

    call gegenbauer_ek_compute ( n, alpha, x, w )
 
    write ( *, '(a)' ) ''
    do i = 1, n
      write ( *, '(10x,2x,g24.16,2x,g24.16)' ) w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine gegenbauer_integral_test ( )

!*****************************************************************************80
!
!! gegenbauer_integral_test() tests gegenbauer_integral().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 June 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  integer n
  real ( kind = rk ) value

  alpha = 0.50D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_INTEGRAL_TEST'
  write ( *, '(a)' ) '  GEGENBAUER_INTEGRAL evaluates'
  write ( *, '(a)' ) '  Integral ( -1 < x < +1 ) x^n * (1-x^2)^(alpha-1/2) dx'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         N         Value'
  write ( *, '(a)' ) ''

  do n = 0, 10
 
    call gegenbauer_integral ( n, alpha, value )

    write ( *, '(2x,i8,2x,g24.16)' ) n, value

  end do
 
  return
end
subroutine gegenbauer_polynomial_value_test ( )

!*****************************************************************************80
!
!! gegenbauer_polynomial_value_test() tests gegenbauer_polynomial_value().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  real ( kind = rk ), allocatable :: c(:,:)
  real ( kind = rk ) fx(1)
  integer m
  integer n
  integer n_data
  real ( kind = rk ) x(1)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_POLYNOMIAL_VALUE_TEST:'
  write ( *, '(a)' ) '  GEGENBAUER_POLYNOMIAL_VALUE evaluates the Gegenbauer polynomial.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       M     ALPHA         X           GPV    GEGENBAUER'
  write ( *, '(a)' ) ''

  n = 1

  n_data = 0

  do

    call gegenbauer_polynomial_values ( n_data, m, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(0:m,1:n) )

    call gegenbauer_polynomial_value ( m, n, alpha, x, c )

    write ( *, '(2x,i6,2x,f8.2,2x,f8.2,2x,f12.4,2x,f12.4)' ) m, alpha, x(1), fx, c(m,1)

    deallocate ( c )

  end do

  return
end
subroutine gegenbauer_ss_compute_test ( )

!*****************************************************************************80
!
!! gegenbauer_ss_compute_test() tests gegenbauer_ss_compute().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  integer i
  integer n
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( : ) :: x

  alpha = 0.50D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'GEGENBAUER_SS_COMPUTE_TEST'
  write ( *, '(a)' ) '  GEGENBAUER_SS_COMPUTE computes a Gauss-Gegenbauer rule;'
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Using parameter ALPHA = ', alpha
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '                 W                         X'
  write ( *, '(a)' ) ''

  do n = 1, 10

    allocate ( w(n) )
    allocate ( x(n) )

    call gegenbauer_ss_compute ( n, alpha, x, w )
 
    write ( *, '(a)' ) ''
    do i = 1, n
      write ( *, '(10x,2x,g24.16,2x,g24.16)' ) w(i), x(i)
    end do

    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine gegenbauer_to_monomial_matrix_test ( )

!*****************************************************************************80
!
!! gegenbauer_to_monomial_matrix_test() tests gegenbauer_to_monomial_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) A(n,n)
  real ( kind = rk ) alpha
  real ( kind = rk ) g(n)
  integer i
  real ( kind = rk ) mono(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'gegenbauer_to_monomial_matrix_test():'
  write ( *, '(a)' ) '  gegenbauer_to_monomial_matrix() evaluates the matrix'
  write ( *, '(a)' ) '  which converts Gegenbauer polyjomial coefficients'
  write ( *, '(a)' ) '  to monomial coefficients.'
  
  alpha = 0.5D+00

  call gegenbauer_to_monomial_matrix ( n, alpha, A )

  call r8mat_print ( n, n, A, '  Gegenbauer to Monomial matrix G:' )    

  do i = 0, n - 1
    g(1:n) = 0.0D+00
    g(i+1) = 1.0D+00
    mono = matmul ( A, g )
    write ( *, '(a)' ) ''
    write ( *, '(a,i1)' ) '  Monomial form of Gegenbauer polynomial ', i
    call r8poly_print ( i + 1, mono, '' )
  end do

  return
end
subroutine imtqlx_test ( )

!*****************************************************************************80
!
!! imtqlx_test() tests imtqlx().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 June 2015
!
!  Author:
!
!    John Burkardt.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) angle
  real ( kind = rk ) d(n)
  real ( kind = rk ) e(n)
  integer i
  real ( kind = rk ) lam(n)
  real ( kind = rk ) lam2(n)
  real ( kind = rk ) qtz(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'IMTQLX_TEST'
  write ( *, '(a)' ) '  IMTQLX takes a symmetric tridiagonal matrix A'
  write ( *, '(a)' ) '  and computes its eigenvalues LAM.'
  write ( *, '(a)' ) '  It also accepts a vector Z and computes Q''*Z,'
  write ( *, '(a)' ) '  where Q is the matrix that diagonalizes A.'

  d(1:n) = 2.0D+00
  e(1:n-1) = -1.0D+00
  e(n) = 0.0D+00
  z(1:n) = 1.0D+00
!
!  On input, LAM is D, and QTZ is Z.
!
  lam(1:n) = d(1:n)
  qtz(1:n) = z(1:n)

  call imtqlx ( n, lam, e, qtz )

  call r8vec_print ( n, lam, '  Computed eigenvalues:' )

  do i = 1, n
    angle = real ( i, kind = rk ) * r8_pi / real ( 2 * ( n + 1 ), kind = rk )
    lam2(i) = 4.0D+00 * ( sin ( angle ) ) ** 2
  end do

  call r8vec_print ( n, lam2, '  Exact eigenvalues:' )

  call r8vec_print ( n, z, '  Vector Z:' )
  call r8vec_print ( n, qtz, '  Vector Q''*Z:' )

  return
end
subroutine monomial_to_gegenbauer_matrix_test ( )

!*****************************************************************************80
!
!! monomial_to_gegenbauer_matrix_test() tests monomial_to_gegenbauer_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) alpha
  real ( kind = rk ) G(0:n-1,0:n-1)
  real ( kind = rk ) Ginv(0:n-1,0:n-1)
  real ( kind = rk ) I(0:n-1,0:n-1)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  monomial_to_gegenbauer_matrix_test():'
  write ( *, '(a)' ) '  monomial_to_gegenbauer_matrix() evaluates the matrix'
  write ( *, '(a)' ) '  which converts monomial polynomial coefficients'
  write ( *, '(a)' ) '  to Gegenbauer coefficients.'
  
  alpha = 0.5D+00

  call gegenbauer_to_monomial_matrix ( n, alpha, G )
  call r8mat_print ( n, n, G, '  Gegenbauer to Monomial matrix G:' )    

  call monomial_to_gegenbauer_matrix ( n, alpha, Ginv )
  call r8mat_print ( n, n, Ginv, '  Monomial to Gegenbauer matrix Ginv:' )    

  I = matmul ( G, Ginv )
  call r8mat_print ( n, n, I, '  I = G * Ginv:' )    

  return
end
subroutine r8_hyper_2f1_test ( )

!*****************************************************************************80
!
!! r8_hyper_2f1_test() tests r8_hyper_2f1().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_HYPER_2F1_TEST:'
  write ( *, '(a)' ) '  R8_HYPER_2F1 evaluates the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B       C       X      ', &
  ' 2F1                       2F1                     DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(tabulated)               (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_values ( n_data, a, b, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call r8_hyper_2f1 ( a, b, c, x, fx2 )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, c, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine r8_psi_test ( )

!*****************************************************************************80
!
!! r8_psi_test() tests r8_psi().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) r8_psi
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_PSI_TEST:'
  write ( *, '(a)' ) '  R8_PSI evaluates the Psi function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Psi(X)                     Psi(X)  ' &
  // '               DIFF'
  write ( * , '(a)' ) '               (Tabulated)                (R8_PSI)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_psi ( x )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine r8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! r8mat_print() prints a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 September 2004
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    real ( kind = rk ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! r8mat_print_some prints some of a real matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    integer ILO, JLO, the first row and column to print.
!
!    integer IHI, JHI, the last row and column to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)'
    return
  end if

  do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i8,6x)' ) j
    end do

    write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) ' '

    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi

      do j2 = 1, inc

        j = j2lo - 1 + j2

        write ( ctemp(j2), '(g14.6)' ) a(i,j)

      end do

      write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

    end do

  end do

  return
end
subroutine r8poly_print ( n, a, title )

!*****************************************************************************80
!
!! r8poly_print() prints a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = a(0) + a(1) * x + ... + a(n-1) * x^(n-1) + a(n) * x^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 July 2015
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the dimension of A.
!
!    real ( kind = rk ) A(0:N), the polynomial coefficients.
!    A(1) is the constant term and
!    A(N) is the coefficient of X^(N-1).
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(1:n)
  integer i
  real ( kind = rk ) mag
  character plus_minus
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ' '

  if ( n < 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( a(n) < 0.0D+00 ) then
    plus_minus = '-'
  else
    plus_minus = ' '
  end if

  mag = abs ( a(n) )

  if ( 3 <= n ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' ) &
      plus_minus, mag, n - 1
  else if ( n == 2 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' ) &
      plus_minus, mag
  else if ( n == 1 ) then
    write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
  end if

  do i = n - 1, 1, -1

    if ( a(i) < 0.0D+00 ) then
      plus_minus = '-'
    else
      plus_minus = '+'
    end if

    mag = abs ( a(i) )

    if ( mag /= 0.0D+00 ) then

      if ( 3 <= i ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' ) &
          plus_minus, mag, i - 1
      else if ( i == 2 ) then
        write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' ) plus_minus, mag
      else if ( i == 1 ) then
        write ( *, ' ( ''         '', a1, g14.6 )' ) plus_minus, mag
      end if
    end if

  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
