program main

!*****************************************************************************80
!
!! jacobi_polynomial_test() tests jacobi_polynomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jacobi_polynomial_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test jacobi_polynomial().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jacobi_polynomial_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests J_POLYNOMIAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) e
  real ( kind = rk ) fx1
  real ( kind = rk ) fx2
  real ( kind = rk ), allocatable :: fx2_vec(:,:)
  integer m
  integer n
  integer n_data
  real ( kind = rk ) x
  real ( kind = rk ) xvec(1)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01:'
  write ( *, '(a)' ) '  J_POLYNOMIAL_VALUES stores values of'
  write ( *, '(a)' ) '  the Jacobi polynomials.'
  write ( *, '(a)' ) '  J_POLYNOMIAL evaluates the polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                                    Tabulated                 Computed'
  write ( *, '(a)' ) '     N     A     B        X           J(N,A,B,X)                    J(N,A,B,X)                     Error'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call j_polynomial_values ( n_data, n, a, b, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( fx2_vec(1,0:n) )

    m = 1
    xvec(1) = x

    call j_polynomial ( m, n, a, b, xvec, fx2_vec )

    fx2 = fx2_vec(1,n)
    e = fx1 - fx2

    write ( *, '(2x,i4,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g8.2)' ) &
      n, a, b, x, fx1, fx2, e

    deallocate ( fx2_vec )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests J_POLYNOMIAL_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 3

  real ( kind = rk ) a
  real ( kind = rk ), dimension ( test_num) :: a_test = (/ &
    0.5D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = rk ) b
  real ( kind = rk ), dimension ( test_num ) :: b_test = (/ &
    0.5D+00, 1.5D+00, 0.5D+00 /)
  integer degree
  real ( kind = rk ), allocatable :: hz(:,:)
  integer test
  character ( len = 80 ) title
  real ( kind = rk ), allocatable :: z(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  J_POLYNOMIAL_ZEROS computes the zeros of J(n,a,b,x);'
  write ( *, '(a)' ) '  Check by calling J_POLYNOMIAL there.'

  do test = 1, 3

    a = a_test(test);
    b = b_test(test);

    do degree = 1, 5

      allocate ( z(1:degree) )
      call j_polynomial_zeros ( degree, a, b, z )
      write ( title, '(a,i1,a,f3.1,a,f3.1,a)' ) 'Zeros for J(', degree, ',', a, ',', b, ')'
      call r8vec_print ( degree, z, title )

      allocate ( hz(degree,degree+1) )
      call j_polynomial ( degree, degree, a, b, z, hz )
      write ( title, '(a,i1,a,f3.1,a,f3.1,a)' ) 'Evaluate J(', degree, ',', a, ',', b, ')'
      call r8vec_print ( degree, hz(1:degree,degree+1), title )

      deallocate ( hz )
      deallocate ( z )

    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests J_QUADRATURE_RULE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer j
  real ( kind = rk ) j_double_product_integral
  real ( kind = rk ), allocatable :: ji(:,:)
  real ( kind = rk ), allocatable :: jj(:,:)
  integer k
  integer n
  real ( kind = rk ) q
  real ( kind = rk ) q_exact
  real ( kind = rk ), allocatable :: w(:)
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  J_QUADRATURE_RULE computes the quadrature rule'
  write ( *, '(a)' ) '  associated with J(n,a,b,x);'

  n = 7
  a = 1.0D+00
  b = 2.5D+00

  allocate ( x(1:n) )
  allocate ( w(1:n) )

  call j_quadrature_rule ( n, a, b, x, w )

  call r8vec2_print ( n, x, w, '      X            W' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Use the quadrature rule to estimate:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J      Q_Estimate         Q_Exact'
  write ( *, '(a)' ) ' '

  do i = 0, 5
    allocate ( ji(n,0:i) )
    call j_polynomial ( n, i, a, b, x, ji )
    do j = i, 5
      allocate ( jj(n,0:j) )
      call j_polynomial ( n, j, a, b, x, jj )
      q = 0.0D+00
      do k = 1, n
        q = q + w(k) * ji(k,i) * jj(k,j)
      end do
      q_exact = j_double_product_integral ( i, j, a, b )
      write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) i, j, q, q_exact
      deallocate ( jj )
    end do
    deallocate ( ji )
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests J_DOUBLE_PRODUCT_INTEGRAL
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  integer j
  real ( kind = rk ) j_double_product_integral
  real ( kind = rk ) q

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04:'
  write ( *, '(a)' ) '  J_DOUBLE_PRODUCT_INTEGRAL returns the value of'
  write ( *, '(a)' ) '  the weighted integral of J(i,a,b,x) * J(j,a,b,x);'

  a = 1.0D+00
  b = 2.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    Q = Integral (-1<x<+1) J(i,a,b,x) J(j,a,b,x) (1-x)^a (1+x)^b dx'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I   J      Q'
  write ( *, '(a)' ) ' '

  do i = 0, 5
    do j = i, 5
      q = j_double_product_integral ( i, j, a, b )
      write ( *, '(2x,i2,2x,i2,2x,g14.6)' ) i, j, q
    end do
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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
!    18 May 2013
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
