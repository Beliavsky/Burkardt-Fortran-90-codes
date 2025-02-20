subroutine monomial_value ( m, n, e, x, value )

!*****************************************************************************80
!
!! MONOMIAL_VALUE evaluates a monomial.
!
!  Discussion:
!
!    This routine evaluates a monomial of the form
!
!      product ( 1 <= i <= m ) x(i)^e(i)
!
!    where the exponents are nonnegative integers.  Note that
!    if the combination 0^0 is encountered, it should be treated
!    as 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 April 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the spatial dimension.
!
!    Input, integer N, the number of points at which the
!    monomial is to be evaluated.
!
!    Input, integer E(M), the exponents.
!
!    Input, real ( kind = rk ) X(M,N), the point coordinates.
!
!    Output, real ( kind = rk ) VALUE(N), the value of the monomial.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer e(m)
  integer i
  real ( kind = rk ) value(n)
  real ( kind = rk ) x(m,n)

  value(1:n) = 1.0D+00

  do i = 1, m
    if ( 0 /= e(i) ) then
      value(1:n) = value(1:n) * x(i,1:n) ** e(i)
    end if
  end do

  return
end
subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 fills an R8MAT with unit pseudorandom numbers.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 August 2004
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns in
!    the array.
!
!    Input/output, integer SEED, the "seed" value, which
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(M,N), the array of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer seed
  real ( kind = rk ) r(m,n)

  do j = 1, n

    do i = 1, m

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + i4_huge
      end if

      r(i,j) = real ( seed, kind = rk ) * 4.656612875D-10

    end do
  end do

  return
end
function square01_area ( )

!*****************************************************************************80
!
!! SQUARE01_AREA: area of the unit square in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) SQUARE01_AREA, the area.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) square01_area

  square01_area = 1.0D+00

  return
end
subroutine square01_monomial_integral ( e, integral )

!*****************************************************************************80
!
!! SQUARE01_MONOMIAL_INTEGRAL: integral over interior of unit square in 2D.
!
!  Discussion:
!
!    The integration region is 
!
!      0 <= X <= 1, 
!      0 <= Y <= 1.
!
!    The monomial is F(X,Y) = X^E(1) * Y^E(2).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer E(2), the exponents.  
!    Each exponent must be nonnegative.
!
!    Output, real ( kind = rk ) INTEGRAL, the integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2

  integer e(m)
  integer i
  real ( kind = rk ) integral

  if ( any ( e(1:m) < 0 ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SQUARE01_MONOMIAL_INTEGRAL - Fatal error!'
    write ( *, '(a)' ) '  All exponents must be nonnegative.'
    stop 1
  end if

  integral = 1.0D+00
  do i = 1, m
    integral = integral / real ( e(i) + 1, kind = rk )
  end do

  return
end
subroutine square01_sample ( n, seed, x )

!*****************************************************************************80
!
!! SQUARE01_SAMPLE samples the interior of the unit square in 2D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 January 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Russell Cheng,
!    Random Variate Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998, pages 168.
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer N, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) X(2,N), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer n

  integer seed
  real ( kind = rk ) x(m,n)

  call r8mat_uniform_01 ( m, n, seed, x )

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
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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
