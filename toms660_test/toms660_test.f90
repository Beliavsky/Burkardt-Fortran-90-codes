program main

!*****************************************************************************80
!
!! toms660_test() tests toms660().
!
!  Discussion:
!
!    A quadratic function is sampled on a 6x6 regular grid of [0,1]x[0,1].
!
!    The interpolant is formed from this data.
!
!    The interpolant is evaluated on a 10x10 regular grid of [0,1]x[0,1]
!    and compared to the original function.
!
!    Since the function and interpolant are quadratic, they should be
!    numerically equal for this test.
!
!  Modified:
!
!    26 January 2012
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 36
  integer, parameter :: nq = 13
  integer, parameter :: nr = 3
  integer, parameter :: nw = 19

  real ( kind = rk ) a(5,n)
  real ( kind = rk ) dx
  real ( kind = rk ) dy
  real ( kind = rk ) eps
  real ( kind = rk ) eq
  real ( kind = rk ) eqx
  real ( kind = rk ) eqy
  real ( kind = rk ) f(n)
  real ( kind = rk ) fq
  real ( kind = rk ) fx
  real ( kind = rk ) fy
  integer i
  integer ier
  integer j
  integer k
  integer lcell(3,3)
  integer lnext(n)
  real ( kind = rk ) p(10)
  real ( kind = rk ) px
  real ( kind = rk ) py
  real ( kind = rk ) q
  real ( kind = rk ) q1
  real ( kind = rk ) qs2val
  real ( kind = rk ) qx
  real ( kind = rk ) qy
  real ( kind = rk ) rmax
  real ( kind = rk ) rq
  real ( kind = rk ) rsq(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) xmin
  real ( kind = rk ) xx
  real ( kind = rk ) y(n)
  real ( kind = rk ) ymin
  real ( kind = rk ) yy
!
!  Quadratic test function and partial derivatives.
!
  fq(xx,yy) = (         ( xx + 2.0D+00 * yy ) / 3.0E+00 )**2
  fx(xx,yy) = 2.0D+00 * ( xx + 2.0E+00 * yy ) / 9.0E+00
  fy(xx,yy) = 4.0D+00 * ( xx + 2.0E+00 * yy ) / 9.0E+00

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms660_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS660().'
!
!  Generate a 6 by 6 grid of nodes in the unit square with the natural ordering.
!
  k = 0
  do j = 5, 0, -1
    do i = 0, 5
      k = k + 1
      x(k) = real ( i, kind = rk ) / 5.0D+00
      y(k) = real ( j, kind = rk ) / 5.0D+00
    end do
  end do
!
!  Compute the data values.
!
  do k = 1, n
    f(k) = fq ( x(k), y(k) )
  end do
!
!  Call QSHEP2 to define the interpolant Q to this data.
!
  call qshep2 ( n, x, y, f, nq, nw, nr, lcell, lnext, xmin, ymin, &
    dx, dy, rmax, rsq, a, ier )

  if ( ier /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'TOMS660_PRB - Error!'
    write ( *, '(a,i8)' ) '  Error in TOMS660, IER = ', ier
    stop
  end if
!
!  Generate a 10 by 10 uniform grid of interpolation points
!  (p(i),p(j)) in the unit square.
!
  do i = 1, 10
    p(i) = real ( i - 1, kind = rk ) / 9.0D+00
  end do
!
!  Compute the machine precision EPS.
!
  eps = epsilon ( eps )
!
!  Compute the interpolation errors.
!
  eq = 0.0D+00
  eqx = 0.0D+00
  eqy = 0.0D+00

  do j = 1, 10

    py = p(j)

    do i = 1, 10

      px = p(i)

      q1 = qs2val ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
        ymin, dx, dy, rmax, rsq, a )

      call qs2grd ( px, py, n, x, y, f, nr, lcell, lnext, xmin, &
        ymin, dx, dy, rmax, rsq, a, q, qx, qy, ier )

      if ( ier /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS660_PRB - Error!'
        write ( *, '(a,i6)' ) '  Error in QS2GRD, IER = ', ier
        stop
      end if

      if ( 3.0D+00 * abs ( q ) * eps < abs ( q1 - q ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS660_PRB - Error!'
        write ( *, '(a)' ) '  The interpolated values Q1 (QS2VAL)'
        write ( *, '(a)' ) '  and Q (QS2GRD) differ.'
        write ( *, '(a,g14.6)' ) '  Q1 = ', q1
        write ( *, '(a,g14.6)' ) '  Q  = ', q
        stop
      end if

      eq = max ( eq, abs ( fq ( px, py ) - q ) )
      eqx = max ( eqx, abs ( fx ( px, py ) - qx ) )
      eqy = max ( eqy, abs ( fy ( px, py ) - qy ) )

    end do

  end do
!
!  Print the maximum errors and the ratio EQ / EPS.
!
  rq = eq / eps

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Maximum absolute errors in the interpolant Q and'
  write ( *, '(a)' ) '  partial derivatives QX and QY relative to machine'
  write ( *, '(a)' ) '  precision EPS.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Function   Max error   Max error/EPS'
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a4,7x,e9.3,5x,f4.2)' ) 'Q   ', eq, rq
  write ( *, '(2x,a4,7x,e9.3)'      ) 'dQdX', eqx
  write ( *, '(2x,a4,7x,e9.3)'      ) 'dQdY', eqy
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS660_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
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
!  Modified:
!
!    06 August 2005
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
