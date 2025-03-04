program main

!*****************************************************************************80
!
!! quad_openmp() estimates an integral, using OpenMP for parallel execution.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) error
  real ( kind = rk ) exact
  real ( kind = rk ), external :: f
  integer i
  integer n
  real ( kind = rk ) total
  real ( kind = rk ) wtime
  real ( kind = rk ) x

  a =  0.0D+00
  b = 10.0D+00
  n = 10000000
  exact = 0.49936338107645674464D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'quad_openmp():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use OpenMP for parallel execution.'
  write ( *, '(a)' ) '  Estimate the integral of f(x) from A to B.'
  write ( *, '(a)' ) '  f(x) = 50 / ( pi * ( 2500 * x * x + 1 ) ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  A        = ', a
  write ( *, '(a,g14.6)' ) '  B        = ', b
  write ( *, '(a,i8)' ) '  N        = ', n
  write ( *, '(a,g14.6)' ) '  Exact    = ', exact

  wtime = omp_get_wtime ( )

  total = 0.0D+00
!$omp parallel shared ( a, b, n ) private ( i, x ) 
!$omp do reduction ( + : total )
  do i = 1, n
    x = ( real ( n - i,     kind = rk ) * a &
        + real (     i - 1, kind = rk ) * b ) &
        / real ( n     - 1, kind = rk )
    total = total + f ( x )
  end do
!$omp end do
!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime

  total = ( b - a ) * total / real ( n, kind = rk )
  error = abs ( total - exact )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Estimate = ', total
  write ( *, '(a,g14.6)' ) '  Error    = ', error
  write ( *, '(a,g14.6)' ) '  Time     = ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'quad_openmp():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function f ( x )

!*****************************************************************************80
!
!! f() evaluates the function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
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
!    Input, real ( kind = rk ) X, the argument.
!
!    Output, real ( kind = rk ) F, the function value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f
  real ( kind = rk ) pi
  real ( kind = rk ) x

  pi = 3.141592653589793D+00
  f = 50.0D+00 / ( pi * ( 2500.0D+00 * x * x + 1.0D+00 ) )

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
