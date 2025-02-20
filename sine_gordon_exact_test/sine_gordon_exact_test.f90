program main

!*****************************************************************************80
!
!! sine_gordon_exact_test() tests sine_gordon_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  integer i
  integer j
  integer, parameter :: nx = 6
  integer, parameter :: ny = 6
  real ( kind = rk ) r
  real ( kind = rk ) u
  real ( kind = rk ) uxy
  real ( kind = rk ), allocatable :: x(:)
  real ( kind = rk ), allocatable :: y(:)

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sine_gordon_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test sine_gordon_exact().'

  a = 1.5D+00

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a,g14.6)' ) '    a     = ', a

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Evaluate solution and residual at selected points (X,Y)'
  
  allocate ( x(1:nx) )
  allocate ( y(1:ny) )

  call r8vec_linspace ( nx, 0.0D+00, 1.0D+00, x )
  call r8vec_linspace ( ny, 0.0D+00, 1.0D+00, y )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      X       Y       U(X,Y)      Resid(X,Y)'
  write ( *, '(a)' ) ''
  do j = 1, ny
    do i = 1, nx
      call sine_gordon_exact ( a, x(i), y(j), u, uxy )
      call sine_gordon_residual ( u, uxy, r )
      write ( *, '(f14.6,f14.6,g14.6,g14.6)' ) x(i), y(j), u, r
    end do
    write ( *, '(a)' ) ''
  end do

  deallocate ( x )
  deallocate ( y )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'sine_gordon_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop ( 0 )
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! r8vec_linspace() creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    real ( kind = rk ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

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
