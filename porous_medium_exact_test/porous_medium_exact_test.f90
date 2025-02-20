program main

!*****************************************************************************80
!
!! porous_medium_exact_test() tests porous_medium_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  interface
    subroutine porous_medium_parameters ( c_in, delta_in, m_in, t0_in, &
      tstop_in, c_out, delta_out, m_out, t0_out, tstop_out )
      integer, parameter :: rk8 = kind ( 1.0D+00 )
      real ( kind = rk8 ), optional :: c_in
      real ( kind = rk8 ), optional :: c_out
      real ( kind = rk8 ), optional :: delta_in
      real ( kind = rk8 ), optional :: delta_out
      real ( kind = rk8 ), optional :: m_in
      real ( kind = rk8 ), optional :: m_out
      real ( kind = rk8 ), optional :: t0_in
      real ( kind = rk8 ), optional :: t0_out
      real ( kind = rk8 ), optional :: tstop_in
      real ( kind = rk8 ), optional :: tstop_out
    end subroutine
  end interface

  real ( kind = rk8 ) c
  real ( kind = rk8 ) delta
  real ( kind = rk8 ) m
  real ( kind = rk8 ) t0
  real ( kind = rk8 ) tstop

  call timestamp ( )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'porous_medium_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test porous_medium_exact(),'
  write ( *, '(a)' ) '  an exact solution of the porous medium equation.'
!
!  Report the current parameter values.
!
  call porous_medium_parameters ( c_out = c, delta_out = delta, m_out = m, &
    t0_out = t0, tstop_out = tstop )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  parameters:'
  write ( *, '(a,g14.6)' ) '    c =     ', c
  write ( *, '(a,g14.6)' ) '    delta = ', delta
  write ( *, '(a,g14.6)' ) '    m =     ', m
  write ( *, '(a,g14.6)' ) '    t0 =    ', t0
  write ( *, '(a,g14.6)' ) '    tstop = ', tstop

  call porous_medium_residual_test ( tstop )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'porous_medium_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  return
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
!    real ( kind = rk8 ) A, B, the first and last entries.
!
!  Output:
!
!    real ( kind = rk8 ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  integer i
  real ( kind = rk8 ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk8 ) * a   &
             + real (     i - 1, kind = rk8 ) * b ) &
             / real ( n     - 1, kind = rk8 )
    end do

  end if

  return
end
subroutine porous_medium_residual_test ( tstop )

!*****************************************************************************80
!
!! porous_medium_residual_test() tests porous_medium_residual().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 11
  integer, parameter :: m = 5

  integer i
  integer j
  real ( kind = rk8 ) r
  real ( kind = rk8 ) t(m)
  real ( kind = rk8 ) tstop
  real ( kind = rk8 ) u
  real ( kind = rk8 ) ut
  real ( kind = rk8 ) ux
  real ( kind = rk8 ) uxx
  real ( kind = rk8 ) x(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'porous_medium_residual_test():'
  write ( *, '(a)' ) '  Evaluate solution and residual at selected points (X,T)'
  
  call r8vec_linspace ( n, -1.0D+00, 1.0D+00, x )
  call r8vec_linspace ( m, 0.0D+00, tstop, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      X       T       U(X,T)      Resid(X,T)'
  write ( *, '(a)' ) ''
  do j = 1, m
    do i = 1, n
      call porous_medium_exact ( x(i), t(j), u, ut, ux, uxx )
      call porous_medium_residual ( x(i), t(j), r )
      write ( *, '(2x,f8.4,2x,f8.4,2x,g10.4,2x,g10.4)' ) x(i), t(j), u, r
    end do
    write ( *, '(a)' ) ''
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

