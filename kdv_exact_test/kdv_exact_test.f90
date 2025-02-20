program main

!*****************************************************************************80
!
!! kdv_exact_test() tests kdv_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'kdv_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test kdv_exact().'

  call kdv_exact_rational_test ( )
  call kdv_exact_sech_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'kdv_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine kdv_exact_rational_test ( )

!*****************************************************************************80
!
!! kdv_exact_rational_test() tests kdv_exact_rational().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 April 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine kdv_parameters ( a_in, v_in, t0_in, tstop_in, &
      a_out, v_out, t0_out, tstop_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: t0_in
      real ( kind = rk ), optional :: t0_out
      real ( kind = rk ), optional :: tstop_in
      real ( kind = rk ), optional :: tstop_out
      real ( kind = rk ), optional :: v_in
      real ( kind = rk ), optional :: v_out
    end subroutine
  end interface

  integer, parameter :: nt = 6
  integer, parameter :: nx = 6

  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) t(nt)
  real ( kind = rk ) t0
  real ( kind = rk ) tstop
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uxxx
  real ( kind = rk ) x(nx)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'kdv_exact_rational_test():'
  write ( *, '(a)' ) '  Test kdv_exact_rational().'

  call kdv_parameters ( t0_out = t0, tstop_out = tstop )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a,g14.6)' ) '    t0    = ', t0
  write ( *, '(a,g14.6)' ) '    tstop = ', tstop

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Evaluate solution and residual at selected points (X,T)'
  
  call r8vec_linspace ( nx, 0.0D+00, 1.0D+00, x )
  call r8vec_linspace ( nt, t0, tstop, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      T       X       U(X,T)      Resid(X,T)'
  write ( *, '(a)' ) ''

  do i = 1, nt
    do j = 1, nx

      call kdv_exact_rational ( x(j), t(i), u, ut, ux, uxx, uxxx )
      call kdv_residual ( u, ut, ux, uxxx, r )

      if ( t(i) == 0.0 .and. x(j) == 0.0 ) then
        write ( *, '(2g14.6,a)' ) t(i), x(j), ' Undefined     Undefined'
      else
        write ( *, '(4g14.6)' ) t(i), x(j), u, r
      end if
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine kdv_exact_sech_test ( )

!*****************************************************************************80
!
!! kdv_exact_sech_test() tests kdv_exact_sech().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine kdv_parameters ( a_in, v_in, t0_in, tstop_in, &
      a_out, v_out, t0_out, tstop_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: t0_in
      real ( kind = rk ), optional :: t0_out
      real ( kind = rk ), optional :: tstop_in
      real ( kind = rk ), optional :: tstop_out
      real ( kind = rk ), optional :: v_in
      real ( kind = rk ), optional :: v_out
    end subroutine
  end interface

  integer, parameter :: nt = 6
  integer, parameter :: nx = 6

  real ( kind = rk ) a
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) t(nt)
  real ( kind = rk ) t0
  real ( kind = rk ) tstop
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) uxxx
  real ( kind = rk ) v
  real ( kind = rk ) x(nx)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'kdv_exact_sech_test():'
  write ( *, '(a)' ) '  Test kdv_exact_sech().'

  call kdv_parameters ( a_out = a, v_out = v, t0_out = t0, tstop_out = tstop)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Parameters:'
  write ( *, '(a,g14.6)' ) '    a     = ', a
  write ( *, '(a,g14.6)' ) '    v     = ', v
  write ( *, '(a,g14.6)' ) '    t0    = ', t0
  write ( *, '(a,g14.6)' ) '    tstop = ', tstop

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Evaluate solution and residual at selected points (X,T)'
  
  call r8vec_linspace ( nx, 0.0D+00, 1.0D+00, x )
  call r8vec_linspace ( nt, t0, tstop, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      T       X       U(X,T)      Resid(X,T)'
  write ( *, '(a)' ) ''

  do i = 1, nt
    do j = 1, nx

      call kdv_exact_sech ( x(j), t(i), u, ut, ux, uxx, uxxx )
      call kdv_residual ( u, ut, ux, uxxx, r )

      if ( t(i) == 0.0 .and. x(j) == 0.0 ) then
        write ( *, '(2g14.6,a)' ) t(i), x(j), ' Undefined      Undefined'
      else
        write ( *, '(4g14.6)' ) t(i), x(j), u, r
      end if
    end do
    write ( *, '(a)' ) ''
  end do

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
