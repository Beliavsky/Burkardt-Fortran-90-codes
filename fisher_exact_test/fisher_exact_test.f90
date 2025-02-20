program main

!*****************************************************************************80
!
!! fisher_exact_test() tests fisher_exact().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: c_in
      real ( kind = rk ), optional :: c_out
      real ( kind = rk ), optional :: k_in
      real ( kind = rk ), optional :: k_out
    end subroutine
  end interface

  real ( kind = rk ) a
  real ( kind = rk ) c
  real ( kind = rk ) k

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'fisher_exact_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test fisher_exact(), for exact solutions of Fisher PDE.'
!
!  Report the current parameter values.
!
  call fisher_parameters ( a_out = a, c_out = c, k_out = k )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  parameters:'
  write ( *, '(a,g14.6)' ) '    a = ', a
  write ( *, '(a,g14.6)' ) '    c = ', c
  write ( *, '(a,g14.6)' ) '    k = ', k

  call f_residual_test ( )

  call fisher_residual_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'fisher_exact_test():'
  write ( *, '(a)' ) '  Normal end of execution.'

  call timestamp ( )

  return
end
subroutine f_residual_test ( )

!*****************************************************************************80
!
!! f_residual_test() tests f_residual().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: c_in
      real ( kind = rk ), optional :: c_out
      real ( kind = rk ), optional :: k_in
      real ( kind = rk ), optional :: k_out
    end subroutine
  end interface

  real ( kind = rk ) c
  real ( kind = rk ) f
  real ( kind = rk ) fz
  real ( kind = rk ) fzz
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) t(6)
  real ( kind = rk ) x(6)
  real ( kind = rk ) z

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'f_residual_test():'
  write ( *, '(a)' ) '  Evaluate F(Z) and residual at selected points Z=X-cT'

  call r8vec_linspace ( 6, 0.0D+00, 1.0D+00, x )
  call r8vec_linspace ( 6, 0.0D+00, 10.0D+00, t )

  call fisher_parameters ( c_out = c )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      Z       F(Z)      Resid(Z)'
  write ( *, '(a)' ) ''

  do j = 1, 6
    do i = 1, 6
      z = x(i) - c * t(j)
      call f_exact ( z, f, fz, fzz )
      call f_residual ( f, fz, fzz, r )
      write ( *, '(3g14.6)' ) z, f, r
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine fisher_residual_test ( )

!*****************************************************************************80
!
!! fisher_residual_test() tests fisher_residual().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 May 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) t(6)
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) x(6)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'fisher_residual_test():'
  write ( *, '(a)' ) '  Evaluate solution and residual at selected points (X,T)'
  
  call r8vec_linspace ( 6, 0.0D+00, 1.0D+00, x )
  call r8vec_linspace ( 6, 0.0D+00, 10.0D+00, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '      X       T       U(X,T)      Resid(X,T)'
  write ( *, '(a)' ) ''

  do j = 1, 6
    do i = 1, 6
      call fisher_exact ( x(i), t(j), u, ut, ux, uxx )
      call fisher_residual ( x(i), t(j), r )
      write ( *, '(4g14.6)' ) x(i), t(j), u, r
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
