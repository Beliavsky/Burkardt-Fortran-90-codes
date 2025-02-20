program main

!*****************************************************************************80
!
!! poisson_2d_test() tests poisson_2d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 October 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny
  external rhs1
  real ( kind = rk8 ) t_start
  real ( kind = rk8 ) t_stop
  real ( kind = rk8 ), external :: u_exact1
  real ( kind = rk8 ), external :: uxxyy_exact1
  real ( kind = rk8 ) wtime

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poisson_2d_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test poisson_2d().'
!
!  Problem 1.
!
  t_start = wtime ( )
  nx = 161
  ny = 161
  call poisson_2d ( nx, ny, rhs1, u_exact1 )
  t_stop = wtime ( )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Wall clock time in seconds is ', t_stop - t_start
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'poisson_2d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine rhs1 ( nx, ny, f )

!*****************************************************************************80
!
!! rhs1() initializes the right hand side "vector" for case #1.
!
!  Discussion:
!
!    It is convenient for us to set up RHS as a 2D array.  However, each
!    entry of RHS is really the right hand side of a linear system of the
!    form
!
!      A * U = F
!
!    In cases where U(I,J) is a boundary value, then the equation is simply
!
!      U(I,J) = F(i,j)
!
!    and F(I,J) holds the boundary data.
!
!    Otherwise, the equation has the form
!
!      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)
!
!    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of
!
!      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer NX, NY, the X and Y grid dimensions.
!
!  Output:
!
!    real ( kind = rk8 ) F(NX,NY), the right hand side data.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer nx
  integer ny

  real ( kind = rk8 ) f(nx,ny)
  real ( kind = rk8 ) fnorm
  integer i
  integer j
  real ( kind = rk8 ) rms_norm
  real ( kind = rk8 ) u_exact1
  real ( kind = rk8 ) uxxyy_exact1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y
!
!  The "boundary" entries of F will store the boundary values of the solution.
!
  x = 0.0D+00
  do j = 1, ny
    y = real ( j - 1, kind = rk8 ) / real ( ny - 1, kind = rk8 )
    f(1,j) = u_exact1 ( x, y )
  end do

  x = 1.0D+00
  do j = 1, ny
    y = real ( j - 1, kind = rk8 ) / real ( ny - 1, kind = rk8 )
    f(nx,j) = u_exact1 ( x, y )
  end do

  y = 0.0D+00
  do i = 1, nx
    x = real ( i - 1, kind = rk8 ) / real ( nx - 1, kind = rk8 )
    f(i,1) = u_exact1 ( x, y )
  end do

  y = 1.0D+00
  do i = 1, nx
    x = real ( i - 1, kind = rk8 ) / real ( nx - 1, kind = rk8 )
    f(i,ny) = u_exact1 ( x, y )
  end do
!
!  The "interior" entries of F store the right hand sides 
!  of the Poisson equation.
!
  do j = 2, ny - 1
    y = real ( j - 1, kind = rk8 ) / real ( ny - 1, kind = rk8 )
    do i = 2, nx - 1
      x = real ( i - 1, kind = rk8 ) / real ( nx - 1, kind = rk8 )
      f(i,j) = - uxxyy_exact1 ( x, y )
    end do
  end do

  fnorm = rms_norm ( nx, ny, f ) 

  write ( *, '(a,g14.6)' ) '  RMS of F = ', fnorm

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
function u_exact1 ( x, y )

!*****************************************************************************80
!
!! u_exact1() evaluates the exact solution for case #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, Y, the coordinates of a point.
!
!  Output:
!
!    real ( kind = rk8 ) U_EXACT1, the value of the exact solution at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) u_exact1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  u_exact1 = sin ( pi * x * y )

  return
end
function uxxyy_exact1 ( x, y )

!*****************************************************************************80
!
!! uxxyy_exact1() evaluates ( d/dx d/dx + d/dy d/dy ) for case #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk8 ) X, Y, the coordinates of a point.
!
!  Output:
!
!    real ( kind = rk8 ) UXXYY_EXACT1, the value of 
!    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) uxxyy_exact1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  uxxyy_exact1 = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y )

  return
end
function wtime ( )

!*****************************************************************************80
!
!! wtime() returns a reading of the wall clock time.
!
!  Discussion:
!
!    To get the elapsed wall clock time, call WTIME before and after a given
!    operation, and subtract the first reading from the second.
!
!    This function is meant to suggest the similar routines:
!
!      "omp_get_wtime ( )" in OpenMP,
!      "MPI_Wtime ( )" in MPI,
!      and "tic" and "toc" in MATLAB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 April 2009
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk ) WTIME, the wall clock reading, in seconds.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer clock_max
  integer clock_rate
  integer clock_reading
  real ( kind = rk ) wtime

  call system_clock ( clock_reading, clock_rate, clock_max )

  wtime = real ( clock_reading, kind = rk ) &
        / real ( clock_rate, kind = rk )

  return
end

