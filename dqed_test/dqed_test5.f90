program main

!*****************************************************************************80
!
!! MAIN is the main program for DQED_PRB5.
!
!  Discussion:
!
!    DQED_PRB5 demonstrates how to call the DBOLS routine directly,
!    bypassing the DQED routine, if the system to be solved is simply
!    a constrained linear system of equations.
!
!  Modified:
!
!    03 September 2007
!
!  Author:
!
!    Steve White wrote the original version of this example.
!    John Burkardt made some modifications.
!
!  Local Parameters:
!
!    Local, integer ( kind = 4 ) IOPT(NI+1), contains a final value of '99', preceded
!    by NI values indicating user options.  When calling DBOLS directly,
!    the option '7' must be specified.
!
!    Local, integer ( kind = 4 ) MDW, the allocated row dimension of W.
!
!    Local, integer ( kind = 4 ) MROWS, the number of rows in W.
!
!    Local, integer ( kind = 4 ) NCOLS, the number of columns in W.
!
!    Local, integer ( kind = 4 ) NI, the number of options in the IOPT array.
!
  implicit none

  integer ( kind = 4 ), parameter :: mrows = 5
  integer ( kind = 4 ), parameter :: ncols = 3
  integer ( kind = 4 ), parameter :: ni = 1
  integer ( kind = 4 ), parameter :: nx = 0

  integer ( kind = 4 ), parameter :: mdw = mrows

  real ( kind = 8 ),dimension ( ncols ) :: bl = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  real ( kind = 8 ),dimension ( ncols ) :: bu = (/ 0.0D+00, 0.0D+00, 0.0D+00 /)
  integer ( kind = 4 ) col
  integer ( kind = 4 ), dimension ( ncols ) :: ind = (/ 1, 1, 1 /)
  integer ( kind = 4 ), dimension ( 1 + ni ) :: iopt = (/ 7, 99 /)
  integer ( kind = 4 ) iw(2*ncols)
  integer ( kind = 4 ) mode
  real ( kind = 8 ) rnorm
  integer ( kind = 4 ) row
  real ( kind = 8 ) rw(5*ncols)
  real ( kind = 8 ), dimension ( mdw, ncols+1 ) :: w = reshape ( (/ &
      4.0D+00,  3.0D+00, -7.0D+00, 7.0D+00,  -1.0D+00, &
     -7.0D+00, -4.0D+00,  3.0D+00,  7.0D+00,  7.0D+00, &
     -7.0D+00,  0.0D+00, -3.0D+00,  1.0D+00,  6.0D+00, &
    -22.0D+00, -6.0D+00, -2.0D+00, 34.0D+00, 13.0D+00 /), (/ mdw, ncols+1 /) )
  real ( kind = 8 ), dimension ( ncols+nx) :: x

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB5'
  write ( *, '(a)' ) '  Demonstrate how to call the DBOLS routine directly'
  write ( *, '(a)' ) '  if the nonlinear constrained system to be solved'
  write ( *, '(a)' ) '  is actually simply a LINEAR constrained system.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The linear system E*X=F is stored as a single '
  write ( *, '(a)' ) '  array W = [ E | F ].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) &
    '  The order of E is ', mrows, ' rows by ', ncols, ' columns.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The E matrix:'
  write ( *, '(a)' ) ' '
  do row = 1, mrows
    write ( *, '(5g12.3)' ) w(row,1:ncols)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The F right hand side:'
  write ( *, '(a)' ) ' '
  do row = 1, mrows
    write ( *, '(g12.3)' ) w(row,ncols+1)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Lower bounds on the solution X:'
  write ( *, '(a)' ) ' '
  do col = 1, ncols
    write ( *, '(g12.3)' ) bl(col)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Call DBOLS for least squares solution minimizing'
  write ( *, '(a)' ) '  the norm of the constrained linear system.'

  call dbols ( w, mdw, mrows, ncols, bl, bu, ind, iopt, x, rnorm, mode, rw, iw )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The solution X:'
  write ( *, '(a)' ) ' '
  do col = 1, ncols
    write ( *, '(g12.3)' ) x(col)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The L2 norm of the residual vector E*X-F = ', rnorm
  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  The error return flag is MODE = ', mode
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DQED_PRB5:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
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

