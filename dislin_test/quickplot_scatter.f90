program main

!*****************************************************************************80
!
!! quickplot_scatter() demonstrates the DISLIN quickplot command QPLOT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 February 2014
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Helmut Michels,
!    The Data Plotting Software DISLIN - version 10.4,
!    Shaker Media GmbH, January 2010,
!    ISBN13: 978-3-86858-517-9.
!
  use dislin

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 100

  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) xvec(n)
  real ( kind = rk ) yvec(n)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_SCATTER:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Demonstrate the DISLIN "quickplot" command QPLSCA'
  write ( *, '(a)' ) '  to make a scatter plot.'
!
!  Generate the data.  
!  We average 4 random values to get data that tends to cluster
!  near (0.5,0.5).
!
  do i = 1, n
    s = 0.0D+00
    do j = 1, 4
      call random_number ( harvest = r )
      s = s + r
    end do
    xvec(i) = s / 4.0D+00
  end do

  do i = 1, n
    s = 0.0D+00
    do j = 1, 4
      call random_number ( harvest = r )
      s = s + r
    end do
    yvec(i) = s / 4.0D+00
  end do
!
!  Specify the format of the output file.
!
  call metafl ( 'png' )
!
!  Indicate that new data overwrites old data.
!
  call filmod ( 'delete' )
!
!  Specify the name of the output graphics file.
!
  call setfil ( 'quickplot_scatter.png' )
!
!  Choose the page size and orientation.
!
  call setpag ( 'usal' )
!
!  For PNG output, reverse the default black background to white.
!
  call scrmod ( 'reverse' )
!
!  Open DISLIN.
!
  call disini ( )
!
!  Label the axes and the plot.
!
  call name ( '<-- X -->', 'X' )
  call name ( '<-- Y -->', 'Y' )
  call titlin ( 'Quick plot by QPLSCA', 2 )
!
!  Draw the curve.
!
  call qplsca ( xvec, yvec, n )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUICKPLOT_SCATTER:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop ( 0 )
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
