program main

!*****************************************************************************80
!
!! random_scatter() uses DISLIN to draw a scatter plot of X Y data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 July 2020
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

  integer, parameter :: n = 500

  integer i
  integer j
  integer pat
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) xvec(n)
  real ( kind = rk ) yvec(n)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'random_scatter():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Use DISLIN to make a scatter plot of random data.'
!
!  Generate the data.  
!  We average 4 random values to get data that tends to cluster
!  near (0.5,0.5).
!
  seed = 123456789

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
  call setfil ( 'random_scatter.png' )
!
!  Choose the page size and orientation.
!  'USA' is 2160 plot units wide and 2790 plot units high.
!  'P' requests PROFILE rather than LANDSCAPE orientation.
!
  call setpag ( 'usap' )
!
!  For PNG output, use reverse the default black background to white.
!
  call scrmod ( 'reverse' )
!
!  Open DISLIN.
!
  call disini ( )
!
!  Plot a border around the page.
!
  call pagera ( )
!
!  Use the COMPLEX font.
!
  call complx ( )
!
!  Define the X and Y sizes of the axis system in plot units.
!
  call axslen ( 1800, 1800 )
!
!  Specify how the lower X, left Y, upper X and right Y axes are labeled.
!
  call setgrf ( 'line', 'line', 'line', 'line' )
!
!  Set the axis origin 180 plot units to the right, and 2610 plot units DOWN.
!
  call axspos ( 180, 2610 )
!
!  Relate the physical coordinates to the axes.
!
  call graf ( 0.0D+00, 1.0D+00, 0.0D+00, 0.1D+00, 0.0D+00, 1.0D+00, 0.0D+00, 0.1D+00 )
!
!  Add a grid, with one grid line for every tick mark in the X and Y axes.
!
  call grid ( 1, 1 )
!
!  Select the shading pattern.
!
  pat = 16
  call shdpat ( pat )
!
!  Set the color to blue.
!
  call color ( "blue" )
!
!  At every data point, draw a circle of radius 0.01.
!
  do i = 1, n
    call rlcirc ( xvec(i), yvec(i), 0.01D+00 )
  end do
!
!  Select character height in plot units.
!
  call height ( 50 )
!
!  We choose "white" for the title, which is actually black
!  because we reversed black and white earlier!
!
  call color ( "white" )
!
!  Define axis system titles.
!
  call titlin ( 'Scatter plot of random data', 1 )
!
!  Draw the title.
!
  call title ( )
!
!  End this plot.
!
  call endgrf ( )
!
!  Close DISLIN.
!
  call disfin ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDOM_SCATTER:'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  stop 0
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp prints the current YMDHMS date as a time stamp.
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

