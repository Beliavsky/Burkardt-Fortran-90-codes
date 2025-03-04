program main

!*****************************************************************************80
!
!! heated_plate() solves the steady heat equation on a rectangle.
!
!  Discussion:
!
!    This code solves the steady state heat equation on a rectangular region.
!
!    The sequential version of this program needs approximately
!    18/eps iterations to complete. 
!
!
!    The physical region, and the boundary conditions, are suggested
!    by this diagram;
!
!                   W = 0
!             +------------------+
!             |                  |
!    W = 100  |                  | W = 100
!             |                  |
!             +------------------+
!                   W = 100
!
!    The region is covered with a grid of M by N nodes, and an N by N
!    array W is used to record the temperature.  The correspondence between
!    array indices and locations in the region is suggested by giving the
!    indices of the four corners:
!
!                  I = 0
!          [0][0]-------------[0][N-1]
!             |                  |
!      J = 0  |                  |  J = N-1
!             |                  |
!        [M-1][0]-----------[M-1][N-1]
!                  I = M-1
!
!    The steady state solution to the discrete heat equation satisfies the
!    following condition at an interior grid point:
!
!      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
!
!    where "Central" is the index of the grid point, "North" is the index
!    of its immediate neighbor to the "north", and so on.
!   
!    Given an approximate solution of the steady state heat equation, a
!    "better" solution is given by replacing each interior point by the
!    average of its 4 neighbors - in other words, by using the condition
!    as an ASSIGNMENT statement:
!
!      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
!
!    If this process is repeated often enough, the difference between 
!    successive estimates of the solution will go to zero.
!
!    This program carries out such an iteration, using a tolerance specified by
!    the user, and writes the final estimate of the solution to a file that can
!    be used for graphic processing.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 July 2008
!
!  Author:
!
!    Original FORTRAN90 version by Michael Quinn.
!    This FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Michael Quinn,
!    Parallel Programming in C with MPI and OpenMP,
!    McGraw-Hill, 2004,
!    ISBN13: 978-0071232654,
!    LC: QA76.73.C15.Q55.
!
!  Parameters:
!
!    Commandline argument 1, real ( kind = rk ) EPS, the error tolerance.  
!
!    Commandline argument 2, character OUTPUT_FILENAME, the name of the file 
!    into which the steady state solution is written when the program 
!    has completed.
!
!  Local:
!
!    real ( kind = rk ) DIFF, the norm of the change in the solution from 
!    one iteration to the next.
!
!    real ( kind = rk ) MEAN, the average of the boundary values, used 
!    to initialize the values of the solution in the interior.
!
!    real ( kind = rk ) U(M,N), the solution at the previous iteration.
!
!    real ( kind = rk ) W(M,N), the solution computed at the latest 
!    iteration.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 500
  integer, parameter :: n = 500

  character ( len = 255 ) :: arg
  real ( kind = rk ) ctime
  real ( kind = rk ) ctime1
  real ( kind = rk ) ctime2
  real ( kind = rk ) diff
  real ( kind = rk ) eps
  integer i
  integer iterations
  integer iterations_print
  real ( kind = rk ) mean
  integer numarg
  character ( len = 255 ) output_filename
  integer output_unit
  real ( kind = rk ) u(m,N)
  real ( kind = rk ) w(m,N)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'heated_plate():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) &
    '  A program to solve for the steady state temperature distribution'
  write ( *, '(a)' ) '  over a rectangular plate.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a)' ) '  Spatial grid of ', m, ' by ', N, ' points.'

  numarg = iargc ( )
!
!  Read EPSILON from the command line or the user.
!
  if ( numarg < 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter EPS, the error tolerance.'
    read ( *, * ) eps
  else
    call getarg ( 1, arg )
    read ( arg, * ) eps
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) &
    '  The iteration will repeat until the change is <= ', eps
  diff = eps
!
!  Read OUTPUT_FILENAME from the command line or the user.
!
  if ( numarg < 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Enter OUTPUT_FILE, the name of the output file.'
    read ( *, '(a)' ) output_filename
  else
    call getarg ( 2, output_filename )
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The steady state solution will be written to "' &
    // trim ( output_filename ) // '".'
!
!  Set the boundary values, which don't change. 
!
  w(2:m-1,1) = 100.0
  w(2:m-1,n) = 100.0
  w(m,1:n) = 100.0
  w(1,1:n) =   0.0
!
!  Average the boundary values, to come up with a reasonable
!  initial value for the interior.
!
  mean = ( &
      sum ( w(2:m-1,1) ) &
    + sum ( w(2:m-1,n) ) &
    + sum ( w(m,1:n)   ) &
    + sum ( w(1,1:n)   ) ) &
    / dble ( 2 * m + 2 * n - 4 )
!
!  Initialize the interior solution to the mean value.
!
  w(2:m-1,2:n-1) = mean
!
!  iterate until the  new solution W differs from the old solution U
!  by no more than EPS.
!
  iterations = 0
  iterations_print = 1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Iteration  Change'
  write ( *, '(a)' ) ' '
  call cpu_time ( ctime1 )

  do while ( eps <= diff )

    u(1:m,1:n) = w(1:m,1:n)

    w(2:m-1,2:n-1) = 0.25 * ( &
        u(1:m-2,2:n-1) &
      + u(3:m,2:n-1) &
      + u(2:m-1,1:n-2) &
      + u(2:m-1,3:n) )

    diff = maxval ( abs ( u(1:m,1:n) - w(1:m,1:n) ) )

    iterations = iterations + 1

    if ( iterations == iterations_print ) then
      write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
      iterations_print = 2 * iterations_print
    end if

  end do

  call cpu_time ( ctime2 )
  ctime = ctime2 - ctime

  write ( *, '(a)' ) ' '
  write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Error tolerance achieved.'
  write ( *, '(a,g14.6)' ) '  CPU time = ', ctime
!
!  Write the solution to the output file.
!
  output_unit = 10

  open ( unit = output_unit, file = output_filename )

  write ( output_unit, * ) m
  write ( output_unit, * ) n
  do i = 1, m
    write ( output_unit, * ) w(i,1:n)
  end do

  close ( unit = output_unit )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Solution written to the output file "' &
    // trim ( output_filename ) // '".'
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HEATED_PLATE:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
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
