program main

!*****************************************************************************80
!
!! hello_mpi() is a simple MPI test program.
!
!  Discussion:
!
!    Each process prints out a "Hello, world!" message.
!
!    The master process also prints out a short message.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    William Gropp, Ewing Lusk, Anthony Skjellum,
!    Using MPI: Portable Parallel Programming with the
!    Message-Passing Interface,
!    Second Edition,
!    MIT Press, 1999,
!    ISBN: 0262571323,
!    LC: QA76.642.G76.
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer error
  integer id
  integer p
  real ( kind = rk ) wtime
!
!  Initialize MPI.
!
  call MPI_Init ( error )
!
!  Get the number of processes.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, error )
!
!  Get the individual process ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, error )
!
!  Print a message.
!
  if ( id == 0 ) then

    wtime = MPI_Wtime ( )

    call timestamp ( )
    write ( *, '(a)' ) ''
    write ( *, '(a,i1,2x,a)' ) 'P', id, 'HELLO_MPI - Master process:'
    write ( *, '(a,i1,2x,a)' ) 'P', id, '  FORTRAN90/MPI version'
    write ( *, '(a,i1,2x,a)' ) 'P', id, '  An MPI test program.'
    write ( *, '(a,i1,2x,a,i8)' ) 'P', id, '  The number of MPI processes is ', p

  end if
!
!  Every MPI process will print this message.
!
  write ( *, '(a,i1,2x,a)' ) 'P', id, '"Hello, world!"'

  if ( id == 0 ) then

    write ( *, '(a)' ) ''
    write ( *, '(a,i1,2x,a)' ) 'P', id, 'HELLO_MPI - Master process:'
    write ( *, '(a,i1,2x,a)' ) 'P', id, '  Normal end of execution: "Goodbye, world!".'

    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a)' ) ''
    write ( *, '(a,i1,2x,a,g14.6,a)' ) &
      'P', id, '  Elapsed wall clock time = ', wtime, ' seconds.'

  end if
!
!  Shut down MPI.
!
  call MPI_Finalize ( error )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a,i1,2x,a)' ) 'P', id, 'HELLO_MPI - Master process:'
    write ( *, '(a,i1,2x,a)' ) 'P', id, '  Normal end of execution.'
    write ( *, '(a)' ) ''
    call timestamp ( )
  end if

  stop
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
!    06 August 2005
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
