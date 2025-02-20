program main

!*****************************************************************************80
!
!! SEARCH_MPI uses MPI to perform a search in parallel.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2012
!
!  Author:
!
!    John Burkardt
!
  use mpi

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer a
  integer b
  integer c
  integer f
  integer, parameter :: i4_huge = 2147483647
  integer id
  integer ierr
  integer j
  integer p
  real ( kind = rk ) wtime

  call MPI_Init ( ierr )

  call MPI_Comm_Rank ( MPI_COMM_WORLD, id, ierr )

  call MPI_Comm_Size ( MPI_COMM_WORLD, p, ierr )

  a = 1
  b = i4_huge
  c = 45

  if ( id == 0 ) then

    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'search_mpi():'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) '  Search the integers from A to B'
    write ( *, '(a)' ) '  for a value J such that F(J) = C.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  A           = ', a
    write ( *, '(a,i12)' ) '  B           = ', b
    write ( *, '(a,i12)' ) '  C           = ', c
  end if

  wtime = MPI_Wtime ( )

  call search ( a, b, c, id, p, j )

  wtime = MPI_Wtime ( )- wtime
!
!  Any process that finds a solution should report it.
!
  if ( j /= -1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a,i1,a,i12)' ) '  Process ', id, ' found     J = ', j
    write ( *, '(a,i12)' ) '  Verify F(J) = ', f ( j )
  end if

  if ( id == 0 ) then
    write ( *, '(a,g14.6)' ) '  Elapsed CPU time is ', wtime
  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'search_mpi():'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine search ( a, b, c, id, p, j )

!*****************************************************************************80
!
!! search() searches integers in [A,B] for a J so that F(J) = C.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer A, B, the search range.
!
!    Input, integer C, the desired function value.
!
!    Input, integer ID, the identifier of this process.
!    0 <= ID < P.
!
!    Input, integer P, the number of processes.
!
!    Output, integer J, the computed solution, or -1
!    if no solution was found.
!
  implicit none

  integer a
  integer b
  integer c
  integer f
  integer fi
  integer i
  integer id
  integer j
  integer p

  j = -1
!
!  May have to be careful here, since I + P could suddenly become
!  negative when we exceed the maximum integer!
!
  do i = a + id, b, p

    fi = f ( i )

    if ( fi == c ) then
      j = i
      return
    end if

  end do

  return
end
function f ( i )

!*****************************************************************************80
!
!! f() is the function we are analyzing.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the argument.
!
!    Input, integer F, the value.
!
  implicit none

  integer f
  integer i
  integer, parameter :: i4_huge = 2147483647
  integer j
  integer k
  integer value

  value = i

  do j = 1, 5

    k = value / 127773

    value = 16807 * ( value - k * 127773 ) - k * 2836

    if ( value < 0 ) then
      value = value + i4_huge
    end if

  end do

  f = value

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
