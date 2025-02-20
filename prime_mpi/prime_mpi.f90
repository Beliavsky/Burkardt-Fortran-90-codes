program main

!*****************************************************************************80
!
!! prime_mpi() counts prime numbers, using MPI for parallel processing.
!
!  Discussion:
!
!    This program calls a version of PRIME_NUMBER that includes
!    MPI calls for parallel processing.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 August 2009
!
!  Author:
!
!    John Burkardt
!
  use mpi

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer id
  integer ierr
  integer n
  integer n_factor
  integer n_hi
  integer n_lo
  integer p
  integer primes
  integer primes_part
  real ( kind = rk ) wtime

  n_lo = 1
  n_hi = 131072
  n_factor = 2
!
!  Initialize MPI.
!
  call MPI_Init ( ierr )
!
!  Get this process's ID.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
!
!  Find out how many processes are available.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

  if ( id == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'prime_mpi():'
    write ( *, '(a)' ) '  FORTRAN90/MPI version'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  An MPI example program to count the number of primes.'
    write ( *, '(a,i8)' ) '  The number of processes is ', p
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '         N        Pi          Time'
    write ( *, '(a)' ) ' '
  end if

  n = n_lo

  do while ( n <= n_hi )

    if ( id == 0 ) then
      wtime = MPI_Wtime ( )
    end if

    call MPI_Bcast ( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )

    call prime_number ( n, id, p, primes_part )

    call MPI_Reduce ( primes_part, primes, 1, MPI_INTEGER, MPI_SUM, 0, &
      MPI_COMM_WORLD, ierr )

    if ( id == 0 ) then
      wtime = MPI_Wtime ( ) - wtime
      write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, wtime
    end if

    n = n * n_factor

  end do 
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'PRIME_MPI():'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '
    call timestamp ( )
  end if

  stop
end
subroutine prime_number ( n, id, p, total )

!*****************************************************************************80
!
!! PRIME_NUMBER returns a part of the number of primes between 1 and N.
!
!  Discussion:
!
!    In order to divide the work up evenly among P processors, processor
!    ID starts at 2+ID and skips by P.
!
!    A naive algorithm is used.
!
!    Mathematica can return the number of primes less than or equal to N
!    by the command PrimePi[N].
!
!                N  PRIME_NUMBER
!
!                1           0
!               10           4
!              100          25
!            1,000         168
!           10,000       1,229
!          100,000       9,592
!        1,000,000      78,498
!       10,000,000     664,579
!      100,000,000   5,761,455
!    1,000,000,000  50,847,534
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 May 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the maximum number to check.
!
!    Input, integer ID, the ID of this process,
!    between 0 and P-1.
!
!    Input, integer P, the number of processes.
!
!    Output, integer TOTAL, the number of prime numbers up to N,
!    starting at 2+ID and skipping by P.
!
  implicit none

  integer i
  integer id
  integer j
  integer n
  integer p
  integer prime
  integer total

  total = 0

  do i = 2+id, n, p

    prime = 1

    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = 0
        exit
      end if
    end do

    total = total + prime

  end do

  return
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
