program main

!*****************************************************************************80
!
!! satisfy_openmp() uses OpenMP to handle a formula satisfaction problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Michael Quinn,
!    Parallel Programming in C with MPI and OpenMP,
!    McGraw-Hill, 2004,
!    ISBN13: 978-0071232654,
!    LC: QA76.73.C15.Q55.
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 23

  integer bvec(n)
  integer circuit_value
  integer i
  integer id
  integer ihi
  integer ihi2
  integer ilo
  integer ilo2
  integer j
  integer proc_num
  integer solution_num_local
  integer solution_num
  integer thread_num
  integer value
  real ( kind = rk ) wtime

  call timestamp ( )

  proc_num = omp_get_num_procs ( )
  thread_num = omp_get_max_threads ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'satisfy_openmp():'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version' 
  write ( *, '(a)' ) '  We have a logical function of N logical arguments.'
  write ( *, '(a)' ) '  We do an exhaustive search of all 2^N possibilities,'
  write ( *, '(a)' ) '  seeking those inputs that make the function TRUE.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of processors available = ', proc_num
  write ( *, '(a,i8)' ) '  The number of threads available    = ', thread_num
!
!  The BIG calculation goes from 0 = ILO <= I < IHI = 2*N.
!  Compute the upper limit.
!
  ilo = 0

  ihi = 2**n

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The number of logical variables is N = ',  n
  write ( *, '(a,i8)' ) '  The number of input vectors to check is ', ihi
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   # Processor       Index    ' // &
    '---------Input Values------------------------'
  write ( *, '(a)' ) ' '
!
!  Processor ID takes the interval ILO2 <= I < IHI2.
!  Using the formulas below yields a set of nonintersecting intervals
!  which cover the original interval [ILO,IHI).
!
  solution_num = 0

  wtime = omp_get_wtime ( )

!$omp parallel &
!$omp shared ( ihi, ilo, thread_num ) &
!$omp private ( bvec, i, id, ilo2, ihi2, j, solution_num_local, value ) &
!$omp reduction ( + : solution_num )

  id = omp_get_thread_num ( )

  ilo2 = ( ( thread_num - id     ) * ilo   &
         + (              id     ) * ihi ) &
         / ( thread_num          )

  ihi2 = ( ( thread_num - id - 1 ) * ilo   &
         + (              id + 1 ) * ihi ) &
         / ( thread_num          )

  solution_num_local = 0

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,a,i8,a,i8)' )  &
    'Processor ', id, ' iterates from ', ilo2, ' <= I < ', ihi2
  write ( *, '(a)' ) ' '
!
!  Check every possible input vector.
!
  do i = ilo2, ihi2 - 1

    call i4_to_bvec ( i, n, bvec )

    value = circuit_value ( n, bvec )

    if ( value == 1 ) then
      solution_num_local = solution_num_local + 1

      write ( *, '(2x,i2,2x,i8,2x,i10,3x,23i2)' )  &
        solution_num_local, id, i, bvec(1:n)

    end if

  end do

  solution_num = solution_num + solution_num_local

!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of solutions found was ', solution_num
  write ( *, '(a,g14.6)' ) '  Elapsed wall clock time (seconds) ', wtime
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SATISFY_OPENMP'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function circuit_value ( n, bvec )

!*****************************************************************************80
!
!! CIRCUIT_VALUE returns the value of a circuit for a given input set.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
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
!    Input, integer N, the length of the input vector.
!
!    Input, integer BVEC(N), the binary inputs.
!
!    Output, integer CIRCUIT_VALUE, the output of the circuit.
!
  implicit none

  integer n

  integer bvec(n)
  integer circuit_value
  logical value

  value = ( bvec(1)  == 1 .or. bvec(2)  == 1 ) &
    .and. ( bvec(2)  == 0 .or. bvec(4)  == 0 ) &
    .and. ( bvec(3)  == 1 .or. bvec(4)  == 1 ) &
    .and. ( bvec(4)  == 0 .or. bvec(5)  == 0 ) &
    .and. ( bvec(5)  == 1 .or. bvec(6)  == 0 ) &
    .and. ( bvec(6)  == 1 .or. bvec(7)  == 0 ) &
    .and. ( bvec(6)  == 1 .or. bvec(7)  == 1 ) &
    .and. ( bvec(7)  == 1 .or. bvec(16) == 0 ) &
    .and. ( bvec(8)  == 1 .or. bvec(9)  == 0 ) &
    .and. ( bvec(8)  == 0 .or. bvec(14) == 0 ) &
    .and. ( bvec(9)  == 1 .or. bvec(10) == 1 ) &
    .and. ( bvec(9)  == 1 .or. bvec(10) == 0 ) &
    .and. ( bvec(10) == 0 .or. bvec(11) == 0 ) &
    .and. ( bvec(10) == 1 .or. bvec(12) == 1 ) &
    .and. ( bvec(11) == 1 .or. bvec(12) == 1 ) &
    .and. ( bvec(13) == 1 .or. bvec(14) == 1 ) &
    .and. ( bvec(14) == 1 .or. bvec(15) == 0 ) &
    .and. ( bvec(15) == 1 .or. bvec(16) == 1 ) &
    .and. ( bvec(15) == 1 .or. bvec(17) == 1 ) &
    .and. ( bvec(18) == 1 .or. bvec(2)  == 1 ) &
    .and. ( bvec(19) == 1 .or. bvec(1)  == 0 ) &
    .and. ( bvec(20) == 1 .or. bvec(2)  == 1 ) &
    .and. ( bvec(20) == 1 .or. bvec(19) == 0 ) &
    .and. ( bvec(20) == 0 .or. bvec(10) == 0 ) &
    .and. ( bvec(1)  == 1 .or. bvec(18) == 1 ) &
    .and. ( bvec(2)  == 0 .or. bvec(21) == 1 ) &
    .and. ( bvec(22) == 0 .or. bvec(21) == 1 ) &
    .and. ( bvec(23) == 0 .or. bvec(21) == 1 ) &
    .and. ( bvec(22) == 0 .or. bvec(21) == 0 ) &
    .and. ( bvec(23) == 1 .or. bvec(21) == 0 )

  if ( value ) then
    circuit_value = 1
  else
    circuit_value = 0
  end if

  return
end
subroutine i4_to_bvec ( i4, n, bvec )

!*****************************************************************************80
!
!! I4_TO_BVEC converts an integer into a binary vector.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I4, the integer.
!
!    Input, integer N, the dimension of the vector.
!
!    Output, integer BVEC(N), the vector of binary remainders.
!
  implicit none

  integer n

  integer bvec(n)
  integer i
  integer i4
  integer i4_copy

  i4_copy = i4

  do i = n, 1, -1
    bvec(i) = mod ( i4_copy, 2 )
    i4_copy = i4_copy / 2
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
