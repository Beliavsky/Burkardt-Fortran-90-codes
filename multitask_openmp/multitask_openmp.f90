program main

!*****************************************************************************80
!
!! multitask_openmp() multitasks, using OpenMP for parallel execution.
!
!  Discussion:
!
!    This program demonstrates how OpenMP can be used for multitasking, that 
!    is, a simple kind of parallel processing in which a certain number of 
!    perhaps quite unrelated tasks must be done.
!
!    The OpenMP SECTIONS directive identifies the portion of the program where
!    the code for these tasks is given.
!
!    The OpenMP SECTION directive is used repeatedly to divide this area of
!    the program into independent tasks.
!
!    The code will get the benefit of parallel processing up to the point where
!    there are as many threads as there are tasks.
!
!    The code will get a substantial speedup if the tasks take roughly the
!    same amount of time.  However, if one task takes substantially more time
!    than the others, this results in a limit to the parallel speedup that is
!    possible.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer prime_num
  integer, allocatable :: primes(:)
  integer sine_num
  real ( kind = rk ), allocatable :: sines(:)
  real ( kind = rk ) wtime
  real ( kind = rk ) wtime1
  real ( kind = rk ) wtime2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MULTITASK_OPENMP():'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'
  write ( *, '(a)' ) '  Demonstrate how OpenMP can "multitask" by using the'
  write ( *, '(a)' ) &
    '  SECTIONS directive to carry out several tasks in parallel.'

  prime_num = 20000
  allocate ( primes(1:prime_num) )
  sine_num = 20000
  allocate ( sines(1:sine_num) )

  wtime = omp_get_wtime ( )

!$omp parallel shared ( prime_num, primes, sine_num, sines )

  !$omp sections

    !$omp section
    wtime1 = omp_get_wtime ( )
    call prime_table ( prime_num, primes )
    wtime1 = omp_get_wtime ( ) - wtime1
    !$omp section
    wtime2 = omp_get_wtime ( )
    call sine_table ( sine_num, sines )
    wtime2 = omp_get_wtime ( ) - wtime2
  !$omp end sections

!$omp end parallel

  wtime = omp_get_wtime ( ) - wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6)' ) '  Number of primes computed was ', prime_num
  write ( *, '(a,i12)' ) '  Last prime was ', primes(prime_num)
  write ( *, '(a,i6)' ) '  Number of sines computed was ', sine_num
  write ( *, '(a,g14.6)' ) '  Last sine computed was ', sines(sine_num)
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Elapsed time = ', wtime
  write ( *, '(a,g14.6)' ) '  Task 1 time = ', wtime1
  write ( *, '(a,g14.6)' ) '  Task 2 time = ', wtime2
!
!  Free memory.
!
  deallocate ( primes )
  deallocate ( sines )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MULTITASK_OPENMP():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine prime_table ( prime_num, primes )

!*****************************************************************************80
!
!! PRIME_TABLE computes a table of the first PRIME_NUM prime numbers.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PRIME_NUM, the number of primes to compute.
!
!    Output, integer PRIMES(PRIME_NUM), the computed primes.
!
  implicit none

  integer prime_num

  integer i
  integer j
  integer p
  logical prime
  integer primes(prime_num)

  i = 2
  p = 0

  do while ( p < prime_num )

    prime = .true.

    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = .false.
        exit
      end if
    end do
      
    if ( prime ) then
      p = p + 1
      primes(p) = i
    end if

    i = i + 1

  end do

  return
end
subroutine sine_table ( sine_num, sines )

!*****************************************************************************80
!
!! SINE_TABLE computes a table of sines.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 October 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer SINE_NUM, the number of sines to compute.
!
!    Output, real ( kind = rk ) SINES(SINE_NUM), the sines.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer sine_num

  real ( kind = rk ) a
  integer i
  integer j
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sines(sine_num)

  do i = 1, sine_num
    sines(i) = 0.0D+00
    do j = 1, i
      a = real ( j - 1, kind = rk ) * r8_pi / real ( sine_num - 1, kind = rk )
      sines(i) = sines(i) + sin ( a )
    end do
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
