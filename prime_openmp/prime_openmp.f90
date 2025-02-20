program main

!*****************************************************************************80
!
!! prime_openmp() counts prime numbers, using OpenMP for parallel execution.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer n_factor
  integer n_hi
  integer n_lo

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP():'
  write ( *, '(a)' ) '  FORTRAN90/OpenMP version'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Number of processors available = ', omp_get_num_procs ( )
  write ( * ,'(a,i8)' ) &
    '  Number of threads =              ', omp_get_max_threads ( )

  n_lo = 1
  n_hi = 131072
  n_factor = 2

  call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

  n_lo = 5
  n_hi = 500000
  n_factor = 10

  call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP():'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  stop 0
end
subroutine prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

!*****************************************************************************80
!
!! PRIME_NUMBER_SWEEP_OPENMP does repeated calls to PRIME_NUMBER.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N_LO, the first value of N.
!
!    Input, integer N_HI, the last value of N.
!
!    Input, integer N_FACTOR, the factor by which to increase N
!    after each iteration.
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_factor
  integer n_hi
  integer n_lo
  integer primes
  real ( kind = rk ) wtime

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_SWEEP_OPENMP():'
  write ( *, '(a)' ) '  Call PRIME_NUMBER() to count the primes from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        Pi          Time'
  write ( *, '(a)' ) ' '

  n = n_lo

  do while ( n <= n_hi )

    wtime = omp_get_wtime ( )

    call prime_number ( n, primes )

    wtime = omp_get_wtime ( ) - wtime

    write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, wtime

    n = n * n_factor

  end do
 
  return
end
subroutine prime_number ( n, total )

!*****************************************************************************80
!
!! PRIME_NUMBER returns the number of primes between 1 and N.
!
!  Discussion:
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
!    Output, integer TOTAL, the number of prime numbers up to N.
!
  implicit none

  integer i
  integer j
  integer n
  integer prime
  integer total

  total = 0

!$omp parallel &
!$omp shared ( n ) &
!$omp private ( i, j, prime )

!$omp do reduction ( + : total )

  do i = 2, n

    prime = 1

    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        prime = 0
        exit
      end if
    end do

    total = total + prime

  end do

!$omp end do

!$omp end parallel

  return
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

