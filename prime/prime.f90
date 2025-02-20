subroutine prime_number ( n, total )

!*****************************************************************************80
!
!! prime_number() returns the number of primes between 1 and N.
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
!    23 April 2009
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

  return
end
subroutine prime_number_sweep ( n_lo, n_hi, n_factor )

!*****************************************************************************80
!
!! prime_number_sweep() does repeated timed calls to PRIME_NUMBER.
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
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_factor
  integer n_hi
  integer n_lo
  integer primes
  real ( kind = rk ) time1
  real ( kind = rk ) time2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_NUMBER_SWEEP():'
  write ( *, '(a)' ) '  Call PRIME_NUMBER to count the primes from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        Pi          Time'
  write ( *, '(a)' ) ' '

  n = n_lo

  do while ( n <= n_hi )

    call cpu_time ( time1 )

    call prime_number ( n, primes )

    call cpu_time ( time2 )

    write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, time2 - time1

    n = n * n_factor

  end do

  return
end

