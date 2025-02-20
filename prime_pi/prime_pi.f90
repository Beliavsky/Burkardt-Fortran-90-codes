function prime_pi1 ( n )

!*****************************************************************************80
!
!! prime_pi1() evaluates prime(n), the number of primes less than or equal to n.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 December 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the value to be tested.
!
!  Output:
!
!    integer value: the value of pi(n).
!
  implicit none

  integer i
  logical is_prime3
  integer n
  integer prime_pi1

  prime_pi1 = 0
  do i = 1, n
    if ( is_prime3 ( i ) ) then
      prime_pi1 = prime_pi1 + 1
    end if
  end do

  return
end
function is_prime3 ( n )

!*****************************************************************************80
!
!! is_prime3() reports whether an integer is prime.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the value to be tested.
!
!  Output:
!
!    logical is_prime3: true if n is a prime.
!
  implicit none

  integer i
  integer ihi
  logical is_prime3
  integer n

  if ( n <= 1 ) then
    is_prime3 = .false.
    return
  end if

  if ( ( n == 2 ) .or. ( n == 3 ) ) then
    is_prime3 = .true.
    return
  end if

  if ( ( mod ( n, 2 ) == 0 ) .or. ( mod ( n, 3 ) == 0 ) ) then
    is_prime3 = .false.
    return
  end if

  ihi = int ( sqrt ( real ( n ) ) )

  do i = 5, ihi, 6
    if ( ( mod ( n, i ) == 0 ) .or. ( mod ( n, i + 2 ) == 0 ) ) then
      is_prime3 = .false.
      return
    end if
  end do

  is_prime3 = .true.

  return
end

