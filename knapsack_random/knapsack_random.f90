subroutine knapsack_random ( n, s )

!*****************************************************************************80
!
!! knapsack_random() returns a random possible solution of a knapsack problem.
!
!  Discussion:
!
!    The subset is represented as a vector of binary digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the total number of items in the set.
!
!  Output:
!
!    integer s[n]: a subset.  Item i is in the subset if s[i] = 1.
!
  implicit none

  integer n

  integer s(n)

  call subset_random ( n, s )

  return
end
subroutine subset_next ( n, s )

!*****************************************************************************80
!
!! subset_next() returns the next subset of n items.
!
!  Discussion:
!
!    The subset is represented as a vector of binary digits.
!    The subsets are listed in numerical order.
!    After the last subset is returned, the sequence begins again at 0.
!
!  Example:
!
!    [ 0, 0, 0 ]
!    [ 0, 0, 1 ]
!    [ 0, 1, 0 ]
!    [ 0, 1, 1 ]
!    [ 1, 0, 0 ]
!    [ 1, 0, 1 ]
!    [ 1, 1, 0 ]
!    [ 1, 1, 1 ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of items in the set.
!
!    integer s[n]: the current subset.
!
!  Output:
!
!    integer s[n]: the next subset.
!
  implicit none

  integer n

  integer i
  integer s(n)
!
!  Starting at the last index, add 1, and if necessary, carry to the left.
!
  do i = n, 1, -1

    if ( s(i) == 0 ) then
      s(i) = 1
      exit
    else
      s(i) = 0
    end if

  end do

  return
end
subroutine subset_random ( n, s )

!*****************************************************************************80
!
!! subset_random() returns a random subset of n items.
!
!  Discussion:
!
!    The subset is represented as a vector of binary digits.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the total number of items in the set.
!
!  Output:
!
!    integer s[n]: a subset.  Item i is in the subset if s[i] = 1.
!
  implicit none

  integer n

  real r(n)
  integer s(n)

  call random_number ( harvest = r(1:n) )
  s = nint ( r )

  return
end
subroutine subset_to_rank ( n, s, index )

!*****************************************************************************80
!
!! subset_to_rank() returns the rank of a subset().
!
!  Discussion:
!
!    The subset is described by a binary vector of n digits.
!    The units digit is the last one.
!    Reading from right to left, we add selected powers of 2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n

  integer i
  integer index
  integer power2
  integer s(n)

  index = 0
  power2 = 1
  do i = n, 1, -1
    index = index + s(i) * power2
    power2 = power2 * 2
  end do

  return
end

