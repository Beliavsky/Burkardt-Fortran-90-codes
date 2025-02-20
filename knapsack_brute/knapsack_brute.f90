subroutine knapsack_brute ( n, v, w, k, vmax, wmax, smax )

!*****************************************************************************80
!
!! knapsack_brute() seeks a solution of the knapsack problem.
!
!  Discussion:
!
!    N valuable items are available, each with given value V and weight W.
!
!    A thief's knapsack can carry no more than K pounds.  
!
!    The thief seeks a selection S of items to carry in the knapsack 
!    of maximum total value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the number of items.
!
!    integer v[n], w[n]: the value and weight of each item.
!
!    integer k: the maximum weight capacity of the knapsack.
!
!  Output:
!
!    integer vmax: the total value of the items in the knapsack.
!
!    integer wmax: the total weight of the items in the knapsack.
!
!    integer smax[n]: a vector of 0's and 1's, indicating which items
!    were selected.
!
  implicit none

  integer n

  integer k
  integer s_test(n)
  integer smax(n)
  integer v(n)
  integer v_test
  integer vmax
  integer w(n)
  integer w_test
  integer wmax

  vmax = 0
  wmax = 0
  smax(1:n) = 0

  s_test(1:n) = 0
  
  do

    w_test = dot_product ( s_test, w )

    if ( w_test <= k ) then

      v_test = dot_product ( s_test, v )

      if ( vmax < v_test ) then
        smax(1:n) = s_test(1:n)
        vmax = v_test
        wmax = w_test
      end if

    end if

    call subset_next ( n, s_test )

    if ( sum ( s_test ) == 0 ) then
      exit
    end if

  end do

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
!    25 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the size of the subset.
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

