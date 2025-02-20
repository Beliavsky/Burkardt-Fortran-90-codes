subroutine i4vec_sort_insert_d ( n, a )

!*****************************************************************************80
!
!! i4vec_sort_insert_d() uses a descending insertion sort on an I4VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Input:
!
!    integer N, the number of items in the vector.
!    N must be positive.
!
!    integer A(N): data to be sorted.
!
!  Output:
!
!    integer A(N): the entries of A have been sorted in descending order.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer j
  integer x

  do i = 2, n

    x = a(i)

    j = i - 1

    do while ( 1 <= j )

      if ( x <= a(j) ) then
        exit
      end if

      a(j+1) = a(j)
      j = j - 1

    end do

    a(j+1) = x

  end do

  return
end
subroutine subset_sum_swap ( n, a, sum_desired, index, sum_achieved )

!*****************************************************************************80
!
!! subset_sum_swap() seeks a solution of the subset sum problem by swapping.
!
!  Discussion:
!
!    Given a collection of N not necessarily distinct positive integers A(I),
!    and a positive integer SUM_DESIRED, select a subset of the values so that
!    their sum is as close as possible to SUM_DESIRED without exceeding it.
!
!  Algorithm:
!
!    Start with no values selected, and SUM_ACHIEVED = 0.
!
!    Consider each element A(I):
!
!      If A(I) is not selected and SUM_ACHIEVED + A(I) <= SUM_DESIRED,
!        select A(I).
!
!      If A(I) is still not selected, and there is a selected A(J)
!      such that SUM_GOT < SUM_ACHIEVED + A(I) - A(J),
!        select A(I) and deselect A(J).
!
!      If no items were selected on this sweep,
!        exit.
!      Otherwise,
!        repeat the search.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 March 2001
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Donald Kreher, Douglas Simpson,
!    Combinatorial Algorithms,
!    CRC Press, 1998,
!    ISBN: 0-8493-3988-X,
!    LC: QA164.K73.
!
!  Input:
!
!    integer N, the number of values.  N must be positive.
!
!    integer A(N), a collection of positive values.
!
!    integer SUM_DESIRED, the desired sum.
!
!  Output:
!
!    integer A(N): A has been sorted into descending order.
!
!    integer INDEX(N); INDEX(I) is 1 if A(I) is part of the
!    sum, and 0 otherwise.
!
!    integer SUM_ACHIEVED, the sum of the selected
!    elements.
!
  implicit none

  integer n

  integer a(n)
  integer i
  integer index(n)
  integer j
  integer nmove
  integer sum_achieved
  integer sum_desired
!
!  Initialize.
!
  sum_achieved = 0
  index(1:n) = 0
!
!  Sort into descending order.
!
  call i4vec_sort_insert_d ( n, a )

  do

    nmove = 0

    do i = 1, n

      if ( index(i) == 0 ) then

        if ( sum_achieved + a(i) <= sum_desired ) then
          index(i) = 1
          sum_achieved = sum_achieved + a(i)
          nmove = nmove + 1
          cycle
        end if

      end if

      if ( index(i) == 0 ) then

        do j = 1, n

          if ( index(j) == 1 ) then

            if ( sum_achieved < sum_achieved + a(i) - a(j) .and. &
              sum_achieved + a(i) - a(j) <= sum_desired ) then
              index(j) = 0
              index(i) = 1
              nmove = nmove + 2
              sum_achieved = sum_achieved + a(i) - a(j)
              exit
            end if

          end if

        end do

      end if

    end do

    if ( nmove <= 0 ) then
      exit
    end if

  end do

  return
end
