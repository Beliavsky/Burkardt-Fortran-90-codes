subroutine knapsack_rational ( n, mass_limit, p, w, x, mass, profit )

!*****************************************************************************80
!
!! knapsack_rational() solves the rational knapsack problem.
!
!  Discussion:
!
!    The rational knapsack problem is a generalization of the 0/1 knapsack
!    problem.  It is mainly used to derive a bounding function for the
!    0/1 knapsack problem.
!
!    The 0/1 knapsack problem is as follows:
!
!      Given:
!        a set of N objects,
!        a profit P(I) and weight W(I) associated with each object,
!        and a weight limit MASS_LIMIT,
!      Determine:
!        a set of choices X(I) which are 0 or 1, that maximizes the profit
!          P = Sum ( 1 <= I <= N ) P(I) * X(I)
!        subject to the constraint
!          Sum ( 1 <= I <= N ) W(I) * X(I) <= MASS_LIMIT.
!
!    By contrast, the rational knapsack problem allows the values X(I)
!    to be any value between 0 and 1.  A solution for the rational knapsack
!    problem is known.  Arrange the objects in order of their "profit density"
!    ratios P(I)/W(I), and then take in order as many of these as you can.
!    If you still have "room" in the weight constraint, then you should
!    take the maximal fraction of the very next object, which will complete
!    your weight limit, and maximize your profit.
!
!    If should be obvious that, given the same data, a solution for
!    the rational knapsack problem will always have a profit that is
!    at least as high as for the 0/1 problem.  Since the rational knapsack
!    maximum profit is easily computed, this makes it a useful bounding
!    function.
!
!    Note that this routine assumes that the objects have already been
!    arranged in order of the "profit density".
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
!    integer N, the number of objects.
!
!    real ( kind = rk8 ) MASS_LIMIT, the weight limit of the chosen objects.
!
!    real ( kind = rk8 ) P(N), the "profit" or value of each object.
!    The entries of P are assumed to be nonnegative.
!
!    real ( kind = rk8 ) W(N), the "weight" or cost of each object.
!    The entries of W are assumed to be nonnegative.
!
!  Output:
!
!    real ( kind = rk8 ) X(N), the choice function for the objects.
!    0.0, the object was not taken.
!    1.0, the object was taken.
!    R, where 0 < R < 1, a fractional amount of the object was taken.
!
!    real ( kind = rk8 ) MASS, the total mass of the objects taken.
!
!    real ( kind = rk8 ) PROFIT, the total profit of the objects taken.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk8 ) mass
  real ( kind = rk8 ) mass_limit
  real ( kind = rk8 ) p(n)
  real ( kind = rk8 ) profit
  real ( kind = rk8 ) w(n)
  real ( kind = rk8 ) x(n)

  mass = 0.0D+00
  profit = 0.0D+00

  do i = 1, n

    if ( mass_limit <= mass ) then
      x(i) = 0.0D+00
    else if ( mass + w(i) <= mass_limit ) then
      x(i) = 1.0D+00
      mass = mass + w(i)
      profit = profit + p(i)
    else
      x(i) = ( mass_limit - mass ) / w(i)
      mass = mass_limit
      profit = profit + p(i) * x(i)
    end if

  end do

  return
end
subroutine knapsack_reorder ( n, p, w )

!*****************************************************************************80
!
!! knapsack_reorder() reorders the knapsack data by "profit density".
!
!  Discussion:
!
!    This routine must be called to rearrange the data before calling
!    routines that handle a knapsack problem.
!
!    The "profit density" for object I is defined as P(I)/W(I).
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
!    integer N, the number of objects.
!
!    real ( kind = rk8 ) P(N), the profit or value of each object.
!
!    real ( kind = rk8 ) W(N), the weight or cost of each object.
!
!  Output:
!
!    real ( kind = rk8 ) P(N), the profit or value of each object after sorting.
!
!    real ( kind = rk8 ) W(N), the weight or cost of each object after sorting.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n

  integer i
  integer j
  real ( kind = rk8 ) p(n)
  real ( kind = rk8 ) t
  real ( kind = rk8 ) w(n)
!
!  Rearrange the objects in order of "profit density".
!
  do i = 1, n
    do j = i + 1, n
      if ( p(i) * w(j) < p(j) * w(i) ) then

        t    = p(i)
        p(i) = p(j)
        p(j) = t

        t    = w(i)
        w(i) = w(j)
        w(j) = t

      end if
    end do
  end do

  return
end
