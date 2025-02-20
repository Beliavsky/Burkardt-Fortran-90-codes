subroutine i4_to_digits_binary ( i, n, c )

!*****************************************************************************80
!
!! i4_to_digits_binary() produces the binary digits of an I4.
!
!  Discussion:
!
!    An I4 is an integer.
!
!  Example:
!
!     I    N     C               Binary
!    --  ---   ---         ------------
!     0    1   0                      0
!     0    2   0, 0                  00
!     1    3   1, 0, 0              100
!     2    3   0, 1, 0              010
!     3    3   1, 1, 0              011
!     4    3   0, 0, 1              100
!     8    3   0, 0, 0           (1)000
!     8    5   0, 0, 0, 1, 0      01000
!    -8    5   0, 0, 0, 1, 0  (-) 01000
!
!     0    3   0, 0, 0
!     1    3   1, 0, 0
!     2    3   0, 1, 0
!     3    3   1, 1, 0
!     4    3   0, 0, 1
!     5    3   1, 0, 1
!     6    3   0, 1, 1
!     7    3   1, 1, 1
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 December 2010
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer I, an integer to be represented.
!
!    integer N, the number of binary digits to produce.
!
!  Output:
!
!    integer C(N), the first N binary digits of I,
!    with C(1) being the units digit.
!
  implicit none

  integer n

  integer c(n)
  integer i
  integer i_copy
  integer j

  i_copy = abs ( i )

  do j = 1, n

    c(j) = mod ( i_copy, 2 )
    i_copy = i_copy / 2

  end do

  return
end
function i4vec_dot_product ( n, x, y )

!*****************************************************************************80
!
!! i4vec_dot_product() computes the dot product of two I4VEC's.
!
!  Discussion:
!
!    An I4VEC is a vector of I4's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 December 2011
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the size of the array.
!
!    integer X(N), Y(N), the arrays.
!
!  Output:
!
!    integer I4VEC_DOT_PRODUCT, the dot product of X and Y.
!
  implicit none

  integer n

  integer i4vec_dot_product
  integer x(n)
  integer y(n)

  i4vec_dot_product = dot_product ( x(1:n), y(1:n) )

  return
end
subroutine subset_sum_brute ( n, weight, target, choice )

!*****************************************************************************72
!
!! subset_sum_brute seeks a subset of a set that has a given sum.
!
!  Discussion:
!
!    This function tries to compute a target value as the sum of
!    a selected subset of a given set of weights.
!
!    This function works by brute force, that is, it tries every
!    possible subset to see if it sums to the desired value.
!
!    Given N weights, every possible selection can be described by 
!    one of the N-digit binary numbers from 0 to 2^N-1.
!
!    It is possible that there may be multiple solutions of the problem.  
!    This function will only return the first solution found.
!
!  Example:
!
!    n = 6
!    target = 22
!    w = (/ 1, 2, 4, 8, 16, 32 /)
!
!    choice = (/ 0, 1, 1, 0, 1, 0 /)
!    w(choice) = 2 + 4 + 16 = 22
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of weights.
!
!    integer WEIGHT(N), the weights.
!
!    integer TARGET, the target value.
!
!  Output:
!
!    integer CHOICE(N), contains a 1 for each
!    weight that is chosen.  If no solution was found, all entries
!    are returned as -1.
!
  implicit none

  integer n

  integer choice(n)
  integer i
  integer i_max
  integer i4vec_dot_product
  integer target
  integer w_sum
  integer weight(n)

  i_max = ( 2 ** n ) - 1

  do i = 0, i_max
!
!  Convert I to a string of binary digits.
!
    call i4_to_digits_binary ( i, n, choice )
!
!  Combine the weights whose binary digit is 1.
!
    w_sum = i4vec_dot_product ( n, choice, weight )
!
!  Return if we matched our target sum.
!
    if ( w_sum == target ) then
      return
    end if
    
  end do

  do i = 1, n
    choice(i) = -1
  end do

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

