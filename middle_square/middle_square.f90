subroutine middle_square_next ( s, d, r )

!*****************************************************************************80
!
!! middle_square_next() computes the next middle square random value.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 September 2022
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Brian Hayes,
!    The Middle of the Square,
!    08 August 2022,
!    http://bit-player.org/
!
!  Input:
!
!    integer S, the current seed, of no more than 2*D digits.
!
!    integer D, HALF the number of digits.  
!    Typical values are 2, 3, or 4.
!
!  Output:
!
!    integer R, the next seed.
!
  implicit none

  integer, parameter :: i20 = selected_int_kind ( 20 )

  integer ( kind = i20 ) d
  integer ( kind = i20 ) r
  integer ( kind = i20 ) s
!
!  Square.
!
  r = s * s
!
!  Drop last two digits.
!
  r = ( r / 10**d )
!
!  Drop all but last four digits.
!
  r = mod ( r, 10**( 2 * d ) )

  return
end

