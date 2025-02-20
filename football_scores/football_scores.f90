subroutine football_scores ( n, s )

!*****************************************************************************80
!
!! football_scores() counts the ways of getting any football score from 0 to 100.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2022
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer n: the highest score to be considered.
!
!  Output:
!
!    integer ( kind = selected_int_kind ( 16 ) ) s(0:n): the number
!    of ways of achieving each score from 0 to n.
!
  implicit none

  integer, parameter :: i16 = selected_int_kind ( 16 )

  integer n

  integer i
  integer ( kind = i16 ) s(0:n)

  do i = 0, n

    s(i) = 0

    if ( i == 0 ) then
      s(i) = 1
    end if
    if ( 1 <= i ) then
      s(i) = s(i) + s(i-1)
    end if
    if ( 2 <= i ) then
      s(i) = s(i) + s(i-2)
    end if
    if ( 3 <= i ) then
      s(i) = s(i) + s(i-3)
    end if
    if ( 6 <= i ) then
      s(i) = s(i) + s(i-6)
    end if
    if ( 7 <= i ) then
      s(i) = s(i) + s(i-7)
    end if
    if ( 8 <= i ) then
      s(i) = s(i) + s(i-8)
    end if

  end do

  return
end

