subroutine counter ( )

!*****************************************************************************80
!
!! counter() simply reports how many times it has been called.
!
!  Discussion:
!
!    If you can get a function like this to work, you have made a good
!    start at implementing persistent data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 May 2021
!
!  Author:
!
!    John Burkardt
!
!  Persistent:
!
!    integer CALLS: the number of times this function was called.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, save :: calls = 0

  calls = calls + 1

  write ( *, '(a,i4,a)' ) 'counter() has been called ', calls, ' times.'

  return
end

