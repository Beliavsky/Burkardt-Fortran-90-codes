program main

!*****************************************************************************80
!
!! c8poly_test() tests c8poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'c8poly_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test c8poly().'

  call c8poly_print_test ( )
  call roots_to_c8poly_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'c8poly_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine c8poly_print_test ( )

!*****************************************************************************80
!
!! c8poly_print_test() tests c8poly_print().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: m = 5

  complex ( kind = ck ), dimension ( 0 : m ) :: c = (/ &
    (  1.0D+00,  2.0D+00 ), &
    ( -3.0D+00,  4.0D+00 ), &
    (  5.0D+00,  6.0D+00 ), &
    (  7.0D+00,  8.0D+00 ), &
    (  9.0D+00, 10.0D+00 ), &
    ( 11.0D+00, 12.0D+00 ) /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'c8poly_print_test():'
  write ( *, '(a)' ) '  c8poly_print() prints a C8POLY.'

  call c8poly_print ( m, c, '  The C8POLY:' )

  return
end
subroutine roots_to_c8poly_test ( )

!*****************************************************************************80
!
!! roots_to_r8poly_test() tests roots_to_c8poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 5

  complex ( kind = ck ) c(0:n)
  complex ( kind = ck ), dimension ( n ) :: x = (/ &
    1.0D+00, -4.0D+00, 3.0D+00, 0.0D+00, 3.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'roots_to_c8poly_test():'
  write ( *, '(a)' ) '  roots_to_c8poly() is given N complex roots,'
  write ( *, '(a)' ) '  and constructs the coefficient vector'
  write ( *, '(a)' ) '  of the corresponding polynomial.'

  call c8vec_print ( n, x, '  Roots:' )

  call roots_to_c8poly ( n, x, c )

  call c8poly_print ( n, c, '  Corresponding polynomial:' )

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! c8vec_print() prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    complex ( kind = ck ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
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
!    15 August 2021
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
