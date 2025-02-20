program main

!*****************************************************************************80
!
!! asa299_test() tests asa299().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa299_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa299().'

  call test01 ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa299_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests simplex_lattice_point_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 4

  integer i
  logical more
  integer, parameter :: t = 4
  integer x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  simplex_lattice_point_next() generates lattice points'
  write ( *, '(a)' ) '  in the simplex'
  write ( *, '(a)' ) '    0 <= X'
  write ( *, '(a)' ) '    sum ( X(1:N) ) <= T'
  write ( *, '(a,i8)' ) '  Here N = ', n
  write ( *, '(a,i8)' ) '  and T =  ', t
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Index        X(1)      X(2)      X(3)      X(4)'
  write ( *, '(a)' ) ' '

  more = .false.

  i = 0

  do

    call simplex_lattice_point_next ( n, t, more, x )

    i = i + 1

    write ( *, '(2x,i8,2x,4(2x,i8))' ) i, x(1:n)

    if ( .not. more )  then
      exit
    end if

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
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

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

