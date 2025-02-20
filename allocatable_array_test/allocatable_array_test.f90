program main

!*****************************************************************************80
!
!! allocatable_array_test() tests the use of allocatable arrays.
!
!  Discussion:
!
!    allocatable_array_test attempts to allocate an allocatable array in a
!    subroutine, and use it in the main program.
!
!    The handling of allocatable arrays as dummy arguments and function values
!    is discussed in the technical report ISO/IEC TR:15581.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none
!
!  In order that the MAIN program can call TEST01 with an allocatable array as
!  a dummy argument, MAIN must have an INTERFACE statement that "explains"
!  what is going on.
!
  interface
    subroutine test01 ( n, a )
      integer, allocatable :: a(:)
      integer n
    end subroutine
  end interface

  integer, allocatable :: a(:)
  integer i
  integer n

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'allocatable_array_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Show how a main program can declare an allocatable array'
  write ( *, '(a)' ) '  which can then be allocated and set in a subroutine,'
  write ( *, '(a)' ) '  and then returned to the main program.'

  call test01 ( n, a )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The array size is N = ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The array contents:'
  write ( *, '(a)' ) ' '
  do i = 1, n
    write ( *, '(2x,i8,2x,i8)' ) i, a(i)
  end do
!
!  Free memory.
!
  deallocate ( a )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'allocatable_array_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( n, a )

!*****************************************************************************80
!
!! test01() allocates and assigns the array.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    integer N: the size of the array.
!
!    integer A(N): an array allocated and evaluated by this routine.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0d+00 )

  integer, allocatable :: a(:)
  integer i
  integer n
  real ( kind = rk ) x

  call random_number ( harvest = x )

  n = int ( 10.0_rk * x ) + 1

  allocate ( a(1:n) )

  do i = 1, n
    a(i) = 101 * i
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
