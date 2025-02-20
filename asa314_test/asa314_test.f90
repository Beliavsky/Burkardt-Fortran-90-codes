program main

!*****************************************************************************80
!
!! asa314_test() tests asa314().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Payne,
!    Inversion of matrices with contents subject to modulo arithmetic,
!    Applied Statistics,
!    Volume 46, Number 2, 1997, pages 295-298.
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ASA314_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ASA314().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa314_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests INVMOD.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Roger Payne,
!    Inversion of matrices with contents subject to modulo arithmetic,
!    Applied Statistics,
!    Volume 46, Number 2, 1997, pages 295-298.
!
  implicit none

  integer, parameter :: nrow = 3

  integer cmod(nrow)
  integer i
  integer ifault
  integer imat(nrow,nrow)
  integer, dimension(nrow,nrow) :: jmat = reshape ( (/ &
    1, 0, 0, 2, 1, 0, 1, 0, 1 /), (/ nrow, nrow /) )
  integer, dimension(nrow,nrow) :: mat = reshape ( (/ &
    1, 0, 0, 1, 1, 0, 2, 0, 1 /), (/ nrow, nrow /) )
  integer rmod(nrow)

  do i = 1, nrow
    cmod(i) = 3
    rmod(i) = 3
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  INVMOD computes the inverse of a matrix'
  write ( *, '(a)' ) '  whose elements are subject to modulo arithmetic.'

  call i4mat_print ( nrow, nrow, mat, '  The matrix to be inverted:' )

  call invmod ( mat, imat, rmod, cmod, nrow, ifault )

  call i4mat_print ( nrow, nrow, imat, '  The computed inverse:' )

  call i4mat_print ( nrow, nrow, jmat, '  The correct inverse:' )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
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

