program main

!*****************************************************************************80
!
!! matrix_chain_dynamic_test() tests matrix_chain_dynamic().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer catalan_number
  integer cost
  integer, allocatable :: dims(:)
  integer i
  integer n
  integer orderings
  integer test

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'matrix_chain_dynamic_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test matrix_chain_dynamic().'

  do test = 1, 10

    if ( test == 1 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 40, 20, 30, 10, 30 /)
    elseif ( test == 2 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 1, 2, 3, 4, 3 /)
    elseif ( test == 3 ) then
      n = 3
      allocate ( dims(1:n) )
      dims = (/ 10, 20, 30 /)
    elseif ( test == 4 ) then
      n = 4
      allocate ( dims(1:n) )
      dims = (/ 10, 30, 5, 60 /)
    elseif ( test == 5 ) then
      n = 2
      allocate ( dims(1:n) )
      dims = (/ 10, 20 /)
    elseif ( test == 6 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 40, 20, 0, 10, 30 /)
    elseif ( test == 7 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 1, 100, 1, 100, 1 /)
    elseif ( test == 8 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 100, 50, 1, 50, 100 /)
    elseif ( test == 9 ) then
      n = 5
      allocate ( dims(1:n) )
      dims = (/ 1, 50, 100, 50, 1 /)
    elseif ( test == 10 ) then
      n = 6
      allocate ( dims(1:n) )
      dims = (/ 4, 10, 3, 12, 20, 7 /)
    end if

    write ( *, '(a)' ) ''
    write ( *, '(a,i2)' ) '  Number of matrices = ', n - 1
    write ( *, '(a)', advance = 'no' ) '  Matrix dimensions'
    do i = 1, n
      write ( *, '(2x,i3)', advance = 'no' ) dims(i)
    end do
    write ( *, '(a)' ) ''
    orderings = catalan_number ( n - 2 )
    write ( *, '(a,i2)' ) '  Number of possible orderings is', orderings
    call matrix_chain_dynamic ( n - 1, dims, cost )
    write ( *, '(a,i8)' ) '  Minimal cost is ', cost

    deallocate ( dims )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'matrix_chain_dynamic_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
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

