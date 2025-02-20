program main

!*****************************************************************************80
!
!! matrix_chain_brute_test() tests matrix_chain_brute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 June 2024
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
  integer i4_factorial
  integer n_dims
  integer n_mats
  integer n_mults
  integer, allocatable :: p(:)
  integer parens
  integer perms
  integer test

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'matrix_chain_brute_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test matrix_chain_brute().'

  do test = 1, 10

    if ( test == 1 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 40, 20, 30, 10, 30 /)
    else if ( test == 2 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 1, 2, 3, 4, 3 /)
    else if ( test == 3 ) then
      n_dims = 3
      allocate ( dims(1:n_dims) )
      dims = (/ 10, 20, 30 /)
    else if ( test == 4 ) then
      n_dims = 4
      allocate ( dims(1:n_dims) )
      dims = (/ 10, 30, 5, 60 /)
    else if ( test == 5 ) then
      n_dims = 2
      allocate ( dims(1:n_dims) )
      dims = (/ 10, 20 /)
    else if ( test == 6 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 40, 20, 0, 10, 30 /)
    else if ( test == 7 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 1, 100, 1, 100, 1 /)
    else if ( test == 8 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 100, 50, 1, 50, 100 /)
    else if ( test == 9 ) then
      n_dims = 5
      allocate ( dims(1:n_dims) )
      dims = (/ 1, 50, 100, 50, 1 /)
    else if ( test == 10 ) then
      n_dims = 6
      allocate ( dims(1:n_dims) )
      dims = (/ 4, 10, 3, 12, 20, 7 /)
    end if

    n_mats = n_dims - 1
    n_mults = n_mats - 1

    allocate ( p(1:n_mults) )

    write ( *, '(a)' ) ''
    write ( *, '(a,i2)' ) '  Test #', test
    write ( *, '(a,i4)' ) '  Number of dimensions =      ', n_dims
    write ( *, '(a,i4)' ) '  Number of matrices =        ', n_mats
    write ( *, '(a,i4)' ) '  Number of multiplications = ', n_mults
    write ( *, '(a)' ) '  Matrix dimensions'
    do i = 1, n_dims
      write ( *, '(2x,i4)', advance = 'no' ) dims(i)
    end do
    write ( *, '(a)' ) ''
    parens = catalan_number ( n_mults )
    write ( *, '(a,i4)' ) &
      '  Number of possible parenthesizations is ', parens
    perms = i4_factorial ( n_mults )
    write ( *, '(a,i4)' ) &
      '  Number of possible permutations is ', perms
    call matrix_chain_brute ( n_mats, dims, cost, p )
    write ( *, '(a,i6)' ) '  Minimal cost is ', cost
    write ( *, '(a)' ) '  Multiplication order:'
    if ( n_mults < 1 ) then
      write ( *, '(a)' ) '  [ Empty ]'
    else
      do i = 1, n_mults
        write ( *, '(2x,i2)', advance = 'no' ) p(i)
      end do
      write ( *, '(a)' ) ''
    end if

    deallocate ( dims )
    deallocate ( p )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'matrix_chain_brute_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

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
