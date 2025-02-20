program main

!*****************************************************************************80
!
!! candy_count_test() tests candy_count().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'candy_count_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test candy_count().'

  call candy_count_vector_test ( )
  call candy_count_matrix_test ( )
  call candy_count_box_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'candy_count_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine candy_count_box_test ( )

!*****************************************************************************80
!
!! candy_count_box_test() tests candy_count_box().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer c
  integer, allocatable :: counts(:)
  integer i
  integer l
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'candy_count_box_test():'
  write ( *, '(a)' ) '  candy_count_box() counts candy types in a 3D box.'
  write ( *, '(a)' ) '  There are LxMxN entries in the box A().'
  write ( *, '(a)' ) '  There are C candy types.'
  write ( *, '(a)' ) '  Candy types are assigned cyclically to matrix entries:'
  write ( *, '(a)' ) '  A(I,J,K) = mod ( i + j + k - 3, c ) + 1'
  write ( *, '(a)' ) '  Count the number of candies of each type.'
!
!  Test 1
!
  c = 4
  l = 7
  m = 10
  n = 13
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Count using candy_count_box()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   L   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_box ( c, l, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, l, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  c = 4
  l = 7
  m = 10
  n = 13
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Repeat calculation using candy_count_box_sum()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   L   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_box_sum ( c, l, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, l, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )
!
!  Test 2
!
  c = 5
  l = 12
  m = 13
  n = 19
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Count using candy_count_box()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   L   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_box ( c, l, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, l, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  c = 5
  l = 12
  m = 13
  n = 19
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Repeat calculation using candy_count_box_sum()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   L   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_box_sum ( c, l, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, l, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  return
end
subroutine candy_count_matrix_test ( )

!*****************************************************************************80
!
!! candy_count_matrix_test() tests candy_count_matrix().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer c
  integer, allocatable :: counts(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'candy_count_matrix_test():'
  write ( *, '(a)' ) '  candy_count_matrix() counts candy types in a matrix.'
  write ( *, '(a)' ) '  There are MxN entries in the matrix A().'
  write ( *, '(a)' ) '  There are C candy types.'
  write ( *, '(a)' ) '  Candy types are assigned cyclically to matrix entries:'
  write ( *, '(a)' ) '  A(I,J) = mod ( i + j - 2, c ) + 1'
  write ( *, '(a)' ) '  Count the number of candies of each type.'
!
!  Test 1
!
  c = 4
  m = 10
  n = 13
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Count using candy_count_matrix()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_matrix ( c, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  c = 4
  m = 10
  n = 13
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Repeat calculation using candy_count_matrix_sum()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_matrix_sum ( c, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )
!
!  Test 2
!
  c = 5
  m = 13
  n = 19
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Count using candy_count_matrix()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_matrix ( c, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  c = 5
  m = 13
  n = 19
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Repeat calculation using candy_count_matrix_sum()'
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   C   M   N  '
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) i
  end do
  write ( *, '(a)' ) ''
  allocate ( counts(1:c) )
  call candy_count_matrix_sum ( c, m, n, counts )
  write ( *, '(2x,i2,2x,i2,2x,i2,2x)', advance = 'no' ) c, m, n
  do i = 1, c
    write ( *, '(i4)', advance = 'no' ) counts(i)
  end do
  write ( *, '(a)' ) ''
  deallocate ( counts )

  return
end
subroutine candy_count_vector_test ( )

!*****************************************************************************80
!
!! candy_count_vector_test() tests candy_count_vector().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer c
  integer, allocatable :: counts(:)
  integer i
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'candy_count_vector_test():'
  write ( *, '(a)' ) '  candy_count_vector() counts candy types in a vector.'
  write ( *, '(a)' ) '  There are N entries in the vector A().'
  write ( *, '(a)' ) '  There are C candy types.'
  write ( *, '(a)' ) '  Candy types are assigned cyclically to vector entries:'
  write ( *, '(a)' ) '  A(I) = mod ( i - 1, c ) + 1'
  write ( *, '(a)' ) '  Count the number of candies of each type.'

  c = 4
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Count using candy_count_vector()'
  write ( *, '(a,i4)' ) '  Fix value of C = ', c
  write ( *, '(a)' ) '  Consider a range of values of N:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N    #1  #2  #3  #4'
  write ( *, '(a)' ) ''

  do n = 3, 10
    allocate ( counts(1:c) )
    call candy_count_vector ( c, n, counts )
    write ( *, '(2x,i2,2x)', advance = 'no' ) n
    do i = 1, c
      write ( *, '(i4)', advance = 'no' ) counts(i)
    end do
    write ( *, '(a)' ) ''
    deallocate ( counts )
  end do

  c = 4
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Repeat calculation, using candy_count_vector_sum()'
  write ( *, '(a,i4)' ) '  Fix value of C = ', c
  write ( *, '(a)' ) '  Consider a range of values of N:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N    #1  #2  #3  #4'
  write ( *, '(a)' ) ''

  do n = 3, 10
    allocate ( counts(1:c) )
    call candy_count_vector_sum ( c, n, counts )
    write ( *, '(2x,i2,2x)', advance = 'no' ) n
    do i = 1, c
      write ( *, '(i4)', advance = 'no' ) counts(i)
    end do
    write ( *, '(a)' ) ''
    deallocate ( counts )
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
