program main

!*****************************************************************************80
!
!! l4lib_test() tests l4lib().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'l4lib_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test l4lib().'

  call i4_to_l4_test ( )
  call i4_to_l4vec_test ( )

  call l4_to_i4_test ( )
  call l4_to_s_test ( )
  call l4_uniform_test ( )
  call l4_xor_test ( )

  call l4mat_print_test ( )
  call l4mat_print_some_test ( )
  call l4mat_transpose_print_test ( )
  call l4mat_transpose_print_some_test ( )
  call l4mat_uniform_test ( )

  call l4vec_next_test ( )
  call l4vec_print_test ( )
  call l4vec_uniform_test ( )

  call s_to_l4_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'l4lib_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine i4_to_l4_test ( )

!*****************************************************************************80
!
!! i4_to_l4_test() tests i4_to_l4(). 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i4
  logical l4

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'I4_TO_L4_TEST'
  write ( *, '(a)' ) '  I4_TO_L4 converts an I4 to an L4.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  I4   L4'
  write ( *, '(a)' ) ''

  do i4 = -5, +5

    call i4_to_l4 ( i4, l4 )
    write ( *, '(2x,i2,2x,l1)' ) i4, l4

  end do

  return
end
subroutine i4_to_l4vec_test ( )

!*****************************************************************************80
!
!! I4_TO_L4VEC_TEST tests I4_TO_L4VEC. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 8

  integer i4
  logical l4vec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'I4_TO_L4VEC_TEST'
  write ( *, '(a)' ) '  I4_TO_L4VEC converts an I4 to an L4VEC.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  I4   L4VEC'
  write ( *, '(a)' ) ''

  do i4 = 0, 10

    call i4_to_l4vec ( i4, n, l4vec )
    write ( *, '(2x,i2,4x,8(1x,1l))' ) i4, l4vec(1:n)

  end do

  return
end
subroutine l4_to_i4_test ( )

!*****************************************************************************80
!
!! L4_TO_I4_TEST tests L4_TO_I4. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i4
  logical l4

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'L4_TO_I4_TEST'
  write ( *, '(a)' ) '  L4_TO_I4 converts an L4 to an I4.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  L4   I4'
  write ( *, '(a)' ) ''

  l4 = .false.
  call l4_to_i4 ( l4, i4 )
  write ( *, '(2x,l1,2x,i1)' ) l4, i4

  l4 = .true.
  call l4_to_i4 ( l4, i4 )
  write ( *, '(2x,l1,2x,i1)' ) l4, i4

  return
end
subroutine l4_to_s_test ( )

!*****************************************************************************80
!
!! L4_TO_S_TEST tests L4_TO_S. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!q
!    01 August 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  logical l4
  character ( len = 5 ) s

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'L4_TO_S_TEST'
  write ( *, '(a)' ) '  L4_TO_S converts an L4 to a string.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  L4   S'
  write ( *, '(a)' ) ''

  l4 = .false.
  call l4_to_s ( l4, s )
  write ( *, '(2x,l1,2x,a)' ) l4, s

  l4 = .true.
  call l4_to_s ( l4, s )
  write ( *, '(2x,l1,2x,a)' ) l4, s

  return
end
subroutine l4_uniform_test ( )

!*****************************************************************************80
!
!! L4_UNIFORM_TEST tests L4_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  logical l4_uniform
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4_UNIFORM computes pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,l1)' ) i, l4_uniform ( seed )
  end do

  return
end
subroutine l4_xor_test ( )

!*****************************************************************************80
!
!! L4_XOR_TEST tests L4_XOR. 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer j
  logical l1
  logical l2
  logical l4
  logical l4_xor

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'L4_XOR_TEST'
  write ( *, '(a)' ) '  L4_XOR computes the exclusive OR of two L4''s'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  L1  L2  L4_XOR(L1,L2)'
  write ( *, '(a)' ) ''

  do j = 0, 1
    l1 = ( j == 1 )
    do i = 0, 1
      l2 = ( i == 1 )
      l4 = l4_xor ( l1, l2 );
      write ( *, '(3x,1l,3x,1l,4x,1l)' ) l1, l2, l4
    end do
  end do

  return
end
subroutine l4mat_print_test ( )

!*****************************************************************************80
!
!! L4MAT_PRINT_TEST tests L4MAT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 November 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 20
  integer, parameter :: n = 50

  logical a(m,n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_PRINT_TEST'
  write ( *, '(a)' ) '  L4MAT_PRINT prints an L4MAT.'

  do i = 1, m
    do j = 1, n
      a(i,j) = ( mod ( i, j ) == 0 )
    end do
  end do

  call l4mat_print ( m, n, a, '  A(I,J) = I is divisible by J' )

  return
end
subroutine l4mat_print_some_test ( )

!*****************************************************************************80
!
!! L4MAT_PRINT_SOME_TEST tests L4MAT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 20
  integer, parameter :: n = 50

  logical a(m,n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  L4MAT_PRINT_SOME prints some of an L4MAT.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Here, our matrix is 20x50, but we print rows 5:15, columns 1:5'

  do i = 1, m
    do j = 1, n
      a(i,j) = ( mod ( i, j ) == 0 )
    end do
  end do

  call l4mat_print_some ( m, n, a, 5, 1, 15, 5, '  A(I,J) = I is divisible by J' )

  return
end
subroutine l4mat_transpose_print_test ( )

!*****************************************************************************80
!
!! L4MAT_TRANSPOSE_PRINT_TEST tests L4MAT_TRANSPOSE_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 20
  integer, parameter :: n = 50

  logical a(m,n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_TRANSPOSE_PRINT_TEST'
  write ( *, '(a)' ) '  L4MAT_TRANSPOSE_PRINT prints the transpose of an L4MAT.'

  do i = 1, m
    do j = 1, n
      a(i,j) = ( mod ( i, j ) == 0 )
    end do
  end do

  call l4mat_transpose_print ( m, n, a, '  A(I,J) = I is divisible by J' )

  return
end
subroutine l4mat_transpose_print_some_test ( )

!*****************************************************************************80
!
!! L4MAT_TRANSPOSE_PRINT_SOME_TEST tests L4MAT_TRANSPOSE_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 20
  integer, parameter :: n = 50

  logical a(m,n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_TRANSPOSE_PRINT_SOME_TEST'
  write ( *, '(a)' ) &
    '  L4MAT_TRANSPOSE_PRINT_SOME prints some of an L4MAT, transposed.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  Here, our matrix is 20x50, but we print rows 5:15, columns 1:5'

  do i = 1, m
    do j = 1, n
      a(i,j) = ( mod ( i, j ) == 0 )
    end do
  end do

  call l4mat_transpose_print_some ( m, n, a, 5, 1, 15, 5, &
    '  A(I,J) = I is divisible by J' )

  return
end
subroutine l4mat_uniform_test ( )

!*****************************************************************************80
!
!! L4MAT_UNIFORM_TEST tests L4MAT_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  logical l(m,n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4MAT_UNIFORM computes a vector of'
  write ( *, '(a)' ) '  pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call l4mat_uniform ( m, n, seed, l )

  call l4mat_print ( m, n, l, '  Uniform L4MAT:' )

  return
end
subroutine l4vec_next_test ( )

!*****************************************************************************80
!
!! L4VEC_NEXT_TEST tests L4VEC_NEXT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 3

  logical l4vec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'L4VEC_NEXT_TEST'
  write ( *, '(a)' ) '  L4VEC_NEXT generates logical vectors.'
  write ( *, '(a)' ) ''
 
  l4vec(1:n) = .false.

  do

    write ( *, '(2x,3l1)' ) l4vec(1:n)

    if ( all ( l4vec(1:n) ) ) then
      exit
    end if

    call l4vec_next ( n, l4vec )
 
  end do

  return
end
subroutine l4vec_print_test ( )

!*****************************************************************************80
!
!! L4VEC_PRINT_TEST tests L4VEC_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 November 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 20

  integer i
  integer j
  logical lvec(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'L4VEC_PRINT_TEST'
  write ( *, '(a)' ) '  L4VEC_PRINT prints an L4VEC.'
  
  lvec(1:n) = .true.
  
  lvec(1) = .false.

  do i = 1, n
    do j = 2, i - 1
      if ( mod ( i, j ) == 0 ) then
        lvec(i) = .false.
        exit
      end if
    end do
  end do

  call l4vec_print ( n, lvec, '  Is I Prime?:' )

  return
end
subroutine l4vec_uniform_test ( )

!*****************************************************************************80
!
!! L4VEC_UNIFORM_TEST tests L4VEC_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  logical l(n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4VEC_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4VEC_UNIFORM computes a vector of'
  write ( *, '(a)' ) '  pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call l4vec_uniform ( n, seed, l )

  call l4vec_print ( n, l, '  Uniform L4VEC:' )

  return
end
subroutine s_to_l4_test ( )

!*****************************************************************************80
!
!! S_TO_L4_TEST tests S_TO_L4.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 10

  logical l
  character ( len = 10 ) s
  logical s_to_l4
  character ( len = 10 ) string(test_num)
  integer test

  string(1) = '0'
  string(2) = 'F'
  string(3) = 'f'
  string(4) = '1'
  string(5) = 'T'
  string(6) = 't'
  string(7) = '  0'
  string(8) = '  1  0'
  string(9) = '  01'
  string(10) = '  Talse'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'S_TO_L4_TEST'
  write ( *, '(a)' ) '  S_TO_L4 reads logical data from a string.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  S   L'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    s = string(test)
    l = s_to_l4 ( s )
    write ( *, '(2x,a10,2x,l1,4x,i2,4x,i2)' ) s, l
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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

