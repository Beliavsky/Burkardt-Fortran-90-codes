program main

!*****************************************************************************80
!
!! bvec_test() tests bvec().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test bvec().'

  call bvec_add_test ( )
  call bvec_complement2_test ( )
  call bvec_mul_test ( )
  call bvec_next_test ( )
  call bvec_next_grlex_test ( )
  call bvec_print_test ( )
  call bvec_sub_test ( )
  call bvec_to_i4_test ( )
  call bvec_uniform_test ( )
  call i4_bclr_test ( )
  call i4_bset_test ( )
  call i4_btest_test ( )
  call i4_to_bvec_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine bvec_add_test ( )

!*****************************************************************************80
!
!! bvec_add_test() tests bvec_add();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer l
  integer test
  integer, parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_add_test():'
  write ( *, '(a)' ) '  bvec_add() adds binary vectors representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J      I + J    BVEC_ADD'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform_ab ( -100, 100 )
    j = i4_uniform_ab ( -100, 100 )

    k = i + j

    call i4_to_bvec ( i, n, bvec1 )
    call i4_to_bvec ( j, n, bvec2 )

    call bvec_add ( n, bvec1, bvec2, bvec3 )
    call bvec_to_i4 ( n, bvec3, l )

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k, l

  end do

  return
end
subroutine bvec_complement2_test ( )

!*****************************************************************************80
!
!! bvec_complement2_test() tests bvec_complement2();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer bvec1(n)
  integer bvec2(n)
  integer i
  integer j
  integer i4_uniform_ab
  integer test
  integer, parameter :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_complement2_test():'
  write ( *, '(a)' ) '  bvec_complement2() returns the two''s complement'
  write ( *, '(a)' ) '  of a (signed) binary vector;'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform_ab ( -100, 100 )

    call i4_to_bvec ( i, n, bvec1 )

    call bvec_complement2 ( n, bvec1, bvec2 )

    call bvec_to_i4 ( n, bvec2, j )

    write ( *, '(a)' ) ' '
    write ( *, '(a,2x,i8)' ) '  I = ', i
    write ( *, '(a,2x,i8)' ) '  J = ', j
    call bvec_print ( n, bvec1, ' ' )
    call bvec_print ( n, bvec2, ' ' )

  end do

  return
end
subroutine bvec_mul_test ( )

!*****************************************************************************80
!
!! bvec_mul_test() tests bvec_mul();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 15

  integer bvec1(n)
  integer bvec2(n)
  integer bvec3(n)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer l
  integer test
  integer, parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_mul_test()'
  write ( *, '(a)' ) '  bvec_mul() multiplies binary vectors '
  write ( *, '(a)' ) '  representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         J        I * J  BVEC_MUL'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform_ab ( -100, 100 )
    j = i4_uniform_ab ( -100, 100 )

    k = i * j

    call i4_to_bvec ( i, n, bvec1 )
    call i4_to_bvec ( j, n, bvec2 )
    call bvec_mul ( n, bvec1, bvec2, bvec3 )
    call bvec_to_i4 ( n, bvec3, l )

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k, l

  end do

  return
end
subroutine bvec_next_test ( )

!*****************************************************************************80
!
!! bvec_next_test() tests bvec_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 4
 
  integer b(n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_next_test():'
  write ( *, '(a)' ) '  bvec_next() computes the "next" BVEC.'
  write ( *, '(a)' ) ''

  b(1:n) = 0

  do i = 0, 16
    call bvec_print ( n, b, '' )
    call bvec_next ( n, b )
  end do

  return
end
subroutine bvec_next_grlex_test ( )

!*****************************************************************************80
!
!! bvec_next_grlex_test() tests bvec_next_grlex().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 4
 
  integer b(n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_next_grlex_test():'
  write ( *, '(a)' ) '  bvec_next_grlex() computes binary vectors in GRLEX order.'
  write ( *, '(a)' ) ''

  b(1:n) = 0

  do i = 0, 16
    write ( *, '(2x,i2,a)', advance = 'no' ) i, ':  '
    do j = 1, n
      write ( *, '(i1)', advance = 'no' ) b(j)
    end do
    write ( *, '(a)' ) ''
    call bvec_next_grlex ( n, b )
  end do

  return
end
subroutine bvec_print_test ( )

!*****************************************************************************80
!
!! bvec_print_test() tests bvec_print().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer, dimension ( n ) :: bvec = (/ &
    1, 0, 0, 1, 0, 1, 1, 1, 0, 0 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_print_test()'
  write ( *, '(a)' ) '  bvec_print() prints a binary vector.'

  call bvec_print ( n, bvec, '  BVEC:' )

  return
end
subroutine bvec_sub_test ( )

!*****************************************************************************80
!
!! bvec_sub_test() tests bvec_sub();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer bvec1(n)
  integer bvec2(n)
  integer bvec4(n)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer l
  integer test
  integer, parameter :: test_num = 10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_sub_test()'
  write ( *, '(a)' ) '  bvec_sub() subtracts binary vectors representing integers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        I        J        I - J    BVEC_SUB'
  write ( *, '(a)' ) ' '

  do test = 1, test_num
    
    i = i4_uniform_ab ( -100, 100 )
    j = i4_uniform_ab ( -100, 100 )

    k = i - j

    call i4_to_bvec ( i, n, bvec1 )
    call i4_to_bvec ( j, n, bvec2 )
    call bvec_sub ( n, bvec1, bvec2, bvec4 )
    call bvec_to_i4 ( n, bvec4, l )

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k, l

  end do

  return
end
subroutine bvec_to_i4_test ( )

!*****************************************************************************80
!
!! bvec_to_i4_test() tests bvec_to_i4();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer bvec(n)
  integer i
  integer i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_to_i4_test():'
  write ( *, '(a)' ) '  bvec_to_i4() converts a signed binary vector'
  write ( *, '(a)' ) '  to an integer;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I --> BVEC  -->  I'
  write ( *, '(a)' ) ' '
  do i = -3, 10
    call i4_to_bvec ( i, n, bvec )
    call bvec_to_i4 ( n, bvec, i2 )
    write ( *, '(2x,i3,2x,10i1,2x,i3)' ) i, bvec(1:n), i2
  end do

  return
end
subroutine bvec_uniform_test ( )

!*****************************************************************************80
!
!! bvec_uniform_test() tests bvec_uniform().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  logical b(n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bvec_uniform_test():'
  write ( *, '(a)' ) '  bvec_uniform() computes a binary vector.'
  write ( *, '(a)' ) ''
  do i = 1, 10
    call bvec_uniform ( n, b )
    call bvec_print ( n, b, '' )
  end do

  return
end
subroutine i4_bclr_test ( )

!*****************************************************************************80
!
!! i4_bclr_test() tests i4_bclr().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 2

  integer i4
  integer, dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  integer i4_bclr
  integer ivec(0:31)
  integer j1
  integer j2
  integer pos
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_bclr_test():'
  write ( *, '(a)' ) '  i4_bclr() sets a given bit to 0.'
  write ( *, '(a)' ) '  ibclr() is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Working on I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit       I4_BCLR         IBCLR'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_bclr ( i4, pos )
      j2 = ibclr ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,i12,2x,i12)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine i4_bset_test ( )

!*****************************************************************************80
!
!! i4_bset_test() tests i4_bset().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 2

  integer i4
  integer, dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  integer i4_bset
  integer ivec(0:31)
  integer j1
  integer j2
  integer pos
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_bset_test():'
  write ( *, '(a)' ) '  i4_bset() sets a given bit to 0.'
  write ( *, '(a)' ) '  ibset() is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Working on I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit       I4_BSET         IBSET'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_bset ( i4, pos )
      j2 = ibset ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,i12,2x,i12)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine i4_btest_test ( )

!*****************************************************************************80
!
!! i4_btest_test() tests i4_btest().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 2

  integer i4
  integer, dimension ( test_num ) :: i4_test = (/ &
    101, -31 /)
  logical i4_btest
  integer ivec(0:31)
  logical j1
  logical j2
  integer pos
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_btest_test():'
  write ( *, '(a)' ) '  i4_btest() reports whether a given bit is 0 or 1.'
  write ( *, '(a)' ) '  btest() is a FORTRAN90 function which does the same.'

  do test = 1, test_num

    i4 = i4_test(test)

    call i4_to_bvec ( i4, 32, ivec )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Analyze the integer I4 = ', i4
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '       Pos     Digit  I4_BTEST     BTEST'
    write ( *, '(a)' ) ' '

    do pos = 0, 31
  
      j1 = i4_btest ( i4, pos )
      j2 = btest ( i4, pos )

      write ( *, '(2x,i8,2x,i8,2x,7x,l1,2x,7x,l1)' ) pos, ivec(pos), j1, j2

    end do

  end do

  return
end
subroutine i4_to_bvec_test ( )

!*****************************************************************************80
!
!! i4_to_bvec_test() tests i4_to_bvec();
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  integer bvec(n)
  integer i
  integer i2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_to_bvec_test():'
  write ( *, '(a)' ) '  i4_to_bvec() converts an integer to a '
  write ( *, '(a)' ) '  signed binary vector;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I --> BVEC  -->  I'
  write ( *, '(a)' ) ' '
  do i = -3, 10
    call i4_to_bvec ( i, n, bvec )
    call bvec_to_i4 ( n, bvec, i2 )
    write ( *, '(2x,i3,2x,10i1,2x,i3)' ) i, bvec(1:n), i2
  end do

  return
end

