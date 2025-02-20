program main

!*****************************************************************************80
!
!! r8bb_test() tests r8bb().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8bb_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8bb().'

  call r8bb_add_test ( )
  call r8bb_dif2_test ( )
  call r8bb_fa_test ( )
  call r8bb_get_test ( )
  call r8bb_indicator_test ( )
  call r8bb_mtv_test ( )
  call r8bb_mv_test ( )
  call r8bb_print_test ( )
  call r8bb_print_some_test ( )
  call r8bb_random_test ( )
  call r8bb_set_test ( )
  call r8bb_sl_test ( )
  call r8bb_to_r8ge_test ( )
  call r8bb_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8bb_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8bb_add_test ( )

!*****************************************************************************80
!
!! R8BB_ADD_TEST tests R8BB_ADD.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 0
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer j
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_ADD_TEST'
  write ( *, '(a)' ) '  R8BB_ADD adds a value to elements of an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Initialize matrix to indicator matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )
!
!  Print initial matrix.
!
  call r8bb_print ( n1, n2, ml, mu, a, '  Matrix before additions:' )
!
!  Add 100 to band diagonal.
!
  do i = 1, n1
    j = i
    value = 100.0D+00
    call r8bb_add ( n1, n2, ml, mu, a, i, j, value )
  end do
!
!  Add 200 to right border.
!
  do i = 1, n1
    do j = n1 + 1, n1 + n2
      value = 200.0D+00
      call r8bb_add ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Add 400 to offdiagonals in lower right dense matrix.
!
  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      if ( i /= j ) then
        value = 400.0D+00
        call r8bb_add ( n1, n2, ml, mu, a, i, j, value )
      end if
    end do
  end do

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB matrix after additions:' )

  return
end
subroutine r8bb_dif2_test ( )

!*****************************************************************************80
!
!! R8BB_DIF2_TEST tests R8BB_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_DIF2_TEST'
  write ( *, '(a)' ) '  R8BB_DIF2 sets up an R8BB second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_dif2 ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB second difference matrix:' )

  return
end
subroutine r8bb_fa_test ( )

!*****************************************************************************80
!
!! R8BB_FA_TEST tests R8BB_FA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  integer info
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_FA_TEST'
  write ( *, '(a)' ) '  R8BB_FA factors an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The border-banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n1 + n2, x )
!
!  Compute the corresponding right hand side.
!
  call r8bb_mv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n1+n2, b, '  The right hand side vector:' )
!
!  Factor the matrix.
!
  call r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

  call r8bb_print ( n1, n2, ml, mu, a, '  The FACTORED border-banded matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8BB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8bb_sl ( n1, n2, ml, mu, a, pivot, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8bb_get_test ( )

!*****************************************************************************80
!
!! R8BB_GET_TEST tests R8BB_GET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 3
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 0
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer seed
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_GET_TEST'
  write ( *, '(a)' ) '  R8BB_GET gets a value of an element of an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set matrix to indicator matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )
!
!  Print matrix.
!
  call r8bb_print ( n1, n2, ml, mu, a, '  The matrix to be queried:' )
!
!  Request random entries.
!
  seed = 123456789

  write ( *, '(a)' ) ''
  do k = 1, 10
    i = i4_uniform_ab ( 1, n1 + n2, seed )
    j = i4_uniform_ab ( 1, n1 + n2, seed )
    call r8bb_get ( n1, n2, ml, mu, a, i, j, value )
    write ( *, '(a,i2,a,i2,a,g14.6)' ) '  A(', i, ',', j, ') = ', value
  end do

  return
end
subroutine r8bb_indicator_test ( )

!*****************************************************************************80
!
!! R8BB_INDICATOR_TEST tests R8BB_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8BB_INDICATOR sets up an R8BB indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB indicator matrix:' )

  return
end
subroutine r8bb_mtv_test ( )

!*****************************************************************************80
!
!! R8BB_MTV_TEST tests R8BB_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_MTV_TEST'
  write ( *, '(a)' ) '  R8BB_MTV computes b=A''*x, where A is an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8bb_mtv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8bb_mv_test ( )

!*****************************************************************************80
!
!! R8BB_MV_TEST tests R8BB_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_MV_TEST'
  write ( *, '(a)' ) '  R8BB_MV computes b=A*x, where A is an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8bb_mv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x:' )

  return
end
subroutine r8bb_print_test ( )

!*****************************************************************************80
!
!! R8BB_PRINT_TEST tests R8BB_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_PRINT_TEST'
  write ( *, '(a)' ) '  R8BB_PRINT prints an R8BB indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB indicator matrix:' )

  return
end
subroutine r8bb_print_some_test ( )

!*****************************************************************************80
!
!! R8BB_PRINT_SOME_TEST tests R8BB_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8BB_PRINT_SOME prints some of an R8BB indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print_some ( n1, n2, ml, mu, a, 7, 7, 8, 8, &
    '  The Lower Right Block:' )

  return
end
subroutine r8bb_random_test ( )

!*****************************************************************************80
!
!! R8BB_RANDOM_TEST tests R8BB_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_RANDOM_TEST'
  write ( *, '(a)' ) '  R8BB_RANDOM randomizes an R8BB matrix;'
  write ( *, '(a)' ) '  R8BB_SL solves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The border-banded matrix:' )

  return
end
subroutine r8bb_set_test ( )

!*****************************************************************************80
!
!! R8BB_SET_TEST tests R8BB_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 4
  integer, parameter :: n2 = 1
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer j
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_SET_TEST'
  write ( *, '(a)' ) '  R8BB_SET sets elements of an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Initialize matrix to zero.
!
  call r8bb_zeros ( n1, n2, ml, mu, a )
!
!  Fill in band matrix.
!
  do i = 1, n1
    do j = 1, n1
      if ( i - ml <= j .and. j <= i + mu ) then
        value = real ( 10 * i + j, kind = rk )
        call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
      end if 
    end do
  end do
!
!  Fill in right border vector.
!
  do i = 1, n1
    do j = n1 + 1, n1 + n2
      value = real ( 10 * i + j, kind = rk )
      call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Fill in lower border vector.
!
  do i = n1 + 1, n1 + n2
    do j = 1, n1
      value = real ( 10 * i + j, kind = rk )
      call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Fill in lower right dense matrix.
!
  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      value = real ( 10 * i + j, kind = rk )
      call r8bb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB matrix:' )

  return
end
subroutine r8bb_sl_test ( )

!*****************************************************************************80
!
!! R8BB_SL_TEST tests R8BB_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  integer info
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_SL_TEST'
  write ( *, '(a)' ) '  R8BB_SL solves a linear system factored by R8BB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_random ( n1, n2, ml, mu, seed, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The border-banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n1+n2, x )
!
!  Compute the corresponding right hand side.
!
  call r8bb_mv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n1+n2, b, '  The right hand side vector:' )
!
!  Factor the matrix.
!
  call r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

  call r8bb_print ( n1, n2, ml, mu, a, '  The FACTORED border-banded matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8BB_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8BB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8bb_sl ( n1, n2, ml, mu, a, pivot, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8bb_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8BB_TO_R8GE_TEST tests R8BB_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8BB_TO_R8GE converts an R8BB matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_indicator ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB matrix:' )

  call r8bb_to_r8ge ( n1, n2, ml, mu, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8bb_zeros_test ( )

!*****************************************************************************80
!
!! R8BB_ZEROS_TEST tests R8BB_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n1 = 6
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( 2 * ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BB_ZEROS_TEST'
  write ( *, '(a)' ) '  R8BB_ZEROS zeros an R8BB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8bb_zeros ( n1, n2, ml, mu, a )

  call r8bb_print ( n1, n2, ml, mu, a, '  The R8BB zero matrix:' )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

