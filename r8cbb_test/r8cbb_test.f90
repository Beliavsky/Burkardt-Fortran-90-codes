program main

!*****************************************************************************80
!
!! r8cbb_test() tests r8cbb().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8cbb_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8cbb().'

  call r8cbb_add_test ( )
  call r8cbb_dif2_test ( )
  call r8cbb_fa_test ( )
  call r8cbb_get_test ( )
  call r8cbb_indicator_test ( )
  call r8cbb_mtv_test ( )
  call r8cbb_mv_test ( )
  call r8cbb_print_test ( )
  call r8cbb_print_some_test ( )
  call r8cbb_random_test ( )
  call r8cbb_set_test ( )
  call r8cbb_sl_test ( )
  call r8cbb_to_r8ge_test ( )
  call r8cbb_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8cbb_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8cbb_add_test ( )

!*****************************************************************************80
!
!! R8CBB_ADD_TEST tests R8CBB_ADD.
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer j
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_ADD_TEST'
  write ( *, '(a)' ) '  R8CBB_ADD adds a value to elements of an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Initialize matrix to indicator matrix.
!
  call r8cbb_indicator ( n1, n2, ml, mu, a )
!
!  Print initial matrix.
!
  call r8cbb_print ( n1, n2, ml, mu, a, '  Matrix before additions:' )
!
!  Add 100 to band diagonal.
!
  do i = 1, n1
    j = i
    value = 100.0D+00
    call r8cbb_add ( n1, n2, ml, mu, a, i, j, value )
  end do
!
!  Add 200 to right border.
!
  do i = 1, n1
    do j = n1 + 1, n1 + n2
      value = 200.0D+00
      call r8cbb_add ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Add 400 to offdiagonals in lower right dense matrix.
!
  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      if ( i /= j ) then
        value = 400.0D+00
        call r8cbb_add ( n1, n2, ml, mu, a, i, j, value )
      end if
    end do
  end do

  call r8cbb_print ( n1, n2, ml, mu, a, '  The matrix after additions:' )

  return
end
subroutine r8cbb_dif2_test ( )

!*****************************************************************************80
!
!! R8CBB_DIF2_TEST tests R8CBB_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_DIF2_TEST'
  write ( *, '(a)' ) '  R8CBB_DIF2 sets up an R8CBB second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cbb_dif2 ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB second difference matrix:' )

  return
end
subroutine r8cbb_fa_test ( )

!*****************************************************************************80
!
!! R8CBB_FA_TEST tests R8CBB_FA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  integer info
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_FA_TEST'
  write ( *, '(a)' ) '  R8CBB_FA factors an R8CBB matrix with no pivoting.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8cbb_mv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call r8cbb_fa ( n1, n2, ml, mu, a, info )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The factored R8CBB matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8CBB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8vec_print ( n, b, '  The right hand side vector:' )

  call r8cbb_sl ( n1, n2, ml, mu, a, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8cbb_get_test ( )

!*****************************************************************************80
!
!! R8CBB_GET_TEST tests R8CBB_GET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer i4_uniform_ab
  integer j
  integer k
  integer seed
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_GET_TEST'
  write ( *, '(a)' ) '  R8CBB_GET gets a value of an element of an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set matrix to indicator matrix.
!
  call r8cbb_indicator ( n1, n2, ml, mu, a )
!
!  Print matrix.
!
  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix to be queried:' )
!
!  Request random entries.
!
  seed = 123456789

  write ( *, '(a)' ) ''
  do k = 1, 10
    i = i4_uniform_ab ( 1, n1 + n2, seed )
    j = i4_uniform_ab ( 1, n1 + n2, seed )
    call r8cbb_get ( n1, n2, ml, mu, a, i, j, value )
    write ( *, '(a,i2,a,i2,a,g14.6)' ) '  A(', i, ',', j, ') = ', value
  end do

  return
end
subroutine r8cbb_indicator_test ( )

!*****************************************************************************80
!
!! R8CBB_INDICATOR_TEST tests R8CBB_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CBB_INDICATOR sets an indicator matrix for R8CBB format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cbb_indicator ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The compact border-banded matrix:' )

  return
end
subroutine r8cbb_mtv_test ( )

!*****************************************************************************80
!
!! R8CBB_MTV_TEST tests R8CBB_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_MTV_TEST'
  write ( *, '(a)' ) '  R8CBB_MTV computes b=A''*x, where A is an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_indicator ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8cbb_mtv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8cbb_mv_test ( )

!*****************************************************************************80
!
!! R8CBB_MV_TEST tests R8CBB_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_MV_TEST'
  write ( *, '(a)' ) '  R8CBB_MV computes b=A*x, where A is an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_indicator ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix A:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8cbb_mv ( n1, n2, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A*x:' )

  return
end
subroutine r8cbb_print_test ( )

!*****************************************************************************80
!
!! R8CBB_PRINT_TEST tests R8CBB_PRINT.
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_PRINT_TEST'
  write ( *, '(a)' ) '  R8CBB_PRINT prints a compressed border banded matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )

  return
end
subroutine r8cbb_print_some_test ( )

!*****************************************************************************80
!
!! R8CBB_PRINT_SOME_TEST tests R8CBB_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CBB_PRINT_SOME prints some of an R8CBB matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print_some ( n1, n2, ml, mu, a, 1, 9, 10, 10, '  Rows 1-10, Cols 9-10' )

  return
end
subroutine r8cbb_random_test ( )

!*****************************************************************************80
!
!! R8CBB_RANDOM_TEST tests R8CBB_RANDOM.
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CBB_RANDOM randomizes a compressed border banded matrix;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )

  return
end
subroutine r8cbb_set_test ( )

!*****************************************************************************80
!
!! R8CBB_SET_TEST tests R8CBB_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  integer i
  integer j
  real ( kind = rk ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_SET_TEST'
  write ( *, '(a)' ) '  R8CBB_SET sets elements of an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Initialize matrix to zero.
!
  call r8cbb_zeros ( n1, n2, ml, mu, a )
!
!  Fill in band matrix.
!
  do i = 1, n1
    do j = 1, n1
      if ( i - ml <= j .and. j <= i + mu ) then
        value = real ( 10 * i + j, kind = rk )
        call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
      end if 
    end do
  end do
!
!  Fill in right border vector.
!
  do i = 1, n1
    do j = n1 + 1, n1 + n2
      value = real ( 10 * i + j, kind = rk )
      call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Fill in lower border vector.
!
  do i = n1 + 1, n1 + n2
    do j = 1, n1
      value = real ( 10 * i + j, kind = rk )
      call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do
!
!  Fill in lower right dense matrix.
!
  do i = n1 + 1, n1 + n2
    do j = n1 + 1, n1 + n2
      value = real ( 10 * i + j, kind = rk )
      call r8cbb_set ( n1, n2, ml, mu, a, i, j, value )
    end do
  end do

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )

  return
end
subroutine r8cbb_sl_test ( )

!*****************************************************************************80
!
!! R8CBB_SL_TEST tests R8CBB_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) b(n)
  integer info
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_SL_TEST'
  write ( *, '(a)' ) '  R8CBB_SL solves a linear system factored by R8CBB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cbb_random ( n1, n2, ml, mu, seed, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8cbb_mv ( n1, n2, ml, mu, a, x, b )
!
!  Factor the matrix
!
  call r8cbb_fa ( n1, n2, ml, mu, a, info )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The factored R8CBB matrix:' )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8CBB_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8CBB_FA claims the matrix is singular.'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the system.
!
  call r8vec_print ( n, b, '  The right hand side vector:' )

  call r8cbb_sl ( n1, n2, ml, mu, a, b )

  call r8vec_print ( n, b, '  Solution:' )

  return
end
subroutine r8cbb_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CBB_TO_R8GE_TEST tests R8CBB_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 July 2016
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
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)
  real ( kind = rk ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CBB_TO_R8GE converts an R8CBB matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cbb_indicator ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB matrix:' )

  call r8cbb_to_r8ge ( n1, n2, ml, mu, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8cbb_zeros_test ( )

!*****************************************************************************80
!
!! R8CBB_ZEROS_TEST tests R8CBB_ZEROS.
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

  integer, parameter :: n1 = 8
  integer, parameter :: n2 = 2
  integer, parameter :: n = n1 + n2
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1
  integer, parameter :: na = ( ml + mu + 1 ) * n1 + 2 * n1 * n2 + n2 * n2

  real ( kind = rk ) a(na)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CBB_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CBB_ZEROS zeros an R8CBB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Matrix suborder N1 = ', n1
  write ( *, '(a,i8)' ) '  Matrix suborder N2 = ', n2
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cbb_zeros ( n1, n2, ml, mu, a )

  call r8cbb_print ( n1, n2, ml, mu, a, '  The R8CBB zero matrix:' )

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

