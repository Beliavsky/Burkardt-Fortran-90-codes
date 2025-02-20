program main

!*****************************************************************************80
!
!! r8cb_test() tests r8cb().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8cb_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8cb().'

  call r8cb_det_test ( )
  call r8cb_dif2_test ( )
  call r8cb_indicator_test ( )
  call r8cb_ml_test ( )
  call r8cb_mtv_test ( )
  call r8cb_mv_test ( )
  call r8cb_np_fa_test ( )
  call r8cb_np_sl_test ( )
  call r8cb_print_test ( )
  call r8cb_print_some_test ( )
  call r8cb_random_test ( )
  call r8cb_to_r8ge_test ( )
  call r8cb_to_r8vec_test ( )
  call r8cb_zeros_test ( )
  call r8vec_to_r8cb_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8cb_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8cb_det_test ( )

!*****************************************************************************80
!
!! R8CB_DET_TEST tests R8CB_DET.
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

  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) a_r8ge(n,n)
  real ( kind = rk ) det
  integer info
  integer pivot(n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_DET_TEST'
  write ( *, '(a)' ) '  R8CB_DET computes the determinant of a matrix factored'
  write ( *, '(a)' ) '  by R8CB_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8cb_random ( n, n, ml, mu, seed, a )

  call r8cb_print ( n, n, ml, mu, a, '  The compact band matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8cb_to_r8ge ( n, n, ml, mu, a, a_r8ge )
!
!  Factor the matrix.
!
  call r8cb_np_fa ( n, ml, mu, a, info )
!
!  Compute the determinant.
!
  call r8cb_det ( n, ml, mu, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8CB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call r8ge_fa ( n, a_r8ge, pivot, info )
!
!  Compute the determinant.
!
  call r8ge_det ( n, a_r8ge, pivot, det )

  write ( *, '(a,g14.6)' ) '  R8GE_DET computes the determinant = ', det

  return
end
subroutine r8cb_dif2_test ( )

!*****************************************************************************80
!
!! R8CB_DIF2_TEST tests R8CB_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_DIF2_TEST'
  write ( *, '(a)' ) '  R8CB_DIF2 computes an R8CB second difference matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_dif2 ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB second difference matrix:' )

  return
end
subroutine r8cb_indicator_test ( )

!*****************************************************************************80
!
!! R8CB_INDICATOR_TEST tests R8CB_INDICATOR.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CB_INDICATOR computes an R8CB indicator matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  return
end
subroutine r8cb_ml_test ( )

!*****************************************************************************80
!
!! R8CB_ML_TEST tests R8CB_ML.
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

  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_ML_TEST'
  write ( *, '(a)' ) '  R8CB_ML computes A*x or A''*x'
  write ( *, '(a)' ) '  for a compact band matrix A which has'
  write ( *, '(a)' ) '  been factored by R8CB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8CB_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8cb_ml ( n, ml, mu, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print ( n, b, b2, '  A*x and PLU*x' )
    else
      call r8vec2_print ( n, b, b2, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8cb_mtv_test ( )

!*****************************************************************************80
!
!! R8CB_MTV_TEST tests R8CB_MTV.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 8
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_MTV_TEST'
  write ( *, '(a)' ) '  R8CB_MTV computes b=A''*x, where A is an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8vec_indicator1 ( m, x )

  call r8vec_print ( m, x, '  The vector x:' )

  call r8cb_mtv ( m, n, ml, mu, a, x, b )

  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8cb_mv_test ( )

!*****************************************************************************80
!
!! R8CB_MV_TEST tests R8CB_MV.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 8
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(m)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_MV_TEST'
  write ( *, '(a)' ) '  R8CB_MV computes b=A*x, where A is an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8vec_indicator1 ( n, x )

  call r8vec_print ( n, x, '  The vector x:' )

  call r8cb_mv ( m, n, ml, mu, a, x, b )

  call r8vec_print ( m, b, '  The product b=A*x:' )

  return
end
subroutine r8cb_np_fa_test ( )

!*****************************************************************************80
!
!! R8CB_NP_FA_TEST tests R8CB_NP_FA.
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

  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_NP_FA_TEST'
  write ( *, '(a)' ) '  R8CB_NP_FA factors an R8CB matrix, no pivoting:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_FA_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8CB_NP_FA claims the matrix is singular.'
      write ( *, '(a,i8)' ) '  The value of info is ', info
      return
    end if
!
!  Solve the system.
!
    call r8cb_np_sl ( n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine r8cb_np_sl_test ( )

!*****************************************************************************80
!
!! R8CB_NP_SL_TEST tests R8CB_NP_SL.
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

  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_NP_SL_TEST'
  write ( *, '(a)' ) '  R8CB_NP_SL solves a linear system factored by R8CB_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    call r8cb_random ( n, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
    if ( job == 0 ) then
      call r8cb_mv ( n, n, ml, mu, a, x, b )
    else
      call r8cb_mtv ( n, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8cb_np_fa ( n, ml, mu, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8CB_NP_SL_TEST Fatal error!'
      write ( *, '(a)' ) '  R8CB_NP_FA claims the matrix is singular.'
      write ( *, '(a,i8)' ) '  The value of info is ', info
      return
    end if
!
!  Solve the system.
!
    call r8cb_np_sl ( n, ml, mu, a, b, job )

    if ( job == 0 ) then
      call r8vec_print ( n, b, '  Solution:' )
    else
      call r8vec_print ( n, b, '  Solution to transposed system:' )
    end if

  end do

  return
end
subroutine r8cb_print_test ( )

!*****************************************************************************80
!
!! R8CB_PRINT_TEST tests R8CB_PRINT.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_PRINT_TEST'
  write ( *, '(a)' ) '  R8CB_PRINT prints an R8CB matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  return
end
subroutine r8cb_print_some_test ( )

!*****************************************************************************80
!
!! R8CB_PRINT_SOME_TEST tests R8CB_PRINT_SOME.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CB_PRINT_SOME prints some of an R8CB matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print_some ( m, n, ml, mu, a, 3, 3, 6, 6, '  Rows 3-6, Cols 3-6:' )

  return
end
subroutine r8cb_random_test ( )

!*****************************************************************************80
!
!! R8CB_RANDOM_TEST tests R8CB_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CB_RANDOM randomizes an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  seed = 123456789
  call r8cb_random ( m, n, ml, mu, seed, a )

  call r8cb_print ( m, n, ml, mu, a, '  The random R8CB matrix:' )

  return
end
subroutine r8cb_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CB_TO_R8GE_TEST tests R8CB_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 8
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(ml+mu+1,n)
  real ( kind = rk ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CB_TO_R8GE converts an R8CB matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8cb_to_r8ge ( m, n, ml, mu, a, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8cb_to_r8vec_test ( )

!*****************************************************************************80
!
!! R8CB_TO_R8VEC_TEST tests R8CB_TO_R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 8
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(ml+mu+1,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x((ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_TO_R8VEC_TEST'
  write ( *, '(a)' ) '  R8CB_TO_R8VEC converts an R8CB matrix to an R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB matrix:' )

  call r8cb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, ml + mu + 1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8cb ( m, n, ml, mu, x, a )

  call r8cb_print ( m, n, ml, mu, a, '  The recovered R8CB matrix:' )

  return
end
subroutine r8cb_zeros_test ( )

!*****************************************************************************80
!
!! R8CB_ZEROS_TEST tests R8CB_ZEROS.
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

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 2
  integer, parameter :: mu = 3

  real ( kind = rk ) a(ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CB_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CB_ZEROS zeros an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_zeros ( n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB zero matrix:' )

  return
end
subroutine r8vec_to_r8cb_test ( )

!*****************************************************************************80
!
!! R8VEC_TO_R8CB_TEST tests R8VEC_TO_R8CB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 8
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(ml+mu+1,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x((ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_TO_R8CB_TEST'
  write ( *, '(a)' ) '  R8VEC_TO_R8CB converts an R8VEC to an R8CB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M      = ', m
  write ( *, '(a,i8)' ) '  Matrix columns N   = ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8cb_indicator ( m, n, ml, mu, a )

  call r8cb_print ( m, n, ml, mu, a, '  The R8CB indicator matrix:' )

  call r8cb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, ml + mu + 1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8cb ( m, n, ml, mu, x, a )

  call r8cb_print ( m, n, ml, mu, a, '  The recovered R8CB indicator matrix:' )

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
