program main

!*****************************************************************************80
!
!! r8gb_test() tests r8gb().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8gb_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8gb().'

  call r8gb_det_test ( )
  call r8gb_dif2_test ( )
  call r8gb_fa_test ( )
  call r8gb_indicator_test ( )
  call r8gb_ml_test ( )
  call r8gb_mtv_test ( )
  call r8gb_mu_test ( )
  call r8gb_mv_test ( )
  call r8gb_nz_num_test ( )
  call r8gb_print_test ( )
  call r8gb_print_some_test ( )
  call r8gb_random_test ( )
  call r8gb_sl_test ( )
  call r8gb_to_r8ge_test ( )
  call r8gb_to_r8vec_test ( )
  call r8gb_trf_test ( )
  call r8gb_trs_test ( )
  call r8gb_zeros_test ( )

  call r8ge_to_r8gb_test ( )

  call r8vec_to_r8gb_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8gb_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8gb_det_test ( )

!*****************************************************************************80
!
!! R8GB_DET_TEST tests R8GB_DET.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 3
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) det
  integer info
  integer pivot(n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_DET_TEST'
  write ( *, '(a)' ) '  R8GB_DET computes the determinant'
  write ( *, '(a)' ) '  for a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  A random R8GB matrix:' )
!
!  Copy the matrix into a general array.
!
  call r8gb_to_r8ge ( m, n, ml, mu, a, a2 )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )
!
!  Compute the determinant.
!
  call r8gb_det ( n, ml, mu, a, pivot, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  R8GB_DET computes the determinant = ', det
!
!  Factor the general matrix.
!
  call r8ge_fa ( n, a2, pivot, info )
!
!  Compute the determinant.
!
  call r8ge_det ( n, a2, pivot, det )

  write ( *, '(a,g14.6)' ) '  R8GE_DET computes the determinant = ', det

  return
end
subroutine r8gb_dif2_test ( )

!*****************************************************************************80
!
!! R8GB_DIF2_TEST tests R8GB_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: ml = 1
  integer, parameter :: mu = 1

  integer, parameter :: row_num = 2 * ml + mu + 1

  real ( kind = rk ) a(row_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_DIF2_TEST'
  write ( *, '(a)' ) '  R8GB_DIF2 returns an R8GB second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu

  call r8gb_dif2 ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB second difference matrix:' )

  return
end
subroutine r8gb_fa_test ( )

!*****************************************************************************80
!
!! R8GB_FA_TEST tests R8GB_FA.
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

  integer, parameter :: m = 5
  integer, parameter :: n = m
  integer, parameter :: ml = 1 
  integer, parameter :: mu = 2 

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_FA_TEST'
  write ( *, '(a)' ) '  R8GB_FA computes the PLU factors'
  write ( *, '(a)' ) '  of a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  The banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GB_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_ml ( n, ml, mu, a, pivot, x, b, job )

  call r8vec_print ( n, b, '  Right hand side of transposed system:' )
!
!  Solve the linear system.
!
  job = 1
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine r8gb_indicator_test ( )

!*****************************************************************************80
!
!! R8GB_INDICATOR_TEST tests R8GB_INDICATOR.
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

  integer, parameter :: m = 10
  integer, parameter :: n = 8
  integer, parameter :: ml = 3
  integer, parameter :: mu = 2

  integer, parameter :: row_num = 2 * ml + mu + 1

  real ( kind = rk ) a(row_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8GB_INDICATOR computes the indicator matrix for a general band matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB indicator matrix:' )

  return
end
subroutine r8gb_ml_test ( )

!*****************************************************************************80
!
!! R8GB_ML_TEST tests R8GB_ML.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer pivot(n)
  integer seed
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_ML_TEST'
  write ( *, '(a)' ) '  R8GB_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  for a general banded matrix A'
  write ( *, '(a)' ) '  where A has been factored by R8GB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    seed = 123456789
    call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8gb_mv ( m, n, ml, mu, a, x, b )
    else
      call r8gb_mtv ( m, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8gb_fa ( n, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GB_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GB_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8gb_ml ( n, ml, mu, a, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print ( n, b, b2, '  A*x and PLU*x' )
    else
      call r8vec2_print ( n, b, b2, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8gb_mtv_test ( )

!*****************************************************************************80
!
!! R8GB_MTV_TEST tests R8GB_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(m)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_MTV_TEST'
  write ( *, '(a)' ) '  R8GB_MTV computes A''*x for an R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  seed = 123456789
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set x.
!
  call r8vec_indicator1 ( m, x )
  call r8vec_print ( m, x, '  x:' );
!
!  Compute b.
!
  call r8gb_mtv ( m, n, ml, mu, a, x, b )
  call r8vec_print ( n, b, '  b=A''*x:' )

  return
end
subroutine r8gb_mu_test ( )

!*****************************************************************************80
!
!! R8GB_MU_TEST tests R8GB_MU.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer pivot(n)
  integer seed
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_MU_TEST'
  write ( *, '(a)' ) '  R8GB_MU computes A*x or A''*X'
  write ( *, '(a)' ) '  for a general banded matrix A'
  write ( *, '(a)' ) '  where A has been factored by R8GB_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  do job = 0, 1
!
!  Set the matrix.
!
    seed = 123456789
    call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8gb_mv ( m, n, ml, mu, a, x, b )
    else
      call r8gb_mtv ( m, n, ml, mu, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8gb_trf ( m, n, ml, mu, a, pivot, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GB_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GB_TRF declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      return
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8gb_mu ( n, ml, mu, a, pivot, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print ( n, b, b2, '  A*x and PLU*x' )
    else
      call r8vec2_print ( n, b, b2, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8gb_mv_test ( )

!*****************************************************************************80
!
!! R8GB_MV_TEST tests R8GB_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(m)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_MV_TEST'
  write ( *, '(a)' ) '  R8GB_MV computes A*x for an R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  seed = 123456789
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set x.
!
  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' );
!
!  Compute b.
!
  call r8gb_mv ( m, n, ml, mu, a, x, b )
  call r8vec_print ( m, b, '  b=A*x:' )

  return
end
subroutine r8gb_nz_num_test ( )

!*****************************************************************************80
!
!! R8GB_NZ_NUM_TEST tests R8GB_NZ_NUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 July 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer diag
  integer j
  integer nz_num
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_NZ_NUM_TEST'
  write ( *, '(a)' ) '  R8GB_NZ_NUM counts the nonzero entries'
  write ( *, '(a)' ) '  for a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Make some zero entries.
!
  do j = 1, n
    do diag = 1, 2*ml+mu+1
      if ( a(diag,j) < 0.3D+00 ) then
        a(diag,j) = 0.0D+00
      end if
    end do
  end do

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  call r8gb_nz_num ( m, n, ml, mu, a, nz_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Nonzero entries = ', nz_num

  return
end
subroutine r8gb_print_test ( )

!*****************************************************************************80
!
!! R8GB_PRINT_TEST tests R8GB_PRINT.
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
  integer, parameter :: ml = 1
  integer, parameter :: mu = 3

  real ( kind = rk ) a(2*ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_PRINT_TEST'
  write ( *, '(a)' ) '  R8GB_PRINT prints a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB matrix:' )

  return
end
subroutine r8gb_print_some_test ( )

!*****************************************************************************80
!
!! R8GB_PRINT_SOME_TEST tests R8GB_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = 10
  integer, parameter :: ml = 1
  integer, parameter :: mu = 3

  real ( kind = rk ) a(2*ml+mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8GB_PRINT_SOME prints some of a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print_some ( m, n, ml, mu, a, 5, 4, 7, 10, '  Rows(5-7), Cols (4-10)' )

  return
end
subroutine r8gb_random_test ( )

!*****************************************************************************80
!
!! R8GB_RANDOM_TEST tests R8GB_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  integer, parameter :: row_num = 2 * ml + mu + 1

  real ( kind = rk ) a(row_num,n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_RANDOM_TEST'
  write ( *, '(a)' ) '  R8GB_RANDOM returns a random R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu

  seed = 123456789
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  The random R8GB matrix:' )

  return
end
subroutine r8gb_sl_test ( )

!*****************************************************************************80
!
!! R8GB_SL_TEST tests R8GB_SL.
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

  integer, parameter :: m = 5
  integer, parameter :: n = m
  integer, parameter :: ml = 1 
  integer, parameter :: mu = 2 

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_SL_TEST'
  write ( *, '(a)' ) '  R8GB_SL solves a linear system.'
  write ( *, '(a)' ) '  that was factored by R8GB_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )

  call r8gb_print ( m, n, ml, mu, a, '  The banded matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_fa ( n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GB_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_ml ( n, ml, mu, a, pivot, x, b, job )

  call r8vec_print ( n, b, '  Right hand side of transposed system:' )
!
!  Solve the linear system.
!
  job = 1
  call r8gb_sl ( n, ml, mu, a, pivot, b, job )
 
  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine r8gb_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8GB_TO_R8GE_TEST tests R8GB_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2007
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

  real ( kind = rk ) a_r8gb(2*ml+mu+1,n)
  real ( kind = rk ) a_r8ge(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8GB_TO_R8GE copies a R8GB matrix to a R8GE matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a_r8gb )

  call r8gb_print ( m, n, ml, mu, a_r8gb, '  The R8GB matrix:' )

  call r8gb_to_r8ge ( m, n, ml, mu, a_r8gb, a_r8ge )

  call r8ge_print ( m, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8gb_to_r8vec_test ( )

!*****************************************************************************80
!
!! R8GB_TO_R8VEC_TEST tests R8GB_TO_R8VEC.
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

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x((2*ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_TO_R8VEC_TEST'
  write ( *, '(a)' ) '  R8GB_TO_R8VEC converts an R8GB matrix to a real vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB indicator matrix:' )

  call r8gb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, 2*ml+mu+1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8gb ( m, n, ml, mu, x, a )

  call r8gb_print ( m, n, ml, mu, a, '  The recovered R8GB indicator matrix:' )

  return
end
subroutine r8gb_trf_test ( )

!*****************************************************************************80
!
!! R8GB_TRF_TEST tests R8GB_TRF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: m = n
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: nrhs = 1

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n,nrhs)
  integer info
  integer job
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_TRF_TEST'
  write ( *, '(a)' ) '  R8GB_TRF computes the PLU factors'
  write ( *, '(a)' ) '  for a general banded matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_trf ( m, n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GB_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'N', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_mu ( n, ml, mu, a, pivot, x, b, job )
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'T', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine r8gb_trs_test ( )

!*****************************************************************************80
!
!! R8GB_TRS_TEST tests R8GB_TRS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: m = n
  integer, parameter :: ml = 1
  integer, parameter :: mu = 2
  integer, parameter :: nrhs = 1

  real ( kind = rk ) a(2*ml+mu+1,n)
  real ( kind = rk ) b(n,nrhs)
  integer info
  integer job
  integer pivot(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_TRS_TEST'
  write ( *, '(a)' ) '  R8GB_TRS solves a general band matrix linear system'
  write ( *, '(a)' ) '  that was factored by R8GB_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8gb_random ( m, n, ml, mu, seed, a )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8gb_mv ( m, n, ml, mu, a, x, b )
!
!  Factor the matrix.
!
  call r8gb_trf ( m, n, ml, mu, a, pivot, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GB_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GB_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'N', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8gb_mu ( n, ml, mu, a, pivot, x, b, job )
!
!  Solve the linear system.
!
  call r8gb_trs ( n, ml, mu, nrhs, 'T', a, pivot, b, info )

  call r8vec_print ( n, b, '  Solution to transposed system:' )

  return
end
subroutine r8gb_zeros_test ( )

!*****************************************************************************80
!
!! R8GB_ZEROS_TEST tests R8GB_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 10
  integer, parameter :: n = 8
  integer, parameter :: ml = 3
  integer, parameter :: mu = 2

  integer, parameter :: row_num = 2 * ml + mu + 1

  real ( kind = rk ) a(row_num,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GB_ZEROS_TEST'
  write ( *, '(a)' ) '  R8GB_ZEROS returns a zero R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =       ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =    ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML  = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU  = ', mu

  call r8gb_zeros ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB zero matrix:' )

  return
end
subroutine r8ge_to_r8gb_test ( )

!*****************************************************************************80
!
!! R8GE_TO_R8GB_TEST tests R8GE_TO_R8GB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 July 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5
  integer, parameter :: ml = 2
  integer, parameter :: mu = 1

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) a_r8gb(2*ml+mu+1,n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_TO_R8GB_TEST'
  write ( *, '(a)' ) '  R8GE_TO_R8GB converts an R8GE matrix to an R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  seed = 123456789
  call r8ge_random ( m, n, seed, a )
  call r8ge_print ( m, n, a, '  The R8GE matrix:' )

  call r8ge_to_r8gb ( m, n, ml, mu, a, a_r8gb )

  call r8gb_print ( m, n, ml, mu, a_r8gb, '  The R8GB matrix:' )

  return
end
subroutine r8vec_to_r8gb_test ( )

!*****************************************************************************80
!
!! R8VEC_TO_R8GB_TEST tests R8VEC_TO_R8GB.
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

  real ( kind = rk ) a(2*ml+mu+1,n)
  integer i
  integer j
  integer k
  real ( kind = rk ) x((2*ml+mu+1)*n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_TO_R8GB_TEST'
  write ( *, '(a)' ) '  R8VEC_TO_R8GB converts a real vector to an R8GB matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =      ', m
  write ( *, '(a,i8)' ) '  Matrix columns N =   ', n
  write ( *, '(a,i8)' ) '  Lower bandwidth ML = ', ml
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8gb_indicator ( m, n, ml, mu, a )

  call r8gb_print ( m, n, ml, mu, a, '  The R8GB indicator matrix:' )

  call r8gb_to_r8vec ( m, n, ml, mu, a, x )

  k = 0
  do j = 1, n
    do i = 1, 2*ml+mu+1
      k = k + 1
      write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
    end do
  end do

  call r8vec_to_r8gb ( m, n, ml, mu, x, a )

  call r8gb_print ( m, n, ml, mu, a, '  The recovered R8GB indicator matrix:' )

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

