program main

!*****************************************************************************80
!
!! r8ge_np_test() tests r8ge_np().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ge_np_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ge_np().'

  call r8ge_np_det_test ( )
  call r8ge_np_fa_test ( )
  call r8ge_np_inverse_test ( )
  call r8ge_np_ml_test ( )
  call r8ge_np_sl_test ( )
  call r8ge_np_trf_test ( )
  call r8ge_np_trm_test ( )
  call r8ge_np_trs_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ge_np_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ge_np_det_test ( )

!*****************************************************************************80
!
!! R8GE_NP_DET_TEST tests R8GE_NP_DET.
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

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) det
  integer info

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_DET_TEST'
  write ( *, '(a)' ) '  R8GE_NP_DET computes the determinant of a matrix'
  write ( *, '(a)' ) '  that was factored by R8GE_NP_FA,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_dif2 ( n, n, a )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_DET_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Get the determinant.
!
  call r8ge_np_det ( n, a, det )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant of -1, 2, -1 matrix is ', det
  write ( *, '(a,g14.6)' ) '  Exact value is ', real ( n + 1, kind = rk )

  return
end
subroutine r8ge_np_fa_test ( )

!*****************************************************************************80
!
!! R8GE_NP_FA_TEST tests R8GE_NP_FA.
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

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_FA_TEST'
  write ( *, '(a)' ) '  R8GE_NP_FA computes the LU factors of a general'
  write ( *, '(a)' ) '  storage matrix without pivoting,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_FA_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_inverse_test ( )

!*****************************************************************************80
!
!! R8GE_NP_INVERSE_TEST tests R8GE_NP_INVERSE.
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

  integer, parameter :: n = 5

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n,n)
  real ( kind = rk ) c(n,n)
  integer info
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_INVERSE'
  write ( *, '(a)' ) '  R8GE_NP_INVERSE computes the inverse of an R8GE matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )

  call r8ge_print ( n, n, a, '  The random matrix:' )
!
!  Factor and invert the matrix.
!
  b(1:n,1:n) = a(1:n,1:n)

  call r8ge_np_fa ( n, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_INVERSE - Warning!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if

  call r8ge_np_inverse ( n, b )

  call r8ge_print ( n, n, b, '  The inverse matrix:' )
!
!  Compute A * B = C.
!
  call r8ge_mm ( n, n, n, a, b, c )

  call r8ge_print ( n, n, c, '  The product:' )

  return
end
subroutine r8ge_np_ml_test ( )

!*****************************************************************************80
!
!! R8GE_NP_ML_TEST tests R8GE_NP_ML.
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

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_ML_TEST'
  write ( *, '(a)' ) '  R8GE_NP_ML computes A*x or A''*X'
  write ( *, '(a)' ) '  where A has been factored by R8GE_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  do job = 0, 1
!
!  Set the matrix.
!
    call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8ge_mv ( n, n, a, x, b )
    else
      call r8ge_mtv ( n, n, a, x, b )
    end if
!
!  Factor the matrix.
!
    call r8ge_np_fa ( n, a, info )

    if ( info /= 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8GE_NP_ML_TEST - Fatal error!'
      write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
      write ( *, '(a,i8)' ) '  The value of INFO is ', info
      cycle
    end if
!
!  Now multiply factored matrix times solution to get right hand side again.
!
    call r8ge_np_ml ( n, a, x, b2, job )

    if ( job == 0 ) then
      call r8vec2_print_some ( n, b, b2, 10, '  A*x and PLU*x' )
    else
      call r8vec2_print_some ( n, b, b2, 10, '  A''*x and (PLU)''*x' )
    end if

  end do

  return
end
subroutine r8ge_np_sl_test ( )

!*****************************************************************************80
!
!! R8GE_NP_SL_TEST tests R8GE_NP_SL.
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

  integer, parameter :: n = 10

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) b(n)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_SL_TEST'
  write ( *, '(a)' ) '  R8GE_NP_SL solves a linear system that was factored'
  write ( *, '(a)' ) '  by R8GE_NP_FA.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( n, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( n, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_fa ( n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_SL_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_FA declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )
 
  call r8vec_print_some ( n, b, 10, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 0
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_ml ( n, a, x, b, job )
!
!  Solve the system
!
  job = 1
  call r8ge_np_sl ( n, a, b, job )

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trf_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRF_TEST tests R8GE_NP_TRF.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: nrhs = 1

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,nrhs)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRF_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRF factors an R8GE matrix without pivoting,'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRF_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trm_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRM_TEST tests R8GE_NP_TRM.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: nrhs = 1

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,nrhs)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRM_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRM computes b=A*x after A has'
  write ( *, '(a)' ) '  been factored by R8GE_NP_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRM_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

  return
end
subroutine r8ge_np_trs_test ( )

!*****************************************************************************80
!
!! R8GE_NP_TRS_TEST tests R8GE_NP_TRS.
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

  integer, parameter :: m = 10
  integer, parameter :: n = m
  integer, parameter :: nrhs = 1

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) b(m,nrhs)
  integer info
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8GE_NP_TRS_TEST'
  write ( *, '(a)' ) '  R8GE_NP_TRS solves a linear system A*x=b'
  write ( *, '(a)' ) '  which has been factored by R8GE_NP_TRF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n
!
!  Set the matrix.
!
  call r8ge_random ( m, n, seed, a )
!
!  Set the desired solution.
!
  x(1:n) = 1.0D+00
!
!  Compute the corresponding right hand side.
!
  call r8ge_mv ( m, n, a, x, b )
!
!  Factor the matrix.
!
  call r8ge_np_trf ( m, n, a, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_NP_TRF declares the matrix is singular!'
    write ( *, '(a,i8)' ) '  The value of INFO is ', info
    return
  end if
!
!  Solve the linear system.
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 0
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system
!
  call r8ge_np_trs ( n, nrhs, 'N', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  job = 1
  call r8ge_np_trm ( m, n, a, x, b, job )
!
!  Solve the system.
!
  call r8ge_np_trs ( n, nrhs, 'T', a, b, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8GE_NP_TRS_TEST - Fatal error!'
    write ( *, '(a)' ) '  R8GE_TRS returned an error condition!'
    write ( *, '(a)' ) '  The value of INFO is ', info
    return
  end if

  call r8vec_print ( n, b, '  Solution of transposed system:' )

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

