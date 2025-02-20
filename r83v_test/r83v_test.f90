program main

!*****************************************************************************80
!
!! r83v_test() tests r83v().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r83v_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test r83v().'

  call r83v_cg_test ( )
  call r83v_copy_test ( )
  call r83v_cr_fa_test ( )
  call r83v_cr_sl_test ( )
  call r83v_cr_sls_test ( )
  call r83v_dif2_test ( )
  call r83v_fs_test ( )
  call r83v_gs_sl_test ( )
  call r83v_indicator_test ( )
  call r83v_jac_sl_test ( )
  call r83v_mtv_test ( )
  call r83v_mv_test ( )
  call r83v_print_test ( )
  call r83v_print_some_test ( )
  call r83v_random_test ( )
  call r83v_res_test ( )
  call r83v_to_r8ge_test ( )
  call r83v_to_r8vec_test ( )
  call r83v_transpose_test ( )
  call r83v_zeros_test ( )

  call r8ge_to_r83v_test ( )

  call r8vec_to_r83v_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r83v_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine r83v_cg_test ( )

!*****************************************************************************80
!
!! R83V_CG_TEST tests R83V_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ) r8vec_norm_affine
  integer seed
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CG_TEST'
  write ( *, '(a)' ) '  R83V_CG applies CG to an R83V matrix.'
!
!  Set A to the second difference matrix.
!
  n = 10

  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  seed = 123456789
  allocate ( x1(1:n) )
  call r8vec_uniform_01 ( n, seed, x1 )
  call r8vec_print ( n, x1, '  Desired solution X1:' )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00

  call r83v_cg ( n, a, b, c, ax, x2 )
  call r8vec_print ( n, x2, '  Computed solution X2' )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83v_res ( n, n, a, b, c, x2, ax, r )

  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_copy_test ( )

!*****************************************************************************80
!
!! R83V_COPY_TEST tests R83V_COPY.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  real ( kind = rk ), allocatable :: b1(:)
  real ( kind = rk ), allocatable :: b2(:)
  real ( kind = rk ), allocatable :: b3(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_COPY_TEST'
  write ( *, '(a)' ) '  R83V_COPY copies an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  A:' )

    allocate ( b1(1:min(m-1,n)) )
    allocate ( b2(1:min(m,n)) )
    allocate ( b3(1:min(m,n-1)) )

    call r83v_copy ( m, n, a1, a2, a3, b1, b2, b3 )
    call r83v_print ( m, n, b1, b2, b3, '  B = copy of A:' )

    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )
    deallocate ( b1 )
    deallocate ( b2 )
    deallocate ( b3 )

  end do

  return
end
subroutine r83v_cr_fa_test ( )

!*****************************************************************************80
!
!! R83V_CR_FA_TEST tests R83V_CR_FA.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: a_cr(:,:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  logical debug
  integer j
  integer n
  real ( kind = rk ), allocatable :: x(:)

  n = 5
  debug = .false.

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_FA_TEST'
  write ( *, '(a)' ) '  R83V_CR_FA factors an R83V matrix'
  write ( *, '(a)' ) '  Once the matrix has been factored, we can call'
  write ( *, '(a)' ) '  R83V_CR_SL to solve a linear system.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  System matrix A:' )
!
!  Factor the matrix once.
!
  allocate ( a_cr(1:3,1:2*n+1) )
  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Print the factor information.
!
  if ( debug ) then
    call r83_print ( n, 2 * n + 1, a_cr, &
      '  Cyclic reduction factor information:' )
  end if

  allocate ( ax(1:n) )
  allocate ( x(1:n) )

  do j = 1, 2

    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Solve linear system number #', j

    if ( j == 1 ) then
      ax(1:n-1) = 0.0D+00
      ax(n) = real ( n + 1, kind = rk )
    else
      ax(1) = 1.0D+00
      ax(2:n-1) = 0.0D+00
      ax(n) = 1.0D+00
    end if
!
!  Solve the linear system.
!
    call r83v_cr_sl ( n, a_cr, ax, x )

    call r8vec_print ( n, x, '  Solution:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_cr_sl_test ( )

!*****************************************************************************80
!
!! R83V_CR_SL_TEST tests R83V_CR_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: a_cr(:,:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer n
  real ( kind = rk ), allocatable :: x(:)

  n = 5

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_SL_TEST'
  write ( *, '(a)' ) '  R83V_CR_SL solves a factored system'
  write ( *, '(a)' ) '  after R83V_CR_FA has factored it.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  a(1:n-1) = 1.0D+00
  b(1:n) = 4.0D+00
  c(1:n-1) = 2.0D+00
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  Input matrix A:' )
!
!  Set the desired solution.
!
  allocate ( x(1:n) )
  call r8vec_indicator1 ( n, x )
!
!  Compute the right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x,ax )

  x(1:n) = 0.0D+00
!
!  Factor the matrix.
!
  allocate ( a_cr(1:3,1:2*n+1) )
  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Solve the linear system.
!
  call r83v_cr_sl ( n, a_cr, ax, x )

  call r8vec_print ( n, x, '  Solution:' )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_cr_sls_test ( )

!*****************************************************************************80
!
!! R83V_CR_SLS_TEST tests R83V_CR_SLS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: a_cr(:,:)
  real ( kind = rk ), allocatable :: ax(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer n
  integer nb
  real ( kind = rk ), allocatable :: x(:,:)

  n = 5

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_CR_SLS_TEST'
  write ( *, '(a)' ) '  R83V_CR_SLS solves multiple systems A*x1:xn=b1:bn'
  write ( *, '(a)' ) '  after R83V_CR_FA has factored it.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
  write ( *, '(a)' ) '  Demonstrate multiple system solution method.'
  write ( *, '(a)' ) ''
!
!  Set the matrix values.
!
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Print the matrix.
!
  call r83v_print ( n, n, a, b, c, '  Input matrix A:' )
!
!  Factor the matrix once.
!
  allocate ( a_cr(3,2*n+1) )

  call r83v_cr_fa ( n, a, b, c, a_cr )
!
!  Set the number of right hand sides.
!
  nb = 2
  allocate ( ax(1:n,1:nb) )

  ax(1:n-1,1) = 0.0D+00
  ax(n,1) = n + 1

  ax(1,2) = 1.0D+00
  ax(2:n-1,2) = 0.0D+00
  ax(n,2) = 1.0D+00

  call r8ge_print ( n, 2, ax, '  Right hand sides b1:b2' )
!
!  Solve the linear systems.
!
  allocate ( x(1:n,1:nb) )

  call r83v_cr_sls ( n, a_cr, nb, ax, x )

  call r8ge_print ( n, nb, x, '  Solutions x1:x2' )
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( a_cr )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x )

  return
end
subroutine r83v_dif2_test ( )

!*****************************************************************************80
!
!! R83V_DIF2_TEST tests R83V_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_DIF2_TEST'
  write ( *, '(a)' ) '  R83V_DIF2 sets an R83V matrix to the second difference.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_dif2 ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Second difference in R83V format:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_fs_test ( )

!*****************************************************************************80
!
!! R83V_FS_TEST tests R83V_FS_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  real ( kind = rk ), allocatable :: b(:)
  integer n
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  n = 10

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_FS_TEST'
  write ( *, '(a)' ) '  R83V_FS factors and solves a linear system'
  write ( *, '(a)' ) '  for an R83V matrix.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  Matrix order N = ', n
!
!  Set the matrix values.
!
  allocate ( a1(1:n-1) )
  allocate ( a2(n) )
  allocate ( a3(n-1) )

  call r83v_dif2 ( n, n, a1, a2, a3 )
!
!  Set the desired solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83v_mv ( n, n, a1, a2, a3, x1, b )

  call r8vec_print  ( n, b, "  The right hand side:" )
!
!  Solve the linear system.
!
  allocate ( x2(1:n) )
  call r83v_fs ( n, a1, a2, a3, b, x2 )

  call r8vec_print ( n, x2, "  Solution:" )
!
!  Free memory.
!
  deallocate ( a1 )
  deallocate ( a2 )
  deallocate ( a3 )
  deallocate ( b )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_gs_sl_test ( )

!*****************************************************************************80
!
!! R83V_GS_SL_TEST tests R83V_GS_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer it_max
  integer n
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_GS_SL_TEST'
  write ( *, '(a)' ) '  R83V_GS_SL applies Gauss-Seidel iteration with an'
  write ( *, '(a)' ) '  R83V matrix to solve a linear system A*x=b.'
!
!  Set A to the second difference matrix.
!
  n = 10
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Set the starting solution.
!
  allocate ( x2(1:n) )
  x2(1:n) = 0.0D+00
!
!  Call the Gauss-Seidel routine.
!
  it_max = 25

  do i = 1, 3

    call r83v_gs_sl ( n, a, b, c, ax, x2, it_max )

    call r8vec_print ( n, x2, '  Current solution estimate:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_indicator_test ( )

!*****************************************************************************80
!
!! R83V_INDICATOR_TEST tests R83V_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_INDICATOR_TEST'
  write ( *, '(a)' ) '  R83V_INDICATOR sets an R83V indicator matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  R83V indicator matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_jac_sl_test ( )

!*****************************************************************************80
!
!! R83V_JAC_SL_TEST tests R83V_JAC_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer it_max
  integer n
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_JAC_SL_TEST'
  write ( *, '(a)' ) '  R83V_JAC_SL applies Jacobi iteration with an R83V'
  write ( *, '(a)' ) '  matrix to solve a linear system A*x=b.'
!
!  Set A to the second difference matrix.
!
  n = 10
  allocate ( a(1:n-1) )
  allocate ( b(1:n  ) )
  allocate ( c(1:n-1) )

  call r83v_dif2 ( n, n, a, b, c )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call r8vec_indicator1 ( n, x1 )
!
!  Compute the corresponding right hand side.
!
  allocate ( ax(1:n) )
  call r83v_mv ( n, n, a, b, c, x1, ax )
!
!  Set the starting vector.
!
  allocate ( x2(1:n) )
  x2(1:n) = 0.0D+00
!
!  Call the Jacobi routine.
!
  it_max = 25

  do i = 1, 3

    call r83v_jac_sl ( n, a, b, c, ax, x2, it_max )

    call r8vec_print ( n, x2, '  Current solution:' )

  end do
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( ax )
  deallocate ( b )
  deallocate ( c )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83v_mtv_test ( )

!*****************************************************************************80
!
!! R83V_MTV_TEST tests R83V_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a_83v(:)
  real ( kind = rk ), allocatable :: a_ge(:,:)
  real ( kind = rk ), allocatable :: ax_83v(:)
  real ( kind = rk ), allocatable :: ax_ge(:)
  real ( kind = rk ), allocatable :: b_83v(:)
  real ( kind = rk ), allocatable :: c_83v(:)
  integer i
  integer m
  integer n
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_MTV_TEST'
  write ( *, '(a)' ) '  R83V_MV computes b=A''*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )

    allocate ( x(1:m) )
    call r8vec_indicator1 ( m, x )

    allocate ( ax_83v(1:n) )
    call r83v_mtv ( m, n, a_83v, b_83v, c_83v, x, ax_83v )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )

    allocate ( ax_ge(1:n) )
    call r8ge_mtv ( m, n, a_ge, x, ax_ge )

    call r8vec2_print ( n, ax_83v, ax_ge, '  Product comparison:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( ax_83v )
    deallocate ( ax_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )
    deallocate ( x )

  end do

  return
end
subroutine r83v_mv_test ( )

!*****************************************************************************80
!
!! R83V_MV_TEST tests R83V_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a_83v(:)
  real ( kind = rk ), allocatable :: a_ge (:,:)
  real ( kind = rk ), allocatable :: ax_83v(:)
  real ( kind = rk ), allocatable :: ax_ge(:)
  real ( kind = rk ), allocatable :: b_83v(:)
  real ( kind = rk ), allocatable :: c_83v(:)
  integer i
  integer m
  integer n
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_MV_TEST'
  write ( *, '(a)' ) '  R83V_MV computes b=A*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )

    allocate ( x(1:n) )
    call r8vec_indicator1 ( n, x )

    allocate ( ax_83v(1:m) )
    call r83v_mv ( m, n, a_83v, b_83v, c_83v, x, ax_83v )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )

    allocate ( ax_ge(1:m) )
    call r8ge_mv ( m, n, a_ge, x, ax_ge )

    call r8vec2_print ( m, ax_83v, ax_ge, '  Product comparison:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( ax_83v )
    deallocate ( ax_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )
    deallocate ( x )

  end do

  return
end
subroutine r83v_print_test ( )

!*****************************************************************************80
!
!! R83V_PRINT_TEST tests R83V_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRINT_TEST'
  write ( *, '(a)' ) '  R83V_PRINT prints an R83V matrix.'

  m = 5
  n = 5

  allocate ( a(1:min(m-1,n)) )
  allocate ( b(1:min(m,n)) )
  allocate ( c(1:min(m,n-1)) )

  call r83v_indicator ( m, n, a, b, c )
  call r83v_print ( m, n, a, b, c, '  R83V matrix:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine r83v_print_some_test ( )

!*****************************************************************************80
!
!! R83V_PRINT_SOME_TEST tests R83V_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R83V_PRINT_SOME prints some of an R83V matrix.'

  m = 5
  n = 5

  allocate ( a(1:min(m-1,n)) )
  allocate ( b(1:min(m,n)) )
  allocate ( c(1:min(m,n-1)) )

  call r83v_indicator ( m, n, a, b, c )
  call r83v_print_some ( m, n, a, b, c, 2, 2, 5, 4, '  Rows 2-5, Cols 2-4:' )

  deallocate ( a )
  deallocate ( b )
  deallocate ( c )

  return
end
subroutine r83v_random_test ( )

!*****************************************************************************80
!
!! R83V_RANDOM_TEST tests R83V_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer m
  integer n
  integer seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_RANDOM_TEST'
  write ( *, '(a)' ) '  R83V_RANDOM randomizes an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Random R83V matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r83v_res_test ( )

!*****************************************************************************80
!
!! R83V_RES_TEST tests R83V_RES.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: ax(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer m
  integer n
  real ( kind = rk ), allocatable :: r(:)
  integer seed
  real ( kind = rk ), allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_RES_TEST'
  write ( *, '(a)' ) '  R83V_RES computes b-A*x, where A is an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    allocate ( ax(1:m) )
    allocate ( r(1:m) )
    allocate ( x(1:n) )

    seed = 123456789
    call r83v_random ( m, n, seed, a, b, c )
    call r8vec_indicator1 ( n, x )
    call r83v_mv ( m, n, a, b, c, x, ax )
    call r83v_res ( m, n, a, b, c, x, ax, r )
    call r8vec_print ( m, r, '  Residual A*x-b:' )

    deallocate ( a )
    deallocate ( ax )
    deallocate ( b )
    deallocate ( c )
    deallocate ( r )
    deallocate ( x )

  end do

  return
end
subroutine r83v_to_r8ge_test ( )

!*****************************************************************************80
!
!! R83V_TO_R8GE_TEST tests R83V_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a_83v(:)
  real ( kind = rk ), allocatable :: a_ge(:,:)
  real ( kind = rk ), allocatable :: b_83v(:)
  real ( kind = rk ), allocatable :: c_83v(:)
  integer i
  integer m
  integer n
  integer seed

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R83V_TO_R8GE converts an R83V matrix to R8GE format.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a_83v(1:min(m-1,n)) )
    allocate ( b_83v(1:min(m,n)) )
    allocate ( c_83v(1:min(m,n-1)) )

    seed = 123456789
    call r83v_random ( m, n, seed, a_83v, b_83v, c_83v )
    call r83v_print ( m, n, a_83v, b_83v, c_83v, '  R83V matrix:' )

    allocate ( a_ge(1:m,1:n) )
    call r83v_to_r8ge ( m, n, a_83v, b_83v, c_83v, a_ge )
    call r8ge_print ( m, n, a_ge, '  R8GE matrix:' )

    deallocate ( a_83v )
    deallocate ( a_ge )
    deallocate ( b_83v )
    deallocate ( c_83v )

  end do

  return
end
subroutine r83v_to_r8vec_test ( )

!*****************************************************************************80
!
!! R83V_TO_R8VEC_TEST tests R83V_TO_R8VEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  integer i
  integer m
  integer n
  integer na

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TO_R8VEC_TEST'
  write ( *, '(a)' ) '  R83V_TO_R8VEC copies an R83V matrix to an R8VEC.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(min(m-1,n)) )
    allocate ( a2(min(m,  n)) )
    allocate ( a3(min(m,n-1)) )
    na = min ( m - 1, n ) + min ( m, n ) + min ( m, n - 1 )
    allocate ( a(na) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, "  R83V matrix A:" )

    call r83v_to_r8vec ( m, n, a1, a2, a3, a )
    call r8vec_print ( na, a, "  Vector version of A:" )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

  end do

  return
end
subroutine r83v_transpose_test ( )

!*****************************************************************************80
!
!! R83V_TRANSPOSE_TEST tests R83V_TRANSPOSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  real ( kind = rk ), allocatable :: b1(:)
  real ( kind = rk ), allocatable :: b2(:)
  real ( kind = rk ), allocatable :: b3(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_TRANSPOSE_TEST'
  write ( *, '(a)' ) &
    '  R83V_TRANSPOSE makes a transposed copy of an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r83v_indicator ( m, n, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  A:' )

    allocate ( b1(1:min(n-1,m)) )
    allocate ( b2(1:min(n,m)) )
    allocate ( b3(1:min(n,m-1)) )

    call r83v_transpose ( m, n, a1, a2, a3, b1, b2, b3 )
    call r83v_print ( n, m, b1, b2, b3, '  B = tranposed copy of A:' )

    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )
    deallocate ( b1 )
    deallocate ( b2 )
    deallocate ( b3 )

  end do

  return
end
subroutine r83v_zeros_test ( )

!*****************************************************************************80
!
!! R83V_ZEROS_TEST tests R83V_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ), allocatable :: c(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83V_ZEROS_TEST'
  write ( *, '(a)' ) '  R83V_ZEROS zeros an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:min(m-1,n)) )
    allocate ( b(1:min(m,n)) )
    allocate ( c(1:min(m,n-1)) )

    call r83v_zeros ( m, n, a, b, c )
    call r83v_print ( m, n, a, b, c, '  Zeroed R83V matrix:' )

    deallocate ( a )
    deallocate ( b )
    deallocate ( c )

  end do

  return
end
subroutine r8ge_to_r83v_test ( )

!*****************************************************************************80
!
!! R8GE_TO_R83V_TEST tests R8GE_TO_R83V.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  integer i
  integer m
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8GE_TO_R83V_TEST'
  write ( *, '(a)' ) '  R8GE_TO_R83V copies an R8GE matrix to an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    allocate ( a(1:m,1:n) )

    call r8ge_indicator ( m, n, a )
    call r8ge_print ( m, n, a, '  R8GE matrix A:' )

    allocate ( a1(1:min(m-1,n)) )
    allocate ( a2(1:min(m,n)) )
    allocate ( a3(1:min(m,n-1)) )

    call r8ge_to_r83v ( m, n, a, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  R83V copy of (some of ) matrix A:' )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

  end do

  return
end
subroutine r8vec_to_r83v_test ( )

!*****************************************************************************80
!
!! R8VEC_TO_R83V_TEST tests R8VEC_TO_R83V.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 February 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: a1(:)
  real ( kind = rk ), allocatable :: a2(:)
  real ( kind = rk ), allocatable :: a3(:)
  integer ahi
  integer bhi
  integer chi
  integer i
  integer m
  integer n
  integer na

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8VEC_TO_R83V_TEST'
  write ( *, '(a)' ) '  R8VEC_TO_R83V copies an R8VEC to an R83V matrix.'
  write ( *, '(a)' ) '  We check three cases, M<N, M=N, M>N.'

  do i = 1, 3

    if ( i == 1 ) then
      m = 3
      n = 5
    else if ( i == 2 ) then
      m = 5
      n = 5
    else if ( i == 3 ) then
      m = 5
      n = 3
    end if

    ahi = min ( m - 1, n )
    bhi = min ( m,     n )
    chi = min ( m,     n - 1 )

    na = ahi + bhi + chi
    allocate ( a(1:na) )

    call r8vec_indicator1 ( na, a )
    call r8vec_print ( na, a, '  R8VEC:' )

    allocate ( a1(1:ahi) )
    allocate ( a2(1:bhi) )
    allocate ( a3(1:chi) )

    call r8vec_to_r83v ( m, n, a, a1, a2, a3 )
    call r83v_print ( m, n, a1, a2, a3, '  R83V matrix:' )

    deallocate ( a )
    deallocate ( a1 )
    deallocate ( a2 )
    deallocate ( a3 )

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

