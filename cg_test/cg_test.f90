program main

!*****************************************************************************80
!
!! cg_test() tests cg().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cg_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test cg().'

  call r8ge_cg_test ( )
  call r83_cg_test ( )
  call r83s_cg_test ( )
  call r83t_cg_test ( )
  call r8pbu_cg_test ( )
  call r8sd_cg_test ( )
  call r8sp_cg_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cg_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ge_cg_test ( )

!*****************************************************************************80
!
!! R8GE_CG_TEST tests R8GE_CG for a full storage matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8GE_CG_TEST'
  write ( *, '(a)' ) '  R8GE_CG applies CG to a full storage matrix.'

  n = 10
!
!  Choose a random symmetric positive definite matrix A.
!
  allocate ( a(1:n,1:n) )
  call spd_random ( n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r8ge_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r8ge_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r8ge_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83_cg_test ( )

!*****************************************************************************80
!
!! R83_CG_TEST tests R83_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83_CG_TEST'
  write ( *, '(a)' ) '  R83_CG applies CG to an R83 matrix.'

  n = 10
!
!  Let A be the -1 2 -1 matrix.
!
  allocate ( a(1:3,1:n) )
  call r83_dif2 ( n, n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r83_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83s_cg_test ( )

!*****************************************************************************80
!
!! R83S_CG_TEST tests R83S_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83S_CG_TEST'
  write ( *, '(a)' ) '  R83S_CG applies CG to an R83S matrix.'

  n = 10
!
!  Let A be the -1 2 -1 matrix.
!
  call r83s_dif2 ( n, n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83s_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r83s_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83s_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r83t_cg_test ( )

!*****************************************************************************80
!
!! R83T_CG_TEST tests R83T_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R83T_CG_TEST'
  write ( *, '(a)' ) '  R83T_CG applies CG to an R83T matrix.'

  n = 10
!
!  Let A be the -1 2 -1 matrix.
!
  allocate ( a(1:3,1:n) )
  call r83t_dif2 ( n, n, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r83t_mv ( n, n, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r83t_cg ( n, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r83t_res ( n, n, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r8pbu_cg_test ( )

!*****************************************************************************80
!
!! R8PBU_CG_TEST tests R8PBU_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer mu
  integer n
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8PBU_CG_TEST'
  write ( *, '(a)' ) '  R8PBU_CG applies CG to an R8PBU matrix.'

  n = 10
  mu = 1
!
!  Let A be the -1 2 -1 matrix.
!
  allocate ( a(1:mu+1,1:n) )
  call r8pbu_dif2 ( n, n, mu, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r8pbu_mv ( n, n, mu, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r8pbu_cg ( n, mu, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r8pbu_res ( n, n, mu, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r8sd_cg_test ( )

!*****************************************************************************80
!
!! R8SD_CG_TEST tests R8SD_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ), allocatable :: b(:)
  real ( kind = rk ) e_norm
  integer n
  integer ndiag
  integer, allocatable :: offset(:)
  real ( kind = rk ), allocatable :: r(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8SD_CG_TEST'
  write ( *, '(a)' ) '  R8SD_CG applies CG to an R8SD matrix.'

  n = 10

  ndiag = 2
  allocate ( offset(1:ndiag) )
  offset(1) = 0
  offset(2) = 1
!
!  Set A to the [-1 2 -1] matrix.
!
  allocate ( a(1:n,1:ndiag) )
  call r8sd_dif2 ( n, n, ndiag, offset, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r8sd_mv ( n, n, ndiag, offset, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r8sd_cg ( n, ndiag, offset, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r8sd_res ( n, n, ndiag, offset, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( offset )
  deallocate ( r )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
subroutine r8sp_cg_test ( )

!*****************************************************************************80
!
!! R8SP_CG_TEST tests R8SP_CG.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 June 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ), allocatable :: b(:)
  integer, allocatable :: col(:)
  real ( kind = rk ) e_norm
  integer n
  integer nz_num
  real ( kind = rk ), allocatable :: r(:)
  integer, allocatable :: row(:)
  real ( kind = rk ) r_norm
  real ( kind = rk ) r8vec_norm_affine
  real ( kind = rk ) r8vec_norm
  real ( kind = rk ), allocatable :: x1(:)
  real ( kind = rk ), allocatable :: x2(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'R8SP_CG_TEST'
  write ( *, '(a)' ) '  R8SP_CG applies CG to an R8SP matrix.'

  n = 10
  nz_num = 3 * n - 2
  allocate ( a(nz_num) )
  allocate ( col(nz_num) )
  allocate ( row(nz_num) )
!
!  Set A to the [-1 2 -1] matrix.
!
  call r8sp_dif2 ( n, n, nz_num, row, col, a )
!
!  Choose a random solution.
!
  allocate ( x1(1:n) )
  call random_number ( harvest = x1(1:n) )
!
!  Compute the corresponding right hand side.
!
  allocate ( b(1:n) )
  call r8sp_mv ( n, n, nz_num, row, col, a, x1, b )
!
!  Call the CG routine.
!
  allocate ( x2(1:n) )
  x2(1:n) = 1.0D+00
  call r8sp_cg ( n, nz_num, row, col, a, b, x2 )
!
!  Compute the residual.
!
  allocate ( r(1:n) )
  call r8sp_res ( n, n, nz_num, row, col, a, x2, b, r )
  r_norm = r8vec_norm ( n, r )
!
!  Compute the error.
!
  e_norm = r8vec_norm_affine ( n, x1, x2 )
!
!  Report.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Number of variables N = ', n
  write ( *, '(a,g14.6)' ) '  Norm of residual ||Ax-b|| = ', r_norm
  write ( *, '(a,g14.6)' ) '  Norm of error ||x1-x2|| = ', e_norm
!
!  Free memory.
!
  deallocate ( a )
  deallocate ( b )
  deallocate ( col )
  deallocate ( r )
  deallocate ( row )
  deallocate ( x1 )
  deallocate ( x2 )

  return
end
