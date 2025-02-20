program main

!*****************************************************************************80
!
!! lcvt_test() tests lcvt().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer sample_function_cvt

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lcvt_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test lcvt().'

  do i = -1, 2
    sample_function_cvt = i
    call test01 ( sample_function_cvt )
  end do

  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LCVT_TEST():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( sample_function_cvt )

!*****************************************************************************80
!
!! TEST01 tests CVT, R8MAT_LATINIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: n = 25

  real ( kind = 8 ) generator(dim_num,n)
  integer i
  integer, parameter :: latin_steps = 3
  integer sample_function_cvt
  integer, parameter :: sample_function_init = 0
  integer, parameter :: sample_num_cvt = 100000
  integer, parameter :: sample_num_steps = 50

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  R8MAT_LATINIZE makes it a Latin Hypersquare.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we vary the sampling used during the'
  write ( *, '(a)' ) '  CVT Latin iteration.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension DIM_NUM =  ', dim_num
  write ( *, '(a,i12)' ) '  Number of generators =       ', n
  write ( *, '(a)' ) ' '

  if ( sample_function_init == -1 ) then
    write ( *, '(a)' ) '  Initialize using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_init == 0 ) then
    write ( *, '(a)' ) '  Initialize using UNIFORM.'
  else if ( sample_function_init == 1 ) then
    write ( *, '(a)' ) '  Initialize using HALTON.'
  else if ( sample_function_init == 2 ) then
    write ( *, '(a)' ) '  Initialize using GRID.'
  else if ( sample_function_init == 3 ) then
    write ( *, '(a)' ) '  USER will initialize data.'
  end if

  if ( sample_function_cvt == -1 ) then
    write ( *, '(a)' ) '  Sample using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_cvt == 0 ) then
    write ( *, '(a)' ) '  Sample using UNIFORM.'
  else if ( sample_function_cvt == 1 ) then
    write ( *, '(a)' ) '  Sample using HALTON.'
  else if ( sample_function_cvt == 2 ) then
    write ( *, '(a)' ) '  Sample using GRID.'
  end if

  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num_cvt
  write ( *, '(a,i8)' ) '  Number of sample steps = ', sample_num_steps

  do i = 1, latin_steps

    call cvt ( dim_num, n, sample_function_init, sample_function_cvt, &
      sample_num_cvt, sample_num_steps, generator )

    call r8mat_transpose_print ( dim_num, n, generator, '  After CVT steps:' )

    call r8mat_latinize ( dim_num, n, generator )

    call r8mat_transpose_print ( dim_num, n, generator, '  After Latin step:' )

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests CVT, R8MAT_LATINIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 2
  integer, parameter :: n = 25

  real ( kind = 8 ) generator(dim_num,n)
  integer i
  integer, parameter :: latin_steps = 3
  integer ngrid
  integer rank
  integer, parameter :: sample_function_cvt = 0
  integer, parameter :: sample_function_init = 3
  integer, parameter :: sample_num_cvt = 100000
  integer, parameter :: sample_num_steps = 50
  integer tuple(dim_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  CVT computes a Centroidal Voronoi Tessellation.'
  write ( *, '(a)' ) '  R8MAT_LATINIZE makes it a Latin Hypersquare.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we initialize the generators to'
  write ( *, '(a)' ) '  grid points; this is an unstable CVT solution.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Spatial dimension DIM_NUM =  ', dim_num
  write ( *, '(a,i12)' ) '  Number of generators =       ', n
  write ( *, '(a)' ) ' '

  if ( sample_function_init == -1 ) then
    write ( *, '(a)' ) '  Initialize using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_init == 0 ) then
    write ( *, '(a)' ) '  Initialize using UNIFORM.'
  else if ( sample_function_init == 1 ) then
    write ( *, '(a)' ) '  Initialize using HALTON.'
  else if ( sample_function_init == 2 ) then
    write ( *, '(a)' ) '  Initialize using GRID.'
  else if ( sample_function_init == 3 ) then
    write ( *, '(a)' ) '  USER will initialize data.'
  end if

  if ( sample_function_init == -1 ) then
    write ( *, '(a)' ) '  Sample using RANDOM_NUMBER (Fortran90 intrinsic).'
  else if ( sample_function_cvt == 0 ) then
    write ( *, '(a)' ) '  Sample using UNIFORM.'
  else if ( sample_function_cvt == 1 ) then
    write ( *, '(a)' ) '  Sample using HALTON.'
  else if ( sample_function_cvt == 2 ) then
    write ( *, '(a)' ) '  Sample using GRID.'
  end if

  write ( *, '(a,i8)' ) '  Number of sample points = ', sample_num_cvt
  write ( *, '(a,i8)' ) '  Number of sample steps =  ', sample_num_steps

  ngrid = 5

  do rank = 0, n-1
    call tuple_next_fast ( ngrid, dim_num, rank, tuple )
    generator(1:dim_num,rank+1) = real ( 2 * tuple(1:dim_num) - 1, kind = 8 ) &
      / real ( 2 * ngrid, kind = 8 )
  end do

  call r8mat_transpose_print ( dim_num, n, generator, &
    '  Initial generators (rows):' )

  do i = 1, latin_steps

    call cvt ( dim_num, n, sample_function_init, sample_function_cvt, &
      sample_num_cvt, sample_num_steps, generator )

    call r8mat_transpose_print ( dim_num, n, generator, &
      '  After CVT steps:' )

    call r8mat_latinize ( dim_num, n, generator )

    call r8mat_transpose_print ( dim_num, n, generator, &
      '  After Latin step:' )

  end do

  return
end
