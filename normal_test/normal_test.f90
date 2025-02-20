program main

!*****************************************************************************80
!
!! normal_test() tests normal().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'normal_test():'
  write ( *, '(a)' ) '  FORTRAN90 version;'
  write ( *, '(a)' ) '  Test normal().'

  call c8_normal_01_test ( )
  call c8vec_normal_01_test ( )

  call i4_normal_ab_test ( )

  call r8_normal_01_test ( )
  call r8_normal_ab_test ( )

  call r8mat_normal_01_test ( )
  call r8mat_normal_ab_test ( )

  call r8vec_normal_01_test ( )
  call r8vec_normal_ab_test ( )
  call r8vec_uniform_01_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'normal_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine c8_normal_01_test ( )

!*****************************************************************************80
!
!! c8_normal_01_test() tests c8_normal_01().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_normal_01
  integer i
  complex ( kind = ck ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'c8_normal_01_test():'
  write ( *, '(a)' ) '  c8_normal_01() computes pseudorandom double precision'
  write ( *, '(a)' ) '  complex values with Normal 01 circular distribution.'
  write ( *, '(a)' ) ''

  do i = 1, 10
    r = c8_normal_01 ( )
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, r
  end do

  return
end
subroutine c8vec_normal_01_test ( )

!*****************************************************************************80
!
!! c8vec_normal_01_test() tests c8vec_normal_01().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 10

  complex ( kind = ck ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'c8vec_normal_01_test():'
  write ( *, '(a)' ) '  c8vec_normal_01() computes a vector of Normal 01 values.'

  call c8vec_normal_01 ( n, x )

  call c8vec_print ( n, x, '  Vector of Normal 01 values:' )
  
  return
end
subroutine r8_normal_01_test ( )

!*****************************************************************************80
!
!! r8_normal_01_test() tests r8_normal_01().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r
  real ( kind = rk ) r8_normal_01

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_NORMAL_01_TEST'
  write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudonormal values '
  write ( *, '(a)' ) '  with mean 0.0 and standard deviation 1.0.'

  do i = 1, 10
    r = r8_normal_01 ( )
    write ( *, '(2x,i8,2x,g14.6)' ) i, r
  end do

  return
end
subroutine r8_normal_ab_test ( )

!*****************************************************************************80
!
!! R8_NORMAL_AB_TEST tests R8_NORMAL_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) mu
  real ( kind = rk ) r
  real ( kind = rk ) r8_normal_ab
  real ( kind = rk ) sigma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_NORMAL_AB_TEST'
  write ( *, '(a)' ) '  R8_NORMAL_AB computes pseudonormal values '
  write ( *, '(a)' ) '  with mean MU and standard deviation SIGMA.'

  mu = 10.0D+00
  sigma = 2.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma
  write ( *, '(a)' ) ' '

  do i = 1, 10
    r = r8_normal_ab ( mu, sigma )
    write ( *, '(2x,i8,2x,f14.8)' ) i, r
  end do

  return
end
subroutine r8mat_normal_01_test ( )

!*****************************************************************************80
!
!! R8MAT_NORMAL_01_TEST tests R8MAT_NORMAL_01.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) r(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8MAT_NORMAL_01_TEST'
  write ( *, '(a)' ) '  R8MAT_NORMAL_01 returns a matrix of Normal 01 values.'

  call r8mat_normal_01 ( m, n, r )

  call r8mat_print ( m, n, r, '  Matrix of Normal 01 values:' )
  
  return
end
subroutine r8mat_normal_ab_test ( )

!*****************************************************************************80
!
!! R8MAT_NORMAL_AB_TEST tests R8MAT_NORMAL_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) mu
  real ( kind = rk ) r(m,n)
  real ( kind = rk ) sigma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8MAT_NORMAL_AB_TEST'
  write ( *, '(a)' ) '  R8MAT_NORMAL_AB returns a matrix of Normal AB values.'

  mu = 100.0D+00
  sigma = 5.0D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma

  call r8mat_normal_ab ( m, n, mu, sigma, r )

  call r8mat_print ( m, n, r, '  Matrix of Normal AB values:' )
  
  return
end
subroutine r8vec_normal_01_test ( )

!*****************************************************************************80
!
!! r8vec_normal_01_test() tests r8vec_normal_01().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) r(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8vec_normal_01_test():'
  write ( *, '(a)' ) '  r8vec_normal_01() computes a vector of Normal 01 values.'

  call r8vec_normal_01 ( n, r )

  call r8vec_print ( n, r, '  Vector of Normal 01 values:' )
  
  return
end
subroutine r8vec_normal_ab_test ( )

!*****************************************************************************80
!
!! R8VEC_NORMAL_AB_TEST tests R8VEC_NORMAL_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 March 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) mu
  real ( kind = rk ) r(n)
  real ( kind = rk ) sigma

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_NORMAL_AB_TEST'
  write ( *, '(a)' ) '  R8VEC_NORMAL_AB computes a vector of Normal AB values.'

  mu = 15.0D+00
  sigma = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  MU = ', mu
  write ( *, '(a,g14.6)' ) '  SIGMA = ', sigma

  call r8vec_normal_ab ( n, mu, sigma, r )

  call r8vec_print ( n, r, '  Vector of Normal AB values:' )
  
  return
end
subroutine r8vec_uniform_01_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01_TEST tests R8VEC_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) r(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01 returns a random R8VEC '
  write ( *, '(a)' ) '  with entries in [0,1].'

  call r8vec_uniform_01 ( n, r )

  call r8vec_print ( n, r, '  Random R8VEC:' )

  return
end
