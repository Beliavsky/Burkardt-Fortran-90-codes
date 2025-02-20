program main

!*****************************************************************************80
!
!! uniform_test() tests uniform().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'uniform_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test uniform().'

  call bvec_uniform_test ( )

  call c8_uniform_01_test ( )
  call c8mat_uniform_01_test ( )
  call c8vec_uniform_01_test ( )

  call ch_uniform_ab_test ( )

  call get_seed_test ( )

  call i4_seed_advance_test ( )

  call i4_uniform_0i_test ( )
  call i4_uniform_ab_test ( )
  call i4mat_uniform_ab_test ( )
  call i4vec_uniform_ab_test ( )

  call l4_uniform_test ( )
  call l4mat_uniform_test ( )
  call l4vec_uniform_test ( )

  call lcrg_anbn_test ( )
  call lcrg_seed_test ( )

  call r8_uniform_01_test ( )
  call r8_uniform_ab_test ( )
  call r8mat_uniform_01_test ( )
  call r8vec_uniform_01_test ( )
  call r8vec_uniform_ab_test ( )

  call r8col_uniform_abvec_test ( )
  call r8row_uniform_abvec_test ( )
  call r8vec_uniform_abvec_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNIFORM_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine bvec_uniform_test ( )

!*****************************************************************************80
!
!! BVEC_UNIFORM_TEST tests BVEC_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  logical b(n)
  integer i
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BVEC_UNIFORM_TEST'
  write ( *, '(a)' ) '  BVEC_UNIFORM computes a binary vector.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ''
  do i = 1, 10
    call bvec_uniform ( n, seed, b )
    call bvec_print ( n, b, '' )
  end do

  return
end
subroutine c8_uniform_01_test ( )

!*****************************************************************************80
!
!! C8_UNIFORM_01_TEST tests C8_UNIFORM_01.
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

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_uniform_01
  integer i
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  C8_UNIFORM_01 computes pseudorandom double precision'
  write ( *, '(a)' ) '  complex values uniformly distributed in the unit'
  write ( *, '(a)' ) '  circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.8,2x,f14.8)' ) i, c8_uniform_01 ( seed )
  end do

  return
end
subroutine c8mat_uniform_01_test ( )

!*****************************************************************************80
!
!! C8MAT_UNIFORM_01_TEST tests C8MAT_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: m = 4
  integer, parameter :: n = 2

  complex ( kind = ck ) c(m,n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8MAT_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  C8MAT_UNIFORM_01 computes pseudorandom complex values '
  write ( *, '(a)' ) '  uniformly distributed in the unit circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call c8mat_uniform_01 ( m, n, seed, c )

  call c8mat_print ( m, n, c, '  Uniform C8MAT:' )

  return
end
subroutine c8vec_uniform_01_test ( )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01_TEST tests C8VEC_UNIFORM_01.
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

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 10

  complex ( kind = ck ) c(n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8VEC_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  C8VEC_UNIFORM_01 computes pseudorandom '
  write ( *, '(a)' ) '  double precision complex values uniformly distributed '
  write ( *, '(a)' ) '  in the unit circle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call c8vec_uniform_01 ( n, seed, c )

  call c8vec_print ( n, c, '  The Uniform C8VEC:' )

  return
end
subroutine ch_uniform_ab_test ( )

!*****************************************************************************80
!
!! CH_UNIFORM_AB_TEST tests CH_UNIFORM_AB.
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

  character ch_uniform_ab
  character chi
  character clo
  integer i
  integer seed

  clo = 'A'
  chi = 'J'
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CH_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  CH_UNIFORM_AB computes pseudorandom characters '
  write ( *, '(a)' ) '  in an interval [CLO,CHI].'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The lower endpoint CLO = "' // clo // '".'
  write ( *, '(a)' ) '  The upper endpoint CHI = "' // chi // '".'
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,a1)' ) i, ch_uniform_ab ( clo, chi, seed )
  end do

  return
end
subroutine get_seed_test ( )

!*****************************************************************************80
!
!! GET_SEED_TEST tests GET_SEED.
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

  real ( kind = rk ) r8_uniform_01
  integer i
  integer j
  integer seed
  integer seed_old

  seed = 12345678
  seed_old = seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GET_SEED_TEST'
  write ( *, '(a)' ) '  GET_SEED picks an initial seed value for UNIFORM.'
  write ( *, '(a)' ) '  The value chosen should vary over time, because'
  write ( *, '(a)' ) '  the seed is based on reading the clock.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  This is just the "calendar" clock, which does'
  write ( *, '(a)' ) '  not change very fast, so calling GET_SEED several'
  write ( *, '(a)' ) '  times in a row may result in the same value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Initial SEED is ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM:'
  write ( *, '(a)' ) ' '

  do j = 1, 3
    write ( *, '(2x,f14.8)' ) r8_uniform_01 ( seed )
  end do

  do i = 1, 4

    do 

      call get_seed ( seed )

      if ( seed /= seed_old ) then
        seed_old = seed
        exit
      end if

    end do

    write ( *, '(a)' ) ' '
    write ( *, '(a,i12)' ) '  New seed from GET_SEED is ', seed
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM_01:'
    write ( *, '(a)' ) ' '

    do j = 1, 3
      write ( *, '(2x,f14.8)' ) r8_uniform_01 ( seed )
    end do

  end do

  return
end
subroutine i4_seed_advance_test ( )

!*****************************************************************************80
!
!! I4_SEED_ADVANCE_TEST tests I4_SEED_ADVANCE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i4_seed_advance
  integer seed
  integer seed_new
  integer step

  seed_new = 12345

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_SEED_ADVANCE_TEST'
  write ( *, '(a)' ) '  I4_SEED_ADVANCE advances the seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Step        SEED input       SEED output'
  write ( *, '(a)' ) ' '

  do step = 1, 10

    seed = seed_new
    seed_new = i4_seed_advance ( seed )

    write ( *, '(2x,i4,2x,i16,2x,i16)' ) step, seed, seed_new

  end do

  return
end
subroutine i4_uniform_0i_test ( )

!*****************************************************************************80
!
!! I4_UNIFORM_0I_TEST tests I4_UNIFORM_0I
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
  integer i4_uniform_0i
  real ( kind = rk ) mean
  integer, parameter :: n = 1000
  integer seed
  real ( kind = rk ) variance
  integer x(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_UNIFORM_0I_TEST'
  write ( *, '(a)' ) '  I4_UNIFORM_0I samples a uniform random'
  write ( *, '(a)' ) '  integer distribution in [0,2**31-1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  Starting with seed = ', seed

  do i = 1, n
    x(i) = i4_uniform_0i ( seed )
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First few values:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,i12)' ) i, x(i)
  end do

  mean = sum ( real ( x(1:n), kind = rk ) / real ( n, kind = rk ) )

  variance = sum ( ( real ( x(1:n), kind = rk ) - mean )**2 ) &
                / real ( n - 1, kind = rk )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
  write ( *, '(a,g14.6)' ) '  Average value was ', mean
  write ( *, '(a,i12)' ) '  Minimum value was ', minval ( x(1:n) )
  write ( *, '(a,i12)' ) '  Maximum value was ', maxval ( x(1:n) )
  write ( *, '(a,g14.6)' ) '  Variance was ', variance

  return
end
subroutine i4_uniform_ab_test ( )

!*****************************************************************************80
!
!! i4_uniform_ab_test() tests i4_uniform_ab().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: a = -100
  integer, parameter :: b = 200
  integer i
  integer i4_uniform_ab
  integer j
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_uniform_ab_test():'
  write ( *, '(a)' ) '  i4_uniform_ab() computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a)' ) ' '

  do i = 1, 20

    j = i4_uniform_ab ( a, b )

    write ( *, '(2x,i8,2x,i8)' ) i, j

  end do

  return
end
subroutine i4mat_uniform_ab_test ( )

!*****************************************************************************80
!
!! I4MAT_UNIFORM_AB_TEST tests I4MAT_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  integer, parameter :: a = -100
  integer, parameter :: b = 200
  integer seed
  integer v(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4MAT_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  I4MAT_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call i4mat_uniform_ab ( m, n, a, b, seed, v )

  call i4mat_print ( m, n, v, '  Uniform I4MAT:' )

  return
end
subroutine i4vec_uniform_ab_test ( )

!*****************************************************************************80
!
!! I4VEC_UNIFORM_AB_TEST tests I4VEC_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 20

  integer, parameter :: a = -100
  integer, parameter :: b = 200
  integer seed
  integer v(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4VEC_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  I4VEC_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
  write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call i4vec_uniform_ab ( n, a, b, seed, v )

  call i4vec_print ( n, v, '  Uniform I4VEC:' )

  return
end
subroutine l4_uniform_test ( )

!*****************************************************************************80
!
!! L4_UNIFORM_TEST tests L4_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 December 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  logical l4_uniform
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4_UNIFORM computes pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,l1)' ) i, l4_uniform ( seed )
  end do

  return
end
subroutine l4mat_uniform_test ( )

!*****************************************************************************80
!
!! L4MAT_UNIFORM_TEST tests L4MAT_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  logical l(m,n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4MAT_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4MAT_UNIFORM computes a vector of'
  write ( *, '(a)' ) '  pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call l4mat_uniform ( m, n, seed, l )

  call l4mat_print ( m, n, l, '  Uniform L4MAT:' )

  return
end
subroutine l4vec_uniform_test ( )

!*****************************************************************************80
!
!! L4VEC_UNIFORM_TEST tests L4VEC_UNIFORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 10

  logical l(n)
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'L4VEC_UNIFORM_TEST'
  write ( *, '(a)' ) '  L4VEC_UNIFORM computes a vector of'
  write ( *, '(a)' ) '  pseudorandom logical values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call l4vec_uniform ( n, seed, l )

  call l4vec_print ( n, l, '  Uniform L4VEC:' )

  return
end
subroutine lcrg_anbn_test ( )

!*****************************************************************************80
!
!! LCRG_ANBN_TEST tests LCRG_ANBN.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer a
  integer an
  integer b
  integer bn
  integer c
  integer j
  integer k
  integer n
  integer u
  integer v
  integer, allocatable, dimension ( : ) :: x
  integer, allocatable, dimension ( : ) :: y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LCRG_ANBN_TEST'
  write ( *, '(a)' ) '  LCRG_ANBN determines a linear congruential random'
  write ( *, '(a)' ) '  number generator equivalent to N steps of a given one.'
!
!  These parameters define the old (1969) IBM 360 random number generator:
!
  a = 16807
  b = 0
  c = 2147483647

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCRG parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  A = ', a
  write ( *, '(a,i12)' ) '  B = ', b
  write ( *, '(a,i12)' ) '  C = ', c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N             A             B'
  write ( *, '(a)' ) ' '

  do n = 0, 10
    call lcrg_anbn ( a, b, c, n, an, bn )
    write ( *, '(2x,i12,2x,i12,2x,i12)' ) n, an, bn
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                           N            In           Out'
  write ( *, '(a)' ) ' '

  k = 0
  u = 12345
  write ( *, '(2x,12x,2x,i12,2x,12x,2x,i12)' ) k, u
  do k = 1, 11
    call lcrg_evaluate ( a, b, c, u, v )
    write ( *, '(2x,12x,2x,i12,2x,i12,2x,i12)' ) k, u, v
    u = v
  end do
!
!  Now try to replicate these results using N procesors.
!
  n = 4
  allocate ( x(1:n) )
  allocate ( y(1:n) )

  call lcrg_anbn ( a, b, c, n, an, bn )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  LCRG parameters:'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  AN = ', an
  write ( *, '(a,i12)' ) '  BN = ', bn
  write ( *, '(a,i12)' ) '  C  = ', c
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             J             N            In           Out'
  write ( *, '(a)' ) ' '

  x(1) = 12345
  do j = 2, n
    call lcrg_evaluate ( a, b, c, x(j-1), x(j) )
  end do

  do j = 1, n
    write ( *, '(2x,i12,2x,i12,2x,12x,2x,i12)' ) j, j-1, x(j)
  end do

  do k = n + 1, 12, n
    do j = 1, n
      call lcrg_evaluate ( an, bn, c, x(j), y(j) )
      write ( *, '(2x,i12,2x,i12,2x,i12,2x,i12)' ) j, k+j-2, x(j), y(j)
      x(j) = y(j)
    end do
  end do

  deallocate ( x )
  deallocate ( y )

  return
end
subroutine lcrg_seed_test ( )

!*****************************************************************************80
!
!! LCRG_SEED_TEST tests LCRG_SEED.
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

  integer a
  integer b
  integer c
  integer i
  real ( kind = rk ) r8_uniform_01
  integer seed
  integer seed_in
  integer seed_lcrg
  integer seed_out
  integer seed_start
  real ( kind = rk ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LCRG_SEED_TEST'
  write ( *, '(a)' ) '  LCRG_SEED directly computes the updated value of a'
  write ( *, '(a)' ) '  seed used by an linear congruential random number'
  write ( *, '(a)' ) '  generator.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I          SEED          SEED          SEED    U'
  write ( *, '(a)' ) '                 Input        Output          LCRG'
  write ( *, '(a)' ) ' '
!
!  These parameters define the old (1969) IBM 360 random number generator:
!
  a = 16807
  b = 0
  c = 2147483647
!
!  This seed value was used in Pierre L'Ecuyer's article.
!
  seed_start = 12345

  seed = seed_start
!
!  Compute 1000 random numbers "the hard way", that is, sequentially.
!  Every now and then, call LCRG_SEED to compute SEED directly.
!
  do i = 1, 1000

    seed_in = seed
    u = r8_uniform_01 ( seed )
    seed_out = seed

    if ( i <= 10 .or. i == 100 .or. i == 1000 ) then

      call lcrg_seed ( a, b, c, i, seed_start, seed_lcrg )

      write ( *, '(2x,i8,2x,i12,2x,i12,2x,i12,2x,g14.6)' ) &
        i, seed_in, seed_out, seed_lcrg, u

    end if

  end do

  return
end
subroutine r8_uniform_01_test ( )

!*****************************************************************************80
!
!! R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
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

  real ( kind = rk ) r8_uniform_01
  integer i
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_01 ( seed )
  end do

  return
end
subroutine r8_uniform_ab_test ( )

!*****************************************************************************80
!
!! R8_UNIFORM_AB_TEST tests R8_UNIFORM_AB.
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

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r8_uniform_ab
  integer i
  integer seed

  a = 5.0D+00
  b = 10.0D+00
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  R8_UNIFORM_AB computes pseudorandom values '
  write ( *, '(a)' ) '  in an interval [A,B].'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
  write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  write ( *, '(a)' ) ' '
  do i = 1, 10
    write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_ab ( a, b, seed )
  end do

  return
end
subroutine r8col_uniform_abvec_test ( )

!*****************************************************************************80
!
!! R8COL_UNIFORM_ABVEC_TEST tests R8COL_UNIFORM_ABVEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ), dimension ( m ) :: a = (/ &
    0.0D+00, 0.20D+00, 10.0D+00, 52.0D+00, -1.0D+00 /)
  real ( kind = rk ), dimension ( m ) :: b = (/ &
    1.0D+00, 0.25D+00, 20.0D+00, 54.0D+00, +1.0D+00 /)
  integer i
  integer seed
  real ( kind = rk ) v(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8COL_UNIFORM_ABVEC_TEST'
  write ( *, '(a)' ) '  R8COL_UNIFORM_ABVEC computes a random R8COL.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8col_uniform_abvec ( m, n, a, b, seed, v )

  write ( *, '(a)' ) ''
  do i = 1, m
    write ( *, '(2x,f8.4,a3,4(2x,f8.4),2x,a3,f8.4)' ) &
      a(i), ':  ', v(i,1:n), '  :', b(i)
  end do
  
  return
end
subroutine r8mat_uniform_01_test ( )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01_TEST tests R8MAT_UNIFORM_01.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  integer seed
  real ( kind = rk ) v(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8MAT_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8MAT_UNIFORM_01 computes a random R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8mat_uniform_01 ( m, n, seed, v )

  call r8mat_print ( m, n, v, '  Uniform R8MAT' )
  
  return
end
subroutine r8row_uniform_abvec_test ( )

!*****************************************************************************80
!
!! R8ROW_UNIFORM_ABVEC_TEST tests R8ROW_UNIFORM_ABVEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 5

  real ( kind = rk ), dimension ( n ) :: a = (/ &
    0.0D+00, 0.20D+00, 10.0D+00, 52.0D+00, -1.0D+00 /)
  real ( kind = rk ), dimension ( n ) :: b = (/ &
    1.0D+00, 0.25D+00, 20.0D+00, 54.0D+00, +1.0D+00 /)
  integer i
  integer seed
  real ( kind = rk ) v(m,n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8ROW_UNIFORM_ABVEC_TEST'
  write ( *, '(a)' ) '  R8ROW_UNIFORM_ABVEC computes a random R8ROW.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8row_uniform_abvec ( m, n, a, b, seed, v )

  write ( *, '(a)' ) ''
  write ( *, '(5(2x,f8.4))' ) b(1:n)
  write ( *, '(a)' ) ''
  do i = 1, m
    write ( *, '(5(2x,f8.4))' ) v(i,1:n)
  end do
  write ( *, '(a)' ) ''
  write ( *, '(5(2x,f8.4))' ) a(1:n)
  
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
!    29 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer seed
  real ( kind = rk ) v(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_UNIFORM_01_TEST'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_01 computes a random R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8vec_uniform_01 ( n, seed, v )

  call r8vec_print ( n, v, '  Uniform R8VEC' )
  
  return
end
subroutine r8vec_uniform_ab_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_AB_TEST tests R8VEC_UNIFORM_AB.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 October 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer seed
  real ( kind = rk ) v(n)

  a = -1.0D+00
  b = +5.0D+00
  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_UNIFORM_AB_TEST'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_AB computes a random R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6)' ) '  ', a, ' <= X <= ', b
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8vec_uniform_ab ( n, a, b, seed, v )

  call r8vec_print ( n, v, '  Uniform R8VEC' )
  
  return
end
subroutine r8vec_uniform_abvec_test ( )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_ABVEC_TEST tests R8VEC_UNIFORM_ABVEC.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ), dimension ( n ) :: a = (/ &
    0.0D+00, 0.20D+00, 10.0D+00, 52.0D+00, -1.0D+00 /)
  real ( kind = rk ), dimension ( n ) :: b = (/ &
    1.0D+00, 0.25D+00, 20.0D+00, 54.0D+00, +1.0D+00 /)
  integer i
  integer seed
  real ( kind = rk ) v(n)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8VEC_UNIFORM_ABVEC_TEST'
  write ( *, '(a)' ) '  R8VEC_UNIFORM_ABVEC computes a random R8VEC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i12)' ) '  The initial seed is ', seed

  call r8vec_uniform_abvec ( n, a, b, seed, v )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I         A         X         B'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i4,2x,f8.4,2x,f8.4,2x,f8.4)' ) i, a(i), v(i), b(i)
  end do
  
  return
end
