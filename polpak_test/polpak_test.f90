program main

!*****************************************************************************80
!
!! polpak_test() tests polpak().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polpak_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test polpak().'

  call agud_test ( )
  call align_enum_test ( )

  call bell_test ( )
  call bell_poly_coef_test ( )
  call benford_test ( )
  call bernoulli_number_test ( )
  call bernoulli_number2_test ( )
  call bernoulli_number3_test ( )
  call bernoulli_poly_test ( )
  call bernoulli_poly2_test ( )
  call bernstein_poly_test ( )
  call bpab_test ( )

  call cardan_poly_test ( )
  call cardan_poly_coef_test ( )
  call cardinal_cos_test ( )
  call cardinal_sin_test ( )
  call catalan_test ( )
  call catalan_row_next_test ( )
  call charlier_test ( )
  call cheby_t_poly_test ( )
  call cheby_t_poly_coef_test ( )
  call cheby_t_poly_zero_test ( )
  call cheby_u_poly_test ( )
  call cheby_u_poly_coef_test ( )
  call cheby_u_poly_zero_test ( )
  call chebyshev_discrete_test ( )
  call collatz_count_test ( )
  call collatz_count_max_test ( )
  call comb_row_next_test ( )
  call commul_test ( )
  call complete_symmetric_poly_test ( )
  call conway_sequence_test ( )
  call cos_power_int_test ( )

  call delannoy_test ( )
  call domino_tiling_num_test ( )

  call euler_mascheroni_test ( )
  call euler_number_test ( )
  call euler_number2_test ( )
  call euler_poly_test ( )
  call eulerian_test ( )

  call f_hofstadter_test ( )
  call fibonacci_direct_test ( )
  call fibonacci_floor_test ( )
  call fibonacci_recursive_test ( )

  call g_hofstadter_test ( )
  call gegenbauer_poly_test ( )
  call gen_hermite_poly_test ( )
  call gen_laguerre_poly_test ( )
  call gud_test ( )

  call h_hofstadter_test ( )
  call harmonic_test ( )
  call harmonic_estimate_test ( )
  call hermite_poly_phys_test ( )
  call hermite_poly_phys_coef_test ( )

  call i4_choose_test ( )
  call i4_factor_test ( )
  call i4_factorial_test ( )
  call i4_factorial2_test ( )
  call i4_is_fibonacci_test ( )
  call i4_is_prime_test ( )
  call i4_is_triangular_test ( )
  call i4_partition_distinct_count_test ( )
  call i4_to_triangle_lower_test ( )

  call jacobi_poly_test ( )
  call jacobi_symbol_test ( )

  call krawtchouk_test ( )

  call laguerre_associated_test ( )
  call laguerre_poly_test ( )
  call laguerre_poly_coef_test ( )
  call lambert_w_test ( )
  call lambert_w_estimate_test ( )
  call legendre_associated_test ( )
  call legendre_associated_normalized_test ( )
  call legendre_function_q_test ( )
  call legendre_poly_test ( )
  call legendre_poly_coef_test ( )
  call legendre_symbol_test ( )
  call lerch_test ( )
  call lock_test ( )

  call meixner_test ( )
  call mertens_test ( )
  call moebius_test ( )
  call motzkin_test ( )

  call normal_01_cdf_inverse_test ( )

  call omega_test ( )

  call pentagon_num_test ( )
  call phi_test ( )
  call pi_estimate_test ( )
  call plane_partition_num_test ( )
  call poly_bernoulli_test ( )
  call poly_coef_count_test ( )
  call prime_test ( )
  call pyramid_num_test ( )
  call pyramid_square_num_test ( )

  call r8_agm_test ( )
  call r8_beta_test ( )
  call r8_choose_test ( )
  call r8_cube_root_test ( )
  call r8_erf_test ( )
  call r8_erf_inverse_test ( )
  call r8_factorial_test ( )
  call r8_factorial_log_test ( )
  call r8_hyper_2f1_test ( )
  call r8_psi_test ( )
  call r8poly_degree_test ( )
  call r8poly_print_test ( )
  call r8poly_value_horner_test ( )

  call sigma_test ( )
  call simplex_num_test ( )
  call sin_power_int_test ( )
  call slices_test ( )
  call spherical_harmonic_test ( )
  call stirling_estimate_test ( )
  call stirling1_table_test ( )
  call stirling2_number_test ( )
  call stirling2_table_test ( )

  call tau_test ( )
  call tetrahedron_num_test ( )
  call triangle_num_test ( )
  call triangle_lower_to_i4_test ( )
  call tribonacci_direct_test ( )
  call tribonacci_recursive_test ( )
  call tribonacci_roots_test ( )
  call trinomial_test ( )

  call v_hofstadter_test ( )
  call vibonacci_test ( )

  call zeckendorf_test ( )
  call zernike_poly_test ( )
  call zernike_poly_coef_test ( )
  call zeta_m1_test ( )
  call zeta_naive_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'polpak_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine agud_test ( )

!*****************************************************************************80
!
!! agud_test() tests agud().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) agud
  real ( kind = rk ) g
  real ( kind = rk ) gud
  integer i
  real ( kind = rk ) x
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'agud_test()'
  write ( *, '(a)' ) '  agud() computes the inverse Gudermannian;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       X              GUD(X)       AGUD(GUD(X))'
  write ( *, '(a)' ) ' '

  do i = 0, 10
    x = 1.0D+00 + real ( i, kind = rk ) / 5.0D+00
    g = gud ( x )
    x2 = agud ( g )
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, g, x2
  end do

  return
end
subroutine align_enum_test ( )

!*****************************************************************************80
!
!! align_enum_test() tests align_enum().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m_max = 10
  integer, parameter :: n_max = 10

  integer align_enum
  integer i
  integer j
  integer table(0:m_max,0:n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'align_enum_test():'
  write ( *, '(a)' ) '  align_enum() counts the number of possible'
  write ( *, '(a)' ) '  alignments of two biological sequences.'

  do i = 0, m_max
    do j = 0, n_max
      table(i,j) = align_enum ( i, j )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Alignment enumeration table:'
  write ( *, '(a)' ) ' '
  write ( *, '(4x,5i5,6i8)' ) ( i, i = 0, n_max )
  write ( *, '(a)' ) ' '
  do i = 0, m_max
    write ( *, '(2x,i2,5i5,6i8)' ) i, table(i,0:n_max)
  end do

  return
end
subroutine bell_test ( )

!*****************************************************************************80
!
!! bell_test() tests bell().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2(0:10)
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bell_test():'
  write ( *, '(a)' ) '  bell() computes Bell numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call bell_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call bell ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2(n)

  end do
 
  return
end
subroutine bell_poly_coef_test ( )

!*****************************************************************************80
!
!! bell_poly_coef_test() tests bell_poly_coef().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 10

  integer, allocatable :: c(:)
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bell_poly_coef_test():'
  write ( *, '(a)' ) '  bell_poly_coef() returns the coefficients of a'
  write ( *, '(a)' ) '  Bell polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Table of polynomial coefficients:'
  write ( *, '(a)' ) ' '

  do n = 0, n_max
    allocate ( c(0:n) )
    call bell_poly_coef ( n, c )
    write ( *, '(2x,i2,a1,11i6)' ) n, ':', c(0:n)
    deallocate ( c )
  end do

  return
end
subroutine benford_test ( )

!*****************************************************************************80
!
!! benford_test() tests benford().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) benford
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'benford_test():'
  write ( *, '(a)' ) '  benford(i) is the Benford probability of the'
  write ( *, '(a)' ) '  initial digit sequence I.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I,  BENFORD(I)'
  write ( *, '(a)' ) ' '

  do i = 1, 9
    write ( *, '(2x,i2,2x,g14.6)' )  i, benford(i)
  end do
 
  return
end
subroutine bernoulli_number_test ( )

!*****************************************************************************80
!
!! bernolli_number_test() tests bernoulli_number().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c0
  real ( kind = rk ) c1(0:30)
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernoulli_number_test():'
  write ( *, '(a)' ) '  bernoulli_number() computes Bernoulli numbers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     Bernoulli'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number ( n, c1 )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, c0, c1(n)

  end do
 
  return
end
subroutine bernoulli_number2_test ( )

!*****************************************************************************80
!
!! bernoulli_number2_test() tests bernoulli_number2().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c0
  real ( kind = rk ) c1(0:30)
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernoulli_number2_test():'
  write ( *, '(a)' ) '  bernoulli_number2() computes Bernoulli numbers;'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     Bernoulli2'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number2 ( n, c1 )
 
    write ( *, '(2x,i4,2g14.6)' ) n, c0, c1(n)

  end do
 
  return
end
subroutine bernoulli_number3_test ( )

!*****************************************************************************80
!
!! bernoulli_number3_test() tests bernoulli_number3().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c0
  real ( kind = rk ) c1
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernoulli_number3_test():'
  write ( *, '(a)' ) '  bernoulli_number3() computes Bernoulli numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   I      Exact     BERNOULLI3'
  write ( *, '(a)' ) ' '
  
  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c0 )

    if ( n_data == 0 ) then
      exit
    end if

    call bernoulli_number3 ( n, c1 )

    write ( *, '(2x,i4,2g14.6)' ) n, c0, c1

  end do
 
  return
end
subroutine bernoulli_poly_test ( )

!*****************************************************************************80
!
!! bernoulli_poly_test() tests bernoulli_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bx
  integer i
  integer, parameter :: n = 15
  real ( kind = rk ) x

  x = 0.2D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernoulli_poly_test():'
  write ( *, '(a)' ) '  bernoulli_poly() evaluates Bernoulli polynomials;'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I          BX'
  write ( *, '(a)' ) ' '
 
  do i = 1, n
    call bernoulli_poly ( i, x, bx )
    write ( *, '(2x,i2,2x,g16.8)' ) i, bx
  end do
 
  return
end
subroutine bernoulli_poly2_test ( )

!*****************************************************************************80
!
!! bernoulli_poly2_test() tests bernoulli_poly2().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bx
  integer i
  integer, parameter :: n = 15
  real ( kind = rk ) x

  x = 0.2D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernoulli_poly2_test():'
  write ( *, '(a)' ) '  bernoulli_poly2() evaluates Bernoulli polynomials. '
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  X = ', x
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  I          BX'
  write ( *, '(a)' ) ' '
 
  do i = 1, n
    call bernoulli_poly2 ( i, x, bx )
    write ( *, '(2x,i2,2x,2g16.8)' ) i, bx
  end do
 
  return
end
subroutine bernstein_poly_test ( )

!*****************************************************************************80
!
!! bernstein_poly_test() tests bernstein_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) b
  real ( kind = rk ) bvec(0:10)
  integer k
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bernstein_poly_test():'
  write ( *, '(a)' ) '  bernstein_poly() evaluates the Bernstein polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N   K   X       Exact        B(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_01_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    call bernstein_poly ( n, x, bvec )

    write ( *, '(2x,i4,i4,f7.4,2g14.6)' ) n, k, x, b, bvec(k)

  end do

  return
end
subroutine bpab_test ( )

!*****************************************************************************80
!
!! bpab_test() tests bpab().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
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
  real ( kind = rk ) bern(0:n)
  integer i
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bpab_test():'
  write ( *, '(a)' ) '  bpab() evaluates Bernstein polynomials.'
  write ( *, '(a)' ) ' '

  x = 0.3D+00
  a = 0.0D+00
  b = 1.0D+00
  call bpab ( n, a, b, x, bern )
 
  write ( *, '(a,i4)' ) '  The Bernstein polynomials of degree ', n
  write ( *, '(a,g14.6)' ) '  based on the interval from ', a
  write ( *, '(a,g14.6)' ) '  to ', b
  write ( *, '(a,g14.6)' ) '  evaluated at X = ', x
  write ( *, '(a)' ) ' '
 
  do i = 0, n
    write ( *, '(2x,i4,2x,g14.6)' )  i, bern(i)
  end do
 
  return
end
subroutine cardan_poly_test ( )

!*****************************************************************************80
!
!! cardan_poly_test() tests cardan_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 10

  real ( kind = rk ) c(0:n_max)
  real ( kind = rk ) cx1
  real ( kind = rk ) cx2(0:n_max)
  integer n
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) s
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cardan_poly_test():'
  write ( *, '(a)' ) '  cardan_poly() evaluates a Cardan polynomial directly.'

  n = n_max
  s = 0.5D+00
  x = 0.25D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Compare CARDAN_POLY_COEF + R8POLY_VAL_HORNER'
  write ( *, '(a)' ) '  versus CARDAN_POLY alone.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Evaluate polynomials at X = ', x
  write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Order, Horner, Direct'
  write ( *, '(a)' ) ' '

  call cardan_poly ( n, x, s, cx2 )

  do n = 0, n_max

    call cardan_poly_coef ( n, s, c )
    cx1 = r8poly_value_horner ( n, c, x )

    write ( *, '(2x,i2,2g14.6)' ) n, cx1, cx2(n)

  end do

  return
end
subroutine cardan_poly_coef_test ( )

!*****************************************************************************80
!
!! CARDAN_POLY_COEF_test() tests CARDAN_POLY_COEF().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 March 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 10

  real ( kind = rk ) c(0:n_max)
  integer n
  real ( kind = rk ) s

  s = 1.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CARDAN_POLY_COEF_test():'
  write ( *, '(a)' ) '  CARDAN_POLY_COEF returns the coefficients of a'
  write ( *, '(a)' ) '  Cardan polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Table of polynomial coefficients:'
  write ( *, '(a)' ) ' '

  do n = 0, n_max
    call cardan_poly_coef ( n, s, c )
    write ( *, '(2x,i2,11f7.0)' ) n, c(0:n)
  end do

  return
end
subroutine cardinal_cos_test ( )

!*****************************************************************************80
!
!! CARDINAL_COS_test() tests CARDINAL_COS().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 May 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 11
  integer, parameter :: n = m + 2

  real ( kind = rk ) c(0:m+1,0:m+1)
  integer i
  integer j
  real ( kind = rk ) :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t(0:m+1)
  real ( kind = rk ) t_hi
  real ( kind = rk ) t_lo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CARDINAL_COS_test():'
  write ( *, '(a)' ) '  CARDINAL_COS evaluates cardinal cosine functions.'
  write ( *, '(a)' ) '  Ci(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).'
  write ( *, '(a)' ) '  A simple check of all pairs should form the identity matrix.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The CARDINAL_COS test matrix:'
  write ( *, '(a)' ) ''

  t_lo = 0.0D+00
  t_hi = r8_pi
  call r8vec_linspace ( n, t_lo, t_hi, t )

  do j = 0, m + 1
    call cardinal_cos ( j, m, n, t, c(0:m+1,j) )
  end do

  do i = 0, m + 1
    write ( *, '(13(2x,f4.1))' ) c(i,0:m+1)
  end do

  return
end
subroutine cardinal_sin_test ( )

!*****************************************************************************80
!
!! CARDINAL_SIN_test() tests CARDINAL_SIN().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 May 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 11

  integer i
  integer j
  real ( kind = rk ) :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) s(0:m+1,0:m+1)
  real ( kind = rk ) t(0:m+1)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'CARDINAL_SIN_test():'
  write ( *, '(a)' ) '  CARDINAL_SIN evaluates cardinal sine functions.'
  write ( *, '(a)' ) '  Si(Tj) = Delta(i,j), where Tj = cos(pi*i/(n+1)).'
  write ( *, '(a)' ) '  A simple check of all pairs should form the identity matrix.'

  call r8vec_linspace ( m + 2, 0.0D+00, r8_pi, t )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  The CARDINAL_SIN test matrix:'
  write ( *, '(a)' ) ''
  do j = 0, m + 1
    call cardinal_sin ( j, m, m + 2, t, s(0:m+1,j) )
  end do

  do i = 0, m + 1
    write ( *, '(13(2x,f4.1))' ) s(i,0:m+1)
  end do

  return
end
subroutine catalan_test ( )

!*****************************************************************************80
!
!! CATALAN_test() tests CATALAN().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2(0:10)
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CATALAN_test():'
  write ( *, '(a)' ) '  CATALAN computes Catalan numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call catalan_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call catalan ( n, c2 )

    write ( *, '(2x,i4,2i8)' ) n, c, c2(n)

  end do
 
  return
end
subroutine catalan_row_next_test ( )

!*****************************************************************************80
!
!! CATALAN_ROW_NEXT_test() tests CATALAN_ROW_NEXT().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer c(0:n)
  integer i
  integer ido

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CATALAN_ROW_NEXT_test():'
  write ( *, '(a)' ) '  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First, compute row 7:'

  ido = 0
  i = 7
  call catalan_row_next ( ido, i, c )
  write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now compute rows one at a time:'
  write ( *, '(a)' ) ' '

  ido = 0
 
  do i = 0, n
    call catalan_row_next ( ido, i, c )
    ido = 1
    write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)
  end do
 
  return
end
subroutine charlier_test ( )

!*****************************************************************************80
!
!! CHARLIER_test() tests CHARLIER().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 5

  real ( kind = rk ) a
  real ( kind = rk ), dimension ( test_num ) :: a_test = (/ &
    0.25D+00, 0.5D+00, 1.0D+00, 2.0D+00, 10.0D+00 /)
  integer i
  integer j
  integer test
  real ( kind = rk ) x
  real ( kind = rk ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHARLIER_TEST:'
  write ( *, '(a)' ) '  CHARLIER evaluates Charlier polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      A         X        P(N,A,X)'
  write ( *, '(a)' ) ' '

  do test = 1, test_num

    a = a_test(test)

    do j = 0, 5

      x = real ( j, kind = rk ) / 2.0D+00

      call charlier ( n, a, x, value )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,g14.6)' ) i, a, x, value(i)

      end do

    end do

  end do

  return
end
subroutine cheby_t_poly_test ( )

!*****************************************************************************80
!
!! CHEBY_T_POLY_test() tests CHEBY_T_POLY().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_T_POLY_test():'
  write ( *, '(a)' ) '  CHEBY_T_POLY evaluates the Chebyshev T polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N      X        Exact F       T(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_t_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call cheby_t_poly ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine cheby_t_poly_coef_test ( )

!*****************************************************************************80
!
!! CHEBY_T_POLY_COEF_test() tests CHEBY_T_POLY_COEF().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_T_POLY_COEF_test():'
  write ( *, '(a)' ) '  CHEBY_T_POLY_COEF determines ' // &
    'the Chebyshev T polynomial coefficients.'

  call cheby_t_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  T(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine cheby_t_poly_zero_test ( )

!*****************************************************************************80
!
!! CHEBY_T_POLY_ZERO_test() tests CHEBY_T_POLY_ZERO().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 4

  real ( kind = rk ), allocatable :: fx(:,:)
  integer i
  integer n
  real ( kind = rk ) z(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_T_POLY_ZERO_TEST:'
  write ( *, '(a)' ) '  CHEBY_T_POLY_ZERO returns zeroes of T(N)(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        T(N)(X)'

  do n = 1, n_max

    call cheby_t_poly_zero ( n, z )

    allocate ( fx(1:n,0:n) )

    call cheby_t_poly ( n, n, z, fx )

    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(i,n)
    end do

    deallocate ( fx )

  end do

  return
end
subroutine cheby_u_poly_test ( )

!*****************************************************************************80
!
!! CHEBY_U_POLY_test() tests CHEBY_U_POLY().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    10 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_U_POLY_test():'
  write ( *, '(a)' ) '  CHEBY_U_POLY evaluates the Chebyshev U polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N      X        Exact F       U(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_u_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call cheby_u_poly ( 1, n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine cheby_u_poly_coef_test ( )

!*****************************************************************************80
!
!! CHEBY_U_POLY_COEF_test() tests CHEBY_U_POLY_COEF().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_U_POLY_COEF_test():'
  write ( *, '(a)' ) '  CHEBY_U_POLY_COEF determines ' // &
    'the Chebyshev U polynomial coefficients.'

  call cheby_u_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  U(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine cheby_u_poly_zero_test ( )

!*****************************************************************************80
!
!! CHEBY_U_POLY_ZERO_test() tests CHEBY_U_POLY_ZERO().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 4

  real ( kind = rk ) fx(0:n_max)
  integer i
  integer n
  real ( kind = rk ) z(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_U_POLY_ZERO_TEST:'
  write ( *, '(a)' ) '  CHEBY_U_POLY_ZERO returns zeroes of the U(N)(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      X        U(N)(X)'
  write ( *, '(a)' ) ' '

  do n = 1, n_max

    call cheby_u_poly_zero ( n, z )

    do i = 1, n

      call cheby_u_poly ( 1, n, z(i), fx )

      write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(n)

    end do

    write ( *, '(a)' ) ' '

  end do

  return
end
subroutine chebyshev_discrete_test ( )

!*****************************************************************************80
!
!! CHEBYSHEV_DISCRETE_test() tests CHEBYSHEV_DISCRETE().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer i
  integer j
  integer m
  real ( kind = rk ) x
  real ( kind = rk ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBYSHEV_DISCRETE_TEST:'
  write ( *, '(a)' ) &
    '  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      M         X        T(N,M,X)'
  write ( *, '(a)' ) ' '

  m = 5

  do j = 0, 5

    x = real ( j, kind = rk ) / 2.0D+00

    call chebyshev_discrete ( n, m, x, value )

    write ( *, '(a)' ) ' '

    do i = 0, n

      write ( *, '(2x,i8,2x,i8,2x,f8.4,2x,g14.6)' ) i, m, x, value(i)

    end do

  end do

  return
end
subroutine collatz_count_test ( )

!*****************************************************************************80
!
!! collatz_count_test() tests collatz_count().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer collatz_count
  integer count
  integer count2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLLATZ_COUNT_TEST:'
  write ( *, '(a)' ) '  COLLATZ_COUNT(N) counts the length of the'
  write ( *, '(a)' ) '  Collatz sequence beginning with N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       COUNT(N)     COUNT(N)'
  write ( *, '(a)' ) '              (computed)    (table)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call collatz_count_values ( n_data, n, count )

    if ( n_data == 0 ) then
      exit
    end if

    count2 = collatz_count ( n )

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, count, count2

  end do

  return
end
subroutine collatz_count_max_test ( )

!*****************************************************************************80
!
!! collatz_count_max_test() tests collatz_count_max().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i_max
  integer j_max
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLLATZ_COUNT_MAX_TEST:'
  write ( *, '(a)' ) '  COLLATZ_COUNT_MAX(N) returns the length of the'
  write ( *, '(a)' ) '  longest Collatz sequence from 1 to N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N     I_MAX     J_MAX'
  write ( *, '(a)' ) ' '

  n = 10

  do while ( n <= 100000 )

    call collatz_count_max ( n, i_max, j_max )

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, i_max, j_max

    n = n * 10

  end do

  return
end
subroutine comb_row_next_test ( )

!*****************************************************************************80
!
!! comb_row_next_test() tests comb_row_next().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 10

  integer c(0:n_max)
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMB_ROW_NEXT_test():'
  write ( *, '(a)' ) '  COMB_ROW computes a row of Pascal''s triangle.'
  write ( *, '(a)' ) ' '
 
  do n = 0, n_max
    call comb_row_next ( n, c )
    write ( *, '(2x,i2,2x,11i5)' ) n, c(0:n)
  end do
 
  return
end
subroutine commul_test ( )

!*****************************************************************************80
!
!! commul_test() tests commul().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer factor(4)
  integer i
  integer ncomb
  integer nfactor

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COMMUL_test():'
  write ( *, '(a)' ) '  COMMUL computes a multinomial coefficient.'
  write ( *, '(a)' ) ' '

  n = 8
  nfactor = 2
  factor(1) = 6
  factor(2) = 2

  call commul ( n, nfactor, factor, ncomb ) 

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  n = 8
  nfactor = 3
  factor(1) = 2
  factor(2) = 2
  factor(3) = 4
  call commul ( n, nfactor, factor, ncomb )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  n = 13
  nfactor = 4
  factor(1) = 5
  factor(2) = 3
  factor(3) = 3
  factor(4) = 2
  call commul ( n, nfactor, factor, ncomb )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  N = ', n
  write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
  do i = 1, nfactor
    write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
  end do
  write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

  return
end
subroutine complete_symmetric_poly_test ( )

!*****************************************************************************80
!
!! complete_symmetric_poly_test() tests complete_symmetric_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  integer nn
  integer rr
  real ( kind = rk ) value
  real ( kind = rk ), dimension ( n ) :: x = (/ &
    1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'COMPLETE_SYMMETRIC_POLY_test():'
  write ( *, '(a)' ) '  COMPLETE_SYMMETRIC_POLY evaluates a complete symmetric.'
  write ( *, '(a)' ) '  polynomial in a given set of variables X.'
 
  call r8vec_print ( n, x, '  Variable vector X:' );

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   N\R    0       1       2       3       4       5'
  write ( *, '(a)' ) ''

  do nn = 0, n
    write ( *, '(2x,i2)', ADVANCE = 'NO' ) nn
    do rr = 0, 5
      call complete_symmetric_poly ( nn, rr, x, value )
      write ( *, '(2x,f6.0)', ADVANCE = 'NO' ) value
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine conway_sequence_test ( )

!*****************************************************************************80
!
!! conway_sequence_test() tests conway_sequence().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 64

  integer A(n)
  integer i

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'conway_sequence_test():'
  write ( *, '(a)' ) '  conway_sequence() returns the first N entries of'
  write ( *, '(a)' ) '  the Conway challenge sequence.'
 
  call conway_sequence ( n, A )

  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(i2,2x,i2)' ) i, A(i)
  end do

  return
end
subroutine cos_power_int_test ( )

!*****************************************************************************80
!
!! cos_power_int_test() tests cos_power_int().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) cos_power_int
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_POWER_INT_TEST:'
  write ( *, '(a)' ) '  COS_POWER_INT returns values of '
  write ( *, '(a)' ) '  the integral of COS(X)^N from A to B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      A         B          N      Exact           Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = cos_power_int ( a, b, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) a, b, n, fx, fx2

  end do

  return
end
subroutine delannoy_test ( )

!*****************************************************************************80
!
!! delannoy_test() tests delannoy().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = 8

  integer a(0:m,0:n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DELANNOY_test():'
  write ( *, '(a)' ) '  DELANNOY computes the Delannoy numbers A(0:M,0:N).'
  write ( *, '(a)' ) '  A(M,N) counts the paths from (0,0) to (M,N).'
  write ( *, '(a)' ) ' '

  call delannoy ( m, n, a )

  do i = 0, m
    write ( *, '(2x,i4,2x,5i4,3i8,i10)' )  i, a(i,0:n)
  end do
 
  return
end
subroutine domino_tiling_num_test ( )

!*****************************************************************************80
!
!! domino_tiling_num_test() tests domino_tiling_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer value

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DOMINO_TILING_NUM_TEST:'
  write ( *, '(a)' ) '  DOMINO_TILING_NUM returns the number of tilings of an'
  write ( *, '(a)' ) '  MxN rectangle by dominoes.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   M   N    Tilings'

  do m = 1, 8
    write ( *, '(a)' ) ''
    do n = 1, m
      call domino_tiling_num ( m, n, value )
      write ( *, '(2x,i2,2x,i2,2x,i10)' ) m, n, value
    end do
  end do

  return
end
subroutine euler_mascheroni_test ( )

!*****************************************************************************80
!
!! euler_mascheroni_test() tests euler_mascheroni().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) euler_mascheroni
  real ( kind = rk ) g
  real ( kind = rk ) g_approx
  integer i
  integer n
  real ( kind = rk ) n_r8
  integer test

  g = euler_mascheroni ( )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'euler_mascheroni_test():'
  write ( *, '(a)' ) '  euler_mascheroni() returns the Euler-Mascheroni constant'
  write ( *, '(a)' ) '  sometimes denoted by "gamma".'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  gamma = limit ( N -> oo ) ( sum ( 1 <= I <= N ) 1 / I ) - log ( N )'
  write ( *, '(a)' ) ''
  write ( *, '(a,g24.16)' ) '  Numerically, g = ', g
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         N      Partial Sum    |gamma - partial sum|'
  write ( *, '(a)' ) ''

  n = 1
  do test = 0, 20
    n_r8 = real ( n, kind = rk )
    g_approx = - log ( n_r8 )    
    do i = 1, n
      g_approx = g_approx + 1.0D+00 / real ( i, kind = rk )
    end do
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, g_approx, abs ( g_approx - g )
    n = n * 2
  end do

  return
end
subroutine euler_number_test ( )

!*****************************************************************************80
!
!! euler_number_test() tests euler_number().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c1
  integer c2(0:12)
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULER_NUMBER_test():'
  write ( *, '(a)' ) '  EULER_NUMBER computes Euler numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N       exact   EULER_NUMBER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c1 )

    if ( n_data == 0 ) then
      exit
    end if

    call euler_number ( n, c2 )

    write ( *, '(2x,i4,2i12,g14.6)' ) n, c1, c2(n)

  end do
 
  return
end
subroutine euler_number2_test ( )

!*****************************************************************************80
!
!! euler_number2_test() tests euler_number2().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c1
  real ( kind = rk ) c2
  real ( kind = rk ) euler_number2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULER_NUMBER2_test():'
  write ( *, '(a)' ) '  EULER_NUMBER2 computes Euler numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N       exact   EULER_NUMBER2'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c1 )

    if ( n_data == 0 ) then
      exit
    end if

    c2 = euler_number2 ( n )

    write ( *, '(2x,i4,i12,g14.6)' ) n, c1, c2

  end do
 
  return
end
subroutine euler_poly_test ( )

!*****************************************************************************80
!
!! euler_poly_test() tests euler_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) euler_poly
  real ( kind = rk ) f
  integer i
  integer, parameter :: n = 15
  real ( kind = rk ) x

  x = 0.5D+00
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULER_POLY_test():'
  write ( *, '(a)' ) '  EULER_POLY evaluates Euler polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N      X             F(X)'
  write ( *, '(a)' ) ' '
   
  do i = 0, n
    f = euler_poly ( i, x )
    write ( *, '(2x,i2,2x,2g14.6)' ) i, x, f
  end do
 
  return
end
subroutine eulerian_test ( )

!*****************************************************************************80
!
!! eulerian_test() tests eulerian().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 7

  integer e(n,n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULERIAN_test():'
  write ( *, '(a)' ) '  EULERIAN evaluates Eulerian numbers.'
  write ( *, '(a)' ) ' '
 
  call eulerian ( n, e )

  do i = 1, n
    write ( *, '(2x,10i6)' )  e(i,1:n)
  end do
 
  return
end
subroutine f_hofstadter_test ( )

!*****************************************************************************80
!
!! f_hofstadter_test() tests f_hofstadter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f
  integer f_hofstadter
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F_HOFSTADTER_test():'
  write ( *, '(a)' ) '  F_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  F function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   F(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    f = f_hofstadter ( i )
    write ( *, '(2x,2i8)' ) i, f
  end do

  return
end
subroutine fibonacci_direct_test ( )

!*****************************************************************************80
!
!! fibonacci_direct_test() tests fibonacci_direct().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f
  integer i
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIBONACCI_DIRECT_test():'
  write ( *, '(a)' ) '  FIBONACCI_DIRECT evalutes a Fibonacci number directly.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I        F(I)'
  write ( *, '(a)' ) ' '

  n = 20
 
  do i = 1, n
    call fibonacci_direct ( i, f )
    write ( *, '(2x,i8,i10)' ) i, f
  end do
 
  return
end
subroutine fibonacci_floor_test ( )

!*****************************************************************************80
!
!! fibonacci_floor_test() tests fibonacci_floor().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f
  integer i
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIBONACCI_FLOOR_test():'
  write ( *, '(a)' ) '  FIBONACCI_FLOOR computes the largest Fibonacci number'
  write ( *, '(a)' ) '  less than or equal to a given positive integer.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N  Fibonacci  Index'
  write ( *, '(a)' ) ' ' 

  do n = 1, 20
    call fibonacci_floor ( n, f, i )
    write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, f, i
  end do
 
  return
end
subroutine fibonacci_recursive_test ( )

!*****************************************************************************80
!
!! fibonacci_recursive_test() tests fibonacci_recursive().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 20

  integer f(n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FIBONACCI_RECURSIVE_test():'
  write ( *, '(a)' ) '  FIBONACCI_RECURSIVE computes the Fibonacci sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       F(N)'
  write ( *, '(a)' ) ' '
 
  call fibonacci_recursive ( n, f )
 
  do i = 1, n
    write ( *, '(2x,i8,i10)' ) i, f(i)
  end do
 
  return
end
subroutine g_hofstadter_test ( )

!*****************************************************************************80
!
!! g_hofstadter_test() tests g_hofstadter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer g
  integer g_hofstadter
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'G_HOFSTADTER_test():'
  write ( *, '(a)' ) '  G_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  G function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   G(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    g = g_hofstadter ( i )
    write ( *, '(2x,2i8)' ) i, g
  end do

  return
end
subroutine gegenbauer_poly_test ( )

!*****************************************************************************80
!
!! gegenbauer_poly_test() tests gegenbauer_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ), allocatable, dimension ( : ) :: c
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEGENBAUER_POLY_TEST:'
  write ( *, '(a)' ) '  GEGENBAUER_POLY computes values of '
  write ( *, '(a)' ) '  the Gegenbauer polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        A           X       GPV      GEGENBAUER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gegenbauer_poly_values ( n_data, n, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(0:n) )

    call gegenbauer_poly ( n, a, x, c )
    fx2 = c(n)

    write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2g14.6)' ) n, a, x, fx, fx2

    deallocate ( c )

  end do

  return
end
subroutine gen_hermite_poly_test ( )

!*****************************************************************************80
!
!! gen_hermite_poly_test() tests gen_hermite_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 6
  integer, parameter :: n = 10

  real ( kind = rk ) c(0:n)
  integer j
  real ( kind = rk ) mu
  real ( kind = rk ), dimension ( test_num ) :: mu_test = (/ &
    0.0D+00, 0.0D+00, 0.1D+00, 0.1D+00, 0.5D+00, 1.0D+00 /)
  integer test
  real ( kind = rk ) x
  real ( kind = rk ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_HERMITE_POLY_test():'
  write ( *, '(a)' ) '  GEN_HERMITE_POLY evaluates the generalized Hermite '
  write ( *, '(a)' ) '  polynomials.'

  do test = 1, test_num

    x = x_test(test)
    mu = mu_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Table of H(N,MU)(X) for'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '    N(max) = ', n
    write ( *, '(a,g14.6)' ) '    MU =  ', mu
    write ( *, '(a,g14.6)' ) '    X =      ', x
    write ( *, '(a)' ) ' '
  
    call gen_hermite_poly ( n, x, mu, c )
 
    do j = 0, n
      write ( *, '(2x,i8,2x,g14.6)' ) j, c(j)
    end do

  end do
 
  return
end
subroutine gen_laguerre_poly_test ( )

!*****************************************************************************80
!
!! gen_laguerre_poly_test() tests gen_laguerre_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 February 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 6
  integer, parameter :: n = 10

  real ( kind = rk ) alpha
  real ( kind = rk ), dimension ( test_num ) :: alpha_test = (/ &
    0.0D+00, 0.0D+00, 0.1D+00, 0.1D+00, 0.5D+00, 1.0D+00 /)
  real ( kind = rk ) c(0:n)
  integer j
  integer test
  real ( kind = rk ) x
  real ( kind = rk ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEN_LAGUERRE_POLY_test():'
  write ( *, '(a)' ) '  GEN_LAGUERRE_POLY evaluates the generalized Laguerre '
  write ( *, '(a)' ) '  polynomials.'

  do test = 1, test_num

    x = x_test(test)
    alpha = alpha_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Table of L(N,ALPHA)(X) for'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '    N(max) = ', n
    write ( *, '(a,g14.6)' ) '    ALPHA =  ', alpha
    write ( *, '(a,g14.6)' ) '    X =      ', x
    write ( *, '(a)' ) ' '
  
    call gen_laguerre_poly ( n, alpha, x, c )
 
    do j = 0, n
      write ( *, '(2x,i8,2x,g14.6)' ) j, c(j)
    end do

  end do
 
  return
end
subroutine gud_test ( )

!*****************************************************************************80
!
!! gud_test() tests gud().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) gud
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GUD_TEST:'
  write ( *, '(a)' ) '  GUD evaluates the Gudermannian function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Exact F       GUD(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gud_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gud ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine h_hofstadter_test ( )

!*****************************************************************************80
!
!! h_hofstadter_test() tests h_hofstadter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer h
  integer h_hofstadter
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'H_HOFSTADTER_test():'
  write ( *, '(a)' ) '  H_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  H function.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   H(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    h = h_hofstadter ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, h
  end do

  return
end
subroutine harmonic_test ( )

!*****************************************************************************80
!
!! harmonic_test() tests harmonic().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h1
  real ( kind = rk ) h2
  real ( kind = rk ) harmonic
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'harmonic_test():'
  write ( *, '(a)' ) '  harmonic() computes Harmonic numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact H(N)    computed H(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call harmonic_values ( n_data, n, h1 )

    if ( n_data == 0 ) then
      exit
    end if

    h2 = harmonic ( n )

    write ( *, '(2x,i8,2x,g20.12,2x,g20.12)' ) n, h1, h2

  end do

  return
end
subroutine harmonic_estimate_test ( )

!*****************************************************************************80
!
!! harmonic_estimate_test() tests harmonic_estimate().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) h1
  real ( kind = rk ) h2
  real ( kind = rk ) h3
  real ( kind = rk ) harmonic
  real ( kind = rk ) harmonic_estimate
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'harmonic_estimate_test():'
  write ( *, '(a)' ) '  harmonic_estimate() estimates Harmonic numbers.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact H(N)    computed H(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call harmonic_values ( n_data, n, h2 )

    if ( n_data == 0 ) then
      exit
    end if

    h1 = harmonic_estimate ( n )

    h3 = harmonic ( n )

    write ( *, '(2x,i8,2x,g20.12,2x,g20.12,2x,g20.12)' ) n, h1, h2, h3

  end do

  return
end
subroutine hermite_poly_phys_test ( )

!*****************************************************************************80
!
!! hermite_poly_phys_test() tests hermite_poly_phys().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLY_PHYS_TEST:'
  write ( *, '(a)' ) '  HERMITE_POLY_PHYS evaluates the Hermite polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       H(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_phys_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call hermite_poly_phys ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine hermite_poly_phys_coef_test ( )

!*****************************************************************************80
!
!! hermite_poly_phys_coef_test() tests hermite_poly_phys_coef().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLY_PHYS_COEF_test():'
  write ( *, '(a)' ) &
    '  HERMITE_POLY_PHYS_COEF: Hermite polynomial coefficients.'

  call hermite_poly_phys_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  H(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  return
end
subroutine i4_choose_test ( )

!*****************************************************************************80
!
!! i4_choose_test() tests i4_choose().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer cnk
  integer i4_choose
  integer k
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_CHOOSE_test():'
  write ( *, '(a)' ) '  I4_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K    CNK'
  write ( *, '(a)' ) ' '
 
  do n = 0, 4
    do k = 0, n
      cnk = i4_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine i4_factor_test ( )

!*****************************************************************************80
!
!! i4_factor_test() tests i4_factor().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: maxfactor = 10

  integer i
  integer j
  integer n
  integer, dimension ( 3 ) :: n_test = (/ &
    60, 664048, 8466763 /)
  integer nfactor
  integer nleft
  integer factor(maxfactor)
  integer power(maxfactor)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_FACTOR_TEST:'
  write ( *, '(a)' ) '  I4_FACTOR tries to factor an I4'

  do i = 1, 3
    n = n_test(i)
    call i4_factor ( n, maxfactor, nfactor, factor, power, nleft )
    write ( *, '(a)' ) ''
    write ( *, '(a,i9)' ) '  Factors of N = ', n
    do j = 1, nfactor
      write ( *, '(i9,a,i4)' ) factor(j), '^', power(j)
    end do
    if ( nleft /= 1 ) then
      write ( *, '(a,i4)' ) '  Unresolved factor NLEFT = ', nleft
    end if
  end do

  return
end
subroutine i4_factorial_test ( )

!*****************************************************************************80
!
!! i4_factorial_test() tests i4_factorial().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fn
  integer fn2
  integer i4_factorial
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_FACTORIAL_TEST:'
  write ( *, '(a)' ) '  I4_FACTORIAL evaluates the factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X       Exact F       I4_FACTORIAL(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = i4_factorial ( n )

    write ( *, '(2x,i4,2i12)' ) n, fn, fn2

  end do

  return
end
subroutine i4_factorial2_test ( )

!*****************************************************************************80
!
!! i4_factorial2_test() tests i4_factorial2().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer fn
  integer fn2
  integer n
  integer n_data
  integer i4_factorial2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_FACTORIAL2_TEST:'
  write ( *, '(a)' ) '  I4_FACTORIAL2 evaluates the double factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   Exact  I4_FACTORIAL2(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_factorial2_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = i4_factorial2 ( n )

    write ( *, '(2x,i4,2i8)' ) n, fn, fn2

  end do

  return
end
subroutine i4_is_fibonacci_test ( )

!*****************************************************************************80
!
!! i4_is_fibonacci_test() tests i4_is_fibonacci().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 February 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer i4
  logical i4_is_fibonacci
  integer :: i4_test(10) = (/- 13, 0, 1, 8, 10, 50, 55, 100, 144, 200 /)
  logical l
  integer test_num

  test_num = 10
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'I4_IS_FIBONACCI_test():'
  write ( *, '(a)' ) '  I4_IS_FIBONACCI returns T or F depending on'
  write ( *, '(a)' ) '  whether I4 is a Fibonacci number.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I4   T/F'
  write ( *, '(a)' ) ''

  do i = 1, test_num

    i4 = i4_test(i)
    l = i4_is_fibonacci ( i4 )
    write ( *, '(2x,i4,4x,l1)' ) i4, l

  end do

  return
end
subroutine i4_is_prime_test ( )

!*****************************************************************************80
!
!! i4_is_prime_test() tests i4_is_prime().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 April 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  logical i4_is_prime

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'I4_IS_PRIME_test():'
  write ( *, '(a)' ) '  I4_IS_PRIME reports whether an integer is prime.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  I     I4_IS_PRIME(I)'
  write ( *, '(a)' ) ''

  do i = -2, 25
    write ( *, '(2x,i6,2x,l1)' ) i, i4_is_prime ( i )
  end do

  return
end
subroutine i4_is_triangular_test ( )

!*****************************************************************************80
!
!! i4_triangular_test() tests i4_is_triangular().
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

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  logical i4_is_triangular
  logical l

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_IS_TRIANGULAR_test():'
  write ( *, '(a)' ) '  I4_IS_TRIANGULAR returns T or F depending on'
  write ( *, '(a)' ) '  whether I is triangular.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I    T/F'
  write ( *, '(a)' ) ' '

  do i = 0, 20

    l = i4_is_triangular ( i )

    write ( *, '(2x,i4,4x,l1)' )  i, l

  end do
 
  return
end
subroutine i4_partition_distinct_count_test ( )

!*****************************************************************************80
!
!! i4_partition_distinct_count_test() tests i4_partition_distinct_count().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data
  integer, parameter :: n_max = 20

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_PARTITION_DISTINCT_COUNT_TEST:'
  write ( *, '(a)' ) '  For the number of partitions of an integer'
  write ( *, '(a)' ) '  into distinct parts,'
  write ( *, '(a)' ) '  I4_PARTITION_DISTINCT_COUNT'
  write ( *, '(a)' ) '  computes any value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           N       Exact F    Q(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_distinct_count_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    if ( n_max < n ) then
      cycle
    end if

    call i4_partition_distinct_count ( n, c2 )

    write ( *, '(2x,3i10)' ) n, c, c2

  end do

  return
end
subroutine i4_to_triangle_lower_test ( )

!*****************************************************************************80
!
!! i4_to_triangle_lower_test() tests i4_to_triangle_lower().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 March 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_TO_TRIANGLE_LOWER_test():'
  write ( *, '(a)' ) '  I4_TO_TRIANGLE_LOWER converts a linear index to a'
  write ( *, '(a)' ) '  lower triangular one.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I  =>   J   K'
  write ( *, '(a)' ) ' '

  do i = 0, 20

    call i4_to_triangle_lower ( i, j, k )

    write ( *, '(2x,i4,4x,i4,i4,4x,i4)' )  i, j, k

  end do
 
  return
end
subroutine jacobi_poly_test ( )

!*****************************************************************************80
!
!! jacobi_poly_test() tests jacobi_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ), allocatable :: c(:)
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JACOBI_POLY_TEST:'
  write ( *, '(a)' ) '  JACOBI_POLY computes values of '
  write ( *, '(a)' ) '  the Jacobi polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       A       B      X       JPV      JACOBI'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_poly_values ( n_data, n, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(n+1) )

    call jacobi_poly ( n, a, b, x, c )
    fx2 = c(n+1)

    write ( *, '(2x,i8,2x,f8.4,2x,f8.4,f10.4,2g14.6)' ) n, a, b, x, fx, fx2

    deallocate ( c )

  end do

  return
end
subroutine jacobi_symbol_test ( )

!*****************************************************************************80
!
!! jacobi_symbol_test() tests jacobi_symbol().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  integer l
  integer p
  integer, dimension ( test_num ) :: p_test = (/ 3, 9, 10, 12 /)
  integer q
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JACOBI_SYMBOL_test():'
  write ( *, '(a)' ) '  JACOBI_SYMBOL computes the Jacobi symbol'
  write ( *, '(a)' ) '  (Q/P), which records if Q is a quadratic '
  write ( *, '(a)' ) '  residue modulo the number P.'

  do test = 1, test_num
    p = p_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Jacobi Symbols for P = ', p
    write ( *, '(a)' ) ' '
    do q = 0, p
      call jacobi_symbol ( q, p, l )
      write ( *, '(2x,3i8)' ) p, q, l
    end do
  end do

  return
end
subroutine krawtchouk_test ( )

!*****************************************************************************80
!
!! krawtchouk_test() tests krawtchouk().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 2

  integer i
  integer j
  integer m
  real ( kind = rk ) p
  real ( kind = rk ), dimension ( test_num ) :: p_test = (/ &
    0.25D+00, 0.5D+00 /)
  integer test
  real ( kind = rk ) x
  real ( kind = rk ) value(0:n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KRAWTCHOUK_TEST:'
  write ( *, '(a)' ) '  KRAWTCHOUK evaluates Krawtchouk polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      P         X           M      K(N,P,X,M)'
  write ( *, '(a)' ) ' '

  m = 5

  do test = 1, test_num

    p = p_test(test)

    do j = 0, 5

      x = real ( j, kind = rk ) / 2.0D+00

      call krawtchouk ( n, p, x, m, value )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,i8,2x,g14.6)' ) &
          i, p, x, m, value(i)

      end do

    end do

  end do

  return
end
subroutine laguerre_associated_test ( )

!*****************************************************************************80
!
!! laguerre_associated_test() tests laguerre_associated().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 6
  integer, parameter :: n = 6

  real ( kind = rk ) c(0:n)
  integer j
  integer m
  integer, dimension ( test_num ) :: m_test = (/ 0, 0, 1, 2, 3, 1 /)
  integer test
  real ( kind = rk ) x
  real ( kind = rk ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_ASSOCIATED_test():'
  write ( *, '(a)' ) '  LAGUERRE_ASSOCIATED evaluates the associated Laguerre'
  write ( *, '(a)' ) '  polynomials.'

  do test = 1, test_num

    m = m_test(test)
    x = x_test(test)

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Table of L(N,M)(X) for'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i4)' ) '  N(max) = ', n
    write ( *, '(a,i4)' ) '  M      = ', m
    write ( *, '(a,g14.6)' ) '  X =      ', x
    write ( *, '(a)' ) ' '
 
    call laguerre_associated ( n, m, x, c )
 
    do j = 0, n
      write ( *, '(2x,i8,g14.6)' ) j, c(j)
    end do
 
  end do

  return
end
subroutine laguerre_poly_test ( )

!*****************************************************************************80
!
!! LAGUERRE_POLY_test() tests LAGUERRE_POLY.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLY_TEST:'
  write ( *, '(a)' ) '  LAGUERRE_POLY evaluates the Laguerre polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       L(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call laguerre_poly ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine laguerre_poly_coef_test ( )

!*****************************************************************************80
!
!! LAGUERRE_POLY_COEF_test() tests LAGUERRE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n,0:n)
  real ( kind = rk ) fact
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLY_COEF_test():'
  write ( *, '(a)' ) &
    '  LAGUERRE_POLY_COEF determines the Laguerre polynomial coefficients.'

  call laguerre_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  L(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do
 
  fact = 1.0D+00

  do i = 0, n

    if ( 0 < i ) then
      fact = fact * real ( i, kind = rk )
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  Factorially scaled L(', i, ')'
    write ( *, '(a)' ) ' '

    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) fact * c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) fact * c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) fact * c(i,j), ' * x**', j
      end if
    end do
    
  end do

  return
end
subroutine lambert_w_test ( )

!*****************************************************************************80
!
!! lambert_w_test() tests lambert_w().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) lambert_w
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lambert_w_test():'
  write ( *, '(a)' ) '  lambert_w() estimates the Lambert W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           X           W(X)        W(X)'
  write ( *, '(a)' ) '                   Tabulated   Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lambert_w ( x )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine lambert_w_estimate_test ( )

!*****************************************************************************80
!
!! lambert_w_estimate_test() tests lambert_w_estimate().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) lambert_w_estimate
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lambert_w_estimate_test():'
  write ( *, '(a)' ) '  lambert_w_estimate() estimates the Lambert W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           X           W(X)        W(X)'
  write ( *, '(a)' ) '                   Tabulated       Estimate'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lambert_w_estimate ( x )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine legendre_associated_test ( )

!*****************************************************************************80
!
!! legendre_associated_test() tests legendre_associated().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) fx2(0:n_max)
  real ( kind = rk ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'legendre_associated_test():'
  write ( *, '(a)' ) &
    '  legendre_associated() evaluates associated Legendre functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M    X     Exact F       PNM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_associated ( n, m, x, fx2 )

    write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

  end do

  return
end
subroutine legendre_associated_normalized_test ( )

!*****************************************************************************80
!
!! legendre_associated_normalized_test() tests legendre_associated_normalized().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) fx2(0:n_max)
  real ( kind = rk ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'legendre_associated_normalized_test():'
  write ( *, '(a)' ) &
    '  legendre_associated_normalized() evaluates associated Legendre functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M    X     Exact F       PNM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_associated_normalized ( n, m, x, fx2 )

    write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

  end do

  return
end
subroutine legendre_function_q_test ( )

!*****************************************************************************80
!
!! legendre_function_q_test() tests legendre_function_q().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'legendre_function_q_test():'
  write ( *, '(a)' ) '  legendre_function_q() evaluates the Legendre Q function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       Q(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_function_q_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_function_q ( n, x, fx2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine legendre_poly_test ( )

!*****************************************************************************80
!
!! legendre_poly_test() tests legendre_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 12

  real ( kind = rk ) fx
  real ( kind = rk ) fp2(0:n_max)
  real ( kind = rk ) fx2(0:n_max)
  integer n
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'legendre_poly_test():'
  write ( *, '(a)' ) '  legendre_poly() evaluates the Legendre PN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N    X      Exact F       P(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call legendre_poly ( n, x, fx2, fp2 )

    write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

  end do

  return
end
subroutine legendre_poly_coef_test ( )

!*****************************************************************************80
!
!! LEGENDRE_POLY_COEF_test() tests LEGENDRE_POLY_COEF.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n,0:n)
  integer i
  integer j

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLY_COEF_test():'
  write ( *, '(a)' ) &
    '  LEGENDRE_POLY_COEF returns Legendre polynomial coefficients.'

  call legendre_poly_coef ( n, c )
 
  do i = 0, n
    write ( *, '(a)' ) ' '
    write ( *, '(a,i2,a)' ) '  P(', i, ')'
    write ( *, '(a)' ) ' '
    do j = i, 0, -1
      if ( j == 0 ) then
        write ( *, '(2x,g14.6)' ) c(i,j)
      else if ( j == 1 ) then
        write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
      else
        write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
      end if
    end do
  end do

  return
end
subroutine legendre_symbol_test ( )

!*****************************************************************************80
!
!! LEGENDRE_SYMBOL_test() tests LEGENDRE_SYMBOL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  integer l
  integer p
  integer, dimension ( test_num ) :: p_test = (/ 7, 11, 13, 17 /)
  integer q
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_SYMBOL_test():'
  write ( *, '(a)' ) '  LEGENDRE_SYMBOL computes the Legendre'
  write ( *, '(a)' ) '  symbol (Q/P) which records whether Q is '
  write ( *, '(a)' ) '  a quadratic residue modulo the prime P.'

  do test = 1, test_num
    p = p_test(test)
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Legendre Symbols for P = ', p
    write ( *, '(a)' ) ' '
    do q = 0, p
      call legendre_symbol ( q, p, l )
      write ( *, '(2x,3i8)' ) p, q, l
    end do
  end do

  return
end
subroutine lerch_test ( )

!*****************************************************************************80
!
!! LERCH_test() tests LERCH.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) lerch
  integer n_data
  integer s
  real ( kind = rk ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LERCH_test():'
  write ( *, '(a)' ) '  LERCH computes the Lerch function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       Z       S       A         Lerch           Lerch'
  write ( *, '(a)' ) '                             Tabulated        Computed'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call lerch_values ( n_data, z, s, a, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = lerch ( z, s, a )

    write ( *, '(2x,f8.4,2x,i4,2x,f8.4,2x,g14.6,2x,g14.6)' ) z, s, a, fx, fx2

  end do
 
  return
end
subroutine lock_test ( )

!*****************************************************************************80
!
!! LOCK_test() tests LOCK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer a(0:n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOCK_test():'
  write ( *, '(a)' ) '  LOCK counts the combinations on a button lock.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      LOCK(I)'
  write ( *, '(a)' ) ' '

  call lock ( n, a )

  do i = 0, n
    write ( *, '(2x,i8,2x,i10)' )  i, a(i)
  end do
 
  return
end
subroutine meixner_test ( )

!*****************************************************************************80
!
!! MEIXNER_test() tests MEIXNER.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5
  integer, parameter :: test_num = 3

  real ( kind = rk ) beta
  real ( kind = rk ) :: beta_test(test_num) = (/ &
    0.5D+00, 1.0D+00, 2.0D+00 /)
  real ( kind = rk ) c
  real ( kind = rk ) :: c_test(test_num) = (/ &
    0.125D+00, 0.25D+00, 0.5D+00 /)
  integer i
  integer j
  integer test
  real ( kind = rk ) v(0:n)
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MEIXNER_TEST:'
  write ( *, '(a)' ) '  MEIXNER evaluates Meixner polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      BETA         C         X        M(N,BETA,C,X)'

  do test = 1, test_num

    beta = beta_test(test)
    c = c_test(test)

    do j = 0, 5

      x = real ( j, kind = rk ) / 2.0D+00

      call meixner ( n, beta, c, x, v )

      write ( *, '(a)' ) ' '

      do i = 0, n

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) &
          i, beta, c, x, v(i)

      end do

    end do

  end do

  return
end
subroutine mertens_test ( )

!*****************************************************************************80
!
!! MERTENS_test() tests MERTENS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 October 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer mertens
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MERTENS_test():'
  write ( *, '(a)' ) '  MERTENS computes the Mertens function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MERTENS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call mertens_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    c2 = mertens ( n )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine moebius_test ( )

!*****************************************************************************80
!
!! MOEBIUS_test() tests MOEBIUS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOEBIUS_test():'
  write ( *, '(a)' ) '  MOEBIUS computes the Moebius function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N     Exact   MOEBIUS(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call moebius_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call moebius ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine motzkin_test ( )

!*****************************************************************************80
!
!! MOTZKIN_test() tests MOTZKIN.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  integer a(0:n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOTZKIN_test():'
  write ( *, '(a)' ) '  MOTZKIN computes the Motzkin numbers A(0:N).'
  write ( *, '(a)' ) '  A(N) counts the paths from (0,0) to (N,0).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I         A(I)'
  write ( *, '(a)' ) ' '

  call motzkin ( n, a )

  do i = 0, n
    write ( *, '(2x,i8,2x,i10)' )  i, a(i)
  end do
 
  return
end
subroutine normal_01_cdf_inverse_test ( )

!*****************************************************************************80
!
!! NORMAL_01_CDF_INVERSE_test() tests NORMAL_01_CDF_INVERSE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  integer n_data
  real ( kind = rk ) x1
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_01_CDF_INVERSE_TEST:'
  write ( *, '(a)' ) '  NORMAL_01_CDF_INVERSE inverts the error function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    FX            X    NORMAL_01_CDF_INVERSE(FX)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x1, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call normal_01_cdf_inverse ( fx, x2 )

    write ( *, '(2x,f8.4,2g14.6)' ) fx, x1, x2

  end do

  return
end
subroutine omega_test ( )

!*****************************************************************************80
!
!! OMEGA_test() tests OMEGA.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OMEGA_test():'
  write ( *, '(a)' ) '  OMEGA counts the distinct prime divisors of an integer N.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '             N      Exact   OMEGA(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call omega_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call omega ( n, c2 )

    write ( *, '(2x,i12,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine pentagon_num_test ( )

!*****************************************************************************80
!
!! PENTAGON_NUM_test() tests PENTAGON_NUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer p

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PENTAGON_NUM_test():'
  write ( *, '(a)' ) '  PENTAGON_NUM computes the pentagonal numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I      Pent(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    call pentagon_num ( n, p )
    write ( *, '(2x,i8,2x,i8)' ) n, p
  end do
 
  return
end
subroutine phi_test ( )

!*****************************************************************************80
!
!! phi_test() tests phi().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PHI_test():'
  write ( *, '(a)' ) '  PHI computes the PHI function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N     Exact     PHI(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call phi_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call phi ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine pi_estimate_test ( )

!*****************************************************************************80
!
!! pi_estimate_test() tests pi_estimate().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_data
  real ( kind = rk ) pi_estimate
  real ( kind = rk ) pi1
  integer pi2
  real ( kind = rk ) r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'pi_estimate_test():'
  write ( *, '(a)' ) '  pi_estimate() estimates Pi(n).'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact Pi(N)    computed Pi(N)    Ratio'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call pi_values ( n_data, n, pi2 )

    if ( n_data == 0 ) then
      exit
    end if

    pi1 = pi_estimate ( n )

    r = pi1 / real ( pi2, kind = rk )

    write ( *, '(2x,i8,2x,g20.12,2x,i20,2x,g20.12)' ) n, pi1, pi2, r

  end do

  return
end
subroutine plane_partition_num_test ( )

!*****************************************************************************80
!
!! plane_partition_num_test() tests plane_partition_num.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer p
  integer plane_partition_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PLANE_PARTITION_NUM_test():'
  write ( *, '(a)' ) '  PLANE_PARTITION_NUM computes the number of.'
  write ( *, '(a)' ) '  plane partitions of the number N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      P(N)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    p = plane_partition_num ( n )
    write ( *, '(2x,i8,2x,i8)' ) n, p
  end do
 
  return
end
subroutine poly_bernoulli_test ( )

!*****************************************************************************80
!
!! poly_bernoulli_test() tests poly_bernoulli().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer b
  integer k
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_BERNOULLI_test():'
  write ( *, '(a)' ) '  POLY_BERNOULLI computes the poly-Bernoulli numbers'
  write ( *, '(a)' ) '  of negative index, B_n^(-k)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N     K    B_N^(-K)'
  write ( *, '(a)' ) ' '

  do k = 0, 6
    write ( *, '(a)' ) ' '
    do n = 0, 6

      call poly_bernoulli ( n, k, b )

      write ( *, '(2x,i4,2x,i4,2x,i12)' ) n, k, b

    end do
  end do

  return
end
subroutine poly_coef_count_test ( )

!*****************************************************************************80
!
!! poly_coef_count_test() tests poly_coef_count().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer degree
  integer dim
  integer poly_coef_count

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLY_COEF_COUNT_test():'
  write ( *, '(a)' ) '  POLY_COEF_COUNT counts the number of coefficients'
  write ( *, '(a)' ) '  in a polynomial of degree DEGREE and dimension DIM'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' Dimension    Degree     Count'

  do dim = 1, 10, 3
    write ( *, '(a)' ) ' '
    do degree = 0, 5
      write ( *, '(2x,i8,2x,i8,2x,i8)' ) &
       dim, degree, poly_coef_count ( dim, degree )
    end do
  end do
 
  return
end
subroutine prime_test ( )

!*****************************************************************************80
!
!! prime_test() tests prime().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer n
  integer prime
  integer prime_max

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'PRIME_test():'
  write ( *, '(a)' ) '  PRIME returns primes from a table.'

  n = -1
  prime_max = prime ( n )
  write ( *, '(a)' ) ''
  write ( *, '(a,i6)' ) '  Number of primes stored is ', prime_max
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     I    Prime(I)'
  write ( *, '(a)' ) ''
  do i = 1, 10
    write ( *, '(4x,i4,2x,i6)' ) i, prime(i)
  end do
  write ( *, '(a)' ) ''
  do i = prime_max - 10, prime_max
    write ( *, '(4x,i4,2x,i6)' ) i, prime(i)
  end do
  
  return
end

subroutine pyramid_num_test ( )

!*****************************************************************************80
!
!! pyramid_num_test() tests pyramid_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer pyramid_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_NUM_test():'
  write ( *, '(a)' ) '  PYRAMID_NUM computes the pyramidal numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    PYR(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, pyramid_num ( n )
  end do
 
  return
end
subroutine pyramid_square_num_test ( )

!*****************************************************************************80
!
!! pyramid_square_num_test() tests pyramid_square_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer pyramid_square_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PYRAMID_SQUARE_NUM_test():'
  write ( *, '(a)' ) '  PYRAMID_SQUARE_NUM computes the pyramidal square numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    PYR(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, pyramid_square_num ( n )
  end do
 
  return
end
subroutine r8_agm_test ( )

!*****************************************************************************80
!
!! r8_agm_test() tests r8_agm().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 December 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) r8_agm

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_AGM_test():'
  write ( *, '(a)' ) '  R8_AGM computes the arithmetic geometric mean.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a,a)' ) '      A           B          ', &
    '   AGM                       AGM                   Diff'
  write ( *, '(a,a)' ) '                             ', &
    '  (Tabulated)             R8_AGM(A,B)'
  write ( *, '(a)' ) ' '
     
  n_data = 0

  do

    call agm_values ( n_data, a, b, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_agm ( a, b )

    write ( *, '(2x,f10.6,2x,f10.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      a, b, fx, fx2, abs ( fx - fx2 )

  end do
     
  return
end
subroutine r8_beta_test ( )

!*****************************************************************************80
!
!! r8_beta_test() tests r8_beta().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fxy
  real ( kind = rk ) fxy2
  integer n_data
  real ( kind = rk ) r8_beta
  real ( kind = rk ) x
  real ( kind = rk ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_BETA_TEST:'
  write ( *, '(a)' ) '  R8_BETA evaluates the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Y        Exact F       R8_BETA(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    fxy2 = r8_beta ( x, y )

    write ( *, '(2x,2f8.4,2g14.6)' ) x, y, fxy, fxy2

  end do

  return
end
subroutine r8_choose_test ( )

!*****************************************************************************80
!
!! r8_choose_test() tests r8_choose().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cnk
  integer k
  integer n
  real ( kind = rk ) r8_choose

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_CHOOSE_test():'
  write ( *, '(a)' ) '  R8_CHOOSE evaluates C(N,K).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N       K      CNK'
  write ( *, '(a)' ) ' '
 
  do n = 0, 4
    do k = 0, n
      cnk = r8_choose ( n, k )
      write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) n, k, cnk
    end do
  end do
 
  return
end
subroutine r8_cube_root_test ( )

!*****************************************************************************80
!
!! r8_cube_root_test() tests r8_cube_root().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 May 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  real ( kind = rk ) r8_cube_root
  real ( kind = rk ) x1
  real ( kind = rk ) y
  real ( kind = rk ) x2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'r8_cube_root_test():'
  write ( *, '(a)' ) '  r8_cube_root() computes the cube root of an R8.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '       X               Y               Y^3'
  write ( *, '(a)' ) ''

  do i = 1, 10
    call random_number ( x1 )
    x1 = -10.0D+00 + 20.0D+00 * x1
    y = r8_cube_root ( x1 )
    x2 = y ** 3
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x1, y, x2
  end do

  return
end
subroutine r8_erf_test ( )

!*****************************************************************************80
!
!! r8_erf_test() tests r8_erf().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_erf
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_ERF_TEST:'
  write ( *, '(a)' ) '  R8_ERF evaluates the error function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X      Exact F       R8_ERF(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_erf ( x )

    write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

  end do

  return
end
subroutine r8_erf_inverse_test ( )

!*****************************************************************************80
!
!! r8_erf_inverse_test() tests r8_erf_inverse().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 August 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r8_erf_inverse
  real ( kind = rk ) fx
  integer n_data
  real ( kind = rk ) x1
  real ( kind = rk ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_ERF_INVERSE_TEST:'
  write ( *, '(a)' ) '  R8_ERF_INVERSE inverts the error function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    FX            X    R8_ERF_INVERSE(FX)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x1, fx )

    if ( n_data == 0 ) then
      exit
    end if

    x2 = r8_erf_inverse ( fx )

    write ( *, '(2x,f8.4,2g14.6)' ) fx, x1, x2

  end do

  return
end
subroutine r8_factorial_test ( )

!*****************************************************************************80
!
!! r8_factorial_test() tests r8_factorial().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fn
  real ( kind = rk ) fn2
  integer n_data
  integer n
  real ( kind = rk ) r8_factorial

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FACTORIAL_TEST:'
  write ( *, '(a)' ) '  R8_FACTORIAL evaluates the factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N       Exact F       R8_FACTORIAL(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = r8_factorial ( n )

    write ( *, '(2x,i4,2g14.6)' ) n, fn, fn2

  end do

  return
end
subroutine r8_factorial_log_test ( )

!*****************************************************************************80
!
!! r8_factorial_log_test() tests r8_factorial_log().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fn
  real ( kind = rk ) fn2
  real ( kind = rk ) r8_factorial_log
  integer n_data
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FACTORIAL_LOG_TEST:'
  write ( *, '(a)' ) '  R8_FACTORIAL_LOG evaluates the logarithm of the '
  write ( *, '(a)' ) '  factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N	   Exact F	 R8_FACTORIAL_LOG(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_log_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    fn2 = r8_factorial_log ( n )

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, fn, fn2

  end do

  return
end
subroutine r8_gamma_test ( )

!*****************************************************************************80
!
!! r8_gamma_test() tests r8_gamma().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 Nvember 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) r8_gamma
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_GAMMA_TEST:'
  write ( *, '(a)' ) '  R8_GAMMA evaluates the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Gamma(X)                   Gamma(X)  ' &
  // '               DIFF'
  write ( * , '(a)' ) '               (Tabulated)              (R8_GAMMA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_gamma ( x )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine r8_hyper_2f1_test ( )

!*****************************************************************************80
!
!! r8_hyper_2f1_test() tests r8_hyper_2f1().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    18 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) r8_hyper_2f1
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_HYPER_2F1_TEST:'
  write ( *, '(a)' ) '  R8_HYPER_2F1 evaluates the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) '      A       B       C       X      ', &
  ' 2F1                       2F1                     DIFF'
  write ( *, '(a,a)' ) '                                     ', &
  '(tabulated)               (computed)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_values ( n_data, a, b, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_hyper_2f1 ( a, b, c, x )

    write ( *, &
    '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, c, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine r8_psi_test ( )

!*****************************************************************************80
!
!! r8_psi_test() tests r8_psi().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    09 February 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n_data
  real ( kind = rk ) r8_psi
  real ( kind = rk ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_PSI_TEST:'
  write ( *, '(a)' ) '  R8_PSI evaluates the Psi function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X         Psi(X)                     Psi(X)  ' &
  // '               DIFF'
  write ( * , '(a)' ) '               (Tabulated)                (R8_PSI)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r8_psi ( x )

    write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine r8poly_degree_test ( )

!*****************************************************************************80
!
!! r8poly_degree_test() tests r8poly_degree().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) :: c1(0:3) = (/ 1.0, 2.0, 3.0, 4.0 /) 
  real ( kind = rk ) :: c2(0:3) = (/ 1.0, 2.0, 3.0, 0.0 /) 
  real ( kind = rk ) :: c3(0:3) = (/ 1.0, 2.0, 0.0, 4.0 /)
  real ( kind = rk ) :: c4(0:3) = (/ 1.0, 0.0, 0.0, 0.0 /)
  real ( kind = rk ) :: c5(0:3) = (/ 0.0, 0.0, 0.0, 0.0 /)
  integer d
  integer m
  integer r8poly_degree

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8POLY_DEGREE_test():'
  write ( *, '(a)' ) '  R8POLY_DEGREE determines the degree of an R8POLY.'

  m = 3

  call r8poly_print ( m, c1, '  The R8POLY:' )
  d = r8poly_degree ( m, c1 )
  write ( *, '(a,i2,a,i2)' ) '  Dimensioned degree = ', m, '  Actual degree = ', d

  call r8poly_print ( m, c2, '  The R8POLY:' )
  d = r8poly_degree ( m, c2 )
  write ( *, '(a,i2,a,i2)' ) '  Dimensioned degree = ', m, '  Actual degree = ', d

  call r8poly_print ( m, c3, '  The R8POLY:' )
  d = r8poly_degree ( m, c3 )
  write ( *, '(a,i2,a,i2)' ) '  Dimensioned degree = ', m, '  Actual degree = ', d

  call r8poly_print ( m, c4, '  The R8POLY:' )
  d = r8poly_degree ( m, c4 )
  write ( *, '(a,i2,a,i2)' ) '  Dimensioned degree = ', m, '  Actual degree = ', d

  call r8poly_print ( m, c5, '  The R8POLY:' )
  d = r8poly_degree ( m, c5 )
  write ( *, '(a,i2,a,i2)' ) '  Dimensioned degree = ', m, '  Actual degree = ', d

  return
end
subroutine r8poly_print_test ( )

!*****************************************************************************80
!
!! r8poly_print_test() tests r8poly_print().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  real ( kind = rk ), dimension ( 0 : m ) :: c = (/ &
    12.0D+00, -3.4D+00, 56.0D+00, 0.0D+00, 0.78D+00, 9.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8POLY_PRINT_test():'
  write ( *, '(a)' ) '  R8POLY_PRINT prints an R8poly.'

  call r8poly_print ( m, c, '  The R8POLY:' )

  return
end
subroutine r8poly_value_horner_test ( )

!*****************************************************************************80
!
!! r8poly_value_horner_test() tests r8poly_value_horner().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 16

  real ( kind = rk ), dimension (0:m) :: c = (/ &
    24.0D+00, -50.0D+00, +35.0D+00, -10.0D+00, 1.0D+00 /)
  integer i
  real ( kind = rk ) p
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) x(n)
  real ( kind = rk ) x_hi
  real ( kind = rk ) x_lo

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8POLY_VALUE_HORNER_test():'
  write ( *, '(a)' ) '  R8POLY_VALUE_HORNER evaluates a polynomial at'
  write ( *, '(a)' ) '  one point, using Horner''s method.'

  call r8poly_print ( m, c, '  The polynomial coefficients:' )

  x_lo = 0.0D+00
  x_hi = 5.0D+00
  call r8vec_linspace ( n, x_lo, x_hi, x )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I    X    P(X)'
  write ( *, '(a)' ) ''

  do i = 1, n
    p = r8poly_value_horner ( m, c, x(i) )
    write ( *, '(2x,i2,2x,f8.4,2x,g14.6)' ) i, x(i), p
  end do

  return
end
subroutine sigma_test ( )

!*****************************************************************************80
!
!! sigma_test() tests sigma().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIGMA_test():'
  write ( *, '(a)' ) '  SIGMA computes the SIGMA function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '     N     Exact   SIGMA(N)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call sigma_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call sigma ( n, c2 )

    write ( *, '(2x,i4,2i10)' ) n, c, c2

  end do
 
  return
end
subroutine simplex_num_test ( )

!*****************************************************************************80
!
!! simplex_num_test() tests simplex_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n
  integer simplex_num
  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_NUM_test():'
  write ( *, '(a)' ) '  SIMPLEX_NUM computes the N-th simplex number'
  write ( *, '(a)' ) '  in M dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      M: 0     1     2     3     4     5'
  write ( *, '(a)' ) '   N'
 
  do n = 0, 10
    write ( *, '(2x,i2)', advance = 'no' ) n
    do m = 0, 5
      value = simplex_num ( m, n )
      write ( *, '(2x,i4)', advance = 'no' ) value
    end do
    write ( *, '(a)' ) ''
  end do
 
  return
end
subroutine sin_power_int_test ( )

!*****************************************************************************80
!
!! sin_power_int_test() tests sin_power_int().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  integer n
  integer n_data
  real ( kind = rk ) sin_power_int

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_POWER_INT_TEST:'
  write ( *, '(a)' ) '  SIN_POWER_INT returns values of '
  write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      A         B          N      Exact           Computed'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = sin_power_int ( a, b, n )

    write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) a, b, n, fx, fx2

  end do

  return
end
subroutine slices_test ( )

!*****************************************************************************80
!
!! slices_test() tests slices().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 August 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_max = 5
  integer, parameter :: slice_max = 8

  integer dim_num
  integer p(dim_max,slice_max)
  integer piece_num
  integer slice_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SLICES_TEST:'
  write ( *, '(a)' ) '  SLICES determines the maximum number of pieces created'
  write ( *, '(a)' ) '  by SLICE_NUM slices in a DIM_NUM space.'

  do dim_num = 1, dim_max
    do slice_num = 1, slice_max
      call slices ( dim_num, slice_num, piece_num )
      p(dim_num,slice_num) = piece_num
    end do
  end do

  call i4mat_print ( dim_max, slice_max, p, '  Slice Array:' )

  return
end
subroutine spherical_harmonic_test ( )

!*****************************************************************************80
!
!! spherical_harmonic_test() tests spherical_harmonic().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 20

  real ( kind = rk ) c(0:n_max)
  integer l
  integer m
  integer n_data
  real ( kind = rk ) phi
  real ( kind = rk ) s(0:n_max)
  real ( kind = rk ) theta
  real ( kind = rk ) yi
  real ( kind = rk ) yi2
  real ( kind = rk ) yr
  real ( kind = rk ) yr2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERICAL_HARMONIC_TEST:'
  write ( *, '(a)' ) '  SPHERICAL_HARMONIC evaluates spherical harmonic'
  write ( *, '(a)' ) '  functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       L       M   THETA    PHI     C              S'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

    if ( n_data == 0 ) then
      exit
    end if

    call spherical_harmonic ( l, m, theta, phi, c, s )

    yr2 = c(l)
    yi2 = s(l)

    write ( *, '(2x,i8,2x,i6,2f8.4,2g14.6)' ) l, m, theta, phi, yr,  yi
    write ( *, '(2x,8x,2x,6x,16x,  2g14.6)' )                   yr2, yi2

  end do

  return
end
subroutine stirling_estimate_test ( )

!*****************************************************************************80
!
!! stirling_estimate_test() tests stirling_estimate().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    28 May 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_data
  real ( kind = rk ) stirling_estimate
  real ( kind = rk ) value1
  real ( kind = rk ) value2

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'stirling_estimate_test():'
  write ( *, '(a)' ) '  stirling_estimate() returns the Stirling estimate for'
  write ( *, '(a)' ) '  log(n!).'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '         N     Estimate    Exact'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call factorial_log_values ( n_data, n, value2 )

    if ( n_data == 0 ) then
      exit
    end if

    value1 = stirling_estimate ( n )

    write ( *, '(2x,i6,2x,g14.6,2x,g14.6)' ) n, value1, value2

  end do

  return
end
subroutine stirling1_table_test ( )

!*****************************************************************************80
!
!! stirling1_table_test() tests stirling1_table().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = m

  integer i
  integer s1(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stirling1_table_test():'
  write ( *, '(a)' ) '  stirling1_table(): Stirling numbers of first kind.'
  write ( *, '(a,i8)' ) '  Get rows 1 through ', m
  write ( *, '(a)' ) ' '
 
  call stirling1_table ( m, n, s1 )
 
  do i = 1, m
    write ( *, '(2x,i8,8i8)' ) i, s1(i,1:n)
  end do
 
  return
end
subroutine stirling2_number_test ( )

!*****************************************************************************80
!
!! stirling2_number_test() tests stirling2_number().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer k
  integer k_max
  integer n
  integer n_max
  integer stirling2_number

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stirling2_number_test():'
  write ( *, '(a)' ) '  stirling2_number() evaluates a Stirling'
  write ( *, '(a)' ) '  number of the second kind.'
  write ( *, '(a)' ) ' '
     
  n_max = 8
  k_max = 8

  do n = 0, n_max
    do k = 0, k_max
      write ( *, '(2x,i4)', advance = 'no' ) stirling2_number ( n, k )
    end do
    write ( *, '(a)' ) ''
  end do
     
  return
end
subroutine stirling2_table_test ( )

!*****************************************************************************80
!
!! stirling2_table_test() tests stirling2_table().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 8
  integer, parameter :: n = m

  integer i
  integer s2(m,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'stirling2_table_test():'
  write ( *, '(a)' ) '  stirling2_table(): Stirling numbers of second kind.'
  write ( *, '(a,i4)' ) '  Get rows 1 through ', m
  write ( *, '(a)' ) ' '
 
  call stirling2_table ( m, n, s2 )
 
  do i = 1, m
    write ( *, '(2x,i8,8i8)' ) i, s2(i,1:n)
  end do
 
  return
end
subroutine tau_test ( )

!*****************************************************************************80
!
!! tau_test() tests tau().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer c
  integer c2
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TAU_test():'
  write ( *, '(a)' ) '  TAU computes the Tau function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '         N  exact C(I)  computed C(I)'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call tau_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    call tau ( n, c2 )

    write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

  end do
 
  return
end
subroutine tetrahedron_num_test ( )

!*****************************************************************************80
!
!! tetrahedron_num_test() tests tetrahedron_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer tetrahedron_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TETRAHEDRON_NUM_test():'
  write ( *, '(a)' ) '  TETRAHEDRON_NUM computes the tetrahedron numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    TETR(I)'
  write ( *, '(a)' ) ' '

  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, tetrahedron_num ( n )
  end do
 
  return
end
subroutine triangle_num_test ( )

!*****************************************************************************80
!
!! triangle_num_test() tests triangle_num().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer triangle_num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_NUM_test():'
  write ( *, '(a)' ) '  TRIANGLE_NUM computes the triangular numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I    TRI(I)'
  write ( *, '(a)' ) ' '
 
  do n = 1, 10
    write ( *, '(2x,i8,2x,i8)' ) n, triangle_num ( n )
  end do
 
  return
end
subroutine triangle_lower_to_i4_test ( )

!*****************************************************************************80
!
!! triangle_lower_to_i4_test() tests triangle_lower_to_i4().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    22 March 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIANGLE_LOWER_TO_I4_test():'
  write ( *, '(a)' ) '  TRIANGLE_LOWER_TO_I4 converts a lower triangular index'
  write ( *, '(a)' ) '  to a linear one.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I, J ==> K'
  write ( *, '(a)' ) ' '

  do i = 1, 4
    do j = 1, i
      call triangle_lower_to_i4 ( i, j, k )
      write ( *, '(2x,i4,i4,4x,i4)' ) i, j, k
    end do
  end do
 
  return
end
subroutine tribonacci_direct_test ( )

!*****************************************************************************80
!
!! tribonacci_direct_test() tests tribonacci_direct().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 May 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer n
  integer t
  integer tribonacci_direct

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tribonacci_direct_test():'
  write ( *, '(a)' ) '  tribonacci_direct() computes the Tribonacci sequence.'
  write ( *, '(a)' ) ''
  
  n = 20
  do i = 1, n
    t = tribonacci_direct ( i )
    write ( *, '(2x,i4,2x,i8)' ) i, t
  end do

  return
end
subroutine tribonacci_recursive_test ( )

!*****************************************************************************80
!
!! tribonacci_recursive_test() tests tribonacci_recursive().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 May 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 22

  integer f(n)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'tribonacci_recursive_test():'
  write ( *, '(a)' ) '  tribonacci_recursive() computes the Tribonacci sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       F(N)'
  write ( *, '(a)' ) ' '
 
  call tribonacci_recursive ( n, f )
 
  do i = 1, n
    write ( *, '(2x,i8,i10)' ) i - 2, f(i)
  end do
 
  return
end
subroutine tribonacci_roots_test ( )

!*****************************************************************************80
!
!! tribonacci_roots_test() tests tribonacci_roots().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 May 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) alpha
  complex ( kind = rk ) beta
  complex ( kind = rk ) gamma
  real ( kind = rk ) p
  complex ( kind = rk ) pc

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'tribonacci_roots_test():'
  write ( *, '(a)' ) '  tribonacci_roots() computes the Tribonacci roots.'
  write ( *, '(a)' ) ''

  call tribonacci_roots ( alpha, beta, gamma )

  p = alpha**3 - alpha**2 - alpha - 1.0D+00
  write ( *, '(a,g14.6,a,g14.6)' ) &
    '  alpha = ', alpha, ', p(alpha) = ', p

  pc = beta**3 - beta**2 - beta - 1.0D+00
  write ( *, '(a,2g14.6,a,2g14.6)' ) &
    '  beta = ', beta, ', p(beta) = ', pc
 
  pc = gamma**3 - gamma**2 - gamma - 1.0D+00
  write ( *, '(a,2g14.6,a,2g14.6)' ) &
    '  gamma = ', gamma, ', p(gamma) = ', pc

  return
end
subroutine trinomial_test ( )

!*****************************************************************************80
!
!! trinomial_test() tests trinomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 April 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k
  integer t
  integer trinomial

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRINOMIAL_test():'
  write ( *, '(a)' ) '  TRINOMIAL evaluates the trinomial coefficient:'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  T(I,J,K) = (I+J+K)! / I! / J! / K!'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     I     J     K    T(I,J,K)'
  write ( *, '(a)' ) ' '
 
  do k = 0, 4
    do j = 0, 4
      do i = 0, 4
        t = trinomial ( i, j, k )
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i8)' ) i, j, k, t
      end do
    end do
  end do
 
  return
end
subroutine v_hofstadter_test ( )

!*****************************************************************************80
!
!! v_hofstadter_test() tests v_hofstadter().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer v
  integer v_hofstadter

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'V_HOFSTADTER_test():'
  write ( *, '(a)' ) '  V_HOFSTADTER evaluates Hofstadter''s recursive'
  write ( *, '(a)' ) '  V function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N   V(N)'
  write ( *, '(a)' ) ' '

  do i = 0, 30
    v = v_hofstadter ( i )
    write ( *, '(2x,i8,2x,i8)' ) i, v
  end do

  return
end
subroutine vibonacci_test ( )

!*****************************************************************************80
!
!! vibonacci_test() tests vibonacci().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 20
  integer, parameter :: n_time = 3

  integer i
  integer j
  integer v(n,n_time)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VIBONACCI_test():'
  write ( *, '(a)' ) '  VIBONACCI computes a Vibonacci sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of times we compute the series: ', n_time
  write ( *, '(a)' ) ' '

  do j = 1, n_time
    call vibonacci ( n, v(1,j) ) 
  end do

  do i = 1, n
    write ( *, '(2x,i8,2x,3i8)' ) i, v(i,1:n_time)
  end do
 
  return
end
subroutine zeckendorf_test ( )

!*****************************************************************************80
!
!! zeckendorf_test() tests zeckendorf().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 October 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m_max = 20

  integer i_list(m_max)
  integer f_list(m_max)
  integer m
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZECKENDORF_test():'
  write ( *, '(a)' ) '  ZECKENDORF computes the Zeckendorf decomposition of'
  write ( *, '(a)' ) '  an integer N into nonconsecutive Fibonacci numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N Sum M Parts'
  write ( *, '(a)' ) ' '

  do n = 1, 100

    call zeckendorf ( n, m_max, m, i_list, f_list )

    write ( *, '(2x,i8,2x,15i4)' ) n, f_list(1:m)

  end do

  return
end
subroutine zernike_poly_test ( )

!*****************************************************************************80
!
!! zernike_poly_test() tests zernike_poly().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 5

  real ( kind = rk ) c(0:n_max)
  integer m
  integer n
  real ( kind = rk ) r8poly_value_horner
  real ( kind = rk ) rho
  real ( kind = rk ) z1
  real ( kind = rk ) z2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZERNIKE_POLY_test():'
  write ( *, '(a)' ) '  ZERNIKE_POLY evaluates a Zernike polynomial directly.'

  rho = 0.987654321D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Z1: Compute polynomial coefficients,'
  write ( *, '(a)' ) '  then evaluate by Horner''s method;'
  write ( *, '(a)' ) '  Z2: Evaluate directly by recursion.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   N   M       Z1              Z2'
  write ( *, '(a)' ) ' '

  do n = 0, 5

    write ( *, '(a)' ) ' '

    do m = 0, n

      call zernike_poly_coef ( m, n, c )
      z1 = r8poly_value_horner ( n, c, rho )

      call zernike_poly ( m, n, rho, z2 )

      write ( *, '(2x,i2,2x,i2,2x,g16.8,2x,g16.8)' ) n, m, z1, z2

    end do

  end do

  return
end
subroutine zernike_poly_coef_test ( )

!*****************************************************************************80
!
!! zernike_poly_coef_test() tests zernike_poly_coef().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    11 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) c(0:n)
  integer m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZERNIKE_POLY_COEF_test():'
  write ( *, '(a)' ) '  ZERNIKE_POLY_COEF determines the Zernike'
  write ( *, '(a)' ) '  polynomial coefficients.'

  do m = 0, n

    call zernike_poly_coef ( m, n, c )
 
    call r8poly_print ( n, c, '  Zernike polynomial' )

  end do

  return
end
subroutine zeta_m1_test ( )

!*****************************************************************************80
!
!! zeta_m1_test() tests zeta_m1().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk ) p
  real ( kind = rk ) tol
  real ( kind = rk ) z1
  real ( kind = rk ) z2
  real ( kind = rk ) zeta_m1

  tol = 1.0D-10
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZETA_M1_test():'
  write ( *, '(a)' ) '  ZETA_M1 computes the Zeta Minus One function.'
  write ( *, '(a,g14.6)' ) '  Requested tolerance = ', tol
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact Zeta    computed Zeta'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call zeta_m1_values ( n_data, p, z1 )

    if ( n_data == 0 ) then
      exit
    end if

    z2 = zeta_m1 ( p, tol )

    write ( *, '(2x,f8.2,2x,g20.12,2x,g20.12)' ) p, z1, z2

  end do

  return
end
subroutine zeta_naive_test ( )

!*****************************************************************************80
!
!! zeta_naive_test() tests zeta_naive().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 June 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer n_data
  real ( kind = rk ) n_real
  real ( kind = rk ) z1
  real ( kind = rk ) z2
  real ( kind = rk ) zeta_naive

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZETA_NAIVE_test():'
  write ( *, '(a)' ) '  ZETA_NAIVE computes the Zeta function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N    exact Zeta    computed Zeta'
  write ( *, '(a)' ) ' '
 
  n_data = 0

  do

    call zeta_values ( n_data, n, z1 )

    if ( n_data == 0 ) then
      exit
    end if

    n_real = real ( n, kind = rk )

    z2 = zeta_naive ( n_real )

    write ( *, '(2x,i8,2x,g20.12,2x,g20.12)' ) n, z1, z2

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
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

