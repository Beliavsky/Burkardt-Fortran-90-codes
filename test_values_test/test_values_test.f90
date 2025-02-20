program main

!*****************************************************************************80
!
!! test_values_test() tests test_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_values_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test test_values().'

  call abram0_values_test ( )
  call abram1_values_test ( )
  call abram2_values_test ( )
  call agm_values_test ( )
  call airy_ai_values_test ( )
  call airy_ai_int_values_test ( )
  call airy_ai_prime_values_test ( )
  call airy_bi_values_test ( )
  call airy_bi_int_values_test ( )
  call airy_bi_prime_values_test ( )
  call airy_cai_values_test ( )
  call airy_cbi_values_test ( )
  call airy_gi_values_test ( )
  call airy_hi_values_test ( )
  call arccos_values_test ( )
  call arccosh_values_test ( )
  call arcsin_values_test ( )
  call arcsinh_values_test ( )
  call arctan_values_test ( )
  call arctan_int_values_test ( )
  call arctan2_values_test ( )
  call arctanh_values_test ( )

  call bei0_values_test ( )
  call bei1_values_test ( )
  call bell_values_test ( )
  call ber0_values_test ( )
  call ber1_values_test ( )
  call bernoulli_number_values_test ( )
  call bernoulli_poly_values_test ( )
  call bernstein_poly_01_values_test ( )
  call bessel_i0_values_test ( )
  call bessel_i0_int_values_test ( )
  call bessel_i0_spherical_values_test ( )
  call bessel_i1_values_test ( )
  call bessel_i1_spherical_values_test ( )
  call bessel_in_values_test ( )
  call bessel_ix_values_test ( )
  call bessel_j_spherical_values_test ( )
  call bessel_j0_values_test ( )
  call bessel_j0_int_values_test ( )
  call bessel_j0_spherical_values_test ( )
  call bessel_j0_zero_values_test ( )
  call bessel_j1_values_test ( )
  call bessel_j1_spherical_values_test ( )
  call bessel_jn_values_test ( )
  call bessel_jx_values_test ( )
  call bessel_k0_values_test ( )
  call bessel_k0_int_values_test ( )
  call bessel_k1_values_test ( )
  call bessel_kn_values_test ( )
  call bessel_kx_values_test ( )
  call bessel_y0_values_test ( )
  call bessel_y0_int_values_test ( )
  call bessel_y0_spherical_values_test ( )
  call bessel_y0_zero_values_test ( )
  call bessel_y1_values_test ( )
  call bessel_y1_spherical_values_test ( )
  call bessel_yn_values_test ( )
  call bessel_yx_values_test ( )
  call beta_values_test ( )
  call beta_cdf_values_test ( )
  call beta_inc_values_test ( )
  call beta_log_values_test ( )
  call beta_noncentral_cdf_values_test ( )
  call beta_pdf_values_test ( )
  call binomial_values_test ( )
  call binomial_cdf_values_test ( )
  call binomial_pdf_values_test ( )
  call bivariate_normal_cdf_values_test ( )

  call c8_log_values_test ( )
  call catalan_values_test ( )
  call cauchy_cdf_values_test ( )
  call cbrt_values_test ( )
  call cheby_t_poly_values_test ( )
  call cheby_t01_poly_values_test ( )
  call cheby_u_poly_values_test ( )
  call cheby_u01_poly_values_test ( )
  call cheby_v_poly_values_test ( )
  call cheby_v01_poly_values_test ( )
  call cheby_w_poly_values_test ( )
  call cheby_w01_poly_values_test ( )
  call chi_values_test ( )
  call chi_square_cdf_values_test ( )
  call chi_square_pdf_values_test ( )
  call chi_square_noncentral_cdf_values_test ( )
  call ci_values_test ( )
  call cin_values_test ( )
  call cinh_values_test ( )
  call clausen_values_test ( )
  call clebsch_gordan_values_test ( )
  call collatz_count_values_test ( )
  call cos_values_test ( )
  call cos_degree_values_test ( )
  call cos_power_int_values_test ( )
  call cosh_values_test ( )
  call cot_values_test ( )
  call cp_values_test ( )

  call datenum_values_test ( )
  call dawson_values_test ( )
  call debye1_values_test ( )
  call debye2_values_test ( )
  call debye3_values_test ( )
  call debye4_values_test ( )
  call dedekind_sum_values_test ( )
  call dielectric_values_test ( )
  call dilogarithm_values_test ( )
  call dixon_elliptic_values_test ( )

  call e1_values_test ( )
  call ei_values_test ( )
  call easter_gregorian_values_test ( )
  call easter_julian_values_test ( )
  call elliptic_ea_values_test ( )
  call elliptic_ek_values_test ( )
  call elliptic_em_values_test ( )
  call elliptic_fa_values_test ( )
  call elliptic_fk_values_test ( )
  call elliptic_fm_values_test ( )
  call elliptic_inc_ea_values_test ( )
  call elliptic_inc_ek_values_test ( )
  call elliptic_inc_em_values_test ( )
  call elliptic_inc_fa_values_test ( )
  call elliptic_inc_fk_values_test ( )
  call elliptic_inc_fm_values_test ( )
  call elliptic_inc_pia_values_test ( )
  call elliptic_inc_pik_values_test ( )
  call elliptic_inc_pim_values_test ( )
  call elliptic_pia_values_test ( )
  call elliptic_pik_values_test ( )
  call elliptic_pim_values_test ( )
  call erf_values_test ( )
  call erfc_values_test ( )
  call euler_number_values_test ( )
  call euler_poly_values_test ( )
  call exp_values_test ( )
  call exp3_int_values_test ( )
  call exponential_01_pdf_values_test ( )
  call exponential_cdf_values_test ( )
  call exponential_pdf_values_test ( )
  call extreme_values_cdf_values_test ( )

  call f_cdf_values_test ( )
  call f_noncentral_cdf_values_test ( )
  call factorial_values_test ( )
  call factorial2_values_test ( )
  call fresnel_cos_values_test ( )
  call fresnel_sin_values_test ( )
  call frobenius_number_data_values_test ( )
  call frobenius_number_order_values_test ( )
  call frobenius_number_order2_values_test ( )

  call gamma_values_test ( )
  call gamma_01_pdf_values_test ( )
  call gamma_cdf_values_test ( )
  call gamma_inc_values_test ( )
  call gamma_inc_p_values_test ( )
  call gamma_inc_q_values_test ( )
  call gamma_inc_tricomi_values_test ( )
  call gamma_log_values_test ( )
  call gamma_pdf_values_test ( )
  call gcd_values_test ( )
  call gegenbauer_poly_values_test ( )
  call geometric_cdf_values_test ( )
  call goodwin_values_test ( )
  call gud_values_test ( )

  call harmonic_values_test ( )
  call hermite_function_values_test ( )
  call hermite_poly_phys_values_test ( )
  call hermite_poly_prob_values_test ( )
  call hyper_1f1_values_test ( )
  call hyper_2f1_values_test ( )
  call hyper_2f1_complex_values_test ( )
  call hypergeometric_cdf_values_test ( )
  call hypergeometric_pdf_values_test ( )
  call hypergeometric_u_values_test ( )

  call i0ml0_values_test ( )
  call i1ml1_values_test ( )
  call i4_fall_values_test ( )
  call i4_gpf_values_test ( )
  call i4_rise_values_test ( )
  call int_values_test ( )
  call inverse_chi_square_pdf_values_test ( )
  call inverse_gamma_pdf_values_test ( )
  call is_prime_values_test ( )

  call jacobi_cn_values_test ( )
  call jacobi_dn_values_test ( )
  call jacobi_poly_values_test ( )
  call jacobi_sn_values_test ( )
  call jed_ce_values_test ( )
  call jed_mjd_values_test ( )
  call jed_rd_values_test ( )
  call jed_weekday_values_test ( )

  call kei0_values_test ( )
  call kei1_values_test ( )
  call ker0_values_test ( )
  call ker1_values_test ( )
  call knapsack_values_test ( )

  call laguerre_associated_values_test ( )
  call laguerre_general_values_test ( )
  call laguerre_polynomial_values_test ( )
  call lambert_w_values_test ( )
  call laplace_cdf_values_test ( )
  call legendre_associated_values_test ( )
  call legendre_associated_normalized_values_test ( )
  call legendre_associated_normalized_sphere_values_test ( )
  call legendre_function_q_values_test ( )
  call legendre_normalized_polynomial_values_test ( )
  call legendre_polynomial_values_test ( )
  call legendre_shifted_polynomial_values_test ( )
  call lerch_values_test ( )
  call lobachevsky_values_test ( )
  call lobatto_polynomial_values_test ( )
  call lobatto_polynomial_derivative_values_test ( )
  call log_values_test ( )
  call log_normal_cdf_values_test ( )
  call log_series_cdf_values_test ( )
  call log10_values_test ( )
  call logarithmic_integral_values_test ( )
  call logistic_cdf_values_test ( )

  call mathieu_even_values_test ( )
  call mathieu_odd_values_test ( )
  call mcnugget_number_values_test ( )
  call mersenne_prime_values_test ( )
  call mertens_values_test ( )
  call mittag_leffler_ea_values_test ( )
  call mittag_leffler_eab_values_test ( )
  call moebius_values_test ( )
  call multinomial_pdf_values_test ( )

  call negative_binomial_cdf_values_test ( )
  call nine_j_values_test ( )
  call normal_01_cdf_values_test ( )
  call normal_01_pdf_values_test ( )
  call normal_cdf_values_test ( )
  call normal_pdf_values_test ( )

  call omega_values_test ( )
  call owen_values_test ( )

  call partition_count_values_test ( )
  call partition_distinct_count_values_test ( )
  call phi_values_test ( )
  call pi_values_test ( )
  call poisson_cdf_values_test ( )
  call polylogarithm_values_test ( )
  call polyomino_chiral_count_values_test ( )
  call polyomino_fixed_count_values_test ( )
  call polyomino_free_count_values_test ( )
  call prandtl_values_test ( )
  call prime_values_test ( )
  call psat_values_test ( )
  call psi_values_test ( )

  call r8_factorial_values_test ( )
  call r8_factorial_log_values_test ( )
  call r8_factorial2_values_test ( )
  call r8_fall_values_test ( )
  call r8_rise_values_test ( )
  call rayleigh_cdf_values_test ( )

  call scaled_inverse_chi_square_pdf_values_test ( )
  call secvir_values_test ( )
  call shi_values_test ( )
  call si_values_test ( )
  call sigma_values_test ( )
  call sin_values_test ( )
  call sin_degree_values_test ( )
  call sin_power_int_values_test ( )
  call sinh_values_test ( )
  call six_j_values_test ( )
  call sound_values_test ( )
  call sphere_unit_area_values_test ( )
  call sphere_unit_volume_values_test ( )
  call spherical_harmonic_values_test ( )
  call sqrt_values_test ( )
  call stirling1_values_test ( )
  call stirling2_values_test ( )
  call stromgen_values_test ( )
  call struve_h0_values_test ( )
  call struve_h1_values_test ( )
  call struve_l0_values_test ( )
  call struve_l1_values_test ( )
  call student_cdf_values_test ( )
  call student_noncentral_cdf_values_test ( )
  call subfactorial_values_test ( )
  call surten_values_test ( )
  call synch1_values_test ( )
  call synch2_values_test ( )

  call tan_values_test ( )
  call tanh_values_test ( )
  call tau_values_test ( )
  call thercon_values_test ( )
  call three_j_values_test ( )
  call tran02_values_test ( )
  call tran03_values_test ( )
  call tran04_values_test ( )
  call tran05_values_test ( )
  call tran06_values_test ( )
  call tran07_values_test ( )
  call tran08_values_test ( )
  call tran09_values_test ( )
  call trigamma_values_test ( )
  call truncated_normal_ab_cdf_values_test ( )
  call truncated_normal_ab_pdf_values_test ( )
  call truncated_normal_a_cdf_values_test ( )
  call truncated_normal_a_pdf_values_test ( )
  call truncated_normal_b_cdf_values_test ( )
  call truncated_normal_b_pdf_values_test ( )
  call tsat_values_test ( )

  call van_der_corput_values_test ( )
  call viscosity_values_test ( )
  call von_mises_cdf_values_test ( )

  call weekday_values_test ( )
  call weibull_cdf_values_test ( )
  call wright_omega_values_test ( )
  call zeta_values_test ( )
  call zeta_m1_values_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_values_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine abram0_values_test ( )

!*****************************************************************************80
!
!! abram0_values_test() tests abram0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'abram0_values_test():'
  write ( *, '(a)' ) '  abram0_values() returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine abram1_values_test ( )

!*****************************************************************************80
!
!! abram1_values_test() tests abram1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ABRAM1_VALUES_TEST:'
  write ( *, '(a)' ) '  ABRAM1_VALUES returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine abram2_values_test ( )

!*****************************************************************************80
!
!! abram2_values_test() tests abram2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ABRAM2_VALUES_TEST:'
  write ( *, '(a)' ) '  ABRAM2_VALUES returns values of '
  write ( *, '(a)' ) '  the Abramowitz function of order 2'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Abram2'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call abram2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine agm_values_test ( )

!*****************************************************************************80
!
!! agm_values_test() tests agm_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AGM_VALUES_TEST:'
  write ( *, '(a)' ) '  AGM_VALUES returns values of '
  write ( *, '(a)' ) '  the arithmetic geometric mean function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          A               B            AGM(A,B)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call agm_values ( n_data, a, b, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) a, b, fx

  end do

  return
end
subroutine airy_ai_values_test ( )

!*****************************************************************************80
!
!! airy_ai_values_test() tests airy_ai_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) ai
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_AI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_AI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Ai'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_values ( n_data, x, ai )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, ai

  end do

  return
end
subroutine airy_ai_int_values_test ( )

!*****************************************************************************80
!
!! airy_ai_int_values_test() tests airy_ai_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_AI_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_AI_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Airy function Ai(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Ai_Int'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine airy_ai_prime_values_test ( )

!*****************************************************************************80
!
!! airy_ai_prime_values_test() tests airy_ai_prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) aip
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_AI_PRIME_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_AI_PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the derivative of the Airy functions A''(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           AiP'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_ai_prime_values ( n_data, x, aip )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, aip

  end do

  return
end
subroutine airy_bi_values_test ( )

!*****************************************************************************80
!
!! airy_bi_values_test() tests airy_bi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bi
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_BI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_BI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Bi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Bi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_values ( n_data, x, bi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, bi

  end do

  return
end
subroutine airy_bi_int_values_test ( )

!*****************************************************************************80
!
!! airy_bi_int_values_test() tests airy_bi_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_BI_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_BI_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Airy function Bi(X)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Bi_Int'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine airy_bi_prime_values_test ( )

!*****************************************************************************80
!
!! airy_bi_prime_values_test() tests airy_bi_prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) bip
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_BI_PRIME_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_BI_PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function derivative B''(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BiP'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_bi_prime_values ( n_data, x, bip )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, bip

  end do

  return
end
subroutine airy_cai_values_test ( )

!*****************************************************************************80
!
!! airy_cai_values_test() tests airy_cai_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) ai
  integer n_data
  complex ( kind = ck ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_CAI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_CAI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Ai(X) with complex argument'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X                         Ai'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_cai_values ( n_data, x, ai )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,(g14.6,g14.6),2x,(g24.16,g24.16))' ) x, ai

  end do

  return
end
subroutine airy_cbi_values_test ( )

!*****************************************************************************80
!
!! airy_cbi_values_test() tests airy_cbi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 April 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) bi
  integer n_data
  complex ( kind = ck ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_CBI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_CBI_VALUES returns values of '
  write ( *, '(a)' ) '  the Airy function Bi(X) with complex argument'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X                         Bi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_cbi_values ( n_data, x, bi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,(g14.6,g14.6),2x,(g24.16,g24.16))' ) x, bi

  end do

  return
end
subroutine airy_gi_values_test ( )

!*****************************************************************************80
!
!! airy_gi_values_test() tests airy_gi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_GI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_GI_VALUES returns values of '
  write ( *, '(a)' ) '  the modified Airy function Gi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Gi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_gi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine airy_hi_values_test ( )

!*****************************************************************************80
!
!! airy_hi_values_test() tests airy_hi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'AIRY_HI_VALUES_TEST:'
  write ( *, '(a)' ) '  AIRY_HI_VALUES returns values of '
  write ( *, '(a)' ) '  the modified Airy function Hi(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Hi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call airy_hi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arccos_values_test ( )

!*****************************************************************************80
!
!! arccos_values_test() tests arccos_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCCOS_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCCOS_VALUES returns values of '
  write ( *, '(a)' ) '  the arc cosine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arccos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arccosh_values_test ( )

!*****************************************************************************80
!
!! arccosh_values_test() tests arccosh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCCOSH_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCCOSH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc cosine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arccosh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arcsin_values_test ( )

!*****************************************************************************80
!
!! arcsin_values_test() tests arcsin_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCSIN_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCSIN_VALUES returns values of '
  write ( *, '(a)' ) '  the arc sine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arcsin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arcsinh_values_test ( )

!*****************************************************************************80
!
!! arcsinh_values_test() tests arcsinh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCSINH_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCSINH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc sine'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arcsinh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arctan_values_test ( )

!*****************************************************************************80
!
!! arctan_values_test() tests arctan_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCTAN_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCTAN_VALUES returns values of '
  write ( *, '(a)' ) '  the arc tangent'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arctan_int_values_test ( )

!*****************************************************************************80
!
!! arctan_int_values_test() tests arctan_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCTAN_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCTAN_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the arctangent integral'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine arctan2_values_test ( )

!*****************************************************************************80
!
!! arctan2_values_test() tests arctan2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    05 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f
  integer n_data
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCTAN2_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCTAN2_VALUES returns values of '
  write ( *, '(a)' ) '  the arc tangent with two arguments.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y          ARCTAN2(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctan2_values ( n_data, x, y, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) x, y, f

  end do

  return
end
subroutine arctanh_values_test ( )

!*****************************************************************************80
!
!! arctanh_values_test() tests arctanh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ARCTANH_VALUES_TEST:'
  write ( *, '(a)' ) '  ARCTANH_VALUES returns values of '
  write ( *, '(a)' ) '  the hyperbolic arc tangent'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call arctanh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bei0_values_test ( )

!*****************************************************************************80
!
!! bei0_values_test() tests bei0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BEI0_VALUES_TEST:'
  write ( *, '(a)' ) '  BEI0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BEI of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bei0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bei1_values_test ( )

!*****************************************************************************80
!
!! bei1_values_test() tests bei1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BEI1_VALUES_TEST:'
  write ( *, '(a)' ) '  BEI1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BEI of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bei1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bell_values_test ( )

!*****************************************************************************80
!
!! bell_values_test() tests bell_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer c
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BELL_VALUES_TEST:'
  write ( *, '(a)' ) '  BELL_VALUES returns values of '
  write ( *, '(a)' ) '  the Bell numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        BELL(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bell_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine ber0_values_test ( )

!*****************************************************************************80
!
!! ber0_values_test() tests ber0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BER0_VALUES_TEST:'
  write ( *, '(a)' ) '  BER0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BER of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ber0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine ber1_values_test ( )

!*****************************************************************************80
!
!! ber1_values_test() tests ber1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BER1_VALUES_TEST:'
  write ( *, '(a)' ) '  BER1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function BER of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           BER1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ber1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bernoulli_number_values_test ( )

!*****************************************************************************80
!
!! bernoulli_number_values_test() tests bernoulli_number_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) c
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BERNOULLI_NUMBER_VALUES_TEST:'
  write ( *, '(a)' ) '  BERNOULLI_NUMBER_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernoulli numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        BERNOULLI_NUMBER(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernoulli_number_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,6x,g24.16)' ) n, c

  end do

  return
end
subroutine bernoulli_poly_values_test ( )

!*****************************************************************************80
!
!! bernoulli_poly_values_test() tests bernoulli_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BERNOULLI_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  BERNOULLI_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernoulli polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            BERNOULLI_POLY(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernoulli_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bernstein_poly_01_values_test ( )

!*****************************************************************************80
!
!! bernstein_poly_01_values_test() tests bernstein_poly_01_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) b
  integer k
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BERNSTEIN_POLY_01_VALUES_TEST:'
  write ( *, '(a)' ) '  BERNSTEIN_POLY_01_VALUES returns values of '
  write ( *, '(a)' ) '  the Bernstein polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         K          X           BERNSTEIN(N,K)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bernstein_poly_01_values ( n_data, n, k, x, b )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) n, k, x, b

  end do

  return
end
subroutine bessel_i0_values_test ( )

!*****************************************************************************80
!
!! bessel_i0_values_test() tests bessel_i0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_I0_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_I0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel I0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_i0_int_values_test ( )

!*****************************************************************************80
!
!! bessel_i0_int_values_test() tests bessel_i0_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_I0_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_I0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel I0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_i0_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_i0_spherical_values_test() tests bessel_i0_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_I0_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_I0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel i0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           i0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_i1_values_test ( )

!*****************************************************************************80
!
!! bessel_i1_values_test() tests bessel_i1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_I1_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_I1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel I1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            I1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_i1_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_i1_spherical_values_test() tests bessel_i1_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_I1_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_I1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel i1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           i1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_i1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_in_values_test ( )

!*****************************************************************************80
!
!! bessel_in_values_test() tests bessel_in_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_IN_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_IN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel In function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X         I(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_in_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_ix_values_test ( )

!*****************************************************************************80
!
!! bessel_ix_values_test() tests bessel_ix_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_IX_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_IX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel In function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         I(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_ix_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_j_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_j_spherical_values_test() tests bessel_j_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel j function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N          X           jn(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j_spherical_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_j0_values_test ( )

!*****************************************************************************80
!
!! bessel_j0_values_test() tests bessel_j0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J0_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel J0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_j0_int_values_test ( )

!*****************************************************************************80
!
!! bessel_j0_int_values_test() tests bessel_j0_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J0_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel J0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_j0_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_j0_spherical_values_test() tests bessel_j0_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J0_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel j0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           j0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_j0_zero_values_test ( )

!*****************************************************************************80
!
!! bessel_j0_zero_values_test() tests bessel_j0_zero_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 January 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer k
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bessel_j0_zero_values_test():'
  write ( *, '(a)' ) '  bessel_j0_zero_values() returns zeros of '
  write ( *, '(a)' ) '  the Bessel J0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          K           X(K)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j0_zero_values ( n_data, k, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i6,2x,g24.16)' ) k, fx

  end do

  return
end
subroutine bessel_j1_values_test ( )

!*****************************************************************************80
!
!! bessel_j1_values_test() tests bessel_j1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J1_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel J1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           J1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_j1_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_j1_spherical_values_test() tests bessel_j1_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_J1_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_J1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel j1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           j1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_j1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_jn_values_test ( )

!*****************************************************************************80
!
!! bessel_jn_values_test() tests bessel_jn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_JN_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_JN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Jn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           J(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_jx_values_test ( )

!*****************************************************************************80
!
!! bessel_jx_values_test() tests bessel_jx_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_JX_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_JX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Jn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         J(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_jx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_k0_values_test ( )

!*****************************************************************************80
!
!! bessel_k0_values_test() tests bessel_k0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_K0_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_K0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel K0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_k0_int_values_test ( )

!*****************************************************************************80
!
!! bessel_k0_int_values_test() tests bessel_k0_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_K0_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_K0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel K0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_k1_values_test ( )

!*****************************************************************************80
!
!! bessel_k1_values_test() tests bessel_k1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_K1_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_K1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel K1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            K1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_k1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_kn_values_test ( )

!*****************************************************************************80
!
!! bessel_kn_values_test() tests bessel_kn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_KN_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_KN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Kn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            K(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_kn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_kx_values_test ( )

!*****************************************************************************80
!
!! bessel_kx_values_test() tests bessel_kx_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_KX_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_KX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Kn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         K(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_kx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_y0_values_test ( )

!*****************************************************************************80
!
!! bessel_y0_values_test() tests bessel_y0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_Y0_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_Y0_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            X            Y0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_y0_int_values_test ( )

!*****************************************************************************80
!
!! bessel_y0_int_values_test() tests bessel_y0_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_Y0_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_Y0_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of the Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Y0_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_y0_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_y0_spherical_values_test() tests bessel_y0_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_Y0_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_Y0_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            y0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_y0_zero_values_test ( )

!*****************************************************************************80
!
!! bessel_y0_zero_values_test() tests bessel_y0_zero_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 March 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer k
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bessel_y0_zero_values_test():'
  write ( *, '(a)' ) '  bessel_y0_zero_values() returns zeros of '
  write ( *, '(a)' ) '  the Bessel Y0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          K           X(K)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y0_zero_values ( n_data, k, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i6,2x,g24.16)' ) k, fx

  end do

  return
end
subroutine bessel_y1_values_test ( )

!*****************************************************************************80
!
!! bessel_y1_values_test() tests bessel_y1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_Y1_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_Y1_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Y1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            Y1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_y1_spherical_values_test ( )

!*****************************************************************************80
!
!! bessel_y1_spherical_values_test() tests bessel_y1_spherical_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_Y1_SPHERICAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_Y1_SPHERICAL_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical Bessel Y1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            y1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_y1_spherical_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine bessel_yn_values_test ( )

!*****************************************************************************80
!
!! bessel_yn_values_test() tests bessel_yn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_YN_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_YN_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Yn function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X              Y(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yn_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine bessel_yx_values_test ( )

!*****************************************************************************80
!
!! bessel_yx_values_test() tests bessel_yx_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    31 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BESSEL_YX_VALUES_TEST:'
  write ( *, '(a)' ) '  BESSEL_YX_VALUES returns values of '
  write ( *, '(a)' ) '  the Bessel Yn function with REAL argument N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '            N             X         Y(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bessel_yx_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,g24.16)' ) n, x, fx

  end do

  return
end
subroutine beta_values_test ( )

!*****************************************************************************80
!
!! beta_values_test() tests beta_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fxy
  integer n_data
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y            BETA(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) x, y, fxy

  end do

  return
end
subroutine beta_cdf_values_test ( )

!*****************************************************************************80
!
!! beta_cdf_values_test() tests beta_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine beta_inc_values_test ( )

!*****************************************************************************80
!
!! beta_inc_values_test() tests beta_inc_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_INC_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_INC_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           BETA_INC(A,B)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_inc_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine beta_log_values_test ( )

!*****************************************************************************80
!
!! beta_log_values_test() tests beta_log_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fxy
  integer n_data
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_LOG_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the Beta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X               Y           BETA_LOG(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_log_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) x, y, fxy

  end do

  return
end
subroutine beta_noncentral_cdf_values_test ( )

!*****************************************************************************80
!
!! beta_noncentral_cdf_values_test() tests beta_noncentral_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_NONCENTRAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '        A         B     LAMBDA          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.2,2x,f8.2,2x,f8.2,2x,g14.6,2x,g24.16)' ) &
    a, b, lambda, x, fx

  end do

  return
end
subroutine beta_pdf_values_test ( )

!*****************************************************************************80
!
!! beta_pdf_values_test() tests beta_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BETA_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Beta PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          ALPHA           BETA            X           PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_pdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) alpha, beta, x, fx

  end do

  return
end
subroutine binomial_values_test ( )

!*****************************************************************************80
!
!! binomial_values_test() tests binomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  integer b
  integer c
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINOMIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  BINOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the binomial numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A         B        C(A,B)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_values ( n_data, a, b, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) a, b, c

  end do

  return
end
subroutine binomial_cdf_values_test ( )

!*****************************************************************************80
!
!! binomial_cdf_values_test() tests binomial_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINOMIAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BINOMIAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Binomial Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A        B            X       CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,i8,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine binomial_pdf_values_test ( )

!*****************************************************************************80
!
!! binomial_pdf_values_test() tests binomial_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BINOMIAL_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BINOMIAL_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Binomial Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A                B                   X       CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call binomial_pdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16,2x,i8,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine bivariate_normal_cdf_values_test ( )

!*****************************************************************************80
!
!! bivariate_normal_cdf_values_test() tests bivariate_normal_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 May 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fxy
  integer n_data
  real ( kind = rk8 ) r
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BIVARIATE_NORMAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  BIVARIATE_NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the bivariate normal CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          X               Y               R           F(R)(X,Y)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call bivariate_normal_cdf_values ( n_data, x, y, r, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) x, y, r, fxy

  end do

  return
end
subroutine c8_log_values_test ( )

!*****************************************************************************80
!
!! c8_log_values_test() tests c8_log_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2019
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) fx
  integer n_data
  complex ( kind = ck ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'C8_LOG_VALUES_TEST:'
  write ( *, '(a)' ) '  C8_LOG_VALUES returns values of the natural'
  write ( *, '(a)' ) '  logarithm of a complex value.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '                X                               LN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call c8_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,a,g14.6,a,g14.6,a,2x,a,g14.6,a,g14.6,a)' ) &
      '(', real ( x ), ',', imag ( x ), ')', &
      '(', real ( fx ), ',', imag ( fx ), ')'

  end do

  return
end
subroutine catalan_values_test ( )

!*****************************************************************************80
!
!! catalan_values_test() tests catalan_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer c
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CATALAN_VALUES_TEST:'
  write ( *, '(a)' ) '  CATALAN_VALUES returns values of '
  write ( *, '(a)' ) '  the Catalan numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          C(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call catalan_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine cauchy_cdf_values_test ( )

!*****************************************************************************80
!
!! cauchy_cdf_values_test() tests cauchy_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CAUCHY_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  CAUCHY_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Cauchy Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       MU              SIGMA           X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cauchy_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine cbrt_values_test ( )

!*****************************************************************************80
!
!! cbrt_values_test() tests cbrt_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CBRT_VALUES_TEST:'
  write ( *, '(a)' ) '  CBRT_VALUES returns values of the cube root.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X         CBRT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cbrt_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cheby_t_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_t_poly_values_test() tests cheby_t_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_T_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_T_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev T polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         T(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_t_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_t01_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_t01_poly_values_test() tests cheby_t01_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_T01_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_T01_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the shifted Chebyshev T polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         T01(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_t01_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_u_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_u_poly_values_test() tests cheby_u_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_U_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_U_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev U polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X          U(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_u_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_u01_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_u01_poly_values_test() tests cheby_u01_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_U01_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_U01_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the shifted Chebyshev U polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X          U01(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_u01_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_v_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_v_poly_values_test() tests cheby_v_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_V_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_V_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev V polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         V(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_v_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_v01_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_v01_poly_values_test() tests cheby_v01_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_V01_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_V01_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the shifted Chebyshev V polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         V01(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_v01_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_w_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_w_poly_values_test() tests cheby_w_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_W_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_W_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Chebyshev W polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         W(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_w_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine cheby_w01_poly_values_test ( )

!*****************************************************************************80
!
!! cheby_w01_poly_values_test() tests cheby_w01_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_W01_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  CHEBY_W01_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the shifted Chebyshev W polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         W01(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cheby_w01_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine chi_values_test ( )

!*****************************************************************************80
!
!! chi_values_test() tests chi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_VALUES_TEST:'
  write ( *, '(a)' ) '  CHI_VALUES returns values of '
  write ( *, '(a)' ) '  the Hyperbolic Cosine Integral function CHI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           CHI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine chi_square_cdf_values_test ( )

!*****************************************************************************80
!
!! chi_square_cdf_values_test() tests chi_square_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_SQUARE_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  CHI_SQUARE_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Chi-Squared Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        X         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine chi_square_pdf_values_test ( )

!*****************************************************************************80
!
!! chi_square_pdf_values_test() tests chi_square_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) df
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_SQUARE_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  CHI_SQUARE_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Chi-Squared Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         DF         X         PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_pdf_values ( n_data, df, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,g24.16)' ) df, x, fx

  end do

  return
end
subroutine chi_square_noncentral_cdf_values_test ( )

!*****************************************************************************80
!
!! chi_square_noncentral_cdf_values_test() tests chi_square_noncentral_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer df
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHI_SQUARE_NONCENTRAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the noncentral Chi-Squared Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        DF        LAMBDA        X         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call chi_square_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i10,2x,g10.4,2x,g10.4,2x,g24.16)' ) df, lambda, x, fx

  end do

  return
end
subroutine ci_values_test ( )

!*****************************************************************************80
!
!! ci_values_test() tests ci_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CI_VALUES_TEST:'
  write ( *, '(a)' ) '  CI_VALUES returns values of '
  write ( *, '(a)' ) '  the Cosine Integral function CI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           CI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ci_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cin_values_test ( )

!*****************************************************************************80
!
!! cin_values_test() tests cin_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CIN_VALUES_TEST:'
  write ( *, '(a)' ) '  CIN_VALUES returns values of '
  write ( *, '(a)' ) '  the Cosine Integral function CIN(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X            CIN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cinh_values_test ( )

!*****************************************************************************80
!
!! cinh_values_test() tests cinh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CINH_VALUES_TEST:'
  write ( *, '(a)' ) '  CINH_VALUES returns values of '
  write ( *, '(a)' ) '  the Hyperbolic Cosine Integral function CINH(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X            CINH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cinh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine clausen_values_test ( )

!*****************************************************************************80
!
!! clausen_values_test() tests clausen_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLAUSEN_VALUES_TEST:'
  write ( *, '(a)' ) '  CLAUSEN_VALUES returns values of '
  write ( *, '(a)' ) '  Clausen''s integral.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            FX'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call clausen_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine clebsch_gordan_values_test ( )

!*****************************************************************************80
!
!! clebsch_gordan_values_test() tests clebsch_gordan_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) j1
  real ( kind = rk8 ) j2
  real ( kind = rk8 ) j3
  real ( kind = rk8 ) m1
  real ( kind = rk8 ) m2
  real ( kind = rk8 ) m3
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CLEBSCH_GORDAN_VALUES_TEST:'
  write ( *, '(a)' ) '  CLEBSCH_GORDAN_VALUES returns values of '
  write ( *, '(a)' ) '  the Clebsch-Gordan coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      M1      M2      M3      CG'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call clebsch_gordan_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, m1, m2, m3, fx

  end do

  return
end
subroutine collatz_count_values_test ( )

!*****************************************************************************80
!
!! collatz_count_values_test() tests collatz_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 March 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer count
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COLLATZ_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  COLLATZ_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the length of the Collatz sequence that'
  write ( *, '(a)' ) '  starts at N.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      COLLATZ_COUNT(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call collatz_count_values ( n_data, n, count )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, count

  end do

  return
end
subroutine cos_values_test ( )

!*****************************************************************************80
!
!! cos_values_test() tests cos_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_VALUES_TEST:'
  write ( *, '(a)' ) '  COS_VALUES returns values of the cosine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COS(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cos_degree_values_test ( )

!*****************************************************************************80
!
!! cos_degree_values_test() tests cos_degree_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_DEGREE_VALUES_TEST:'
  write ( *, '(a)' ) '  COS_DEGREE_VALUES returns values of the cosine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COS_DEGREE(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_degree_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cos_power_int_values_test ( )

!*****************************************************************************80
!
!! cos_power_int_values_test() tests cos_power_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COS_POWER_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  COS_POWER_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of COS(X)^N from A to B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B              N    COS_POWER_INT(A,B,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cos_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,i8,2x,g14.6)' ) a, b, n, fx

  end do

  return
end
subroutine cosh_values_test ( )

!*****************************************************************************80
!
!! cosh_values_test() tests cosh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COSH_VALUES_TEST:'
  write ( *, '(a)' ) &
    '  COSH_VALUES returns values of the hyperbolic cosine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COSH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cosh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cot_values_test ( )

!*****************************************************************************80
!
!! cot_values_test() tests cot_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'COT_VALUES_TEST:'
  write ( *, '(a)' ) '  COT_VALUES returns values of the cotangent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           COT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cot_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine cp_values_test ( )

!*****************************************************************************80
!
!! cp_values_test() tests cp_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) cp
  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cp_values_test():'
  write ( *, '(a)' ) '  cp_values() returns values of '
  write ( *, '(a)' ) '  the specific heat CP '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T               P          CP(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call cp_values ( n_data, tc, p, cp )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, cp

  end do

  return
end
subroutine datenum_values_test ( )

!*****************************************************************************80
!
!! datenum_values_test() tests datenum_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 December 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  real ( kind = rk8 ) date_num
  integer m
  integer n_data
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'datenum_values_test():'
  write ( *, '(a)' ) '  datenum_values() returns values of '
  write ( *, '(a)' ) '  the MATLAB DATENUM for a given Y/M/D date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Y     M     D     DATENUM'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call datenum_values ( n_data, y, m, d, date_num )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,i4,2x,f6.0)' ) y, m, d, date_num

  end do

  return
end
subroutine dawson_values_test ( )

!*****************************************************************************80
!
!! dawson_values_test() tests dawson_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dawson_values_test():'
  write ( *, '(a)' ) '  dawson_values() returns values of '
  write ( *, '(a)' ) '  Dawson''s Integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            DAWSON(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dawson_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine debye1_values_test ( )

!*****************************************************************************80
!
!! debye1_values_test() tests debye1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'debye1_values_test():'
  write ( *, '(a)' ) '  debye1_values() returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine debye2_values_test ( )

!*****************************************************************************80
!
!! debye2_values_test() tests debye2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'debye2_values_test():'
  write ( *, '(a)' ) '  debye2_values() returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine debye3_values_test ( )

!*****************************************************************************80
!
!! debye3_values_test() tests debye3_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'debye3_values_test():'
  write ( *, '(a)' ) '  debye3_values() returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 3.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE3(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye3_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine debye4_values_test ( )

!*****************************************************************************80
!
!! debye4_values_test() tests debye4_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'debye4_values_test():'
  write ( *, '(a)' ) '  debye4_values() returns values of '
  write ( *, '(a)' ) '  Debye''s function of order 4.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           DEBYE4(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call debye4_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine dedekind_sum_values_test ( )

!*****************************************************************************80
!
!! dedekind_sum_values_test() tests dedekind_sum_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 July 2009
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  integer n
  integer n_data
  integer p
  integer q

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dedekind_sum_values_test():'
  write ( *, '(a)' ) '  dedekind_sum_values() returns values of the'
  write ( *, '(a)' ) '  Dedekind sum function: (N/D) = Dedekind_Sum ( P, Q ).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     P     Q     N     D'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dedekind_sum_values ( n_data, p, q, n, d )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) p, q, n, d

  end do

  return
end
subroutine dielectric_values_test ( )

!*****************************************************************************80
!
!! dielectric_values_test() tests dielectric_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) eps
  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dielectric_values_test():'
  write ( *, '(a)' ) '  dielectric_values() returns values of '
  write ( *, '(a)' ) '  the dielectric function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T            P            EPS(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dielectric_values ( n_data, tc, p, eps )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, eps

  end do

  return
end
subroutine dilogarithm_values_test ( )

!*****************************************************************************80
!
!! dilogarithm_values_test() tests dilogarithm_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dilogarithm_values_test():'
  write ( *, '(a)' ) '  dilogarithm_values() returns values of'
  write ( *, '(a)' ) '  the dilogarithm function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          DILOGARITHM(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dilogarithm_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine dixon_elliptic_values_test ( )

!*****************************************************************************80
!
!! dixon_elliptic_values_test() tests dixon_elliptic_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) c
  integer n_data
  real ( kind = rk8 ) s
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'dixon_elliptic_test():'
  write ( *, '(a)' ) '  dixon_elliptic() returns values of the'
  write ( *, '(a)' ) '  Dixon elliptic functions cm(x) and sm(x).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      x     cm(x)     sm(x)  cm(x)^3+sm(x)^3'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call dixon_elliptic_values ( n_data, x, c, s )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) x, c, s, c**3+s**3

  end do

  return
end
subroutine e1_values_test ( )

!*****************************************************************************80
!
!! e1_values_test() tests e1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'e1_values_test():'
  write ( *, '(a)' ) '  e1_values() returns values of'
  write ( *, '(a)' ) '  the exponential integral function E1(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          E1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call e1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine ei_values_test ( )

!*****************************************************************************80
!
!! ei_values_test() tests ei_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ei_values_test():'
  write ( *, '(a)' ) '  ei_values() returns values of'
  write ( *, '(a)' ) '  the exponential integral function EI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          EI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ei_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine easter_gregorian_values_test ( )

!*****************************************************************************80
!
!! easter_gregorian_values_test() tests easter_gregorian_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  integer m
  integer n_data
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'easter_gregorian_values_test():'
  write ( *, '(a)' ) '  easter_gregorian_values() returns values of'
  write ( *, '(a)' ) '  the occurrence of Easter in the Gregorian calendar.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   D   M     Y'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call easter_gregorian_values ( n_data, d, m, y )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i2,2x,i2,2x,i4)' ) d, m, y

  end do

  return
end
subroutine easter_julian_values_test ( )

!*****************************************************************************80
!
!! easter_julian_values_test() tests easter_julian_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  integer m
  integer n_data
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'easter_julian_values_test():'
  write ( *, '(a)' ) '  easter_julian_values() returns values of'
  write ( *, '(a)' ) '  the occurrence of Easter in the Julian calendar.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   D   M     Y'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call easter_julian_values ( n_data, d, m, y )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i2,2x,i2,2x,i4)' ) d, m, y

  end do

  return
end
subroutine elliptic_ea_values_test ( )

!*****************************************************************************80
!
!! elliptic_ea_values_test() tests elliptic_ea_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_ea_values_test():'
  write ( *, '(a)' ) '  elliptic_ea_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameter angle ALPHA in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          ALPHA        E(ALPHA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_ea_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_ek_values_test ( )

!*****************************************************************************80
!
!! elliptic_ek_values_test() tests elliptic_ek_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_em_values_test():'
  write ( *, '(a)' ) '  elliptic_em_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameter K.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          K            E(K)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_ek_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_em_values_test ( )

!*****************************************************************************80
!
!! elliptic_em_values_test() tests elliptic_em_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_em_values_test():'
  write ( *, '(a)' ) '  elliptic_em_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameter modulus M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          M            E(M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_em_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_fa_values_test ( )

!*****************************************************************************80
!
!! elliptic_fa_values_test() tests elliptic_fa_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_fa_values_test():'
  write ( *, '(a)' ) '  elliptic_fa_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameter angle ALPHA in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          ALPHA        F(ALPHA)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_fa_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_fk_values_test ( )

!*****************************************************************************80
!
!! elliptic_fk_values_test() tests elliptic_fk_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_fk_values_test():'
  write ( *, '(a)' ) '  elliptic_fk_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameter K.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          K            F(K)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_fk_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_fm_values_test ( )

!*****************************************************************************80
!
!! elliptic_fm_values_test() tests elliptic_fm_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'elliptic_fm_values_test():'
  write ( *, '(a)' ) '  elliptic_fm_values() returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameter modulus M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          M            F(M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_fm_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine elliptic_inc_ea_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_ea_values_test() tests elliptic_inc_ea_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) ea
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_EA_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_EA_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameters PHI, A.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        A            E(PHI,A)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_ea_values ( n_data, phi, a, ea )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, a, ea

  end do

  return
end
subroutine elliptic_inc_ek_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_ek_values_test() tests elliptic_inc_ek_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) ek
  real ( kind = rk8 ) k
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_EK_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_EK_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameters PHI, K.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        K            E(PHI,K)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_ek_values ( n_data, phi, k, ek )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, k, ek

  end do

  return
end
subroutine elliptic_inc_em_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_em_values_test() tests elliptic_inc_em_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) em
  real ( kind = rk8 ) m
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_EM_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_EM_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the second'
  write ( *, '(a)' ) '  kind, with parameters PHI, M.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        M            E(PHI,M)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_em_values ( n_data, phi, m, em )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, m, em

  end do

  return
end
subroutine elliptic_inc_fa_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_fa_values_test() tests elliptic_inc_fa_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fa
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_FA_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_FA_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameters PHI, A.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        A            F(PHI,A)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_fa_values ( n_data, phi, a, fa )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, a, fa

  end do

  return
end
subroutine elliptic_inc_fk_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_fk_values_test() tests elliptic_inc_fk_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fk
  real ( kind = rk8 ) k
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_FK_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_FK_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameters PHI, K.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        K            F(PHI,K)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_fk_values ( n_data, phi, k, fk )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, k, fk

  end do

  return
end
subroutine elliptic_inc_fm_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_fm_values_test() tests elliptic_inc_fm_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fm
  real ( kind = rk8 ) m
  integer n_data
  real ( kind = rk8 ) phi

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_FM_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_FM_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the first'
  write ( *, '(a)' ) '  kind, with parameters PHI, M.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI        M            F(PHI,M)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_fm_values ( n_data, phi, m, fm )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, m, fm

  end do

  return
end
subroutine elliptic_inc_pia_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_pia_values_test() tests elliptic_inc_pia_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) phi
  real ( kind = rk8 ) pia

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_PIA_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_PIA_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameters PHI, N and A.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI           N             A            Pi(PHI,N,A)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_pia_values ( n_data, phi, n, a, pia )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, n, a, pia

  end do

  return
end
subroutine elliptic_inc_pik_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_pik_values_test() tests elliptic_inc_pik_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) k
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) phi
  real ( kind = rk8 ) pik

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_PIK_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_PIK_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameters PHI, N and K.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI           N             K            Pi(PHI,N,K)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_pik_values ( n_data, phi, n, k, pik )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, n, k, pik

  end do

  return
end
subroutine elliptic_inc_pim_values_test ( )

!*****************************************************************************80
!
!! elliptic_inc_pim_values_test() tests elliptic_inc_pim_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) m
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) phi
  real ( kind = rk8 ) pim

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ELLIPTIC_INC_PIM_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_INC_PIM_VALUES stores values of'
  write ( *, '(a)' ) '  the incomplete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameters PHI, N and M.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    PHI           N             M            Pi(PHI,N,M)'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call elliptic_inc_pim_values ( n_data, phi, n, m, pim )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.6,2x,f12.6,2x,f12.6,2x,f24.16)' ) phi, n, m, pim

  end do

  return
end
subroutine elliptic_pia_values_test ( )

!*****************************************************************************80
!
!! elliptic_pia_values_test() tests elliptic_pia_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) pia

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPTIC_PIA_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_PIA_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameter angle A in degrees.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          N         A        Pi(N,A)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_pia_values ( n_data, n, a, pia )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) n, a, pia

  end do

  return
end
subroutine elliptic_pik_values_test ( )

!*****************************************************************************80
!
!! elliptic_pik_values_test() tests elliptic_pik_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) k
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) pik

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPTIC_PIK_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_PIK_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameter K.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          N         K        Pi(N,K)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_pik_values ( n_data, n, k, pik )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) n, k, pik

  end do

  return
end
subroutine elliptic_pim_values_test ( )

!*****************************************************************************80
!
!! elliptic_pim_values_test() tests elliptic_pim_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 May 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) m
  real ( kind = rk8 ) n
  integer n_data
  real ( kind = rk8 ) pim

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ELLIPTIC_PIM_VALUES_TEST:'
  write ( *, '(a)' ) '  ELLIPTIC_PIM_VALUES returns values of '
  write ( *, '(a)' ) '  the complete elliptic integral of the third'
  write ( *, '(a)' ) '  kind, with parameter M.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          N         M        Pi(N,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call elliptic_pim_values ( n_data, n, m, pim )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) n, m, pim

  end do

  return
end
subroutine erf_values_test ( )

!*****************************************************************************80
!
!! erf_values_test() tests erf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERF_VALUES_TEST:'
  write ( *, '(a)' ) '  ERF_VALUES returns values of'
  write ( *, '(a)' ) '  the error function ERF(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           ERF(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine erfc_values_test ( )

!*****************************************************************************80
!
!! erfc_values_test() tests erfc_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ERFC_VALUES_TEST:'
  write ( *, '(a)' ) '  ERFC_VALUES returns values of'
  write ( *, '(a)' ) '  the complementary error function ERFC(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           ERFC(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erfc_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine euler_number_values_test ( )

!*****************************************************************************80
!
!! euler_number_values_test() tests euler_number_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer c
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'euler_number_values_test()::'
  write ( *, '(a)' ) '  euler_number_values() returns values of '
  write ( *, '(a)' ) '  the Euler numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N        EULER_NUMBER(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_number_values ( n_data, n, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, c

  end do

  return
end
subroutine euler_poly_values_test ( )

!*****************************************************************************80
!
!! euler_poly_values_test() tests euler_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EULER_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  EULER_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Euler polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X             EULER_POLY(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call euler_poly_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine exp_values_test ( )

!*****************************************************************************80
!
!! exp_values_test() tests exp_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXP_VALUES_TEST:'
  write ( *, '(a)' ) '  EXP_VALUES returns values of the exponential function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           EXP(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine exp3_int_values_test ( )

!*****************************************************************************80
!
!! exp3_int_values_test() tests exp3_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXP3_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  EXP3_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the EXP3 Integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X          EXP3_INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exp3_int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine exponential_01_pdf_values_test ( )

!*****************************************************************************80
!
!! exponential_01_pdf_values_test() tests exponential_01_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENTIAL_01_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  EXPONENTIAL_01_PDF_VALUES returns values of'
  write ( *, '(a)' ) '  the standard Exponential PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exponential_01_pdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine exponential_cdf_values_test ( )

!*****************************************************************************80
!
!! exponential_cdf_values_test() tests exponential_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENTIAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  EXPONENTIAL_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the Exponential CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          LAMBDA          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exponential_cdf_values ( n_data, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) lambda, x, fx

  end do

  return
end
subroutine exponential_pdf_values_test ( )

!*****************************************************************************80
!
!! exponential_pdf_values_test() tests exponential_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXPONENTIAL_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  EXPONENTIAL_PDF_VALUES returns values of'
  write ( *, '(a)' ) '  the Exponential PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          LAMBDA          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call exponential_pdf_values ( n_data, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) lambda, x, fx

  end do

  return
end
subroutine extreme_values_cdf_values_test ( )

!*****************************************************************************80
!
!! extreme_values_cdf_values_test() tests extreme_values_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EXTREME_VALUES_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  EXTREME_VALUES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Extreme Values Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       ALPHA          BETA             X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call extreme_values_cdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, beta, x, fx

  end do

  return
end
subroutine f_cdf_values_test ( )

!*****************************************************************************80
!
!! f_cdf_values_test() tests f_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  integer b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  F_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the F cumulative density function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A         B          X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine f_noncentral_cdf_values_test ( )

!*****************************************************************************80
!
!! f_concentral_cdf_values_test() tests f_noncentral_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer a
  integer b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'F_NONCENTRAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  F_NONCENTRAL_CDF_VALUES returns values of'
  write ( *, '(a)' ) '  the F cumulative density function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A         B      LAMBDA       X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call f_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,f10.6,2x,g14.6,2x,g14.6)' ) a, b, lambda, x, fx

  end do

  return
end
subroutine factorial_values_test ( )

!*****************************************************************************80
!
!! factorial_values_test() tests factorial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'factorial_values_test():'
  write ( *, '(a)' ) '  factorial_values() returns values of '
  write ( *, '(a)' ) '  the factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N            N!'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine factorial2_values_test ( )

!*****************************************************************************80
!
!! factorial2_values_test() tests factorial2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'factorial2_values_test():'
  write ( *, '(a)' ) '  factorial2_values() returns values of '
  write ( *, '(a)' ) '  the double factorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N           N!!'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call factorial2_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine fresnel_cos_values_test ( )

!*****************************************************************************80
!
!! fresnel_cos_values_test() tests fresnel_cos_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRESNEL_COS_VALUES_TEST:'
  write ( *, '(a)' ) '  FRESNEL_COS_VALUES returns values of'
  write ( *, '(a)' ) '  the Fresnel cosine integral C(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           C(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call fresnel_cos_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine fresnel_sin_values_test ( )

!*****************************************************************************80
!
!! fresnel_sin_values_test() tests fresnel_sin_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FRESNEL_SIN_VALUES_TEST:'
  write ( *, '(a)' ) '  FRESNEL_SIN_VALUES returns values of'
  write ( *, '(a)' ) '  the Fresnel sine integral S(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           S(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call fresnel_sin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine frobenius_number_data_values_test ( )

!*****************************************************************************80
!
!! frobenius_number_data_values_test() tests frobenius_number_data_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, allocatable, dimension ( : ) :: c
  integer f
  integer n_data
  integer order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FROBENIUS_NUMBER_DATA_VALUES_TEST:'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_DATA_VALUES returns the corresponding'
  write ( *, '(a)' ) '  coin denominations for a Frobenius problem.'

  n_data = 0

  do

    call frobenius_number_order_values ( n_data, order )

    if ( n_data == 0 ) then
      exit
    end if

    allocate ( c(1:order) )

    call frobenius_number_data_values ( n_data, order, c, f )

    write ( *, '(a)'        ) ' '
    write ( *, '(a,i8)'    ) '  Order = ', order
    write ( *, '(10(2x,i6))' ) c(1:order)
    write ( *, '(a,i12)'   ) '  Frobenius number = ', f

    deallocate ( c )

  end do

  return
end
subroutine frobenius_number_order_values_test ( )

!*****************************************************************************80
!
!! frobenius_number_order_values_test() tests frobenius_number_order_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  integer order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FROBENIUS_NUMBER_ORDER_VALUES_TEST:'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_ORDER_VALUES returns the order for'
  write ( *, '(a)' ) '  a Frobenius problem;'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   Problem   Order'
  write ( *, '(a)' ) ''
  n_data = 0

  do

    call frobenius_number_order_values ( n_data, order )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4)'    ) n_data, order

  end do

  return
end
subroutine frobenius_number_order2_values_test ( )

!*****************************************************************************80
!
!! frobenius_number_order2_values_test() tests frobenius_number_order2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer c1
  integer c2
  integer f
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FROBENIUS_NUMBER_ORDER2_VALUES_TEST:'
  write ( *, '(a)' ) '  FROBENIUS_NUMBER_ORDER2_VALUES returns values of '
  write ( *, '(a)' ) '  the Frobenius number of order 2.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        C1        C2     F(C1,C2)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call frobenius_number_order2_values ( n_data, c1, c2, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) c1, c2, f

  end do

  return
end
subroutine gamma_values_test ( )

!*****************************************************************************80
!
!! gamma_values_test() tests gamma_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_VALUES returns values of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GAMMA(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine gamma_01_pdf_values_test ( )

!*****************************************************************************80
!
!! gamma_01_pdf_values_test() tests gamma_01_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_01_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_01_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the standard Gamma Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      ALPHA    X                  PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_01_pdf_values ( n_data, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, x, fx

  end do

  return
end
subroutine gamma_cdf_values_test ( )

!*****************************************************************************80
!
!! gamma_cdf_values_test() tests gamma_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the GAMMA Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       MU             SIGMA            X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine gamma_inc_values_test ( )

!*****************************************************************************80
!
!! gamma_inc_values_test() tests gamma_inc_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_INC_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               X            GAMMA_INC(A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine gamma_inc_p_values_test ( )

!*****************************************************************************80
!
!! gamma_inc_p_values_test() tests gamma_inc_p_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_P_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_INC_P_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               X            GAMMA_INC_P(A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_p_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine gamma_inc_q_values_test ( )

!*****************************************************************************80
!
!! gamma_inc_q_values_test() tests gamma_inc_q_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_Q_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_INC_Q_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               X            GAMMA_INC_Q(A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_q_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine gamma_inc_tricomi_values_test ( )

!*****************************************************************************80
!
!! gamma_inc_tricomi_values_test() tests gamma_inc_tricomi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_INC_TRICOMI_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_INC_TRICOMI_VALUES returns values of '
  write ( *, '(a)' ) '  the incomplete Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               X            GAMMA_INC_TRICOMI(A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_tricomi_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine gamma_log_values_test ( )

!*****************************************************************************80
!
!! gamma_log_values_test() tests gamma_log_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_LOG_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the Gamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GAMMA_LOG(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine gamma_pdf_values_test ( )

!*****************************************************************************80
!
!! gamma_pdf_values_test() tests gamma_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GAMMA_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  GAMMA_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  a Gamma Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       BETA            ALPHA          X               PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_pdf_values ( n_data, beta, alpha, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) beta, alpha, x, fx

  end do

  return
end
subroutine gcd_values_test ( )

!*****************************************************************************80
!
!! gcd_values_test() tests gcd_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer f
  integer m
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'gcd_values_test():'
  write ( *, '(a)' ) '  gcd_values() returns values of the greatest common divisor function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         m         n      gcd(m,n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gcd_values ( n_data, m, n, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8)' ) m, n, f

  end do

  return
end
subroutine gegenbauer_poly_values_test ( )

!*****************************************************************************80
!
!! gegenbauer_poly_values_test() tests gegenbauer_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEGENBAUER_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  GEGENBAUER_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Gegenbauer polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       A            X     G(N,A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gegenbauer_poly_values ( n_data, n, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g10.4,2x,g14.6)' ) n, a, x, fx

  end do

  return
end
subroutine geometric_cdf_values_test ( )

!*****************************************************************************80
!
!! geometric_cdf_values_test() tests geometric_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) cdf
  integer n_data
  real ( kind = rk8 ) p
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GEOMETRIC_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  GEOMETRIC_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Geometric Probability Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X        P         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call geometric_cdf_values ( n_data, x, p, cdf )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) x, p, cdf

  end do

  return
end
subroutine goodwin_values_test ( )

!*****************************************************************************80
!
!! goodwin_values_test() tests goodwin_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GOODWIN_VALUES_TEST:'
  write ( *, '(a)' ) '  GOODWIN_VALUES returns values of'
  write ( *, '(a)' ) '  the Goodwin and Staton function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call goodwin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine gud_values_test ( )

!*****************************************************************************80
!
!! gud_values_test() tests gud_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'GUD_VALUES_TEST:'
  write ( *, '(a)' ) '  GUD_VALUES returns values of '
  write ( *, '(a)' ) '  the Gudermannian function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            GUD(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gud_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine harmonic_values_test ( )

!*****************************************************************************80
!
!! harmonic_values_test() tests harmonic_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) h
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'harmonic_values_test():'
  write ( *, '(a)' ) '  harmonic_values() returns values of the Harmonic number sequence.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          H(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call harmonic_values ( n_data, n, h )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, h

  end do

  return
end
subroutine hermite_function_values_test ( )

!*****************************************************************************80
!
!! hermite_function_values_test() tests hermite_function_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_FUNCTION_VALUES_TEST:'
  write ( *, '(a)' ) '  HERMITE_FUNCTION_VALUES returns values of '
  write ( *, '(a)' ) '  the Hermite functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               Hf(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_function_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine hermite_poly_phys_values_test ( )

!*****************************************************************************80
!
!! hermite_poly_phys_values_test() tests hermite_poly_phys_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLY_PHYS_VALUES_TEST:'
  write ( *, '(a)' ) '  HERMITE_POLY_PHYS_VALUES returns values of '
  write ( *, '(a)' ) '  the physicist''s Hermite polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               H(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_phys_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine hermite_poly_prob_values_test ( )

!*****************************************************************************80
!
!! hermite_poly_prob_values_test() tests hermite_poly_prob_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 February 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HERMITE_POLY_PROB_VALUES_TEST:'
  write ( *, '(a)' ) '  HERMITE_POLY_PROB_VALUES returns values of '
  write ( *, '(a)' ) '  the probabilist''s Hermite polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       X               He(N,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hermite_poly_prob_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine hyper_1f1_values_test ( )

!*****************************************************************************80
!
!! hyper_1f1_values_test() tests hyper_1f1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPER_1F1_VALUES_TEST:'
  write ( *, '(a)' ) '  HYPER_1F1_VALUES returns values of '
  write ( *, '(a)' ) '  the hypergeometric 1F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A           B           X       Hyper_1F1(A,B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_1f1_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine hyper_2f1_values_test ( )

!*****************************************************************************80
!
!! hyper_2f1_values_test() tests hyper_2f1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPER_2F1_VALUES_TEST:'
  write ( *, '(a)' ) '  HYPER_2F1_VALUES returns values of '
  write ( *, '(a)' ) '  the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A           B           C            X       Hyper_2F1(A,B,C,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_values ( n_data, a, b, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,f10.6,2x,g24.16)' ) a, b, c, x, fx

  end do

  return
end
subroutine hyper_2f1_complex_values_test ( )

!*****************************************************************************80
!
!! hyper_2f1_complex_values_test() tests hyper_2f1_complex_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 December 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  complex ( kind = ck ) fz
  integer n_data
  complex ( kind = ck ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hyper_2f1_complex_values_test():'
  write ( *, '(a)' ) '  hyper_2f1_complex_values returns values of '
  write ( *, '(a)' ) '  the hypergeometric 2F1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A           B           C            Z       Hyper_2F1(A,B,C;Z)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hyper_2f1_complex_values ( n_data, a, b, c, z, fz )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,(f10.6,f10.6),2x,(g14.6,g14.6))' ) &
      a, b, c, z, fz

  end do

  return
end
subroutine hypergeometric_cdf_values_test ( )

!*****************************************************************************80
!
!! hypergeometric_cdf_values_test() tests hypergeometric_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  integer pop
  integer sam
  integer suc
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'hypergeometric_cdf_values_test:'
  write ( *, '(a)' ) '  hypergeometric_cdf_values returns values of '
  write ( *, '(a)' ) '  the Hypergeometric Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       SAM       SUC       POP         X      HyperCDF(S,S,P)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_cdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) sam, suc, pop, x, fx

  end do

  return
end
subroutine hypergeometric_pdf_values_test ( )

!*****************************************************************************80
!
!! hypergeometric_pdf_values_test() tests hypergeometric_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  integer pop
  integer sam
  integer suc
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERGEOMETRIC_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Hypergeometric Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '       SAM       SUC       POP         X      HyperPDF(S,S,P)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_pdf_values ( n_data, sam, suc, pop, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8,2x,g24.16)' ) sam, suc, pop, x, fx

  end do

  return
end
subroutine hypergeometric_u_values_test ( )

!*****************************************************************************80
!
!! hypergeometric_u_values_test() tests hypergeometric_u_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'HYPERGEOMETRIC_U_VALUES_TEST:'
  write ( *, '(a)' ) '  HYPERGEOMETRIC_U_VALUES returns values of '
  write ( *, '(a)' ) '  the hypergeometric U function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '         A           B           X       HyperU(A,B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call hypergeometric_u_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f10.6,2x,f10.6,2x,f10.6,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine i0ml0_values_test ( )

!*****************************************************************************80
!
!! i0ml0_values_test() tests i0ml0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I0ML0_VALUES_TEST:'
  write ( *, '(a)' ) '  I0ML0_VALUES returns values of'
  write ( *, '(a)' ) '  the I0ML0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i0ml0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine i1ml1_values_test ( )

!*****************************************************************************80
!
!! i1ml1_values_test() tests i1ml1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I1ML1_VALUES_TEST:'
  write ( *, '(a)' ) '  I1ML1_VALUES returns values of'
  write ( *, '(a)' ) '  the I1ML1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i1ml1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine i4_fall_values_test ( )

!*****************************************************************************80
!
!! i4_fall_values_test() tests i4_fall_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fmn
  integer m
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_fall_values_test():'
  write ( *, '(a)' ) '  i4_fall_values() returns some exact values'
  write ( *, '(a)' ) '  of the falling factorial function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M         N      I4_FALL(M,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_fall_values ( n_data, m, n, fmn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) m, n, fmn

  end do

  return
end
subroutine i4_gpf_values_test ( )

!*****************************************************************************80
!
!! i4_gpf_values_test() tests i4_gpf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 November 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer gpf
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'i4_gpf_values_test():'
  write ( *, '(a)' ) '  i4_gpf_values() returns some exact values'
  write ( *, '(a)' ) '  of the greatest prime factor function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N      i4_gpf(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_gpf_values ( n_data, n, gpf )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, gpf

  end do

  return
end
subroutine i4_rise_values_test ( )

!*****************************************************************************80
!
!! i4_rise_values_test() tests i4_rise_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fmn
  integer m
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'I4_RISE_VALUES_TEST:'
  write ( *, '(a)' ) '  I4_RISE_VALUES returns some exact values'
  write ( *, '(a)' ) '  of the rising factorial function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         M         N      I4_RISE(M,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call i4_rise_values ( n_data, m, n, fmn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) m, n, fmn

  end do

  return
end
subroutine int_values_test ( )

!*****************************************************************************80
!
!! int_values_test() tests int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INT_VALUES_TEST:'
  write ( *, '(a)' ) '  INT_VALUES returns values of the INT function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X      INT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call int_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.4,2x,f14.4)' ) x, fx

  end do

  return
end
subroutine inverse_chi_square_pdf_values_test ( )

!*****************************************************************************80
!
!! inverse_chi_square_pdf_values_test() tests inverse_chi_square_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) df
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INVERSE_CHI_SQUARE_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  INVERSE_CHI_SQUARE_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the inverse chi square Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          DF              X            PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call inverse_chi_square_pdf_values ( n_data, df, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) df, x, fx

  end do

  return
end
subroutine inverse_gamma_pdf_values_test ( )

!*****************************************************************************80
!
!! inverse_gamma_pdf_values_test() tests inverse_gamma_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'INVERSE_GAMMA_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  INVERSE_GAMMA_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the inverse gamma Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          ALPHA           BETA             X            PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call inverse_gamma_pdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, beta, x, fx

  end do

  return
end
subroutine is_prime_values_test ( )

!*****************************************************************************80
!
!! is_prime_values_test() tests is_prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 December 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer n
  integer n_data
  logical tf

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'is_prime_values_test():'
  write ( *, '(a)' ) '  is_prime_values() returns values of the is_prime() function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         n      is_prime(n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call is_prime_values ( n_data, n, tf )

    if ( n_data == 0 ) then
      exit
    end if

    if ( tf ) then
      write ( *, '(2x,i8,2x,a)' ) n, 'true'
    else
      write ( *, '(2x,i8,2x,a)' ) n, 'false'
    end if

  end do

  return
end
subroutine jacobi_cn_values_test ( )

!*****************************************************************************80
!
!! jacobi_cn_values_test() tests jacobi_cn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2020
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) k
  real ( kind = rk8 ) m
  integer n_data
  real ( kind = rk8 ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jacobi_cn_values_test:'
  write ( *, '(a)' ) '  jacobi_cn_values_test() returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic CN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      U         M       CN(U,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_cn_values ( n_data, u, a, k, m, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,g24.16)' ) u, m, fx

  end do

  return
end
subroutine jacobi_dn_values_test ( )

!*****************************************************************************80
!
!! jacobi_dn_values_test() tests jacobi_dn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 November 2020
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) k
  real ( kind = rk8 ) m
  integer n_data
  real ( kind = rk8 ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jacobi_dn_values_test:'
  write ( *, '(a)' ) '  jacobi_dn_values() returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic DN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      U         M       DN(U,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_dn_values ( n_data, u, a, k, m, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,g24.16)' ) u, m, fx

  end do

  return
end
subroutine jacobi_poly_values_test ( )

!*****************************************************************************80
!
!! jacobi_poly_values_test() tests jacobi_poly_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 April 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JACOBI_POLY_VALUES_TEST:'
  write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns values of '
  write ( *, '(a)' ) '  the Jacobi polynomial.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         A         B      X       J(N,A,B)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_poly_values ( n_data, n, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,f24.16,2x,g24.16)' ) n, a, b, x, fx

  end do

  return
end
subroutine jacobi_sn_values_test ( )

!*****************************************************************************80
!
!! jacobi_sn_values_test() tests jacobi_sn_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 June 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) k
  real ( kind = rk8 ) m
  integer n_data
  real ( kind = rk8 ) u

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jacobi_sn_values_test:'
  write ( *, '(a)' ) '  jacobi_sn_values() returns values of '
  write ( *, '(a)' ) '  the Jacobi elliptic SN function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      U         M       SN(U,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jacobi_sn_values ( n_data, u, a, k, m, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,g14.6)' ) u, m, fx

  end do

  return
end
subroutine jed_ce_values_test ( )

!*****************************************************************************80
!
!! jed_ce_values_test() tests jed_ce_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  real ( kind = rk8 ) f
  real ( kind = rk8 ) jed
  integer n_data
  integer m
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JED_CE_VALUES_TEST:'
  write ( *, '(a)' ) '  JED_CE_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  YMDF, the corresponding year, month, day, fraction.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED          Y   M   D    F'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_ce_values ( n_data, jed, y, m, d, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,i6,2x,i2,2x,i2,2x,f6.4)' ) jed, y, m, d, f

  end do

  return
end
subroutine jed_mjd_values_test ( )

!*****************************************************************************80
!
!! jed_mjd_values_test() tests jed_mjd_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) jed
  integer n_data
  real ( kind = rk8 ) mjd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JED_MJD_VALUES_TEST:'
  write ( *, '(a)' ) '  JED_MJD_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  MJD, the corresponding Modified Julian Day count.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED           MJD'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_mjd_values ( n_data, jed, mjd )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,f12.2)' ) jed, mjd

  end do

  return
end
subroutine jed_rd_values_test ( )

!*****************************************************************************80
!
!! jed_rd_values_test() tests jed_rd_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) jed
  integer n_data
  real ( kind = rk8 ) rd

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JED_RD_VALUES_TEST:'
  write ( *, '(a)' ) '  JED_RD_VALUES returns:'
  write ( *, '(a)' ) '  JED, a Julian Ephemeris Date, and'
  write ( *, '(a)' ) '  RD, the corresponding Reingold Dershowitz Day count.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED            RD'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_rd_values ( n_data, jed, rd )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,f12.2)' ) jed, rd

  end do

  return
end
subroutine jed_weekday_values_test ( )

!*****************************************************************************80
!
!! jed_weekday_values_test() tests jed_weekday_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) jed
  integer n_data
  integer weekday
  character ( len = 9 ), dimension ( 7 ) :: weekday_name = (/ &
    'Sunday   ', 'Monday   ', 'Tuesday  ', 'Wednesday', 'Thursday ', &
    'Friday   ', 'Saturday ' /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'JED_WEEKDAY_VALUES_TEST:'
  write ( *, '(a)' ) '  JED_WEEKDAY_VALUES returns Julian Ephemeris Dates '
  write ( *, '(a)' ) '  (JED) and the corresponding weekday'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        JED     #  Weekday'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call jed_weekday_values ( n_data, jed, weekday )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.2,2x,i1,2x,a9)' ) jed, weekday, weekday_name(weekday)

  end do

  return
end
subroutine kei0_values_test ( )

!*****************************************************************************80
!
!! kei0_values_test() tests kei0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KEI0_VALUES_TEST:'
  write ( *, '(a)' ) '  KEI0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KEI of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KEI0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call kei0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine kei1_values_test ( )

!*****************************************************************************80
!
!! kei1_values_test() tests kei1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KEI1_VALUES_TEST:'
  write ( *, '(a)' ) '  KEI1_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KEI of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KEI1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call kei1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine ker0_values_test ( )

!*****************************************************************************80
!
!! ker0_values_test() tests ker0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'KER0_VALUES_TEST:'
  write ( *, '(a)' ) '  KER0_VALUES returns values of '
  write ( *, '(a)' ) '  the Kelvin function KER of order 0'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KER0'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ker0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine ker1_values_test ( )

!*****************************************************************************80
!
!! ker1_values_test() tests ker1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 June 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ker1_values_test():'
  write ( *, '(a)' ) '  ker1_values() returns values of '
  write ( *, '(a)' ) '  the Kelvin function KER of order 1'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           KER1'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call ker1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine knapsack_values_test ( )

!*****************************************************************************80
!
!! knapsack_values_test() tests knapsack_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 November 2024
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer k
  integer n
  integer n_data
  real ( kind = rk8 ) r
  integer, allocatable :: s(:)
  integer, allocatable :: v(:)
  integer, allocatable :: w(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'knapsack_values_test():'
  write ( *, '(a)' ) '  Test knapsack_values().'

  n_data = 0

  do
!
!  First call returns value of n.
!
    n = 0
    call knapsack_values ( n_data, n, v, w, s, k )

    if ( n == 0 ) then
      exit
    end if

    allocate ( v(1:n) )
    allocate ( w(1:n) )
    allocate ( s(1:n) )
!
!  Second call returns data, and increments n_data.
!
    call knapsack_values ( n_data, n, v, w, s, k )

    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Problem #', n_data
    write ( *, '(a,i4)' ) '  Number of items is ', n
    write ( *, '(a,i6)' ) '  Knapsack weight limit is ', k
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '   Item 0/1  Value  Weight  Value/Weight'
    write ( *, '(a)' ) ''
    do i = 1, n
      r = real ( v(i), kind = rk8 ) / real ( w(i), kind = rk8 )
      write ( *, '(2x,i5,2x,i2,2x,i8,2x,i8,2x,f7.2)' ) &
        i, s(i), v(i), w(i), r
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a,2x,i2,2x,i8,2x,i8,2x,f7.2)' ) &
      '  Taken  ', &
      sum ( s(1:n) ), &
      dot_product ( s, v ), &
      dot_product ( s, w ), &
      real ( dot_product ( s, v ), kind = rk8 ) / &
      real ( dot_product ( s, w ), kind = rk8 )

    deallocate ( s )
    deallocate ( w )
    deallocate ( v )

  end do

  return
end
subroutine laguerre_associated_values_test ( )

!*****************************************************************************80
!
!! laguerre_associated_values_test() tests laguerre_associated_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_ASSOCIATED_VALUES_TEST:'
  write ( *, '(a)' ) '  LAGUERRE_ASSOCIATED_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Laguerre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M      X            L(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine laguerre_general_values_test ( )

!*****************************************************************************80
!
!! laguerre_general_values_test() tests laguerre_general_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_GENERAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LAGUERRE_GENERAL_VALUES returns values of '
  write ( *, '(a)' ) '  the generalized Laguerre function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       A      X            L(N,A)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_general_values ( n_data, n, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g24.16)' ) n, a, x, fx

  end do

  return
end
subroutine laguerre_polynomial_values_test ( )

!*****************************************************************************80
!
!! laguerre_polynomial_values_test() tests laguerre_polynomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAGUERRE_POLYNOMIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LAGUERRE_POLYNOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the Laguerre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            L(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laguerre_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine lambert_w_values_test ( )

!*****************************************************************************80
!
!! lambert_w_values_test() tests lambert_w_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 June 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lambert_w_values_test():'
  write ( *, '(a)' ) '  lambert_w_values() returns values of '
  write ( *, '(a)' ) '  the Lambert W function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           Branch       W(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lambert_w_values ( n_data, x, fx, b )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,i2,2x,g24.16)' ) x, b, fx

  end do

  return
end
subroutine laplace_cdf_values_test ( )

!*****************************************************************************80
!
!! laplace_cdf_values_test() tests laplace_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LAPLACE_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  LAPLACE_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Laplace CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          MU              BETA            X           CDF(MU,BETA,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call laplace_cdf_values ( n_data, mu, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, beta, x, fx

  end do

  return
end
subroutine legendre_associated_values_test ( )

!*****************************************************************************80
!
!! legendre_associated_values_test() tests legendre_associated_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Legendre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine legendre_associated_normalized_values_test ( )

!*****************************************************************************80
!
!! legendre_associated_normalized_values_test() tests legendre_associated_normalized_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns values of '
  write ( *, '(a)' ) '  the normalized associated Legendre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine legendre_associated_normalized_sphere_values_test ( )

!*****************************************************************************80
!
!! legendre_associated_normalized_sphere_values_test() tests legendre_associated_normalized_sphere_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer m
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES returns values of '
  write ( *, '(a)' ) '  the associated Legendre polynomials normalized for a sphere.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M          X            P(N,M)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_associated_normalized_sphere_values ( n_data, n, m, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g24.16)' ) n, m, x, fx

  end do

  return
end
subroutine legendre_function_q_values_test ( )

!*****************************************************************************80
!
!! legendre_function_q_values_test() tests legendre_function_q_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_FUNCTION_Q_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_FUNCTION_Q_VALUES returns values of '
  write ( *, '(a)' ) '  the Legendre QN polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           Q(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_function_q_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, x, fx

  end do

  return
end
subroutine legendre_normalized_polynomial_values_test ( )

!*****************************************************************************80
!
!! legendre_normalized_polynomial_values_test() tests legendre_normalized_polynomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_NORMALIZED_POLYNOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the normalized Legendre polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            Pn(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_normalized_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine legendre_polynomial_values_test ( )

!*****************************************************************************80
!
!! legendre_polynomial_values_test() tests legendre_polynomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_POLYNOMIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_POLYNOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the Legendre PN polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            P(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine legendre_shifted_polynomial_values_test ( )

!*****************************************************************************80
!
!! legendre_shifted_polynomial_values_test() tests legendre_shifted_polynomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LEGENDRE_SHIFTED_POLYNOMIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LEGENDRE_SHIFTED_POLYNOMIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the shifted Legendre PN polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X            P(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call legendre_shifted_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine lerch_values_test ( )

!*****************************************************************************80
!
!! lerch_values_test() tests lerch_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  integer s
  real ( kind = rk8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LERCH_VALUES_TEST:'
  write ( *, '(a)' ) '  LERCH_VALUES returns values of '
  write ( *, '(a)' ) '  the Lerch transcendent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Z            S        A        CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lerch_values ( n_data, z, s, a, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,i8,2x,g10.4,2x,g24.16)' ) z, s, a, fx

  end do

  return
end
subroutine lobachevsky_values_test ( )

!*****************************************************************************80
!
!! lobachevsky_values_test() tests lobachevsky_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOBACHEVSKY_VALUES_TEST:'
  write ( *, '(a)' ) '  LOBACHEVSKY_VALUES returns values of'
  write ( *, '(a)' ) '  the Lobachevsky function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lobachevsky_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine lobatto_polynomial_values_test ( )

!*****************************************************************************80
!
!! lobatto_polynomial_values_test() tests lobatto_polynomial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lobatto_polynomial_values_test():'
  write ( *, '(a)' ) '  lobatto_polynomial_values() returns values of '
  write ( *, '(a)' ) '  the complete Lobatto polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           Lo(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lobatto_polynomial_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine lobatto_polynomial_derivative_values_test ( )

!*****************************************************************************80
!
!! lobatto_polynomial_derivative_values_test() tests lobatto_polynomial_derivative_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 May 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'lobatto_polynomial_derivative_values_test():'
  write ( *, '(a)' ) '  lobatto_polynomial_derivative_values() returns '
  write ( *, '(a)' ) '  derivatives of the complete Lobatto polynomials.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N          X           Lo''(N)(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call lobatto_polynomial_derivative_values ( n_data, n, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g14.6,2x,g24.16)' ) n, x, fx

  end do

  return
end
subroutine log_values_test ( )

!*****************************************************************************80
!
!! log_values_test() tests log_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 September 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG_VALUES_TEST:'
  write ( *, '(a)' ) &
    '  LOG_VALUES returns values of the natural logarithm function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            LN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine log_normal_cdf_values_test ( )

!*****************************************************************************80
!
!! log_normal_cdf_values_test() tests log_normal_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG_NORMAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  LOG_NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Log Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_normal_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine log_series_cdf_values_test ( )

!*****************************************************************************80
!
!! log_series_cdf_values_test() tests log_series_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG_SERIES_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  LOG_SERIES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Log Series Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     T          N   CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log_series_cdf_values ( n_data, t, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,i8,2x,g24.16)' ) t, n, fx

  end do

  return
end
subroutine log10_values_test ( )

!*****************************************************************************80
!
!! log10_values_test() tests log10_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOG10_VALUES_TEST:'
  write ( *, '(a)' ) &
    '  LOG10_VALUES returns values of the logarithm base 10 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X            LOG10(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call log10_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine logarithmic_integral_values_test ( )

!*****************************************************************************80
!
!! logarithmic_integral_values_test() tests logarithmic_integral_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOGARITHMIC_INTEGRAL_VALUES_TEST:'
  write ( *, '(a)' ) '  LOGARITHMIC_INTEGAL_VALUES returns values of'
  write ( *, '(a)' ) '  the logarithmic integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            LI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call logarithmic_integral_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine logistic_cdf_values_test ( )

!*****************************************************************************80
!
!! logistic_cdf_values_test() tests logistic_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'LOGISTIC_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  LOGISTIC_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Logistic Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     BETA     X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call logistic_cdf_values ( n_data, mu, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, beta, x, fx

  end do

  return
end
subroutine mathieu_even_values_test ( )

!*****************************************************************************80
!
!! mathieu_even_values_test() tests mathieu_even_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  integer n_data
  integer q
  integer r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATHIEU_EVEN_VALUES_TEST:'
  write ( *, '(a)' ) '  MATHIEU_EVEN_VALUES returns values of '
  write ( *, '(a)' ) '  the eigenvalues of the Mathieu differential'
  write ( *, '(a)' ) '  equation associated with even periodic solutions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         Q       A(R,Q)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mathieu_even_values ( n_data, r, q, a )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) r, q, a

  end do

  return
end
subroutine mathieu_odd_values_test ( )

!*****************************************************************************80
!
!! mathieu_odd_values_test() tests mathieu_odd_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    19 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) b
  integer n_data
  integer q
  integer r

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MATHIEU_ODD_VALUES_TEST:'
  write ( *, '(a)' ) '  MATHIEU_ODD_VALUES returns values of '
  write ( *, '(a)' ) '  the eigenvalues of the Mathieu differential'
  write ( *, '(a)' ) '  equation associated with odd periodic solutions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         R         Q       B(R,Q)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mathieu_odd_values ( n_data, r, q, b )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) r, q, b

  end do

  return
end
subroutine mcnugget_number_values_test ( )

!*****************************************************************************80
!
!! mcnugget_number_values_test() tests mcnugget_number_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 September 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer m
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mcnugget_number_values_test()::'
  write ( *, '(a)' ) '  mcnugget_number_values() returns values of '
  write ( *, '(a)' ) '  the McNugget numbers.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    N  M(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mcnugget_number_values ( n_data, n, m )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i3,2x,i3)' ) n, m

  end do

  return
end
subroutine mersenne_prime_values_test ( )

!*****************************************************************************80
!
!! mersenne_prime_values_test() tests mersenne_prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 January 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mersenne_prime_values_test():'
  write ( *, '(a)' ) '  mersenne_prime_values() returns indices of Mersenne primes.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       n         mersenne_prime(n)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mersenne_prime_values ( n_data, n )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n_data, n

  end do

  return
end
subroutine mertens_values_test ( )

!*****************************************************************************80
!
!! mertens_values_test() tests mertens_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mertens_values_test():'
  write ( *, '(a)' ) '  mertens_values() returns values of the Mertens function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N         MERTENS(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mertens_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine mittag_leffler_ea_values_test ( )

!*****************************************************************************80
!
!! mittag_leffler_ea_values_test() tests mittag_leffler_ea_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 February 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MITTAG_LEFFLER_EA_VALUES_TEST:'
  write ( *, '(a)' ) '  MITTAG_LEFFLER_EA_VALUES returns values of '
  write ( *, '(a)' ) '  the one-parameter Mittag-Leffler function E(A;X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A        X         E(A;X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mittag_leffler_ea_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.4,2x,f12.4,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine mittag_leffler_eab_values_test ( )

!*****************************************************************************80
!
!! mittag_leffler_eab_values_test() tests mittag_leffler_eab_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 February 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MITTAG_LEFFLER_EAB_VALUES_TEST:'
  write ( *, '(a)' ) '  MITTAG_LEFFLER_EAB_VALUES returns values of '
  write ( *, '(a)' ) '  the two-parameter Mittag-Leffler function E(A;X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         A        B             X         EB(A;B,X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call mittag_leffler_eab_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f12.4,2x,f12.4,2x,f12.4,2x,g24.16)' ) a, b, x, fx

  end do

  return
end
subroutine moebius_values_test ( )

!*****************************************************************************80
!
!! moebius_values_test() tests moebius_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MOEBIUS_VALUES_TEST:'
  write ( *, '(a)' ) '  MOEBIUS_VALUES returns values of '
  write ( *, '(a)' ) '  the Moebius function.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '       N         MU(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call moebius_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine multinomial_pdf_values_test ( )

!*****************************************************************************80
!
!! multinomial_pdf_values_test() tests multinomial_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer m
  integer n
  integer n_data1
  integer n_data2
  real ( kind = rk8 ), allocatable :: p(:)
  real ( kind = rk8 ) pdf
  integer, allocatable :: x(:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MULTINOMIAL_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  MULTINOMIAL_PDF_VALUES stores values of the Multinomial PDF.'
  write ( *, '(a)' ) '  Given M possible outcomes on a single trial,'
  write ( *, '(a)' ) '  with each outcome having probability P,'
  write ( *, '(a)' ) '  PDF is the probability that after N trials,'
  write ( *, '(a)' ) '  outcome I occurred X(I) times.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     N     M     I      P        X        PDF()'

  n_data1 = 0
  n_data2 = 0

  do

    call multinomial_pdf_sizes ( n_data1, m )

    if ( n_data1 == 0 ) then
      exit
    end if

    allocate ( p(1:m) )
    allocate ( x(1:m) )

    call multinomial_pdf_values ( n_data2, m, n, p, x, pdf )

    write ( *, '(a)' ) ''
    do i = 1, m
      write ( *, '(14x,i4,2x,f8.4,2x,i4)' ) i, p(i), x(i)
    end do
    write ( *, '(2x,i4,2x,i4,24x,g14.6)' ) n, m, pdf

    deallocate ( p )
    deallocate ( x )

  end do

  return
end
subroutine negative_binomial_cdf_values_test ( )

!*****************************************************************************80
!
!! negative_binomial_cdf_values_test() tests negative_binomial_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) cdf
  integer f
  integer n_data
  real ( kind = rk8 ) p
  integer s

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NEGATIVE_BINOMIAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  NEGATIVE_BINOMIAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Negative Binomial Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       F       S         P         CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call negative_binomial_cdf_values ( n_data, f, s, p, cdf )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) f, s, p, cdf

  end do

  return
end
subroutine nine_j_values_test ( )

!*****************************************************************************80
!
!! nine_j_values_test() tests nine_j_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) j1
  real ( kind = rk8 ) j2
  real ( kind = rk8 ) j3
  real ( kind = rk8 ) j4
  real ( kind = rk8 ) j5
  real ( kind = rk8 ) j6
  real ( kind = rk8 ) j7
  real ( kind = rk8 ) j8
  real ( kind = rk8 ) j9
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NINE_J_VALUES_TEST:'
  write ( *, '(a)' ) '  NINE_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 9J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      J4      J5      J6' // &
    '      J7      J8      J9        NINE_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call nine_j_values ( n_data, j1, j2, j3, j4, j5, j6, j7, j8, j9, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,' // &
      'f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, j4, j5, j6, j7, j8, j9, fx

  end do

  return
end
subroutine normal_01_cdf_values_test ( )

!*****************************************************************************80
!
!! normal_01_cdf_values_test() tests normal_01_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_01_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  NORMAL_01_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_cdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine normal_01_pdf_values_test ( )

!*****************************************************************************80
!
!! normal_01_pdf_values_test() tests normal_01_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_01_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  NORMAL_01_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X                  PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_01_pdf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine normal_cdf_values_test ( )

!*****************************************************************************80
!
!! normal_cdf_values_test() tests normal_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  NORMAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     SIGMA    X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_cdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine normal_pdf_values_test ( )

!*****************************************************************************80
!
!! normal_pdf_values_test() tests normal_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 July 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NORMAL_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  NORMAL_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Normal Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      MU     SIGMA    X                  PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call normal_pdf_values ( n_data, mu, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) mu, sigma, x, fx

  end do

  return
end
subroutine omega_values_test ( )

!*****************************************************************************80
!
!! omega_values_test() tests omega_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OMEGA_VALUES_TEST:'
  write ( *, '(a)' ) '  OMEGA_VALUES returns values of '
  write ( *, '(a)' ) '  the Omega function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     N           OMEGA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call omega_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, fn

  end do

  return
end
subroutine owen_values_test ( )

!*****************************************************************************80
!
!! owen_values_test() tests owen_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) h
  integer n_data
  real ( kind = rk8 ) t

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'OWEN_VALUES_TEST:'
  write ( *, '(a)' ) '  OWEN_VALUES returns values of '
  write ( *, '(a)' ) '  the Owen T function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      H       A       T'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call owen_values ( n_data, h, a, t )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) h, a, t

  end do

  return
end
subroutine partition_count_values_test ( )

!*****************************************************************************80
!
!! partition_count_values_test() tests partition_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTITION_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  PARTITION_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the integer partition count function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         P(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_count_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine partition_distinct_count_values_test ( )

!*****************************************************************************80
!
!! partition_distinct_count_values_test() tests partition_distinct_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PARTITION_DISTINCT_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  PARTITION_DISTINCT_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the integer distinct partition count function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         Q(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call partition_distinct_count_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine phi_values_test ( )

!*****************************************************************************80
!
!! phi_values_test() tests phi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PHI_VALUES_TEST:'
  write ( *, '(a)' ) '  PHI_VALUES returns values of '
  write ( *, '(a)' ) '  the PHI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         PHI(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call phi_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine pi_values_test ( )

!*****************************************************************************80
!
!! pi_values_test() tests pi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PI_VALUES_TEST:'
  write ( *, '(a)' ) '  PI_VALUES returns values of '
  write ( *, '(a)' ) '  the PI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '             N         PI(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call pi_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, fn

  end do

  return
end
subroutine poisson_cdf_values_test ( )

!*****************************************************************************80
!
!! poisson_cdf_values_test() tests poisson_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  integer n_data
  integer x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POISSON_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  POISSON_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Poisson Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A          X    CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call poisson_cdf_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,i8,2x,g24.16)' ) a, x, fx

  end do

  return
end
subroutine polylogarithm_values_test ( )

!*****************************************************************************80
!
!! polylogarithm_values_test() tests polylogarithm_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n
  integer n_data
  real ( kind = rk8 ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYLOGARITHM_VALUES_TEST:'
  write ( *, '(a)' ) '  POLYLOGARITHM_VALUES returns values of '
  write ( *, '(a)' ) '  the polylogarithm.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N      Z        FX'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call polylogarithm_values ( n_data, n, z, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g10.4,2x,g24.16)' ) n, z, fx

  end do

  return
end
subroutine polyomino_chiral_count_values_test ( )

!*****************************************************************************80
!
!! polyomino_chiral_count_values_test() tests polyomino_chiral_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )
  integer, parameter :: ik = selected_int_kind ( 18 )

  integer n_data
  integer ( kind = ik ) number
  integer order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYOMINO_CHIRAL_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  POLYOMINO_CHIRAL_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the number of chiral polyominoes of given order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   ORDER         NUMBER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call polyomino_chiral_count_values ( n_data, order, number )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i24)' ) order, number

  end do

  return
end
subroutine polyomino_fixed_count_values_test ( )

!*****************************************************************************80
!
!! polyomino_fixed_count_values_test() tests polyomino_fixed_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 April 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )
  integer, parameter :: ik = selected_int_kind ( 18 )

  integer n_data
  integer ( kind = ik ) number
  integer order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYOMINO_FIXED_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  POLYOMINO_FIXED_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the number of fixed polyominoes of given order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   ORDER         NUMBER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call polyomino_fixed_count_values ( n_data, order, number )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i24)' ) order, number

  end do

  return
end
subroutine polyomino_free_count_values_test ( )

!*****************************************************************************80
!
!! polyomino_free_count_values_test() tests polyomino_free_count_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 April 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )
  integer, parameter :: ik = selected_int_kind ( 18 )

  integer n_data
  integer ( kind = ik ) number
  integer order

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'POLYOMINO_FREE_COUNT_VALUES_TEST:'
  write ( *, '(a)' ) '  POLYOMINO_FREE_COUNT_VALUES returns values of '
  write ( *, '(a)' ) '  the number of free polyominoes of given order.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   ORDER         NUMBER'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call polyomino_free_count_values ( n_data, order, number )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i24)' ) order, number

  end do

  return
end
subroutine prandtl_values_test ( )

!*****************************************************************************80
!
!! prandtl_values_test() tests prandtl_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) pr
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRANDTL_VALUES_TEST:'
  write ( *, '(a)' ) '  PRANDTL_VALUES returns values of '
  write ( *, '(a)' ) '  the Prandtl number of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            Pr(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call prandtl_values ( n_data, tc, p, pr )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, pr

  end do

  return
end
subroutine prime_values_test ( )

!*****************************************************************************80
!
!! prime_values_test() tests prime_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer n_data
  integer p

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PRIME_VALUES_TEST:'
  write ( *, '(a)' ) '  PRIME_VALUES returns values of '
  write ( *, '(a)' ) '  the prime function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '           N          P[N]'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call prime_values ( n_data, n, p )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i12,2x,i12)' ) n, p

  end do

  return
end
subroutine psat_values_test ( )

!*****************************************************************************80
!
!! psat_values_test() tests psat_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) psat
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSAT_VALUES_TEST:'
  write ( *, '(a)' ) '  PSAT_VALUES returns values of '
  write ( *, '(a)' ) '  the saturation pressure of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            PSAT(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psat_values ( n_data, tc, psat )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) tc, psat

  end do

  return
end
subroutine psi_values_test ( )

!*****************************************************************************80
!
!! psi_values_test() tests psi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PSI_VALUES_TEST:'
  write ( *, '(a)' ) '  PSI_VALUES returns values of '
  write ( *, '(a)' ) '  the PSI function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            PSI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call psi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine r8_factorial_values_test ( )

!*****************************************************************************80
!
!! r8_factorial_values_test() tests r8_factorial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FACTORIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  R8_FACTORIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the factorial function (using real arithmetic).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        N       Factorial(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, fn

  end do

  return
end
subroutine r8_factorial_log_values_test ( )

!*****************************************************************************80
!
!! r8_factorial_log_values_test() tests r8_factorial_log_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FACTORIAL_LOG_VALUES_TEST:'
  write ( *, '(a)' ) '  R8_FACTORIAL_LOG_VALUES returns values of '
  write ( *, '(a)' ) '  the logarithm of the factorial function '
  write ( *, '(a)' ) '  (using real arithmetic).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        N       Log(Factorial(N))'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial_log_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, fn

  end do

  return
end
subroutine r8_factorial2_values_test ( )

!*****************************************************************************80
!
!! r8_factorial2_values_test() tests r8_factorial2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FACTORIAL2_VALUES_TEST:'
  write ( *, '(a)' ) '  R8_FACTORIAL2_VALUES returns values of '
  write ( *, '(a)' ) '  the double factorial function (using real arithmetic).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        N       Factorial2(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_factorial2_values ( n_data, n, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, f

  end do

  return
end
subroutine r8_fall_values_test ( )

!*****************************************************************************80
!
!! r8_fall_values_test() tests r8_fall_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_FALL_VALUES_TEST:'
  write ( *, '(a)' ) '  R8_FALL_VALUES returns some exact values'
  write ( *, '(a)' ) '  of the falling factorial function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X         N      R8_FALL(X,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_fall_values ( n_data, x, n, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.4,2x,i8,2x,f12.4)' ) x, n, f

  end do

  return
end
subroutine r8_rise_values_test ( )

!*****************************************************************************80
!
!! r8_rise_values_test() tests r8_rise_values().
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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) f
  integer n
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8_RISE_VALUES_TEST:'
  write ( *, '(a)' ) '  R8_RISE_VALUES returns some exact values'
  write ( *, '(a)' ) '  of the rising factorial function:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X         N      R8_RISE(X,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call r8_rise_values ( n_data, x, n, f )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.4,2x,i8,2x,g14.6)' ) x, n, f

  end do

  return
end
subroutine rayleigh_cdf_values_test ( )

!*****************************************************************************80
!
!! rayleigh_cdf_values_test() tests rayleigh_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RAYLEIGH_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  RAYLEIGH_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Rayleigh Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      SIGMA    X                  CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call rayleigh_cdf_values ( n_data, sigma, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g24.16)' ) sigma, x, fx

  end do

  return
end
subroutine scaled_inverse_chi_square_pdf_values_test ( )

!*****************************************************************************80
!
!! scaled_inverse_chi_square_pdf_values_test() tests scaled_inverse_chi_square_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 August 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) df
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x
  real ( kind = rk8 ) xi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SCALED_INVERSE_CHI_SQUARE_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  SCALED_INVERSE_CHI_SQUARE_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the scaled inverse chi square Probability Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          DF              XI             X            PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call scaled_inverse_chi_square_pdf_values ( n_data, df, xi, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) df, xi, x, fx

  end do

  return
end
subroutine secvir_values_test ( )

!*****************************************************************************80
!
!! secvir_values_test() tests secvir_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) tc
  real ( kind = rk8 ) vir

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SECVIR_VALUES_TEST:'
  write ( *, '(a)' ) '  SECVIR_VALUES returns values of '
  write ( *, '(a)' ) '  the second virial coefficient of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            VIR(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call secvir_values ( n_data, tc, vir )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) tc, vir

  end do

  return
end
subroutine shi_values_test ( )

!*****************************************************************************80
!
!! shi_values_test() tests shi_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SHI_VALUES_TEST:'
  write ( *, '(a)' ) '  SHI_VALUES returns values of '
  write ( *, '(a)' ) '  the Hyperbolic Sine Integral function SHI(X).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SHI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call shi_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine si_values_test ( )

!*****************************************************************************80
!
!! si_values_test() tests si_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SI_VALUES_TEST:'
  write ( *, '(a)' ) '  SI_VALUES returns values of '
  write ( *, '(a)' ) '  the sine integral function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            SI(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call si_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine sigma_values_test ( )

!*****************************************************************************80
!
!! sigma_values_test() tests sigma_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIGMA_VALUES_TEST:'
  write ( *, '(a)' ) '  SIGMA_VALUES returns values of '
  write ( *, '(a)' ) '  the SIGMA function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         SIGMA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sigma_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine sin_values_test ( )

!*****************************************************************************80
!
!! sin_values_test() tests sin_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_VALUES_TEST:'
  write ( *, '(a)' ) '  SIN_VALUES returns values of the sine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SIN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine sin_degree_values_test ( )

!*****************************************************************************80
!
!! sin_degree_values_test() tests sin_degree_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_DEGREE_VALUES_TEST:'
  write ( *, '(a)' ) '  SIN_DEGREE_VALUES returns values of the sine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SIN_DEGREE(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_degree_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine sin_power_int_values_test ( )

!*****************************************************************************80
!
!! sin_power_int_values_test() tests sin_power_int_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIN_POWER_INT_VALUES_TEST:'
  write ( *, '(a)' ) '  SIN_POWER_INT_VALUES returns values of '
  write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B              N    SIN_POWER_INT(A,B,N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sin_power_int_values ( n_data, a, b, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,i8,2x,g14.6)' ) a, b, n, fx

  end do

  return
end
subroutine sinh_values_test ( )

!*****************************************************************************80
!
!! sinh_values_test() tests sinh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SINH_VALUES_TEST:'
  write ( *, '(a)' ) &
    '  SINH_VALUES returns values of the hyperbolic sine function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           SINH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sinh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine six_j_values_test ( )

!*****************************************************************************80
!
!! six_j_values_test() tests six_j_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) j1
  real ( kind = rk8 ) j2
  real ( kind = rk8 ) j3
  real ( kind = rk8 ) j4
  real ( kind = rk8 ) j5
  real ( kind = rk8 ) j6
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIX_J_VALUES_TEST:'
  write ( *, '(a)' ) '  SIX_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 6J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      J4      J5      J6        SIX_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, j4, j5, j6, fx

  end do

  return
end
subroutine sound_values_test ( )

!*****************************************************************************80
!
!! sound_values_test() tests sound_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) c
  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SOUND_VALUES_TEST:'
  write ( *, '(a)' ) '  SOUND_VALUES returns values of '
  write ( *, '(a)' ) '  the spead of sound in water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            C(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sound_values ( n_data, tc, p, c )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, c

  end do

  return
end
subroutine sphere_unit_area_values_test ( )

!*****************************************************************************80
!
!! sphere_unit_area_values_test() tests sphere_unit_area_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_UNIT_AREA_VALUES_TEST:'
  write ( *, '(a)' ) '  SPHERE_UNIT_AREA_VALUES returns values of '
  write ( *, '(a)' ) '  the area of the unit sphere in various dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      N            Area'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_area_values ( n_data, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,g24.16)' ) n, fx

  end do

  return
end
subroutine sphere_unit_volume_values_test ( )

!*****************************************************************************80
!
!! sphere_unit_volume_values_test() tests sphere_unit_volume_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 May 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  integer n

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERE_UNIT_VOLUME_VALUES_TEST:'
  write ( *, '(a)' ) '  SPHERE_UNIT_VOLUME_VALUES returns values of '
  write ( *, '(a)' ) '  the volume of the unit sphere in various dimensions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      N            Volume'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sphere_unit_volume_values ( n_data, n, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,g24.16)' ) n, fx

  end do

  return
end
subroutine spherical_harmonic_values_test ( )

!*****************************************************************************80
!
!! spherical_harmonic_values_test() tests spherical_harmonic_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer l
  integer m
  integer n_data
  real ( kind = rk8 ) phi
  real ( kind = rk8 ) theta
  real ( kind = rk8 ) yi
  real ( kind = rk8 ) yr

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SPHERICAL_HARMONIC_VALUES_TEST:'
  write ( *, '(a)' ) '  SPHERICAL_HARMONIC_VALUES returns values of '
  write ( *, '(a)' ) '  the spherical harmonic functions.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '   L   M    THETA       PHI       Yr                         Yi'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call spherical_harmonic_values ( n_data, l, m, theta, phi, yr, yi )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i2,2x,i2,2x,f8.4,2x,f8.4,2x,g24.16,2x,g24.16)' ) &
      l, m, theta, phi, yr, yi

  end do

  return
end
subroutine sqrt_values_test ( )

!*****************************************************************************80
!
!! sqrt_values_test() tests sqrt_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    24 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SQRT_VALUES_TEST:'
  write ( *, '(a)' ) '  SQRT evaluates the square root function.'
  write ( *, '(a)' ) '  SQRT_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         X      SQRT(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call sqrt_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f14.4,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine stirling1_values_test ( )

!*****************************************************************************80
!
!! stirling1_values_test() tests stirling1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer n_data
  integer m
  integer s1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STIRLING1_VALUES_TEST:'
  write ( *, '(a)' ) '  STIRLING1_VALUES returns values of '
  write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N         M        S1(N,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stirling1_values ( n_data, n, m, s1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s1

  end do

  return
end
subroutine stirling2_values_test ( )

!*****************************************************************************80
!
!! stirling2_values_test() tests stirling2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer n_data
  integer m
  integer s2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STIRLING2_VALUES_TEST:'
  write ( *, '(a)' ) '  STIRLING2_VALUES returns values of '
  write ( *, '(a)' ) '  the Stirling numbers of the first kind.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N       M        S2(N,M)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stirling2_values ( n_data, n, m, s2 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,i12)' ) n, m, s2

  end do

  return
end
subroutine stromgen_values_test ( )

!*****************************************************************************80
!
!! stromgen_values_test() tests stromgen_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STROMGEN_VALUES_TEST:'
  write ( *, '(a)' ) '  STROMGEN_VALUES returns values of'
  write ( *, '(a)' ) '  the Stromgen function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X          F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call stromgen_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine struve_h0_values_test ( )

!*****************************************************************************80
!
!! struve_h0_values_test() tests struve_h0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRUVE_H0_VALUES_TEST:'
  write ( *, '(a)' ) '  STRUVE_H0_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve H0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            H0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine struve_h1_values_test ( )

!*****************************************************************************80
!
!! struve_h1_values_test() tests struve_h1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRUVE_H1_VALUES_TEST:'
  write ( *, '(a)' ) '  STRUVE_H1_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve H1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            H1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_h1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine struve_l0_values_test ( )

!*****************************************************************************80
!
!! struve_l0_values_test() tests struve_l0_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRUVE_L0_VALUES_TEST:'
  write ( *, '(a)' ) '  STRUVE_L0_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve L0 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            L0(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l0_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine struve_l1_values_test ( )

!*****************************************************************************80
!
!! struve_l1_values_test() tests struve_l1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STRUVE_L1_VALUES_TEST:'
  write ( *, '(a)' ) '  STRUVE_L1_VALUES returns values of '
  write ( *, '(a)' ) '  the Struve L1 function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            L1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call struve_l1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine student_cdf_values_test ( )

!*****************************************************************************80
!
!! student_cdf_values_test() tests student_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 November 2005
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) c
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STUDENT_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  STUDENT_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Student T Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      C     X       CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_cdf_values ( n_data, c, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,g24.16)' ) c, x, fx

  end do

  return
end
subroutine student_noncentral_cdf_values_test ( )

!*****************************************************************************80
!
!! student_noncentral_cdf_values_test() tests student_noncentral_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer df
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'STUDENT_NONCENTRAL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  STUDENT_NONCENTRAL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the noncentral Student T Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      LAMBDA       DF     CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call student_noncentral_cdf_values ( n_data, df, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g10.4,2x,g10.4,2x,i8,2x,g14.6)' ) x, lambda, df, fx

  end do

  return
end
subroutine subfactorial_values_test ( )

!*****************************************************************************80
!
!! subfactorial_values_test() tests subfactorial_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SUBFACTORIAL_VALUES_TEST:'
  write ( *, '(a)' ) '  SUBFACTORIAL_VALUES returns values of '
  write ( *, '(a)' ) '  the subfactorial function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N     Subfactorial[N]'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call subfactorial_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine surten_values_test ( )

!*****************************************************************************80
!
!! surten_values_test() tests surten_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SURTEN_VALUES_TEST:'
  write ( *, '(a)' ) '  SURTEN_VALUES returns values of '
  write ( *, '(a)' ) '  the surface tension of water '
  write ( *, '(a)' ) '  as a function of temperature.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            SIGMA(T)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call surten_values ( n_data, tc, sigma )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) tc, sigma

  end do

  return
end
subroutine synch1_values_test ( )

!*****************************************************************************80
!
!! synch1_values_test() tests synch1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SYNCH1_VALUES_TEST:'
  write ( *, '(a)' ) '  SYNCH1_VALUES returns values of '
  write ( *, '(a)' ) '  the synchrotron radiation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            S1(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch1_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine synch2_values_test ( )

!*****************************************************************************80
!
!! synch2_values_test() tests synch2_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SYNCH2_VALUES_TEST:'
  write ( *, '(a)' ) '  SYNCH2_VALUES returns values of '
  write ( *, '(a)' ) '  the synchrotron radiation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            S2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call synch2_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tan_values_test ( )

!*****************************************************************************80
!
!! tan_values_test() tests tan_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'tan_values_test():'
  write ( *, '(a)' ) '  tan_values() returns values of the tangent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           TAN(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tan_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine tanh_values_test ( )

!*****************************************************************************80
!
!! tanh_values_test() tests tanh_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    23 June 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TANH_VALUES_TEST:'
  write ( *, '(a)' ) &
    '  TANH_VALUES returns values of the hyperbolic tangent function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          X           TANH(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tanh_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine tau_values_test ( )

!*****************************************************************************80
!
!! tau_values_test() tests tau_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer fn
  integer n
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TAU_VALUES_TEST:'
  write ( *, '(a)' ) '  TAU_VALUES returns values of '
  write ( *, '(a)' ) '  the TAU function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       N         TAU(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tau_values ( n_data, n, fn )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i12)' ) n, fn

  end do

  return
end
subroutine thercon_values_test ( )

!*****************************************************************************80
!
!! thercon_values_test() tests thercon_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 March 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) lambda
  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'THERCON_VALUES_TEST:'
  write ( *, '(a)' ) '  THERCON_VALUES returns values of '
  write ( *, '(a)' ) '  the thermal conductivity of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      T            P            LAMBDA(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call thercon_values ( n_data, tc, p, lambda )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, lambda

  end do

  return
end
subroutine three_j_values_test ( )

!*****************************************************************************80
!
!! three_j_values_test() tests three_j_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  real ( kind = rk8 ) j1
  real ( kind = rk8 ) j2
  real ( kind = rk8 ) j3
  real ( kind = rk8 ) m1
  real ( kind = rk8 ) m2
  real ( kind = rk8 ) m3
  integer n_data

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'THREE_J_VALUES_TEST:'
  write ( *, '(a)' ) '  THREE_J_VALUES returns values of '
  write ( *, '(a)' ) '  the Wigner 3J coefficient.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      J1      J2      J3      M1      M2      M3        THREE_J'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' ) &
      j1, j2, j3, m1, m2, m3, fx

  end do

  return
end
subroutine tran02_values_test ( )

!*****************************************************************************80
!
!! tran02_values_test() tests tran02_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN02_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN02_VALUES returns values of '
  write ( *, '(a)' ) '  the order 2 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T2(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran02_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran03_values_test ( )

!*****************************************************************************80
!
!! tran03_values_test() tests tran03_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN03_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN03_VALUES returns values of '
  write ( *, '(a)' ) '  the order 3 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T3(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran03_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran04_values_test ( )

!*****************************************************************************80
!
!! tran04_values_test() tests tran04_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN04_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN04_VALUES returns values of '
  write ( *, '(a)' ) '  the order 4 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T4(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran04_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran05_values_test ( )

!*****************************************************************************80
!
!! tran05_values_test() tests tran05_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN05_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN05_VALUES returns values of '
  write ( *, '(a)' ) '  the order 5 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T5(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran05_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran06_values_test ( )

!*****************************************************************************80
!
!! tran06_values_test() tests tran06_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN06_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN06_VALUES returns values of '
  write ( *, '(a)' ) '  the order 6 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T6(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran06_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran07_values_test ( )

!*****************************************************************************80
!
!! tran07_values_test() tests tran07_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN07_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN07_VALUES returns values of '
  write ( *, '(a)' ) '  the order 7 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T7(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran07_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran08_values_test ( )

!*****************************************************************************80
!
!! tran08_values_test() tests tran08_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN08_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN08_VALUES returns values of '
  write ( *, '(a)' ) '  the order 8 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T8(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran08_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine tran09_values_test ( )

!*****************************************************************************80
!
!! tran09_values_test() tests tran09_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRAN09_VALUES_TEST:'
  write ( *, '(a)' ) '  TRAN09_VALUES returns values of '
  write ( *, '(a)' ) '  the order 9 transportation function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            T9(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tran09_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) x, fx

  end do

  return
end
subroutine trigamma_values_test ( )

!*****************************************************************************80
!
!! trigamma_values_test() tests trigamma_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRIGAMMA_VALUES_TEST:'
  write ( *, '(a)' ) '  TRIGAMMA_VALUES returns values of '
  write ( *, '(a)' ) '  the TriGamma function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X            F(X)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call trigamma_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g24.16)' ) x, fx

  end do

  return
end
subroutine truncated_normal_ab_cdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_ab_cdf_values_test() tests truncated_normal_ab_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_AB_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the truncated Normal CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       A         B         X        CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_ab_cdf_values ( n_data, mu, sigma, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, a, b, x, fx

  end do

  return
end
subroutine truncated_normal_ab_pdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_ab_pdf_values_test() tests truncated_normal_ab_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_AB_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_AB_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the truncated Normal PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       A         B         X        PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_ab_pdf_values ( n_data, mu, sigma, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, a, b, x, fx

  end do

  return
end
subroutine truncated_normal_a_cdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_a_cdf_values_test() tests truncated_normal_a_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_A_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_A_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the lower truncated Normal CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       A         X        CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_a_cdf_values ( n_data, mu, sigma, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, a, x, fx

  end do

  return
end
subroutine truncated_normal_a_pdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_a_pdf_values_test() tests truncated_normal_a_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_A_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_A_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the lower Truncated Normal PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       A         X        PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_a_pdf_values ( n_data, mu, sigma, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, a, x, fx

  end do

  return
end
subroutine truncated_normal_b_cdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_b_cdf_values_test() tests truncated_normal_b_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_B_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the upper truncated Normal CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       B         X        CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_b_cdf_values ( n_data, mu, sigma, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, b, x, fx

  end do

  return
end
subroutine truncated_normal_b_pdf_values_test ( )

!*****************************************************************************80
!
!! truncated_normal_b_pdf_values_test() tests truncated_normal_b_pdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  real ( kind = rk8 ) mu
  integer n_data
  real ( kind = rk8 ) sigma
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TRUNCATED_NORMAL_B_PDF_VALUES_TEST:'
  write ( *, '(a)' ) '  TRUNCATED_NORMAL_B_PDF_VALUES returns values of '
  write ( *, '(a)' ) '  the upper truncated Normal PDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '      MU       SIGMA       B         X        PDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call truncated_normal_b_pdf_values ( n_data, mu, sigma, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.1,2x,f8.1,2x,f8.1,2x,f8.1,2x,g24.16)' ) &
      mu, sigma, b, x, fx

  end do

  return
end
subroutine tsat_values_test ( )

!*****************************************************************************80
!
!! tsat_values_test() tests tsat_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TSAT_VALUES_TEST:'
  write ( *, '(a)' ) '  TSAT_VALUES returns values of '
  write ( *, '(a)' ) '  the saturation temperature '
  write ( *, '(a)' ) '  as a function of pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      P           Tsat(P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call tsat_values ( n_data, p, tc )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6)' ) p, tc

  end do

  return
end
subroutine van_der_corput_values_test ( )

!*****************************************************************************80
!
!! van_der_corput_values_test() tests van_der_corput_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer base
  integer n_data
  integer seed
  real ( kind = rk8 ) value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VAN_DER_CORPUT_VALUES_TEST:'
  write ( *, '(a)' ) '  VAN_DER_CORPUT_VALUES returns values of '
  write ( *, '(a)' ) '  the van der Corput sequence in a given base.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      BASE      SEED    VDC(BASE,SEED)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call van_der_corput_values ( n_data, base, seed, value )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,i8,2x,g16.8)' ) base, seed, value

  end do

  return
end
subroutine viscosity_values_test ( )

!*****************************************************************************80
!
!! viscosity_values_test() tests viscosity_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) eta
  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) tc

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VISCOSITY_VALUES_TEST:'
  write ( *, '(a)' ) '  VISCOSITY_VALUES returns values of '
  write ( *, '(a)' ) '  the viscosity of water '
  write ( *, '(a)' ) '  as a function of temperature and pressure.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          T               P          ETA(T,P)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call viscosity_values ( n_data, tc, p, eta )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) tc, p, eta

  end do

  return
end
subroutine von_mises_cdf_values_test ( )

!*****************************************************************************80
!
!! von_mises_cdf_values_test() tests von_mises_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'VON_MISES_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  VON_MISES_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the von Mises CDF.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '          A               B               X           CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call von_mises_cdf_values ( n_data, a, b, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) a, b, x, fx

  end do

  return
end
subroutine weekday_values_test ( )

!*****************************************************************************80
!
!! weekday_values_test() tests weekday_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 January 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer d
  integer m
  integer n_data
  integer w
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEEKDAY_VALUES_TEST:'
  write ( *, '(a)' ) '  WEEKDAY_VALUES returns values of '
  write ( *, '(a)' ) '  the weekday for a given Y/M/D date.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     Y     M     D     W'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call weekday_values ( n_data, y, m, d, w )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) y, m, d, w

  end do

  return
end
subroutine weibull_cdf_values_test ( )

!*****************************************************************************80
!
!! weibull_cdf_values_test() tests weibull_cdf_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) alpha
  real ( kind = rk8 ) beta
  real ( kind = rk8 ) fx
  integer n_data
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'WEIBULL_CDF_VALUES_TEST:'
  write ( *, '(a)' ) '  WEIBULL_CDF_VALUES returns values of '
  write ( *, '(a)' ) '  the Weibull Cumulative Density Function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       ALPHA          BETA             X              CDF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call weibull_cdf_values ( n_data, alpha, beta, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g24.16)' ) alpha, beta, x, fx

  end do

  return
end
subroutine wright_omega_values_test ( )

!*****************************************************************************80
!
!! wright_omega_values_test() tests wright_omega_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) fz
  integer n_data
  complex ( kind = ck ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'wright_omega_values_test():'
  write ( *, '(a)' ) '  wright_omega_values() returns values of '
  write ( *, '(a)' ) '  the Wright Omega function with complex argument'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          Z                         FZ'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call wright_omega_values ( n_data, z, fz )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,(g14.6,g14.6),2x,(g24.16,g24.16))' ) z, fz

  end do

  return
end
subroutine zeta_values_test ( )

!*****************************************************************************80
!
!! zeta_values_test() tests zeta_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 February 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n
  integer n_data
  real ( kind = rk8 ) zeta

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zeta_values_test():'
  write ( *, '(a)' ) '  zeta_values() returns values of the Riemann Zeta function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       ZETA(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call zeta_values ( n_data, n, zeta )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,i8,2x,g24.16)' ) n, zeta

  end do

  return
end
subroutine zeta_m1_values_test ( )

!*****************************************************************************80
!
!! zeta_m1_values_test() tests zeta_m1_values().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer n_data
  real ( kind = rk8 ) p
  real ( kind = rk8 ) zeta_m1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ZETA_M1_VALUES_TEST:'
  write ( *, '(a)' ) '  ZETA_M1_VALUES returns values of '
  write ( *, '(a)' ) '  the Riemann Zeta Minus One function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N       ZETA_M1(N)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call zeta_m1_values ( n_data, p, zeta_m1 )

    if ( n_data == 0 ) then
      exit
    end if

    write ( *, '(2x,f8.2,2x,g24.16)' ) p, zeta_m1

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

  integer, parameter :: rk8 = kind ( 1.0D+00 )

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
