program main

!*****************************************************************************80
!
!! bdmlib_test() tests bdmlib().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!   02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: acid_num = 20
  integer, parameter :: comp_max = 9
  integer, parameter :: iunit = 1

  character acid_sym(acid_num)
  real ( kind = rk ) beta(acid_num,comp_max)
  real ( kind = rk ) beta_sum(comp_max)
  integer comp_label(comp_max)
  integer comp_num
  real ( kind = rk ) comp_weight(comp_max)
  integer ierror
  character ( len = 30 ) mixture_file_name

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bdmlib_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test bdmlib().'
!
!  Read information about the mixture.
!
  mixture_file_name = 'mixture.txt'

  open ( unit = iunit, file = mixture_file_name, form = 'formatted' )

  call mixture_read ( acid_num, acid_sym, beta, beta_sum, &
    comp_label, comp_max, comp_num, comp_weight, ierror, iunit )

  close ( unit = iunit )
!
!  Print the amino acid parameters.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Numeric key for amino acid abbreviations.'
  write ( *, '(a)' ) ' '

  call amino_print ( acid_num, acid_sym )
!
!  Print the component parameters.
!
  call comp_param_print ( acid_num, acid_sym, comp_max, comp_num, &
    beta, beta_sum, comp_weight )
!
!  Check the parameters.
!
  call dirichlet_mix_check ( comp_num, acid_num, acid_num, beta, comp_weight )
!
!  Perform the simple test of generating 10 isoleucines in a row.
!
  call test01 ( acid_num, beta, comp_num, comp_weight )
!
!  Now test a random sampling.
!
  call test02 ( acid_num, beta, comp_num, comp_weight )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'bdmlib_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( acid_num, beta, comp_num, comp_weight )

!*****************************************************************************80
!
!! test01() generates 10 isoleucine events in a row.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer acid_num
  integer comp_num

  real ( kind = rk ) alpha(comp_num)
  real ( kind = rk ) alpha_0(comp_num)
  real ( kind = rk ) beta(acid_num,comp_num)
  real ( kind = rk ) comp_weight(comp_num)
  real ( kind = rk ) comp_weight_est(comp_num)
  integer event_i
  integer event_num
  real ( kind = rk ) p(comp_num)
  real ( kind = rk ) p_hat(comp_num)
  integer site_num
  integer x_sample(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Generate a (nonrandom) sequence of'
  write ( *, '(a)' ) '  10 isoleucine results in a row.'
  write ( *, '(a)' ) ' '
!
!  Initialize information about the Bayesian process.
!
  alpha_0(1:comp_num) = 1.0D+00

  alpha(1:comp_num) = alpha_0(1:comp_num)

  p(1:comp_num) = 1.0D+00 / real ( comp_num, kind = rk )
  p_hat(1:comp_num) = 1.0D+00 / real ( comp_num, kind = rk )

  event_num = 10

  call r8vec_print ( comp_num, comp_weight, 'Exact component weights:' )

  call r8vec_print ( comp_num, alpha, 'Initial ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
  call dirichlet_mean ( comp_num, alpha, comp_weight_est )

  call r8vec_print ( comp_num, comp_weight_est, &
    'Initial estimated component weights:' )

  site_num = 1

  do event_i = 1, event_num
!
!  Observe a single isoleucine.
!
    x_sample(1:acid_num) = 0
    x_sample(8) = site_num
!
!  Update ALPHA, the estimated weight parameters, based on X.
!
    call event_process ( acid_num, alpha, beta, comp_num, p, p_hat, site_num, &
      x_sample )

    call r8vec_print ( comp_num, alpha, 'Current ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
    call dirichlet_mean ( comp_num, alpha, comp_weight_est )

    call r8vec_print ( comp_num, comp_weight_est, &
      'Estimated component weights:' )

  end do

  return
end
subroutine test02 ( acid_num, beta, comp_num, comp_weight )

!*****************************************************************************80
!
!! test02() generates random events.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer acid_num
  integer comp_num

  real ( kind = rk ) alpha(comp_num)
  real ( kind = rk ) alpha_0(comp_num)
  real ( kind = rk ) beta(acid_num,comp_num)
  integer comp_sample
  real ( kind = rk ) comp_weight(comp_num)
  real ( kind = rk ) comp_weight_est(comp_num)
  integer event_i
  integer event_num
  real ( kind = rk ) p(comp_num)
  real ( kind = rk ) p_hat(comp_num)
  real ( kind = rk ) p_sample(acid_num)
  integer site_num
  integer x_sample(acid_num)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  Generate many random events.'
  write ( *, '(a)' ) '  We should be able to approximate the'
  write ( *, '(a)' ) '  exact component weights.'
  write ( *, '(a)' ) ' '
!
!  Initialize information about the Bayesian process.
!
  alpha_0(1:comp_num) = 1.0D+00

  alpha(1:comp_num) = alpha_0(1:comp_num)

  p(1:comp_num) = 1.0D+00 / real ( comp_num, kind = rk )
  p_hat(1:comp_num) = 1.0D+00 / real ( comp_num, kind = rk )

  event_num = 1000

  call r8vec_print ( comp_num, comp_weight, 'Exact component weights:' )

  call r8vec_print ( comp_num, alpha, 'Initial ALPHA:' )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
  call dirichlet_mean ( comp_num, alpha, comp_weight_est )

  call r8vec_print ( comp_num, comp_weight_est, &
    'Initial estimated component weights:' )

  site_num = 10

  do event_i = 1, event_num
!
!  Randomly choose COMP_SAMPLE, the component PDF to sample.
!
    call discrete_sample ( comp_num, comp_weight, comp_sample )
!
!  Now generate the probabilities P_SAMPLE for a multinomial PDF by sampling
!  the COMP_SAMPLE-th Dirichlet distribution.
!
    call dirichlet_sample ( acid_num, beta(1,comp_sample), p_sample )
!
!  Now generate the event X_SAMPLE by sampling the multinomial PDF with
!  the given P_SAMPLE parameters.
!
    call multinomial_sample ( site_num, acid_num, p_sample, x_sample )
!
!  Update ALPHA, the estimated weight parameters, based on X.
!
    call event_process ( acid_num, alpha, beta, comp_num, p, p_hat, site_num, &
      x_sample )
!
!  Based on the current ALPHA's, compute the mean/expected value/estimate 
!  for the weights.
!
    call dirichlet_mean ( comp_num, alpha, comp_weight_est )

    if ( event_i <= 10 .or. mod ( event_i, event_num/10 ) == 0 ) then

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 'Event ', event_i

      call r8vec_print ( comp_num, alpha, 'Current ALPHA:' )

      call r8vec_print ( comp_num, comp_weight_est, &
        'Estimated component weights:' )

    end if

  end do

  return
end
 
