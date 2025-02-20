subroutine problem_size ( chain_num, cr_num, gen_num, pair_num, par_num )

!*****************************************************************************80
!
!! PROBLEM_SIZE sets information having to do with dimensions.
!
!  Discussion:
!
!    Although this problem involves multivariate normal distributions for
!    the initial sample and the density function, in both cases the covariance
!    matrix is diagonal, so we don't have to work very hard to set up the
!    calculations involving the distributions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jasper Vrugt, CJF ter Braak, CGH Diks, Bruce Robinson, James Hyman, 
!    Dave Higdon,
!    Accelerating Markov Chain Monte Carlo Simulation by Differential 
!    Evolution with Self-Adaptive Randomized Subspace Sampling,
!    International Journal of Nonlinear Sciences and Numerical Simulation,
!    Volume 10, Number 3, March 2009, pages 271-288.
!
!  Parameters:
!
!    Output, integer CHAIN_NUM, the total number of chains.
!    3 <= CHAIN_NUM.
!
!    Output, integer CR_NUM, the total number of CR values.
!    1 <= CR_NUM.
!
!    Output, integer GEN_NUM, the total number of generations.
!    2 <= GEN_NUM.
!
!    Output, integer PAIR_NUM, the number of pairs of 
!    crossover chains.
!    0 <= PAIR_NUM.
!
!    Output, integer PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
  implicit none

  integer chain_num
  integer cr_num
  integer gen_num
  integer pair_num
  integer par_num

  chain_num = 10
  cr_num = 3
  gen_num = 10
  pair_num = 3
  par_num = 10

  return
end
subroutine problem_value ( chain_filename, gr_filename, gr_threshold, &
  jumpstep, limits, par_num, printstep, restart_read_filename, &
  restart_write_filename  )

!*****************************************************************************80
!
!! PROBLEM_VALUE sets information, including numeric data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, character ( len = * ) CHAIN_FILENAME, the "base" filename
!    to be used for the chain files.  If this is the empty string '',
!    then the chain files will not be written.  This name should 
!    include a string of 0's which will be replaced by the chain 
!    indices.  For example, "chain000.txt" would work as long as the
!    number of chains was 1000 or less.
!
!    Output, character ( len = * ) GR_FILENAME, the name of the file
!    in which values of the Gelman-Rubin statistic will be recorded,
!    or '' if this file is not to be written.
!
!    Output, real ( kind = rk ) GR_THRESHOLD, the convergence tolerance for
!    the Gelman-Rubin statistic.
!
!    Output, integer JUMPSTEP, forces a "long jump" every
!    JUMPSTEP generations.
!
!    Output, real ( kind = rk ) LIMITS(2,PAR_NUM), lower and upper bounds
!    for each parameter.
!
!    Input, integer PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
!    Output, integer PRINTSTEP, the interval between generations 
!    on which the Gelman-Rubin statistic will be computed and written to a file.
!
!    Output, character ( len = * ) RESTART_READ_FILENAME, the name of the file
!    containing restart information, or '' if this is not a restart run.
!
!    Output, character ( len = * ) RESTART_WRITE_FILENAME, the name of the file
!    to be written, containing restart information, or '' if a restart file 
!    is not to be written.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer par_num

  character ( len = * ) chain_filename
  character ( len = * ) gr_filename
  real ( kind = rk ) gr_threshold
  integer jumpstep
  real ( kind = rk ) limits(2,par_num)
  integer printstep
  character ( len = * ) restart_read_filename
  character ( len = * ) restart_write_filename

  chain_filename = 'problem2_chain00.txt'
  gr_filename = 'problem2_gr.txt'
  gr_threshold = 1.2D+00
  jumpstep = 5
  limits(1,1:par_num) = -10.0D+00
  limits(2,1:par_num) = +10.0D+00
  printstep = 10
  restart_read_filename = ''
  restart_write_filename = 'problem2_restart.txt'

  return
end
function prior_density ( par_num, zp )

!*****************************************************************************80
!
!! PRIOR_DENSITY evaluates the prior density function.
!
!  Discussion:
!
!    "The initial sample was generated from a normal distribution with
!    variance-covariance matrix 5 * Id."
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2013
!
!  Author:
!
!    John Burkardt.
!
!  Parameters:
!
!    Input, integer PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
!    Input, real ( kind = rk ) ZP(PAR_NUM), the argument of the density
!    function.
!
!    Output, real ( kind = rk ) PRIOR_DENSITY, the value of the prior
!    density function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer par_num

  integer i
  real ( kind = rk ) prior_density
  real ( kind = rk ) r8_normal_pdf
  real ( kind = rk ) zp(par_num)
  real ( kind = rk ), parameter :: zp_mean_1d = 0.0D+00
  real ( kind = rk ), parameter :: zp_sd_1d = sqrt ( 5.0D+00 )

  prior_density = 1.0D+00
  do i = 1, par_num
    prior_density = prior_density &
      * r8_normal_pdf ( zp_mean_1d, zp_sd_1d, zp(i) )
  end do

  return
end
subroutine prior_sample ( par_num, zp )

!*****************************************************************************80
!
!! PRIOR_SAMPLE samples from the prior distribution.
!
!  Discussion:
!
!    "The initial sample was generated from a normal distribution with
!    variance-covariance matrix 5 * Id."
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
!    Output, real ( kind = rk ) ZP(PAR_NUM), the sample from the distribution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer par_num

  integer i
  real ( kind = rk ) r8_normal_sample
  real ( kind = rk ) zp(par_num)
  real ( kind = rk ), parameter :: zp_mean_1d = 0.0D+00
  real ( kind = rk ), parameter :: zp_sd_1d = sqrt ( 5.0D+00 )

  do i = 1, par_num
    zp(i) = r8_normal_sample ( zp_mean_1d, zp_sd_1d )
  end do

  return
end
function sample_likelihood ( par_num, zp )

!*****************************************************************************80
!
!! SAMPLE_LIKELIHOOD computes the log likelihood function.
!
!  Discussion:
!
!    "A 10-dimensional twisted Gaussian density function, given by the
!    unnormalized density
!      pi_b(x) proportional to pi(phi_b(x)),
!    with
!      phi_b(x) = (x1, x2+bx1^2-100b,x3,...,x10).
!    Here, pi signifies the density of a multivariate normal distribution,
!    Nd(0,Sigma), with Sigma=diag(100,1,...,1), and phi_b is a function that
!    is used to transform pi to a twisted distribution."
!
!    The value b = 0.01 corresponds to a mildly nonlinear problem,
!    and b = 0.1 to a highly nonlinear problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PAR_NUM, the total number of parameters.
!    1 <= PAR_NUM.
!
!    Input, real ( kind = rk ) ZP(PAR_NUM), a sample.
!
!    Output, real ( kind = rk ) SAMPLE_LIKELIHOOD, the log likelihood function 
!    for the sample.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer par_num

  real ( kind = rk ), parameter :: b = 0.01D+00
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) sample_likelihood
  real ( kind = rk ) xcx
  real ( kind = rk ) y2
  real ( kind = rk ) zp(par_num)

  y2 = zp(2) + b * zp(1) ** 2 - 100.0D+00 * b

  xcx = zp(1) ** 2 / 100.0D+00 + y2 ** 2
  do i = 3, par_num
    xcx = xcx + zp(i) ** 2
  end do

  sample_likelihood = &
    - 0.5D+00 * real ( par_num, kind = rk ) * log ( 2.0D+00 * r8_pi ) &
    - 0.5D+00 * log ( 100.0D+00 ) &
    - 0.5D+00 * xcx

  return
end

