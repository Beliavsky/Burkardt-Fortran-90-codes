program main

!*****************************************************************************80
!
!! asa189_test() tests asa189().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a_est
  real ( kind = rk ) b
  real ( kind = rk ) b_est
  integer c
  real ( kind = rk ) mu
  real ( kind = rk ) mu_est
  integer sample_log
  integer sample_num
  real ( kind = rk ) theta
  real ( kind = rk ) theta_est

  call timestamp ( )

  a = 2.0D+00
  b = 3.0D+00
  c = 4

  call beta_binomial_check ( a, b, c )

  mu = a / ( a + b )
  theta = 1.0D+00 / ( a + b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa189_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa189().'
  write ( *, '(a)' ) '  Applied Statistics Algorithm 189'
  write ( *, '(a)' ) '  Estimate the parameters A and B of a '
  write ( *, '(a)' ) '  beta binomial distribution.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Exact values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '          A             B                MU         THETA'
  write ( *, '(a)' ) ' '
  write ( *, '(6x,4g14.6)' ) a, b, mu, theta
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Estimated values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Samples   A_est         B_est        MU_est     THETA_est'
  write ( *, '(a)' ) ' '

  do sample_log = 2, 5

    sample_num = 10**sample_log

    call analyze ( sample_num, a, b, c, mu_est, theta_est )
!
!  Convert the ASA189 "THETA" and "MU" parameters to "A" and "B".
!
    a_est = mu_est / theta_est
    b_est = ( 1.0D+00 - mu_est ) / theta_est

    write ( *, '(i6,4g14.6)' ) sample_num, a_est, b_est, mu_est, theta_est

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa189_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine analyze ( sample_num, a, b, c, mu_est, theta_est )

!*****************************************************************************80
!
!! analyze() generates data and analyzes it with ASA189.
!
!  Discussion:
!
!    The routine to generate the samples uses parameter A, B and C,
!    while ASA189 prefers a related form MU, THETA.  The calling routine
!    has to figure out how these are related.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 January 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer SAMPLE_NUM, the number of samples to generate.
!
!    Input, real ( kind = rk ) A, real ( kind = rk ) B, integer C,
!    the parameters to use in the beta binomial distribution.
!
!    Output, real  ( kind = rk ) MU_EST, THETA_EST, the estimates of MU and THETA
!    produced by ASA189.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mrl = 4
  integer sample_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer c
  real ( kind = rk ) ccrit
  integer ifault
  integer in(sample_num)
  integer iter
  real ( kind = rk ) mu_est
  real ( kind = rk ) mu_se
  integer rl(mrl,3)
  real ( kind = rk ) rnl
  integer sample_i
  real ( kind = rk ) theta_est
  real ( kind = rk ) theta_se
  integer x(sample_num)
!
!  Generate the sample data using the exact parameters A, B.
!
  do sample_i = 1, sample_num
    call beta_binomial_sample ( a, b, c, x(sample_i) )
  end do
!
!  Analyze the sample data, trying to estimate the parameters
!  in the form "MU", "THETA".  Note that ASA189 expects to be told
!  the value of C (although C could vary from one sample to the next).
!
  in(1:sample_num) = c

  iter = 10
  ccrit = 0.001D+00

  call bbml ( sample_num, x, in, rl, mrl, iter, ccrit, mu_est, &
    theta_est, mu_se, theta_se, rnl, ifault )

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

