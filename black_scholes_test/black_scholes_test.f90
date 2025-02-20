program main

!*****************************************************************************80
!
!! black_scholes_test() tests black_scholes().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( );
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'black_scholes_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test black_scholes().'

  call asset_path_test ( )
  call binomial_test ( )
  call bsf_test ( )
  call forward_test ( )
  call mc_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'black_scholes_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine asset_path_test ( )

!*****************************************************************************80
!
!! asset_path_test() tests asset_path().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 100

  character ( len = 255 ) command_filename
  integer command_unit
  character ( len = 255 ) data_filename
  integer data_unit
  integer i
  real ( kind = rk ) mu
  real ( kind = rk ) s(0:n)
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asset_path_test():'
  write ( *, '(a)' ) '  asset_path() demonstrates the simulated of an asset price path.'

  s0 = 2.0D+00
  mu = 0.1D+00
  sigma = 0.3D+00
  t1 = 1.0D+00

  write ( *, '(a,g14.6)' ) ' '
  write ( *, '(a,g14.6)' ) '  The asset price at time 0      S0    = ', s0
  write ( *, '(a,g14.6)' ) '  The asset expected growth rate MU    = ', mu
  write ( *, '(a,g14.6)' ) '  The asset volatility           SIGMA = ', sigma
  write ( *, '(a,g14.6)' ) '  The expiry date                T1    = ', t1
  write ( *, '(a,i6)' )    '  The number of time steps       N     = ', n

  call asset_path ( s0, mu, sigma, t1, n, s )

  call r8vec_print_part ( n + 1, s, 10, '  Partial results:' )

  data_filename = 'asset_path_data.txt'
  call get_unit ( data_unit )
  open ( unit = data_unit, file = data_filename, status = 'replace' )
  do i = 0, n
    write ( data_unit, '(i6,2x,g14.6)' ) i, s(i)
  end do
  close ( unit = data_unit )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Created graphics data file "' // trim ( data_filename ) // '".'

  command_filename = 'asset_path_commands.txt'
  call get_unit ( command_unit )
  open ( unit = command_unit, file = command_filename, status = 'replace' )
  write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) '# Usage:'
  write ( command_unit, '(a)' ) '#  gnuplot < ' // trim ( command_filename )
  write ( command_unit, '(a)' ) '#'
  write ( command_unit, '(a)' ) 'set term png'
  write ( command_unit, '(a)' ) 'set output "asset_path.png"'
  write ( command_unit, '(a)' ) 'set xlabel "Time"'
  write ( command_unit, '(a)' ) 'set ylabel "Value"'
  write ( command_unit, '(a)' ) 'set title "Asset value over time"'
  write ( command_unit, '(a)' ) 'set grid'
  write ( command_unit, '(a)' ) 'set style data lines'
  write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) // &
    '" using 1:2 lw 3 linecolor rgb "blue"'
  write ( command_unit, '(a)' ) 'quit'
  close ( unit = command_unit )
  write ( *, '(a)' ) &
    '  Created command file "' // trim ( command_filename ) // '".'

  return
end
subroutine binomial_test ( )

!*****************************************************************************80
!
!! binomial_test() tests binomial().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) e
  integer m
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'binomial_test():'
  write ( *, '(a)' ) '  binomial() demonstrates the binomial method'
  write ( *, '(a)' ) '  for option valuation.'

  s0 = 2.0D+00
  e = 1.0D+00
  r = 0.05D+00
  sigma = 0.25D+00
  t1 = 3.0D+00
  m = 256

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The asset price at time 0 S0    = ', s0
  write ( *, '(a,g14.6)' ) '  The exercise price        E     = ', e
  write ( *, '(a,g14.6)' ) '  The interest rate         R     = ', r
  write ( *, '(a,g14.6)' ) '  The asset volatility      SIGMA = ', sigma
  write ( *, '(a,g14.6)' ) '  The expiry date           T1    = ', t1
  write ( *, '(a,i8)' )    '  The number of intervals   M     = ', m

  call binomial ( s0, e, r, sigma, t1, m, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The option value is ', c

  return
end
subroutine bsf_test ( )

!*****************************************************************************80
!
!! bsf_test() tests bsf().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ) e
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t0
  real ( kind = rk ) t1

  write ( *, '(a)' )
  write ( *, '(a)' ) 'bsf_test():'
  write ( *, '(a)' ) '  bsf() demonstrates the Black-Scholes formula'
  write ( *, '(a)' ) '  for option valuation.'

  s0 = 2.0D+00
  t0 = 0.0D+00
  e = 1.0D+00
  r = 0.05D+00
  sigma = 0.25D+00
  t1 = 3.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The asset price at time T0 S0    = ', s0
  write ( *, '(a,g14.6)' ) '  The time                   T0    = ', t0
  write ( *, '(a,g14.6)' ) '  The exercise price         E     = ', e
  write ( *, '(a,g14.6)' ) '  The interest rate          R     = ', r
  write ( *, '(a,g14.6)' ) '  The asset volatility       SIGMA = ', sigma
  write ( *, '(a,g14.6)' ) '  The expiry date            T1    = ', t1

  call bsf ( s0, t0, e, r, sigma, t1, c )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The option value C = ', c

  return
end
subroutine forward_test ( )

!*****************************************************************************80
!
!! forward_test() tests forward().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e
  integer i
  integer nt
  integer nx
  real ( kind = rk ) r
  real ( kind = rk ) s
  real ( kind = rk ) sigma
  real ( kind = rk ) smax
  real ( kind = rk ) smin
  real ( kind = rk ) t1
  real ( kind = rk ), allocatable :: u(:,:)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'forward_test():'
  write ( *, '(a)' ) '  forward() demonstrates the forward difference method'
  write ( *, '(a)' ) '  for option valuation.'

  e = 4.0D+00
  r = 0.03D+00
  sigma = 0.50D+00
  t1 = 1.0D+00
  nx = 11
  nt = 29
  smax = 10.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  The exercise price        E =     ', e
  write ( *, '(a,g14.6)' ) '  The interest rate         R =     ', r
  write ( *, '(a,g14.6)' ) '  The asset volatility      SIGMA = ', sigma;
  write ( *, '(a,g14.6)' ) '  The expiry date           T1 =    ', t1
  write ( *, '(a,i8)' ) '  The number of space steps NX =    ', nx
  write ( *, '(a,i8)' ) '  The number of time steps  NT =    ', nt
  write ( *, '(a,g14.6)' ) '  The value of              SMAX =  ', smax

  allocate ( u(nx-1,nt+1) )

  call forward ( e, r, sigma, t1, nx, nt, smax, u )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Initial   Option'
  write ( *, '(a)' ) '  Value     Value' 
  write ( *, '(a)' ) ' '

  smin = 0.0D+00
  do i = 1, nx - 1 
    s = ( ( nx - i - 1 ) * smin + i * smax ) / real ( nx - 1, kind = rk )
    write ( *, '(2x,g14.6,2x,g14.6)' ) s, u(i,nt+1)
  end do

  deallocate ( u )

  return
end
subroutine mc_test ( )

!*****************************************************************************80
!
!! mc_test() tests mc().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) conf(2)
  real ( kind = rk ) e
  integer m
  real ( kind = rk ) r
  real ( kind = rk ) s0
  real ( kind = rk ) sigma
  real ( kind = rk ) t1

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'mc_test():'
  write ( *, '(a)' ) '  mc() demonstrates the Monte Carlo method'
  write ( *, '(a)' ) '  for option valuation.'

  s0 = 2.0D+00
  e = 1.0D+00
  r = 0.05D+00
  sigma = 0.25D+00
  t1 = 3.0D+00
  m = 1000000

  write ( *, '(a)' ) ' '
  write ( *, '(a, g14.6)' ) '  The asset price at time 0, S0    = ', s0
  write ( *, '(a, g14.6)' ) '  The exercise price         E     = ', e
  write ( *, '(a, g14.6)' ) '  The interest rate          R     = ', r
  write ( *, '(a, g14.6)' ) '  The asset volatility       SIGMA = ', sigma
  write ( *, '(a, g14.6)' ) '  The expiry date            T1    = ', t1
  write ( *, '(a, i8)' )    '  The number of simulations  M     = ', m

  call mc ( s0, e, r, sigma, t1, m, conf )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6,a,g14.6,a)' ) &
    '  The confidence interval is [', conf(1), ',', conf(2), '].'

  return
end

