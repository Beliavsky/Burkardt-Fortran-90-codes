program main

!*****************************************************************************80
!
!! fftw_test() tests fftw().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fftw_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test fftw().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'fftw_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01: complex 1D data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: ik = selected_int_kind ( 15 ) 
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 100

  include "fftw3.f90"

  integer i
  complex ( kind = ck ) in(n)
  complex ( kind = ck ) in2(n)
  complex ( kind = ck ) out(n)
  integer ( kind = ik ) plan_backward
  integer ( kind = ik ) plan_forward
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Demonstrate FFTW on a single vector '
  write ( *, '(a)' ) '  of complex data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
  write ( *, '(a)' ) '  the data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'
!
!  Compute the data, a complex vector of length N.
!
  call c8vec_uniform_01 ( n, seed, in )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,2g14.6)' ) i, in(i)
  end do
!
!  Make a plan for the FFT, and forward transform the data.
!
  call dfftw_plan_dft_1d_ ( plan_forward, n, in, out, FFTW_FORWARD, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_forward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
  end do
!
!  Make a plan for the backward FFT, and recover the original data.
!
  call dfftw_plan_dft_1d_ ( plan_backward, n, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_backward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divided by N:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,2g14.6)' ) i, in2(i) / real ( n, kind = rk )
  end do
!
!  Discard the information associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02: real 1D data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: ik = selected_int_kind ( 15 ) 
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 100
  integer, parameter :: nc = 51

  include "fftw3.f90"

  integer i
  real ( kind = rk ) in(n)
  real ( kind = rk ) in2(n)
  complex ( kind = ck ) out(nc)
  integer ( kind = ik ) plan_backward
  integer ( kind = ik ) plan_forward
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Demonstrate FFTW on a single vector'
  write ( *, '(a)' ) '  of real data.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover '
  write ( *, '(a)' ) '  data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'
!
!  Set up the input data, a real vector of length N.
!
  call r8vec_uniform_01 ( n, seed, in )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, in(i)
  end do
!
!  Set up a plan, and execute the plan to transform the IN data to
!  the OUT FFT coefficients.
!
  call dfftw_plan_dft_r2c_1d_ ( plan_forward, n, in, out, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_forward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, nc
    write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
  end do
!
!  Set up a plan, and execute the plan to backtransform the
!  complex FFT coefficients in OUT to real data.
!
  call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, out, in2, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_backward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divide by N:'
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i4,2x,g14.6)' ) i, in2(i) / real ( n, kind = rk )
  end do
!
!  Release the memory associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03: complex 2D data.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: ik = selected_int_kind ( 15 ) 
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 8
  integer, parameter :: ny = 10

  include "fftw3.f90"

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) r8_uniform_01
  integer i
  complex ( kind = ck ) in(nx,ny)
  complex ( kind = ck ) in2(nx,ny)
  integer j
  complex ( kind = ck ) out(nx,ny)
  integer ( kind = ik ) plan_backward
  integer ( kind = ik ) plan_forward
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Demonstrate FFTW on a 2D complex array'
  write ( *, '(a,i8)' ) '  NX = ', nx
  write ( *, '(a,i8)' ) '  NY = ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
  write ( *, '(a)' ) '  the data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'
!
!  Compute the data.
!
  do j = 1, ny
    do i = 1, nx
      a = r8_uniform_01 ( seed )
      b = r8_uniform_01 ( seed )
      in(i,j) = cmplx ( a, b, kind = ck )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, in(i,j)
    end do
  end do
!
!  Make a plan for the FFT, and forward transform the data.
!
  call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, in, out, FFTW_FORWARD, &
    FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_forward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
    end do
  end do
!
!  Make a plan for the backward FFT, and recover the original data.
!
  call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, out, in2, FFTW_BACKWARD, &
    FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_backward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) &
        i, j, in2(i,j) / real ( nx * ny, kind = rk )
    end do
  end do
!
!  Discard the information associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04: real 2D data.
!
!  Discussion:
!
!    In contrast to the C example, in FORTRAN it is the FIRST dimension
!    of the complex coefficient array that is "half" the size.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    30 July 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: ik = selected_int_kind ( 15 ) 
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nx = 8
  integer, parameter :: ny = 10
  integer, parameter :: nh = ( nx / 2 ) + 1

  include "fftw3.f90"

  real ( kind = rk ) r8_uniform_01
  integer i
  real ( kind = rk ) in(nx,ny)
  real ( kind = rk ) in2(nx,ny)
  integer j
  complex ( kind = ck ) out(nh,ny)
  integer ( kind = ik ) plan_backward
  integer ( kind = ik ) plan_forward
  integer seed

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Demonstrate FFTW on a 2D real array'
  write ( *, '(a,i8)' ) '  NX = ', nx
  write ( *, '(a,i8)' ) '  NY = ', ny
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Transform data to FFT coefficients.'
  write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
  write ( *, '(a)' ) '  the data.'
  write ( *, '(a)' ) '  Compare recovered data to original data.'
!
!  Compute the data.
!
  do j = 1, ny
    do i = 1, nx
      in(i,j) = r8_uniform_01 ( seed )
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Input Data:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, j, in(i,j)
    end do
  end do
!
!  Make a plan for the FFT, and forward transform the data.
!
  call dfftw_plan_dft_r2c_2d_ ( plan_forward, nx, ny, in, out, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_forward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Output FFT Coefficients:'
  write ( *, '(a)' ) ' '

  do i = 1, nh
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
    end do
  end do
!
!  Make a plan for the backward FFT, and recover the original data.
!
  call dfftw_plan_dft_c2r_2d_ ( plan_backward, nx, ny, out, in2, FFTW_ESTIMATE )

  call dfftw_execute_ ( plan_backward )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
  write ( *, '(a)' ) ' '

  do i = 1, nx
    do j = 1, ny
      write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) &
        i, j, in2(i,j) / real ( nx * ny, kind = rk )
    end do
  end do
!
!  Discard the information associated with the plans.
!
  call dfftw_destroy_plan_ ( plan_forward )
  call dfftw_destroy_plan_ ( plan_backward )

  return
end
subroutine c8vec_uniform_01 ( n, seed, c )

!*****************************************************************************80
!
!! C8VEC_UNIFORM_01 returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of values to compute.
!
!    Input/output, integer SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, complex ( kind = ck ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) c(n)
  integer i
  integer, parameter :: i4_huge = 2147483647
  real ( kind = rk ) r
  integer k
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  integer seed
  real ( kind = rk ) theta

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r = sqrt ( real ( seed, kind = rk ) * 4.656612875D-10 )

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    theta = 2.0D+00 * pi * ( real ( seed, kind = rk ) * 4.656612875D-10 )

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = ck )

  end do

  return
end
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    An R8 is a real ( kind = rk ) value.
!
!    For now, the input quantity SEED is an integer variable.
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2**31 - 1 )
!      r8_uniform_01 = seed / ( 2**31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley Interscience, page 95, 1998.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input/output, integer SEED, the "seed" value, which should
!    NOT be 0. On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: i4_huge = 2147483647
  integer k
  real ( kind = rk ) r8_uniform_01
  integer seed

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + i4_huge
  end if

  r8_uniform_01 = real ( seed, kind = rk ) * 4.656612875D-10

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
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
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

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
