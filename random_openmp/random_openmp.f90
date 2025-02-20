program main

!*****************************************************************************80
!
!! random_openmp() demonstrates the use of random numbers in an OpenMP program.
!
!  Discussion:
!
!    This program simply explores one issue in the generation of random
!    numbers in a parallel program.  If the random number generator uses
!    an integer seed to determine the next entry, then it is not easy for
!    a parallel program to reproduce the same exact sequence.
!
!    But what is worse is that it might not be clear how the separate
!    OpenMP threads should handle the SEED value - as a shared or private
!    variable?  It seems clear that each thread should have a private
!    seed that is initialized to a distinct value at the beginning of
!    the computation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2012
!
!  Author:
!
!    John Burkardt
!
  use omp_lib

  implicit none

  integer n
  integer seed

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'random_openmp():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  An OpenMP program using random numbers.'
  write ( *, '(a)' ) '  The random numbers depend on a seed.'
  write ( *, '(a)' ) '  We need to insure that each OpenMP thread'
  write ( *, '(a)' ) '  starts with a different seed.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) &
    '  Number of processors available = ', omp_get_num_procs ( )
  write ( * ,'(a,i8)' ) &
    '  Number of threads =              ', omp_get_max_threads ( )

  n = 100
  seed = 123456789
  call monte_carlo ( n, seed )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'random_openmp():'
  write ( *, '(a)' ) '  Normal end of execution.'

  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine monte_carlo ( n, seed )

!*****************************************************************************80
!
!! MONTE_CARLO carries out a Monte Carlo calculation with random values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of values to generate.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
  use omp_lib

  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer my_id
  integer my_seed
  integer seed
  real ( kind = rk ) x(n)

!$omp master
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Thread   Seed  I   X(I)'
  write ( *, '(a)' ) ' '
!$omp end master

!$omp parallel private ( i, my_id, my_seed ) shared ( n, x )
  my_id = omp_get_thread_num ( )
  my_seed = seed + my_id
  write ( *, '(2x,i6,2x,i12)' ) my_id, my_seed
!$omp do
  do i = 1, n
    call random_value ( my_seed, x(i) )
    write ( * , '(2x,i6,2x,i12,2x,i6,2x,g14.6)' ) my_id, my_seed, i, x(i)
  end do
!$omp end do

!$omp end parallel

  return
end
subroutine random_value ( seed, r )

!*****************************************************************************80
!
!! RANDOM_VALUE generates a random value R.
!
!  Discussion:
!
!    This is not a good random number generator.  It is a SIMPLE one.
!    It illustrates a model which works by accepting an integer seed value
!    as input, performing some simple operation on the seed, and then
!    producing a "random" real value using some simple transformation.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    03 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) R, the random value.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  integer seed

  seed = mod ( seed, 65536 )
  seed = mod ( ( 3125 * seed ), 65536 )
  r = real ( seed, kind = rk ) / 65536.0D+00

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
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
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
