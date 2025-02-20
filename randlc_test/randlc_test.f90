program main

!*****************************************************************************80
!
!! randlc_test() tests randlc().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'randlc_test():'
  write ( *, '(a)' ) '  Fortran90 version:'
  write ( *, '(a)' ) '  Test randlc().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'randlc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  real ( kind = rk8 ) randlc
  real ( kind = rk8 ) seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  randlc() computes pseudorandom values '
  write ( *, '(a)' ) '  in the interval [0,1].'

  seed = 123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  The initial seed is ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         I          RANDLC'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    write ( *, '(2x,i8,2x,f14.6)' ) i, randlc ( seed )
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: n = 1000

  integer i
  real ( kind = rk8 ) randlc
  real ( kind = rk8 ) seed
  real ( kind = rk8 ) seed_in
  real ( kind = rk8 ) seed_out
  real ( kind = rk8 ) u(n)
  real ( kind = rk8 ) u_avg
  real ( kind = rk8 ) u_var

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  RANDLC computes a sequence of uniformly'
  write ( *, '(a)' ) '  distributed pseudorandom numbers.'

  seed = 123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Initial SEED = ', seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First 10 values:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    u(i) = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) &
      i, seed_in, seed_out, u(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,i6,a)' ) '  Now call RANDLC ', n, ' times.'

  u_avg = 0.0D+00
  do i = 1, n
    u(i) = randlc ( seed )
    u_avg = u_avg + u(i)
  end do

  u_avg = u_avg / real ( n, kind = rk8 )

  u_var = 0.0D+00
  do i = 1, n
    u_var = u_var + ( u(i) - u_avg )**2
  end do
  u_var = u_var / real ( n - 1, kind = rk8 )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Average value = ', u_avg
  write ( *, '(a,g14.6)' ) '  Expecting       ', 0.5D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Variance =      ', u_var
  write ( *, '(a,g14.6)' ) '  Expecting       ', 1.0D+00 / 12.0D+00

  return
  end
  subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests RANDLC.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  real ( kind = rk8 ) randlc
  real ( kind = rk8 ) seed
  real ( kind = rk8 ) seed_in
  real ( kind = rk8 ) seed_out
  real ( kind = rk8 ) seed_save
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  RANDLC computes a sequence of pseudorandom numbers'
  write ( *, '(a)' ) '  but all computations depend on the seed value.'
  write ( *, '(a)' ) '  In this test, we show how a sequence of "random"'
  write ( *, '(a)' ) '  values can be manipulated by accessing the seed.'

  seed = 1066.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Set SEED to ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call RANDLC 10 times, and watch SEED.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10

    seed_in = seed

    if ( i == 5 ) then
      seed_save = seed
    end if
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = seed_save

  write ( *, '(a)' ) ' '
  write ( *, '(a,f15.0)' ) '  Reset SEED to its value at step 5, = ', seed
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Now call RANDLC 10 times, and watch how SEED'
  write ( *, '(a)' ) '  and RANDLC restart themselves.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = 0.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  What happens with an initial zero SEED?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  seed = -123456789.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  What happens with an initial negative SEED?'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       I            Input           Output      RANDLC'
  write ( *, '(a)' ) '                     SEED             SEED'
  write ( *, '(a)' ) ' '

  do i = 1, 10
    seed_in = seed
    x = randlc ( seed )
    seed_out = seed
    write ( *, '(2x,i6,2x,f15.0,2x,f15.0,2x,f10.6)' ) i, seed_in, seed_out, x
  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! RANDLC_TEST04 tests RANDLC_JUMP.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer i
  integer k
  integer klog
  real ( kind = rk8 ) randlc
  real ( kind = rk8 ) randlc_jump
  real ( kind = rk8 ) seed
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RANDLC_TEST04'
  write ( *, '(a)' ) '  RANDLC_JUMP jumps directly to the K-th value'
  write ( *, '(a)' ) '  returned by RANDLC.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         K X(hard way)     X(jump)'
  write ( *, '(a)' ) ' '

  k = 1

  do klog = 1, 10

    seed = 123456789.0D+00
    do i = 1, k
      x1 = randlc ( seed )
    end do

    seed = 123456789.0D+00
    x2 = randlc_jump ( seed, k )

    write ( *, '(2x,i8,2x,f10.6,2x,f10.6)' ) k, x1, x2

    k = k * 2

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

