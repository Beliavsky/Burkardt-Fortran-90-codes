program main

!*****************************************************************************80
!
!! ZIGGURAT_TEST() tests ZIGGURAT().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: sample_num = 1000000

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ziggurat_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ziggurat().'
!
!  Make sure that SEED controls the sequence and can restart it.
!
  call test001 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Measure the time it takes to generate a lot of variables.
!
  call test05 ( sample_num )
  call test06 ( sample_num )
  call test07 ( sample_num )
  call test08 ( sample_num )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ziggurat_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test001 ( )

!*****************************************************************************80
!
!! TEST001 tests SHR3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer j
  integer seed
  integer shr3
  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST001'
  write ( *, '(a)' ) '  SHR3 returns pseudorandom uniformly distributed'
  write ( *, '(a)' ) '  integers.'

  do j = 0, 2

    if ( mod ( j, 2 ) == 0 ) then
      seed = 123456789
    else
      seed = 987654321
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6,2x,i12)' ) 0, seed
    write ( *, '(a)' ) ' '

    do i = 1, 10
      value = shr3 ( seed )
      write ( *, '(2x,i6,2x,i12,2x,i12)' ) i, seed, value
    end do

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests R4_UNI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer i
  integer j
  real r4_uni
  integer seed
  real value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  R4_UNI returns pseudorandom uniformly distributed'
  write ( *, '(a)' ) '  real numbers between 0 and 1.'

  do j = 0, 2

    if ( mod ( j, 2 ) == 0 ) then
      seed = 123456789
    else
      seed = 987654321
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6,2x,i12)' ) 0, seed
    write ( *, '(a)' ) ' '

    do i = 1, 10
      value = r4_uni ( seed )
      write ( *, '(2x,i6,2x,i12,2x,g14.6)' ) i, seed, value
    end do

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests R4_NOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fn(128)
  integer i
  integer j
  integer kn(128)
  real r4_nor
  integer seed
  real value
  real wn(128)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  R4_NOR returns pseudorandom normally distributed'
  write ( *, '(a)' ) '  real numbers between 0 and 1.'

  call r4_nor_setup ( kn, fn, wn )

  do j = 0, 2

    if ( mod ( j, 2 ) == 0 ) then
      seed = 123456789
    else
      seed = 987654321
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6,2x,i12)' ) 0, seed
    write ( *, '(a)' ) ' '

    do i = 1, 10
      value = r4_nor ( seed, kn, fn, wn )
      write ( *, '(2x,i6,2x,i12,2x,g14.6)' ) i, seed, value
    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests R4_EXP.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fe(256)
  integer i
  integer j
  integer ke(256)
  real r4_exp
  integer seed
  real value
  real we(256)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  R4_EXP returns pseudorandom exponentially distributed'
  write ( *, '(a)' ) '  real numbers between 0 and 1.'

  call r4_exp_setup ( ke, fe, we )

  do j = 0, 2

    if ( mod ( j, 2 ) == 0 ) then
      seed = 123456789
    else
      seed = 987654321
    end if

    write ( *, '(a)' ) ' '
    write ( *, '(2x,i6,2x,i12)' ) 0, seed
    write ( *, '(a)' ) ' '

    do i = 1, 10
      value = r4_exp ( seed, ke, fe, we )
      write ( *, '(2x,i6,2x,i12,2x,g14.6)' ) i, seed, value
    end do

  end do

  return
end
subroutine test05 ( sample_num )

!*****************************************************************************80
!
!! TEST05 times SHR3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer sample
  integer sample_num
  integer seed
  integer shr3
  real ( kind = rk ) time1
  real ( kind = rk ) time2
  integer value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Measure the time it takes SHR3 to generate'
  write ( *, '(i8,a)' ) sample_num, ' integers.'

  seed = 123456789

  call cpu_time ( time1 )

  do sample = 1, sample_num
    value = shr3 ( seed )
  end do

  call cpu_time ( time2 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a)' ) time2 - time1, ' seconds'

  return
end
subroutine test06 ( sample_num )

!*****************************************************************************80
!
!! TEST06 times R4_UNI.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real r4_uni
  integer sample
  integer sample_num
  integer seed
  real ( kind = rk ) time1
  real ( kind = rk ) time2
  real value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Measure the time it takes R4_UNI to generate'
  write ( *, '(i8,a)' ) sample_num, ' uniform deviates.'

  seed = 123456789

  call cpu_time ( time1 )

  do sample = 1, sample_num
    value = r4_uni ( seed )
  end do

  call cpu_time ( time2 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a)' ) time2 - time1, ' seconds'

  return
end
subroutine test07 ( sample_num )

!*****************************************************************************80
!
!! TEST07 times R4_NOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real fn(128)
  integer kn(128)
  real r4_nor
  integer sample
  integer sample_num
  integer seed
  real ( kind = rk ) time1
  real ( kind = rk ) time2
  real value
  real wn(128)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Measure the time it takes R4_NOR to generate'
  write ( *, '(i8,a)' ) sample_num, ' normal deviates.'

  call r4_nor_setup ( kn, fn, wn )

  seed = 123456789

  call cpu_time ( time1 )

  do sample = 1, sample_num
    value = r4_nor ( seed, kn, fn, wn )
  end do

  call cpu_time ( time2 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a)' ) time2 - time1, ' seconds'

  return
end
subroutine test08 ( sample_num )

!*****************************************************************************80
!
!! TEST08 times R4_EXP.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real fe(256)
  integer ke(256)
  real r4_exp
  integer sample
  integer sample_num
  integer seed
  real ( kind = rk ) time1
  real ( kind = rk ) time2
  real value
  real we(256)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Measure the time it takes R4_EXP to generate'
  write ( *, '(i8,a)' ) sample_num, ' normal deviates.'

  call r4_exp_setup ( ke, fe, we )

  seed = 123456789

  call cpu_time ( time1 )

  do sample = 1, sample_num
    value = r4_exp ( seed, ke, fe, we )
  end do

  call cpu_time ( time2 )

  write ( *, '(a)' ) ' '
  write ( *, '(2x,g14.6,a)' ) time2 - time1, ' seconds'

  return
end

