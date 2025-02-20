program main

!*****************************************************************************80
!
!! r8ci_test() tests r8ci().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ci_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ci().'

  call r8ci_det_test ( )
  call r8ci_dif2_test ( )
  call r8ci_eval_test ( )
  call r8ci_indicator_test ( )
  call r8ci_mtv_test ( )
  call r8ci_mv_test ( )
  call r8ci_print_test ( )
  call r8ci_print_some_test ( )
  call r8ci_random_test ( )
  call r8ci_sl_test ( )
  call r8ci_to_r8ge_test ( )
  call r8ci_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ci_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ci_det_test ( )

!*****************************************************************************80
!
!! R8CI_DET_TEST tests R8CI_DET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) det
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_DET_TEST'
  write ( *, '(a)' ) '  R8CI_DET finds the determinant of '
  write ( *, '(a)' ) '  a real circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  call r8ci_det ( n, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Computed determinant = ', det

  return
end
subroutine r8ci_dif2_test ( )

!*****************************************************************************80
!
!! R8CI_DIF2_TEST tests R8CI_DIF2.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_DIF2_TEST'
  write ( *, '(a)' ) '  R8CI_DIF2 sets up an R8CI periodic second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_dif2 ( n, a )

  call r8ci_print ( n, a, '  The R8CI second difference matrix:' )

  return
end
subroutine r8ci_eval_test ( )

!*****************************************************************************80
!
!! R8CI_EVAL_TEST tests R8CI_EVAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  complex ( kind = ck ) lambda(n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_EVAL_TEST'
  write ( *, '(a)' ) '  R8CI_EVAL finds the eigenvalues of '
  write ( *, '(a)' ) '  a real circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  call r8ci_eval ( n, a, lambda )

  call c8vec_print ( n, lambda, '  The eigenvalues:' )

  return
end
subroutine r8ci_indicator_test ( )

!*****************************************************************************80
!
!! R8CI_INDICATOR_TEST tests R8CI_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8CI_INDICATOR sets up an R8CI indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_mtv_test ( )

!*****************************************************************************80
!
!! R8CI_MTV_TEST tests R8CI_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_MTV_TEST'
  write ( *, '(a)' ) '  R8CI_MTV computes b=A''*x, where A is an R8CI matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )
  call r8ci_print ( n, a, '  The circulant matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector x:' );

  call r8ci_mtv ( n, a, x, b )
  call r8vec_print ( n, b, '  The product b=A''*x:' )

  return
end
subroutine r8ci_mv_test ( )

!*****************************************************************************80
!
!! R8CI_MV_TEST tests R8CI_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_MV_TEST'
  write ( *, '(a)' ) '  R8CI_MV computes b=A*x, where A is an R8CI matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )
  call r8ci_print ( n, a, '  The circulant matrix A:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  The vector x:' );

  call r8ci_mv ( n, a, x, b )
  call r8vec_print ( n, b, '  The product b=A*x:' )

  return
end
subroutine r8ci_print_test ( )

!*****************************************************************************80
!
!! R8CI_PRINT_TEST tests R8CI_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_PRINT_TEST'
  write ( *, '(a)' ) '  R8CI_PRINT prints an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_print_some_test ( )

!*****************************************************************************80
!
!! R8CI_PRINT_SOME_TEST tests R8CI_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8CI_PRINT_SOME prints some of an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print_some ( n, a, 2, 3, 6, 5, '  Rows 2-6, Cols 3-5:' )

  return
end
subroutine r8ci_random_test ( )

!*****************************************************************************80
!
!! R8CI_RANDOM_TEST tests R8CI_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_RANDOM_TEST'
  write ( *, '(a)' ) '  R8CI_RANDOM sets a random R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  return
end
subroutine r8ci_sl_test ( )

!*****************************************************************************80
!
!! R8CI_SL_TEST tests R8CI_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10

  real ( kind = rk ) a(n)
  real ( kind = rk ) b(n)
  integer job
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_SL_TEST'
  write ( *, '(a)' ) '  R8CI_SL solves a circulant system.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ci_random ( n, seed, a )

  call r8ci_print ( n, a, '  The circulant matrix:' )

  do job = 0, 1
!
!  Set the desired solution.
!
    call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
    if ( job == 0 ) then
      call r8ci_mv ( n, a, x, b )
    else
      call r8ci_mtv ( n, a, x, b )
    end if
!
!  Solve the linear system.
!
    call r8ci_sl ( n, a, b, x, job )

    if ( job == 0 ) then
      call r8vec_print ( n, x, '  Solution:' )
    else
      call r8vec_print ( n, x, '  Solution to transposed system:' )
    end if

  end do
 
  return
end
subroutine r8ci_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8CI_TO_R8GE_TEST tests R8CI_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)
  real ( kind = rk ) a_r8ge(n,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8CI_TO_R8GE converts an R8CI matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_indicator ( n, a )

  call r8ci_print ( n, a, '  The R8CI matrix:' )

  call r8ci_to_r8ge ( n, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix:' )

  return
end
subroutine r8ci_zeros_test ( )

!*****************************************************************************80
!
!! R8CI_ZEROS_TEST tests R8CI_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8CI_ZEROS_TEST'
  write ( *, '(a)' ) '  R8CI_ZEROS zeros an R8CI matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ci_zeros ( n, a )

  call r8ci_print ( n, a, '  The zero R8CI matrix:' )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

