program main

!*****************************************************************************80
!
!! r8but_test() tests r8but().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8but_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8but().'

  call r8but_det_test ( )
  call r8but_indicator_test ( )
  call r8but_mtv_test ( )
  call r8but_mv_test ( )
  call r8but_print_test ( )
  call r8but_print_some_test ( )
  call r8but_random_test ( )
  call r8but_sl_test ( )
  call r8but_slt_test ( )
  call r8but_to_r8ge_test ( )
  call r8but_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8but_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8but_det_test ( )

!*****************************************************************************80
!
!! R8BUT_DET_TEST tests R8BUT_DET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 10

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) det
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_DET_TEST'
  write ( *, '(a)' ) '  R8BUT_DET computes the determinant of an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8but_det ( n, mu, a, det )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Determinant = ', det

  return
end
subroutine r8but_indicator_test ( )

!*****************************************************************************80
!
!! R8BUT_INDICATOR_TEST tests R8BUT_INDICATOR.
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

  integer, parameter :: n = 6
  integer, parameter :: mu = 2

  real ( kind = rk ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8BUT_INDICATOR sets up an R8BUT indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N     = ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu
!
!  Set the matrix.
!
  call r8but_indicator ( n, mu, a )

  call r8but_print ( n, mu, a, '  The R8BUT indicator matrix:' )

  return
end
subroutine r8but_mtv_test ( )

!*****************************************************************************80
!
!! R8BUT_MTV_TEST tests R8BUT_MTV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_MTV_TEST'
  write ( *, '(a)' ) '  R8BUT_MTV computes b=A''*x, where A is an R8BUT matrix.'

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8but_mtv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8but_mv_test ( )

!*****************************************************************************80
!
!! R8BUT_MV_TEST tests R8BUT_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    29 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_MV_TEST'
  write ( *, '(a)' ) '  R8BUT_MV computes b=A*x, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  call r8vec_indicator1 ( n, x )
  call r8vec_print ( n, x, '  x:' )

  call r8but_mv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b=A*x:' )

  return
end
subroutine r8but_print_test ( )

!*****************************************************************************80
!
!! R8BUT_PRINT_TEST tests R8BUT_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a(mu+1,n)
  integer :: seed = 123456789
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRINT_TEST'
  write ( *, '(a)' ) '  R8BUT_PRINT prints an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  return
end
subroutine r8but_print_some_test ( )

!*****************************************************************************80
!
!! R8BUT_PRINT_SOME_TEST tests R8BUT_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 10

  real ( kind = rk ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8BUT_PRINT_SOME prints some of an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_indicator ( n, mu, a )

  call r8but_print_some ( n, mu, a, 1, 2, 4, 4, '  Rows 1:4, Cols 2:4:' )

  return
end
subroutine r8but_random_test ( )

!*****************************************************************************80
!
!! R8BUT_RANDOM_TEST tests R8BUT_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    24 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a(mu+1,n)
  integer :: seed = 123456789
 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_RANDOM_TEST'
  write ( *, '(a)' ) '  R8BUT_RANDOM randomizes an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

  return
end
subroutine r8but_sl_test ( )

!*****************************************************************************80
!
!! R8BUT_SL_TEST tests R8BUT_SL.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 10

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_SL_TEST'
  write ( *, '(a)' ) '  R8BUT_SL solves A*x=b, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8but_mv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b:' )
!
!  Solve the linear system.
!
  call r8but_sl ( n, mu, a, b )
  call r8vec_print ( n, b, '  x:' )

  return
end
subroutine r8but_slt_test ( )

!*****************************************************************************80
!
!! R8BUT_SLT_TEST tests R8BUT_SLT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 10

  real ( kind = rk ) a(mu+1,n)
  real ( kind = rk ) b(n)
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_SLT_TEST'
  write ( *, '(a)' ) '  R8BUT_SLT solves A''*x=b, where A is an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a )
  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )
!
!  Set the desired solution.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the corresponding right hand side.
!
  call r8but_mtv ( n, mu, a, x, b )
  call r8vec_print ( n, b, '  b:' )
!
!  Solve the linear system.
!
  call r8but_slt ( n, mu, a, b )
  call r8vec_print ( n, b, '  x:' )

  return
end
subroutine r8but_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8BUT_TO_R8GE_TEST tests R8BUT_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a_r8but(mu+1,n)
  real ( kind = rk ) a_r8ge(n,n)
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8BUT_TO_R8GE converts a matrix from R8BUT to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_random ( n, mu, seed, a_r8but )

  call r8but_print ( n, mu, a_r8but, '  The R8BUT matrix:' )

  call r8but_to_r8ge ( n, mu, a_r8but, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix' )

  return
end
subroutine r8but_zeros_test ( )

!*****************************************************************************80
!
!! R8BUT_ZEROS_TEST tests R8BUT_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 October 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: mu = 3
  integer, parameter :: n = 5

  real ( kind = rk ) a(mu+1,n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8BUT_ZEROS_TEST'
  write ( *, '(a)' ) '  R8BUT_ZEROS zeros an R8BUT matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N =     ', n
  write ( *, '(a,i8)' ) '  Upper bandwidth MU = ', mu

  call r8but_zeros ( n, mu, a )

  call r8but_print ( n, mu, a, '  The R8BUT matrix:' )

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

