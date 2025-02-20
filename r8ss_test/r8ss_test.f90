program main

!*****************************************************************************80
!
!! r8ss_test() tests r8ss().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    27 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ss_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8ss().'

  call r8ss_dif2_test ( )
  call r8ss_indicator_test ( )
  call r8ss_mv_test ( )
  call r8ss_print_test ( )
  call r8ss_print_some_test ( )
  call r8ss_random_test ( )
  call r8ss_to_r8ge_test ( )
  call r8ss_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ss_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ss_dif2_test ( )

!*****************************************************************************80
!
!! R8SS_DIF2_TEST tests R8SS_DIF2.
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

  real ( kind = rk ) a(2*n-1)
  integer diag(n)
  integer na

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_DIF2_TEST'
  write ( *, '(a)' ) '  R8SS_DIF2 sets the R8SS second difference matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ss_dif2 ( n, na, diag, a )

  call r8ss_print ( n, na, diag, a, '  The R8SS second difference matrix:' )

  return
end
subroutine r8ss_indicator_test ( )

!*****************************************************************************80
!
!! R8SS_INDICATOR_TEST tests R8SS_INDICATOR.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  integer diag(n)
  integer na

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_INDICATOR_TEST'
  write ( *, '(a)' ) '  R8SS_INDICATOR computes an R8SS indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ss_indicator ( n, na, diag, a )

  call r8ss_print ( n, na, diag, a, '  The R8SS indicator matrix:' )

  return
end
subroutine r8ss_mv_test ( )

!*****************************************************************************80
!
!! R8SS_MV_TEST tests R8SS_MV.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n)
  real ( kind = rk ) b2(n)
  integer diag(n)
  integer i
  integer ij
  integer ilo
  integer j
  integer na
  integer :: seed = 123456789
  real ( kind = rk ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_MV_TEST'
  write ( *, '(a)' ) '  R8SS_MV computes A*x for an R8SS matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ss_random ( n, seed, na, diag, a )
!
!  Replace the random entries by marker values.
!
  ij = 0
  do j = 1, n

    if ( j == 1 ) then
      ilo = 1
    else
      ilo = diag(j-1) - diag(j) + j + 1
    end if

    do i = ilo, j
      ij = ij + 1
      a(ij) = real ( 10 * i + j, kind = rk )
    end do

  end do
!
!  Print information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero entries stored is ', na

  call i4vec_print ( n, diag, '  Diagonal storage indices:' )

  call r8ss_print ( n, na, diag, a, '  The R8SS matrix:' )
!
!  Copy the matrix into a general matrix.
!
  call r8ss_to_r8ge ( n, na, diag, a, a2 )
!
!  Set the vector X.
!
  call r8vec_indicator1 ( n, x )
!
!  Compute the product.
!
  call r8ss_mv ( n, na, diag, a, x, b )
!
!  Compute the product using the general matrix.
!
  call r8ge_mv ( n, n, a2, x, b2 )
!
!  Compare the results.
!
  call r8vec2_print ( n, b, b2, '  R8SS_MV verse R8GE_MV' )

  return
end
subroutine r8ss_print_test ( )

!*****************************************************************************80
!
!! R8SS_PRINT_TEST tests R8SS_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    20 September 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  integer diag(n)
  integer i
  integer ij
  integer ilo
  integer j
  integer na
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_PRINT_TEST'
  write ( *, '(a)' ) '  R8SS_PRINT prints an R8SS matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ss_random ( n, seed, na, diag, a )
!
!  Replace the random entries by marker values.
!
  ij = 0
  do j = 1, n

    if ( j == 1 ) then
      ilo = 1
    else
      ilo = diag(j-1) - diag(j) + j + 1
    end if

    do i = ilo, j
      ij = ij + 1
      a(ij) = real ( 10 * i + j, kind = rk )
    end do

  end do
!
!  Print information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero entries stored is ', na

  call i4vec_print ( n, diag, '  Diagonal storage indices:' )

  call r8ss_print ( n, na, diag, a, '  The R8SS matrix:' )

  return
end
subroutine r8ss_print_some_test ( )

!*****************************************************************************80
!
!! R8SS_PRINT_SOME_TEST tests R8SS_PRINT_SOME.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 9

  real ( kind = rk ) a((n*(n+1))/2)
  integer diag(n)
  integer i
  integer ij
  integer ilo
  integer j
  integer na
  integer :: seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_PRINT_SOME_TEST'
  write ( *, '(a)' ) '  R8SS_PRINT_SOME prints some of an R8SS matrix'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
  call r8ss_random ( n, seed, na, diag, a )
!
!  Replace the random entries by marker values.
!
  ij = 0
  do j = 1, n

    if ( j == 1 ) then
      ilo = 1
    else
      ilo = diag(j-1) - diag(j) + j + 1
    end if

    do i = ilo, j
      ij = ij + 1
      a(ij) = real ( 10 * i + j, kind = rk )
    end do

  end do
!
!  Print information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero entries stored is ', na

  call i4vec_print ( n, diag, '  Diagonal storage indices:' )

  call r8ss_print_some ( n, na, diag, a, 2, 1, 8, 5, '  Rows 2-8, Cols 1:5' )

  return
end
subroutine r8ss_random_test ( )

!*****************************************************************************80
!
!! R8SS_RANDOM_TEST tests R8SS_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  integer diag(n)
  integer na
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_RANDOM_TEST'
  write ( *, '(a)' ) '  R8SS_RANDOM randomizes an R8SS matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ss_random ( n, seed, na, diag, a )

  call r8ss_print ( n, na, diag, a, '  The random R8SS matrix:' )

  return
end
subroutine r8ss_to_r8ge_test ( )

!*****************************************************************************80
!
!! R8SS_TO_R8GE_TEST tests R8SS_TO_R8GE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  real ( kind = rk ) a_r8ge(n,n)
  integer diag(n)
  integer na
  integer seed

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_TO_R8GE_TEST'
  write ( *, '(a)' ) '  R8SS_TO_R8GE converts an R8SS matrix to R8GE format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  seed = 123456789
  call r8ss_random ( n, seed, na, diag, a )

  call r8ss_print ( n, na, diag, a, '  The R8SS matrix:' )

  call r8ss_to_r8ge ( n, na, diag, a, a_r8ge )

  call r8ge_print ( n, n, a_r8ge, '  The R8GE matrix' )

  return
end
subroutine r8ss_zeros_test ( )

!*****************************************************************************80
!
!! R8SS_ZEROS_TEST tests R8SS_ZEROS.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 June 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 5

  real ( kind = rk ) a((n*(n+1))/2)
  integer diag(n)
  integer na

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'R8SS_ZEROS_TEST'
  write ( *, '(a)' ) '  R8SS_INDICATOR zeros an R8SS matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix order N = ', n

  call r8ss_zeros ( n, na, diag, a )
!
!  Print information.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Number of nonzero entries stored is ', na

  call i4vec_print ( n, diag, '  Diagonal storage indices:' )

  call r8ss_print ( n, na, diag, a, '  The R8SS zero matrix:' )

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

