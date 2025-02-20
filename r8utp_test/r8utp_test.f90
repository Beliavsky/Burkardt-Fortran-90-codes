program main

!*****************************************************************************80
!
!! r8utp_test() tests r8utp().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_test():'
  write ( *, '(a)' ) '  FORTRAN90 version:'
  write ( *, '(a)' ) '  Test r8utp().'

  call r8ge_to_r8utp_test ( )

  call r8utp_det_test ( )
  call r8utp_indicator_test ( )
  call r8utp_print_test ( )
  call r8utp_print_some_test ( )
  call r8utp_random_test ( )
  call r8utp_size_test ( )
  call r8utp_to_r8ge_test ( )
  call r8utp_zeros_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine r8ge_to_r8utp_test ( )

!*****************************************************************************80
!
!! r8ge_to_r8utp_test() tests r8ge_to_r8utp().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a_ge(m,n)
  real ( kind = rk ), allocatable :: a_utp(:)
  integer i
  integer j
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8ge_to_r8utp_test():'
  write ( *, '(a)' ) '  r8ge_to_r8utp() converts an r8ge matrix to r8utp format.'

  do j = 1, n
    do i = 1, m
      a_ge(i,j) = 10 * i + j
    end do
  end do

  call r8ge_print ( m, n, a_ge, '  The random r8ge matrix:' )

  mn = r8utp_size ( m, n )
  allocate ( a_utp(mn) )

  call r8ge_to_r8utp ( m, n, a_ge, a_utp )

  call r8utp_print ( m, n, a_utp, '  The r8utp matrix' )

  deallocate ( a_utp )

  return
end
subroutine r8utp_det_test ( )

!*****************************************************************************80
!
!! r8utp_det_test() tests r8utp_det().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 5

  real ( kind = rk ), allocatable :: a(:)
  real ( kind = rk ) det
  integer i
  integer j
  integer k
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_det_test():'
  write ( *, '(a)' ) '  r8utp_det() computes the determinant of an r8utp matrix.'

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  k = 1
  do j = 1, n
    do i = 1, min ( j, m )
      a(k) = k
      k = k + 1
    end do
  end do

  call r8utp_print ( m, n, a, '  The matrix A:' )
!
!  Compute the determinant.
!
  call r8utp_det ( m, n, a, det )

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Determinant is ', det

  deallocate ( a )

  return
end
subroutine r8utp_indicator_test ( )

!*****************************************************************************80
!
!! r8utp_indicator_test() tests r8utp_indicator().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ), allocatable :: a(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_indicator_test():'
  write ( *, '(a)' ) '  r8utp_indicator() sets up an indicator matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  call r8utp_indicator ( m, n, a )

  call r8utp_print ( m, n, a, '  The r8utp indicator matrix:' )

  deallocate ( a )

  return
end
subroutine r8utp_print_test ( )

!*****************************************************************************80
!
!! r8utp_print_test() tests r8utp_print().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 6
  integer, parameter :: n = 4

  real ( kind = rk ), allocatable :: a(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_print_test()'
  write ( *, '(a)' ) '  r8utp_print() prints an r8utp matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  call r8utp_indicator ( m, n, a )

  call r8utp_print ( m, n, a, '  The r8utp matrix:' )

  deallocate ( a )

  return
end
subroutine r8utp_print_some_test ( )

!*****************************************************************************80
!
!! r8utp_print_some_test() tests r8utp_print_some().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 4
  integer, parameter :: n = 6

  real ( kind = rk ), allocatable :: a(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_print_some_test():'
  write ( *, '(a)' ) '  r8utp_print_some() prints some of an r8utp matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
  write ( *, '(a,i8)' ) '  Matrix columns N = ', n

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  call r8utp_indicator ( m, n, a )

  call r8utp_print_some ( m, n, a, 1, 4, 3, 6, '  Some of the matrix:' )

  deallocate ( a )

  return
end
subroutine r8utp_random_test ( )

!*****************************************************************************80
!
!! r8utp_random_test() tests r8utp_random().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ), allocatable :: a(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_random_test():'
  write ( *, '(a)' ) '  r8utp_random() randomizes an r8utp matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,i8)' ) '  Matrix order M, N = ', m, n

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  call r8utp_random ( m, n, a )

  call r8utp_print ( m, n, a, '  The matrix:' )
 
  deallocate ( a )

  return
end
subroutine r8utp_size_test ( )

!*****************************************************************************80
!
!! r8utp_size_test() tests r8utp_size().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer j
  integer k
  integer m
  integer mn
  integer mn2
  integer n
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_size_test():'
  write ( *, '(a)' ) '  r8utp_size() determines storage for an R8UTP matrix.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   M   N   Size  Size(check)'
  write ( *, '(a)' ) ''

  m = 4

  do k = 1, 3

    if ( k == 1 ) then
      n = 3
    else if ( k == 2 ) then
      n = 4
    else
      n = 6
    end if

    mn = r8utp_size ( m, n )

    mn2 = 0
    do j = 1, n
      do i = 1, min ( j, m )
        mn2 = mn2 + 1;
      end do
    end do

    write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) m, n, mn, mn2

  end do

  return
end
subroutine r8utp_to_r8ge_test ( )

!*****************************************************************************80
!
!! r8utp_to_r8ge_test() tests r8utp_to_r8ge().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ) a_ge(m,n)
  real ( kind = rk ), allocatable :: a_ut(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_to_r8ge_test():'
  write ( *, '(a)' ) '  r8utp_to_r8ge() converts an r8utp matrix to r8ge format.'

  mn = r8utp_size ( m, n )

  allocate ( a_ut(1:mn) )

  call r8utp_random ( m, n, a_ut )

  call r8utp_print ( m, n, a_ut, '  The random r8utp matrix:' )

  call r8utp_to_r8ge ( m, n, a_ut, a_ge )

  call r8ge_print ( m, n, a_ge, '  The r8ge matrix' )

  deallocate ( a_ut )

  return
end
subroutine r8utp_zeros_test ( )

!*****************************************************************************80
!
!! r8utp_zeros_test() tests r8utp_zeros().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 August 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5
  integer, parameter :: n = 4

  real ( kind = rk ), allocatable :: a(:)
  integer mn
  integer r8utp_size

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'r8utp_zeros_test():'
  write ( *, '(a)' ) '  r8utp_zeros() zeros out a matrix in r8utp format.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8,i8)' ) '  Matrix order M, N = ', m, n

  mn = r8utp_size ( m, n )

  allocate ( a(1:mn) )

  call r8utp_zeros ( m, n, a )

  call r8utp_print ( m, n, a, '  The matrix:' )
 
  deallocate ( a )

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
!    15 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

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

