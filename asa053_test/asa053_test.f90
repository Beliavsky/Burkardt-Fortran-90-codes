program main

!*****************************************************************************80
!
!! asa053_test() tests asa053().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'asa053_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa053().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'asa053_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests wshrt().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: np = 3

  real ( kind = rk ) :: d((np*(np+1))/2) = (/ &
    3.0D+00, &
    2.0D+00, 4.0D+00, &
    1.0D+00, 2.0D+00, 5.0D+00 /)
  integer n
  real ( kind = rk ) sa((np*(np+1))/2)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01()'
  write ( *, '(a)' ) '  wshrt() generates a single Wishart deviate.'

  n = 1

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  The number of variables is ', np
  write ( *, '(a,i4)' ) '  The number of degrees of freedom is ', n

  call r8utp_print ( np, d, '  The upper Cholesky factor:' )

  call wshrt ( d, n, np, sa )

  call r8pp_print ( np, sa, '  The sample matrix:' )

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests wshrt().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: np = 3

  integer, parameter :: npp = ( np * ( np + 1 ) ) / 2

  real ( kind = rk ) :: d(npp) = (/ &
    3.0D+00, &
    2.0D+00, 4.0D+00, &
    1.0D+00, 2.0D+00, 5.0D+00 /)
  integer i
  integer j
  integer k
  integer ki
  integer kj
  integer n
  real ( kind = rk ) s_average(npp)
  real ( kind = rk ) sa(npp)
  real ( kind = rk ) sigma(np,np)
  integer test_num

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  wshrt() generates many Wishart samples.'
  write ( *, '(a)' ) '  Compare to D'' * D * np / n.'

  n = 2

  write ( *, '(a)' ) ''
  write ( *, '(a,i4)' ) '  The number of variables is ', np
  write ( *, '(a,i4)' ) '  The number of degrees of freedom is ', n

  call r8utp_print ( np, d, '  The upper Cholesky factor:' )

  s_average(1:npp) = 0.0D+00

  test_num = 100000
  do i = 1, test_num
    call wshrt ( d, n, np, sa )
    s_average(1:npp) = s_average(1:npp) + sa(1:npp)
  end do

  s_average(1:npp) = s_average(1:npp) / dble ( test_num )

  call r8pp_print ( np, s_average, '  The averaged matrix:' )
!
!  Compare the result to ( D' * D ) * np / n.
!
  sigma(1:np,1:np) = 0.0D+00
  
  do i = 1, np
    do j = 1, np
      do k = 1, min ( i, j )
        ki = k + ( i * ( i - 1 ) ) / 2
        kj = k + ( j * ( j - 1 ) ) / 2
        sigma(i,j) = sigma(i,j) + d(ki) * d(kj)
      end do
      sigma(i,j) = sigma(i,j) * dble ( np ) / dble ( n )
    end do
  end do

  call r8mat_print ( np, np, sigma, '  Expected result:' )

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
!    27 August 2021
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

