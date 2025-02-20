program main

!*****************************************************************************80
!
!! walsh_transform_test() tests walsh_transform().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 April 2023
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'walsh_transform_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test walsh_transform().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'walsh_transform_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests FWT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 16

  integer i
  integer j
  real ( kind = rk ) w(n)
  real ( kind = rk ) work(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  FWT computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      call random_number ( harvest = w(1:n) )
    else
      do i = 1, n
        w(i) = real ( i, kind = rk )
      end do
    end if

    x(1:n) = w(1:n)
    call fwt ( n, w, work )
    y(1:n) = w(1:n) / real ( n, kind = rk )
    call fwt ( n, w, work )
    z(1:n) = w(1:n) / real ( n, kind = rk )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), z(i)
    end do

  end do
    
  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests WALSH.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 16

  integer i
  integer j
  real ( kind = rk ) w(n)
  real ( kind = rk ) work(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  WALSH computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      call random_number ( harvest = w(1:n) )
    else
      do i = 1, n
        w(i) = real ( i, kind = rk )
      end do
    end if

    x(1:n) = w(1:n)
    call walsh ( n, w, work )
    y(1:n) = w(1:n) / real ( n, kind = rk )
    call walsh ( n, w, work )
    z(1:n) = w(1:n) / real ( n, kind = rk )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) i, x(i), y(i), z(i)
    end do

  end do
    
  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests HAAR, HAARIN and HNORM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 16

  integer i
  integer j
  real ( kind = rk ) w(n)
  real ( kind = rk ) work(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  HAAR computes a Haar transform.'
  write ( *, '(a)' ) '  HNORM normalizes the transformed data.'
  write ( *, '(a)' ) '  HAARIN computes an inverse Haar transform.'

  do j = 1, 2

    if ( j == 1 ) then
      call random_number ( harvest = w(1:n) )
    else
      do i = 1, n
        w(i) = real ( i, kind = rk )
      end do
    end if

    x(1:n) = w(1:n)

    call haar ( n, w, work )

    y(1:n) = w(1:n)

    call hnorm ( n, w )

    z(1:n) = w(1:n)

    call haarin ( n, w, work )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
        i, x(i), y(i), z(i), w(i)
    end do

  end do
    
  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests FFWT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 March 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 16

  integer i
  integer j
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  FFWT computes a fast Walsh transform.'

  do j = 1, 2

    if ( j == 1 ) then
      call random_number ( harvest = w(1:n) )
    else
      do i = 1, n
        w(i) = real ( i, kind = rk )
      end do
    end if

    x(1:n) = w(1:n)
    call ffwt ( n, w )
    y(1:n) = w(1:n) / real ( n, kind = rk )
    call ffwt ( n, w )
    z(1:n) = w(1:n) / real ( n, kind = rk )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N'
    write ( *, '(a)' ) ' '
    do i = 1, n
      write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) &
        i, x(i), y(i), z(i)
    end do

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

