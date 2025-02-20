program main

!*****************************************************************************80
!
!! hypergeometric_test() tests hypergeometric().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypergeometric_test():'
  write ( *, '(a)' ) '  hypergeometric() evaluates the generalized'
  write ( *, '(a)' ) '  hypergeometric function pFq(z).'

  call hypergeometric_test01 ( )
  call hypergeometric_test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypergeometric_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  call timestamp ( )

  stop
end
subroutine hypergeometric_test01 ( )

!*****************************************************************************80
!
!! hypergeometric_test01() tests hypergeometric().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ), allocatable :: a(:)
  complex ( kind = ck ), allocatable :: b(:)
  complex ( kind = ck ) exact
  integer ip
  integer iq
  integer ix
  integer lnpfq
  integer nsigfig
  complex ( kind = ck ) pfq
  complex ( kind = ck ) value
  complex ( kind = ck ) z

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypergeometric_test01():'
  write ( *, '(a)' ) '  Test hypergeometric() on the example in the documentation.'
  write ( *, '(a)' ) ''

  ip = 2
  allocate ( a(1:ip) )
  iq = 3
  allocate ( b(1:iq) )
  
  a(1) = cmplx ( 1.0D+00,  1.0D+00, kind = ck )
  a(2) = cmplx ( 1.0D+00,  0.0D+00, kind = ck )
  b(1) = cmplx ( 2.0D+00, -1.0D+00, kind = ck )
  b(2) = cmplx ( 3.0D+00,  0.0D+00, kind = ck )
  b(3) = cmplx ( 3.0D+00,  0.0D+00, kind = ck )

  z = cmplx ( 1.5D+00, 0.0D+00, kind = ck )

  lnpfq = 0
  ix = 0
  nsigfig = 10

  value = pfq ( a, b, ip, iq, z, lnpfq, ix, nsigfig )

  exact = cmplx ( 1.02992154295955, 0.106416425916656, kind = ck )

  write ( *, '(2g24.16)' ) value
  write ( *, '(2g24.16)' ) exact

  deallocate ( a )
  deallocate ( b )

  return
end
subroutine hypergeometric_test02 ( )

!*****************************************************************************80
!
!! hypergeometric_test02() tests hypergeometric().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ), allocatable :: a(:)
  real ( kind = rk ) arg
  complex ( kind = ck ), allocatable :: b(:)
  integer ip
  integer iq
  integer ix
  integer lnpfq
  integer n
  integer nsigfig
  complex ( kind = ck ) pfq
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  complex ( kind = ck ) value
  complex ( kind = ck ) z

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'hypergeometric_test02():'
  write ( *, '(a)' ) '  Test hypergeometric() on a sequence pFq(z) with'
  write ( *, '(a)' ) '  q=p+1.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   n   p   q     z                pFq(z)' 

  do n = 1, 5

    arg = - ( n * pi / 2.0D+00 )**2
    write ( *, '(a)' ) ''

    do ip = 2, 9

      allocate ( a(1:ip) )
      iq = ip + 1
      allocate ( b(1:iq) )
  
      a(1:ip) = cmplx ( 1.25D+00,  0.0D+00, kind = ck )
      b(1) = cmplx ( 0.5D+00, 0.0D+00, kind = ck )
      b(2:iq) = cmplx ( 2.25D+00,  0.0D+00, kind = ck )

      z = cmplx ( arg, 0.0D+00, kind = ck )

      lnpfq = 0
      ix = 0
      nsigfig = 10

      value = pfq ( a, b, ip, iq, z, lnpfq, ix, nsigfig )

      write ( *, '(2x,i2,2x,i2,2x,i2,2x,g14.8,2x,g24.16)' ) &
        n, ip, iq, real ( z, kind = rk ), real ( value, kind = rk )

      deallocate ( a )
      deallocate ( b )

    end do

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

