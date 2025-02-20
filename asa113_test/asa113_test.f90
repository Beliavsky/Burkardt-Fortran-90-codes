program main

!*****************************************************************************80
!
!! asa113_test() tests asa113().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa113_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa113().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa113_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests asa113().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: k = 5
  integer, parameter :: m = 100
  integer, parameter :: n = 2

  real ( kind = rk ) a(m,n)
  integer c(m)
  real ( kind = rk ) c_center(k,n)
  integer c_size(k)
  integer ci
  real ( kind = rk ) critvl
  integer i
  integer ifault
  integer j
  integer ntrans1
  integer ntrans2
  real ( kind = rk ) wss(k)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Test asa113() classification algorithm.'
!
!  Read the data
!
  open ( unit = 1, file = 'points_100.txt', status = 'old' )

  do i = 1, m
    read ( 1, * ) ( a(i,j), j = 1, n )
  end do

  close ( unit = 1 )
!
!  Print first five points.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  First five points:'
  write ( *, '(a)' ) ' '
  do i = 1, 5
    write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
  end do
!
!  Assign points randomly to classes.
!
  do i = 1, m
    c(i) = mod ( i, k ) + 1
  end do
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = rk )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl = 0.0D+00
  do i = 1, k
    critvl = critvl + wss(i)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '        Initial CRITVL = ', critvl
!
!  Compute the clusters.
!
  ntrans1 = -1
  ntrans2 = -1

  do

    call trnsfr ( a, c, c_size, m, k, n, critvl, ntrans1, ifault )

    if ( ntrans1 == 0 .and. ntrans2 == 0 ) then
      exit
    end if

    write ( *, '(a,g14.6)' ) '  After TRNSFR, CRITVL = ', critvl

    call swap ( a, c, c_size, m, k, n, critvl, ntrans2, ifault )

    if ( ntrans1 == 0 .and. ntrans2 == 0 ) then
      exit
    end if

    write ( *, '(a,g14.6)' ) '    After SWAP, CRITVL = ', critvl

  end do
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = rk )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Cluster  Population  Energy'
  write ( *, '(a)' ) ' '

  do i = 1, k
    write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, c_size(i), wss(i)
  end do
  write ( *, '(a)' ) ' '
  write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', m, critvl

  return
end
subroutine crswap ( a, c, c_size, m, k, n, critvl, i1, i2, c1, c2, iswitch, &
  inc )

!*****************************************************************************80
!
!! crswap() determines the effect of swapping two objects.
!
!  Discussion:
!
!    This computation is very inefficient.  It is only set up so that we
!    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
!    ASA 136.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Input:
!
!    real ( kind = rk ) A(M,N), the data values.  There are M objects,
!    each having spatial dimension N.
!
!    integer C(M), the classification of each object.
!
!    integer C_SIZE(K), the number of objects in each class.
!
!    integer M, the number of objects.
!
!    integer K, the number of classes.
!
!    integer N, the number of spatial dimensions, or variates,
!    of the objects.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!    integer I1, I2, the objects to be swapped.
!
!    integer C1, C2, the current classes of objects I1 and I2.
!
!    integer ISWITCH:
!    1, indicates that I1 and I2 should be temporarily swapped, the
!       change in CRITVL should be computed, and then I1 and I2 restored.
!    2, indicates that I1 and I2 will be swapped.
!
!  Output:
!
!    real ( kind = rk ) INC, the change to CRITVL that would occur if I1 and
!    I2 were swapped.  This is only computed for ISWITCH = 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer c(m)
  real ( kind = rk ) c_center(k,n)
  integer c_size(k)
  integer c1
  integer c2
  integer ci
  real ( kind = rk ) critvl
  real ( kind = rk ) critvl_new
  integer i
  integer i1
  integer i2
  real ( kind = rk ) inc
  integer iswitch
  integer j
  real ( kind = rk ) wss(k)

  if ( iswitch == 2 ) then
    return
  end if
!
!  Move object I1 from class C1 to class C2.
!  Move object I2 from class C2 to class C1.
!
  c(i1) = c2
  c(i2) = c1
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = rk )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl_new = 0.0D+00
  do i = 1, k
    critvl_new = critvl_new + wss(i)
  end do

  inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!  Move object I2 from class C1 to class C2.
!
  c(i1) = c1
  c(i2) = c2

  return
end
subroutine crtran ( a, c, c_size, m, k, n, critvl, i1, c1, c2, iswitch, inc )

!*****************************************************************************80
!
!! crtran() determines the effect of moving an object to another class.
!
!  Discussion:
!
!    This computation is very inefficient.  It is only set up so that we
!    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
!    ASA 136.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Parameters:
!
!    real ( kind = rk ) A(M,N), the data values.  There are M objects,
!    each having spatial dimension N.
!
!    integer C(M), the classification of each object.
!
!    integer C_SIZE(K), the number of objects in each class.
!
!    integer M, the number of objects.
!
!    integer K, the number of classes.
!
!    integer N, the number of spatial dimensions, or variates,
!    of the objects.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!    integer I1, the object to be transferred.
!
!    integer C1, C2, the current class of object I1, and the
!    class to which it may be transferred.
!
!    integer ISWITCH:
!    1, indicates that I1 should be temporarily transferred, the
!       change in CRITVL should be computed, and then I1 restored.
!    2, indicates that I1 will be permanently transferred.
!
!  Output:
!
!    real ( kind = rk ) INC, the change to CRITVL that would occur if I1 were
!    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  integer c(m)
  real ( kind = rk ) c_center(k,n)
  integer c_size(k)
  integer c1
  integer c2
  integer ci
  real ( kind = rk ) critvl
  real ( kind = rk ) critvl_new
  integer i
  integer i1
  real ( kind = rk ) inc
  integer iswitch
  integer j
  real ( kind = rk ) wss(k)

  if ( iswitch == 2 ) then
    return
  end if
!
!  Move object I from class C1 to class C2.
!
  c(i1) = c2
  c_size(c1) = c_size(c1) - 1
  c_size(c2) = c_size(c2) + 1
!
!  Define the critical value as the sum of the squares of the distances
!  of the points to their cluster center.
!
  do i = 1, k
    c_size(i) = 0
    do j = 1, n
      c_center(i,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    ci = c(i)
    c_size(ci) = c_size(ci) + 1
    do j = 1, n
      c_center(ci,j) = c_center(ci,j) + a(i,j)
    end do
  end do

  do i = 1, k
    do j = 1, n
      c_center(i,j) = c_center(i,j) / real ( c_size(i), kind = rk )
    end do
  end do

  do i = 1, k
    wss(i) = 0.0D+00
  end do

  do i = 1, m
    ci = c(i)
    do j = 1, n
      wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
    end do
  end do

  critvl_new = 0.0D+00
  do i = 1, k
    critvl_new = critvl_new + wss(i)
  end do

  inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!
  c(i1) = c1
  c_size(c1) = c_size(c1) + 1
  c_size(c2) = c_size(c2) - 1

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
!    28 August 2021
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

