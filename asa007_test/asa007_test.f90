program main

!*****************************************************************************80
!
!! asa007_test() tests asa007().
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

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa007_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa007().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa007_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests syminv().
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

  integer, parameter :: n_max = 15

  real ( kind = rk ) a((n_max*(n_max+1))/2)
  real ( kind = rk ) afull(n_max,n_max)
  real ( kind = rk ) c((n_max*(n_max+1))/2)
  real ( kind = rk ) cfull(n_max,n_max)
  real ( kind = rk ) cta
  real ( kind = rk ) diff
  integer i
  integer ifault
  integer j
  integer k
  integer n
  integer nullty
  real ( kind = rk ) w(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  syminv() computes the inverse of a symmetric positive'
  write ( *, '(a)' ) '  definite matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the matrix A which is'
  write ( *, '(a)' ) '  N+1 on the diagonal and'
  write ( *, '(a)' ) '  N   on the off diagonals.'

  do n = 1, n_max
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
    k = 0
    do i = 1, n
      do j = 1, i - 1
        k = k + 1
        a(k) = real ( n, kind = rk )
      end do
      k = k + 1
      a(k) = real ( n + 1, kind = rk )
    end do

    call syminv ( a, n, c, w, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        cfull(i,j) = c(k)
        cfull(j,i) = c(k)
      end do
      k = k + 1
      cfull(j,j) = c(k)
    end do

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        afull(i,j) = a(k)
        afull(j,i) = a(k)
      end do
      k = k + 1
      afull(j,j) = a(k)
    end do
!
!  Compute C * A - I.
!
    diff = 0.0D+00
    do i = 1, n
      do j = 1, n
        cta = 0.0D+00
        do k = 1, n
          cta = cta + cfull(i,k) * afull(k,j)
        end do
        if ( i .eq. j ) then
          diff = diff + ( 1.0D+00 - cta )**2
        else
          diff = diff + cta**2
        end if
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests syminv().
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

  integer, parameter :: n_max = 15

  real ( kind = rk ) a((n_max*(n_max+1))/2)
  real ( kind = rk ) afull(n_max,n_max)
  real ( kind = rk ) c((n_max*(n_max+1))/2)
  real ( kind = rk ) cfull(n_max,n_max)
  real ( kind = rk ) cta
  real ( kind = rk ) diff
  integer i
  integer ifault
  integer j
  integer k
  integer n
  integer nullty
  real ( kind = rk ) w(n_max)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  syminv() computes the inverse of a symmetric positive'
  write ( *, '(a)' ) '  definite matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'

  do n = 1, n_max
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
    k = 0
    do i = 1, n
      do j = 1, i
        k = k + 1
        a(k) = 1.0D+00 / real ( i + j - 1, kind = rk )
      end do
    end do

    call syminv ( a, n, c, w, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        cfull(i,j) = c(k)
        cfull(j,i) = c(k)
      end do
      k = k + 1
      cfull(j,j) = c(k)
    end do

    k = 0
    do j = 1, n
      do i = 1, j - 1
        k = k + 1
        afull(i,j) = a(k)
        afull(j,i) = a(k)
      end do
      k = k + 1
      afull(j,j) = a(k)
    end do
!
!  Compute C * A - I.
!
    diff = 0.0D+00
    do i = 1, n
      do j = 1, n
        cta = 0.0D+00
        do k = 1, n
          cta = cta + cfull(i,k) * afull(k,j)
        end do
        if ( i .eq. j ) then
          diff = diff + ( 1.0D+00 - cta )**2
        else
          diff = diff + cta**2
        end if
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

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
 
