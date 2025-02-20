program main

!*****************************************************************************80
!
!! asa006_test() tests asa006().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa006_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test asa006().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa006_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests cholesky().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 15

  real ( kind = rk ) a((n_max*(n_max+1))/2)
  real ( kind = rk ) diff
  integer i
  integer ifault
  integer j
  integer k
  integer l
  integer n
  integer nn
  integer nullty
  real ( kind = rk ) u((n_max*(n_max+1))/2)
  real ( kind = rk ) ufull(n_max,n_max)
  real ( kind = rk ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  cholesky() computes the Cholesky factorization'
  write ( *, '(a)' ) '  of a symmetric positive definite matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the matrix A which is'
  write ( *, '(a)' ) '  N+1 on the diagonal and'
  write ( *, '(a)' ) '  N   on the off diagonals.'

  do n = 1, n_max

    nn = ( n * ( n + 1 ) ) / 2
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

    call cholesky ( a, n, nn, u, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu ) * ( a(k) - utu )
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests cholesky().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 15

  real ( kind = rk ) a((n_max*(n_max+1))/2)
  real ( kind = rk ) diff
  integer i
  integer ifault
  integer j
  integer k
  integer l
  integer n
  integer nn
  integer nullty
  real ( kind = rk ) u((n_max*(n_max+1))/2)
  real ( kind = rk ) ufull(n_max,n_max)
  real ( kind = rk ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  cholesky() computes the Cholesky factorization'
  write ( *, '(a)' ) '  of a symmetric positive definite matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'

  do n = 1, n_max

    nn = ( n * ( n + 1 ) ) / 2
!
!  Set A to the Hilbert matrix.
!
    k = 0
    do i = 1, n
      do j = 1, i
        k = k + 1
        a(k) = 1.0D+00 / real ( i + j - 1, kind = rk )
      end do
    end do

    call cholesky ( a, n, nn, u, nullty, ifault )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu ) * ( a(k) - utu )
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() tests subchl().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n_max = 15
  integer, parameter :: nn_max = (n_max*(n_max+1))/2

  real ( kind = rk ) a(nn_max)
  integer b(n_max)
  real ( kind = rk ) det
  real ( kind = rk ) diff
  integer i
  integer ifault
  integer j
  integer k
  integer l
  integer n
  integer nullty
  real ( kind = rk ) u(nn_max)
  real ( kind = rk ) ufull(n_max,n_max)
  real ( kind = rk ) utu

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  subchl() computes the Cholesky factor of a submatrix '
  write ( *, '(a)' ) '  of a symmetric positive definite matrix.'
  write ( *, '(a)' ) '  A compressed storage format is used.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
  write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  For this particular matrix, we expect the'
  write ( *, '(a)' ) '  errors to grow rapidly.'
!
!  Set A to the N_MAX order Hilbert matrix.
!
  k = 0
  do i = 1, n_max
    do j = 1, i
      k = k + 1
      a(k) = 1.0D+00 / real ( i + j - 1, kind = rk )
    end do
  end do

  do n = 1, n_max

    do i = 1, n
      b(i) = i
    end do

    call subchl ( a, b, n, u, nullty, ifault, nn_max, det )

    write ( *, '(a)' ) ' '
    write ( *, '(a,i8)' ) '  Matrix order N = ', n
    write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty
    write ( *, '(a,g14.6)' ) '  Matrix determinant DET = ', det

    k = 0
    do j = 1, n
      do i = 1, j
        k = k + 1
        ufull(i,j) = u(k)
      end do
      do i = j + 1, n
        ufull(i,j) = 0.0D+00
      end do
    end do
!
!  Compute U' * U and compare to A.
!
    k = 0
    diff = 0.0D+00
    do i = 1, n
      do j = 1, i
        k = k + 1
        utu = 0.0D+00
        do l = 1, n
          utu = utu + ufull(l,i) * ufull(l,j)
        end do
        diff = diff + ( a(k) - utu )**2
      end do
    end do

    diff = sqrt ( diff )

    write ( *, '(a,g14.6)' ) '  RMS ( A - U''*U ) = ', diff
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
!    26 August 2021
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

