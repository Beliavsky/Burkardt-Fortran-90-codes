subroutine c8poly_print ( n, c, title )

!*****************************************************************************80
!
!! c8poly_print() prints out a polynomial.
!
!  Discussion:
!
!    The power sum form is:
!
!      p(x) = c(0) + c(1) * z + ... + c(n-1) * z^(n-1) + c(n) * z^(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the degree of the polynomial.
!
!    complex ( kind = ck ) c(0:N), the polynomial coefficients.
!    c(0) is the constant term and
!    c(N) is the coefficient of Z^N.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) c(0:n)
  integer i
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
  end if

  write ( *, '(a)' ) ''

  if ( n < 0 ) then
    write ( *, '( ''  p(x) = 0'' )' )
    return
  end if

  if ( 2 <= n ) then
    write ( *, '( ''  p(x) = ('',g14.6,g14.6,'') * z ^ '', i3 )' ) &
      c(n), n
  else if ( n == 1 ) then
    write ( *, '( ''  p(x) = ('',g14.6,g14.6,'') * z'' )' ) &
      c(n)
  else if ( n == 0 ) then
    write ( *, '( ''  p(x) = ('', g14.6, g14.6, '')'' )' ) c(n)
  end if

  do i = n - 1, 0, -1

    if ( 2 <= i ) then
      write ( *, ' ( ''         ('',g14.6,g14.6,'') * z ^ '', i3 )' ) c(i), i
    else if ( i == 1 ) then
      write ( *, ' ( ''         ('',g14.6,g14.6,'') * z'' )' ) c(i)
    else if ( i == 0 ) then
      write ( *, ' ( ''         ('',g14.6,g14.6, '')'' )' ) c(i)
    end if

  end do

  return
end
subroutine roots_to_c8poly ( n, x, c )

!*****************************************************************************80
!
!! roots_to_c8poly() converts polynomial roots to polynomial coefficients.
!
!  Discussion:
!
!    p(x) = c(0) x^n + x(1) x^(n-1) + ... + c(n)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of roots specified.
!
!    complex ( kind = ck ) X(N), the roots.
!
!  Output:
!
!    complex ( kind = ck ) C(0:N), the coefficients of the polynomial.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) c(0:n)
  integer i
  integer j
  complex ( kind = ck ) x(n)
!
!  Initialize C to (0, 0, ..., 0, 1).
!  Essentially, we are setting up a divided difference table.
!
  c(0:n-1) = 0.0D+00
  c(n) = 1.0D+00
!
!  Convert to standard polynomial form by shifting the abscissas
!  of the divided difference table to 0.
!
  do j = 1, n
    do i = 1, n + 1 - j
      c(n-i) = c(n-i) - x(n+1-i-j+1) * c(n-i+1)
    end do
  end do

  return
end
