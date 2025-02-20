subroutine ffwt ( n, x )

!*****************************************************************************80
!
!! ffwt() performs an in-place fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Walsh transform on an input series X
!    leaving the transformed results in X. 
!    X is dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) hold
  integer i
  integer i4_log_2
  integer ii
  integer j
  integer j2
  integer js
  integer k
  integer l
  integer m
  integer mw
  integer mw1
  integer nw
  integer nz
  integer nz2
  integer nzi
  integer nzn
  integer two_power(24)
  real ( kind = rk ) x(n)
  real ( kind = rk ) z

  m = i4_log_2 ( n )

  do i = 1, m
    two_power(i) = 2 ** ( m - i )
  end do

  do l = 1, m

    nz = 2 ** ( l - 1 )
    nzi = 2 * nz
    nzn = n / nzi
    nz2 = nz / 2
    if ( nz2 == 0 ) then
      nz2 = 1
    end if

    do i = 1, nzn

      js = ( i - 1 ) * nzi
      z = 1.0D+00

      do ii = 1, 2

        do j = 1, nz2
          js = js + 1
          j2 = js + nz
          hold = x(js) + z * x(j2)
          z = - z
          x(j2) = x(js) + z * x(j2)
          x(js) = hold
          z = - z
        end do

        if ( l == 1 ) then
          exit
        end if

        z = - 1.0D+00

      end do

    end do
  end do
!
!  Bit reversal section.
!
  nw = 0
  do k = 1, n
!
!  Choose correct index and switch elements if not already switched.
!
    if ( k < nw + 1 ) then
      hold = x(nw+1)
      x(nw+1) = x(k)
      x(k) = hold
    end if
!
!  Bump up series by 1.
!
    do i = 1, m

      ii = i
      if ( nw < two_power(i) ) then
        exit
      end if
      mw = nw / two_power(i)
      mw1 = mw / 2
      if ( mw <= 2 * mw1 ) then
        exit
      end if

      nw = nw - two_power(i)

    end do

    nw = nw + two_power(ii)

  end do

  return
end
subroutine fwt ( n, x, y )

!*****************************************************************************80
!
!! fwt() performs a fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Walsh transform on an input series X
!    leaving the transformed results in X. 
!    X is dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
!
!    Note that the program text in the reference included the line
!      y(jd) = abs ( x(j) - x(j2) )
!    which has been corrected to:
!      y(jd) = x(j) - x(j2)
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
!    Workspace, real ( kind = rk ) Y(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_log_2
  integer j
  integer j2
  integer jd
  integer js
  integer l
  integer m
  integer n2
  integer nx
  integer ny
  integer nz
  integer nzi
  integer nzn
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  n2 = n / 2
  m = i4_log_2 ( n )

  do l = 1, m

    ny = 0
    nz = 2 ** ( l - 1 )
    nzi = 2 * nz
    nzn = n / nzi

    do i = 1, nzn

      nx = ny + 1
      ny = ny + nz
      js = ( i - 1 ) * nzi
      jd = js + nzi + 1

      do j = nx, ny
        js = js + 1
        j2 = j + n2
        y(js) = x(j) + x(j2)
        jd = jd - 1
        y(jd) = x(j) - x(j2)
      end do

    end do

    x(1:n) = y(1:n)

  end do

  return
end
subroutine haar ( n, x, y )

!*****************************************************************************80
!
!! haar() performs a Haar transform.
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
!    Workspace, real ( kind = rk ) Y(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_log_2
  integer j
  integer jj
  integer k
  integer l
  integer l2
  integer l3
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  k = i4_log_2 ( n )

  do i = 1, k

    l = k + 1 - i
    l2 = 2 ** ( l - 1 )

    y(1:2*l2) = x(1:2*l2)

    do j = 1, l2
       l3 = l2 + j
       jj = 2 * j - 1
       x(j) = y(jj) + y(jj+1)
       x(l3) = y(jj) - y(jj+1)
    end do

  end do

  return
end
subroutine haarin ( n, x, y )

!*****************************************************************************80
!
!! HAARIN inverts a Haar transform.
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
!    Workspace, real ( kind = rk ) Y(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_log_2
  integer j
  integer jj
  integer jj1
  integer k
  integer l
  integer lj
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  k = i4_log_2 ( n )

  do i = 1, k

    l = 2 ** ( i - 1 )

    y(1:2*l) = x(1:2*l)

    do j = 1, l
      lj = l + j
      jj = 2 * j
      jj1 = jj - 1
      x(jj) = y(j) - y(lj)
      x(jj1) = y(j) + y(lj)
    end do

  end do

  return
end
subroutine hnorm ( n, x )

!*****************************************************************************80
!
!! HNORM computes normalization factors for a forward or inverse Haar transform.
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_log_2
  integer ii
  integer jmax
  integer jmin
  integer k
  real ( kind = rk ) wlk
  real ( kind = rk ) x(n)

  k = i4_log_2 ( n )

  x(1) = x(1) / 2.0D+00 ** k

  if ( 1 <= k ) then
    x(2) = x(2) / 2.0D+00 ** k
  end if

  do ii = 2, k

    i = ii - 1
    wlk = 1.0D+00 / 2.0D+00 ** ( k - i )
    jmin = 2 ** i + 1
    jmax = 2 ** ii

    x(jmin:jmax) = x(jmin:jmax) * wlk

  end do

  return
end
function i4_log_2 ( i )

!*****************************************************************************80
!
!! I4_LOG_2 returns the integer part of the logarithm base 2 of an I4.
!
!  Discussion:
!
!    For positive I4_LOG_2(I), it should be true that
!      2^I4_LOG_2(X) <= |I| < 2^(I4_LOG_2(I)+1).
!    The special case of I4_LOG_2(0) returns -HUGE().
!
!    An I4 is an integer value.
!
!  Example:
!
!     I  I4_LOG_2
!
!     0  -1
!     1,  0
!     2,  1
!     3,  1
!     4,  2
!     5,  2
!     6,  2
!     7,  2
!     8,  3
!     9,  3
!    10,  3
!   127,  6
!   128,  7
!   129,  7
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    13 January 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number whose logarithm base 2
!    is desired.
!
!    Output, integer I4_LOG_2, the integer part of the
!    logarithm base 2 of the absolute value of I.
!
  implicit none

  integer i
  integer i_abs
  integer i4_log_2
  integer, parameter :: i4_huge = 2147483647

  if ( i == 0 ) then

    i4_log_2 = - i4_huge

  else

    i4_log_2 = 0

    i_abs = abs ( i )

    do while ( 2 <= i_abs )
      i_abs = i_abs / 2
      i4_log_2 = i4_log_2 + 1
    end do

  end if

  return
end
function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of I4 division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!    An I4 is an integer value.
!
!  Example:
!
!        I     J     MOD I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer I, the number to be divided.
!
!    Input, integer J, the number that divides I.
!
!    Output, integer I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer i
  integer i4_modp
  integer j
  integer value

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
    stop 1
  end if

  value = mod ( i, j )

  if ( value < 0 ) then
    value = value + abs ( j )
  end if

  i4_modp = value

  return
end
function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Discussion:
!
!    An I4 is an integer value.
!
!    There appears to be a bug in the GFORTRAN compiler which can lead to
!    erroneous results when the first argument of I4_WRAP is an expression.
!    In particular:
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i4_wrap ( i + 1, 1, 3 )
!      end if
!    end do
!
!    was, when I = 3, returning I4 = 3.  So I had to replace this with
!
!    do i = 1, 3
!      if ( test ) then
!        i4 = i + 1
!        i4 = i4_wrap ( i4, 1, 3 )
!      end if
!    end do
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  Value
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 September 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer IVAL, a value.
!
!    Input, integer ILO, IHI, the desired bounds.
!
!    Output, integer I4_WRAP, a "wrapped" version of the value.
!
  implicit none

  integer i4_modp
  integer i4_wrap
  integer ihi
  integer ilo
  integer ival
  integer jhi
  integer jlo
  integer value
  integer wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    value = jlo
  else
    value = jlo + i4_modp ( ival - jlo, wide )
  end if

  i4_wrap = value

  return
end
subroutine r8vec_shift_circular ( shift, n, x )

!*****************************************************************************80
!
!! R8VEC_SHIFT_CIRCULAR performs a circular shift on an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    11 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer SHIFT, the amount by which each entry is to
!    be shifted.
!
!    Input, integer N, the length of the vector.
!
!    Input/output, real ( kind = rk ) X(N), the vector to be shifted.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer i4_wrap
  integer j
  integer shift
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)

  y(1:n) = x(1:n)

  do i = 1, n
    j = i4_wrap ( i - shift, 1, n )
    x(i) = y(j)
  end do

  return
end
subroutine walsh ( n, x, y )

!*****************************************************************************80
!
!! walsh() performs a fast Walsh transform.
!
!  Discussion:
!
!    This routine performs a fast Wash transform on an input series X
!    leaving the transformed results in X.  The array Y is working space.
!    X and Y are dimensioned N, which must be a power of 2.
!    The results of this Walsh transform are in sequency order.
!
!    The output sequence could be normalized by dividing by N.
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
!    Ken Beauchamp
!
!  Reference:
!
!    Ken Beauchamp,
!    Walsh functions and their applications,
!    Academic Press, 1975,
!    ISBN: 0-12-084050-2,
!    LC: QA404.5.B33.
!
!  Parameters:
!
!    Input, integer N, the number of items in X.
!    N must be a power of 2.
!
!    Input/output, real ( kind = rk ) X(N), the data to be transformed.
!
!    Workspace, real ( kind = rk ) Y(N).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  integer i
  integer i1
  integer i4_log_2
  integer is
  integer j
  integer j1
  integer l
  integer m
  integer n1
  integer n2
  real ( kind = rk ) w
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n/2)
  real ( kind = rk ) z

  n2 = n / 2
  m = i4_log_2 ( n )
  z = - 1.0D+00

  do j = 1, m

    n1 = 2**( m - j + 1 )
    j1 = 2**( j - 1 )

    do l = 1, j1

      is = ( l - 1 ) * n1 + 1
      i1 = 0
      w = z

      do i = is, is + n1 - 1, 2
        a = x(i)
        x(is+i1) = a + x(i+1)
        i1 = i1 + 1
        y(i1) = ( x(i+1) - a ) * w
        w = w * z
      end do

      x(n1/2+is:n1+is-1) = y(1:n1/2)

    end do

  end do

  return
end
