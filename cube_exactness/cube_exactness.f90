subroutine legendre_3d_exactness ( a, b, n, x, y, z, w, t )

!*****************************************************************************80
!
!! LEGENDRE_3D_EXACTNESS: monomial exactness for the 3D Legendre integral.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(3), the lower limits of integration.
!
!    Input, real ( kind = rk ) B(3), the upper limits of integration.
!
!    Input, integer N, the number of points in the rule.
!
!    Input, real ( kind = rk ) X(N), Y(N), Z(N), the quadrature points.
!
!    Input, real ( kind = rk ) W(N), the quadrature weights.
!
!    Input, integer T, the maximum total degree.
!    0 <= T.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(3)
  real ( kind = rk ) e
  integer i
  integer j
  integer k
  integer p(3)
  real ( kind = rk ) q
  real ( kind = rk ) s
  integer t
  integer tt
  real ( kind = rk ) v(n)
  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)
  real ( kind = rk ) y(n)
  real ( kind = rk ) z(n)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Quadrature rule for the 3D Legendre integral.'
  write ( *, '(a,i3)' ) '  Number of points in rule is ', n
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   D   I       J       K          Relative Error'

  do tt = 0, t

    write ( *, '(2x,i2)' )  tt

    do k = 0, tt

      do j = 0, tt - k

        i = tt - j - k

        p(1) = i
        p(2) = j
        p(3) = k

        call legendre_3d_monomial_integral ( a, b, p, s )

        v(1:n) = x(1:n) ** p(1) * y(1:n) ** p(2) * z(1:n) ** p(3)

        q = dot_product ( w, v )

        if ( s == 0.0D+00 ) then
          e = abs ( q )
        else
          e = abs ( q - s ) / abs ( s )
        end if

        write ( *, '(2x,i6,2x,i6,2x,i6,2x,f24.16)' ) p(1:3), e

      end do

    end do

  end do

  return
end
subroutine legendre_3d_monomial_integral ( a, b, p, value )

!*****************************************************************************80
!
!! LEGENDRE_3D_MONOMIAL_INTEGRAL the Legendre integral of a monomial.
!
!  Discussion:
!
!    The Legendre integral to be evaluated has the form
!
!      I(f) = integral ( z1 <= z <= z2 )
!             integral ( y1 <= y <= y2 ) 
!             integral ( x1 <= x <= x2 ) x^i y^j z^k dx dy dz
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    15 August 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(3), the lower limits of integration.
!
!    Input, real ( kind = rk ) B(3), the upper limits of integration.
!
!    Input, integer P(3), the exponents of X and Y.
!
!    Output, real ( kind = rk ) VALUE, the value of the exact integral.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(3)
  real ( kind = rk ) b(3)
  integer p(3)
  real ( kind = rk ) value

  value = ( b(1) ** ( p(1) + 1 ) - a(1) ** ( p(1) + 1 ) ) &
        / real ( p(1) + 1, kind = rk ) &
        * ( b(2) ** ( p(2) + 1 ) - a(2) ** ( p(2) + 1 ) ) &
        / real ( p(2) + 1, kind = rk ) &
        * ( b(3) ** ( p(3) + 1 ) - a(3) ** ( p(3) + 1 ) ) &
        / real ( p(3) + 1, kind = rk )

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
!  Parameters:
!
!    None
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
