subroutine biharmonic_exact_r1 ( X, Y, a, b, c, d, e, f, g, R )

!*****************************************************************************80
!
!! biharmonic_exact_r1() evaluates exact biharmonic residual for W(X,Y) #1.
!
!  Discussion:
!
!   Note the formula for W:
!
!   W = ( a      * cosh ( g * X ) &
!      + b      * sinh ( g * X ) &
!      + c * X * cosh ( g * X ) &
!      + d * X * sinh ( g * X ) ) & 
!      * &
!      ( e * cos ( g * Y ) &
!      + f * sin ( g * Y ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f, g: parameters.
!
!  Output:
!
!    real R: the residual evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) R
  real ( kind = rk8 ) Wxxxx
  real ( kind = rk8 ) Wxxyy
  real ( kind = rk8 ) Wyyyy
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  Wxxxx =  &
    g**3 * ( e * cos ( g * Y ) + f * sin ( g * Y ) ) &
    * &
    (     a * g     * cosh ( g * X ) &
    +     b * g     * sinh ( g * X ) &
    +     c * g * X * cosh ( g * X ) &
    + 4 * c         * sinh ( g * X ) &
    +     d * g * X * sinh ( g * X ) &
    + 4 * d         * cosh ( g * X ) )

  Wxxyy = &
    -g**3 * ( e * cos ( g * Y ) + f * sin ( g * Y ) ) &
    * &
    (     a * g     * cosh ( g * X ) &
    +     b * g     * sinh ( g * X ) & 
    +     c * g * X * cosh ( g * X ) & 
    + 2 * c         * sinh ( g * X )  &
    +     d * g * X * sinh ( g * X ) & 
    + 2 * d         * cosh ( g * X ) )

  Wyyyy = &
    g**4 * ( e * cos ( g * Y ) + f * sin ( g * Y ) ) &
    * &
    ( a     * cosh ( g * X ) &
    + b     * sinh ( g * X ) &
    + c * X * cosh ( g * X ) &
    + d * X * sinh ( g * X ) )

  R = Wxxxx + 2.0 * Wxxyy + Wyyyy

  return
end
subroutine biharmonic_exact_r2 ( X, Y, a, b, c, d, e, f, g, R )

!*****************************************************************************80
!
!! biharmonic_exact_r2() evaluates exact biharmonic residual for W(X,Y) #2.
!
!  Discussion:
!
!   Note the formula for W:
!
!  W = ( a      * cos ( g * X ) &
!      + b      * sin ( g * X ) &
!      + c * X * cos ( g * X ) &
!      + d * X * sin ( g * X ) ) & 
!      *
!      ( e * cosh ( g * Y ) &
!      + f * sinh ( g * Y ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f, g: parameters.
!
!  Output:
!
!    real R: the residual evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) R
  real ( kind = rk8 ) Wxxxx
  real ( kind = rk8 ) Wxxyy
  real ( kind = rk8 ) Wyyyy
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  Wxxxx =  g**3 * ( e * cosh ( g * Y ) + f * sinh ( g * Y ) ) &
    * &
    (     a * g     * cos ( g * X ) &
    +     b * g     * sin ( g * X ) &
    +     c * g * X * cos ( g * X ) &
    + 4 * c         * sin ( g * X ) &
    +     d * g * X * sin ( g * X ) &
    - 4 * d         * cos ( g * X ) )

  Wxxyy =  - g**3 * ( e * cosh ( g * Y ) + f * sinh ( g * Y ) ) &
    * &
    (     a * g     * cos ( g * X ) &
    +     b * g     * sin ( g * X ) &
    +     c * g * X * cos ( g * X ) &
    + 2 * c         * sin ( g * X ) &
    +     d * g * X * sin ( g * X ) &
    - 2 * d         * cos ( g * X ) )

  Wyyyy =  g**4 * ( e * cosh ( g * Y ) + f * sinh ( g * Y ) ) &
    * &
    ( a     * cos ( g * X ) &
    + b     * sin ( g * X ) &
    + c * X * cos ( g * X ) &
    + d * X * sin ( g * X ) )

  R = Wxxxx + 2.0 * Wxxyy + Wyyyy

  return
end
subroutine biharmonic_exact_r3 ( X, Y, a, b, c, d, e, f, R )

!*****************************************************************************80
!
!! biharmonic_exact_r3() evaluates exact biharmonic residual for W(X,Y) #3.
!
!  Discussion:
!
!   Note the formula for W:
!
!   R = sqrt ( ( X - e )**2 + ( Y - f )**2 )
!
!   W =   a     * R**2 * log ( R ) &
!       + b     * R**2 &
!       + c * log ( R ) &
!       + d
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f: parameters.
!
!  Output:
!
!    real R: the residual evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) R
  real ( kind = rk8 ) Wxxxx
  real ( kind = rk8 ) Wxxyy
  real ( kind = rk8 ) Wyyyy
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  Wxxxx =  2 * ( 8 * a * ( e - X )**4 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - 12 * a * ( e - X )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 ) + 3 * a - 24 * c * ( e - X )**4 &
    / ( ( e - X )**2 + ( f - Y )**2 )**3 + 24 * c * ( e - X )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - 3 * c &
    / ( ( e - X )**2 + ( f - Y )**2 ) ) &
    / ( ( e - X )**2 + ( f - Y )**2 )

  Wxxyy =  2 * ( 8 * a * ( e - X )**2 * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - 2 * a * ( e - X )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 ) - 2 * a * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 ) + a - 24 * c * ( e - X )**2 * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**3 + 4 * c * ( e - X )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 + 4 * c * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - c &
    / ( ( e - X )**2 + ( f - Y )**2 ) ) &
    / ( ( e - X )**2 + ( f - Y )**2 )

  Wyyyy =  2 * ( 8 * a * ( f - Y )**4 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - 12 * a * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 ) + 3 * a - 24 * c * ( f - Y )**4 &
    / ( ( e - X )**2 + ( f - Y )**2 )**3 + 24 * c * ( f - Y )**2 &
    / ( ( e - X )**2 + ( f - Y )**2 )**2 - 3 * c &
    / ( ( e - X )**2 + ( f - Y )**2 ) ) &
    / ( ( e - X )**2 + ( f - Y )**2 )

  R = Wxxxx + 2.0 * Wxxyy + Wyyyy

  return
end
subroutine biharmonic_exact_w1 ( X, Y, a, b, c, d, e, f, g, W )

!*****************************************************************************80
!
!! biharmonic_exact1() evaluates exact biharmonic solution W(X,Y) #1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f, g: parameters.
!
!  Output:
!
!    real W: the solution evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) W
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  W = ( a     * cosh ( g * X ) &
      + b     * sinh ( g * X ) &
      + c * X * cosh ( g * X ) &
      + d * X * sinh ( g * X ) ) & 
      * &
      ( e * cos ( g * Y ) &
      + f * sin ( g * Y ) )

  return
end
subroutine biharmonic_exact_w2 ( X, Y, a, b, c, d, e, f, g, W )

!*****************************************************************************80
!
!! biharmonic_exact_w2() evaluates exact biharmonic solution W(X,Y) #2.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f, g: parameters.
!
!  Output:
!
!    real W: the solution evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) g
  real ( kind = rk8 ) W
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  W = ( a     * cos ( g * X ) &
      + b     * sin ( g * X ) &
      + c * X * cos ( g * X ) &
      + d * X * sin ( g * X ) ) & 
      * &
      ( e * cosh ( g * Y ) &
      + f * sinh ( g * Y ) )

  return
end
subroutine biharmonic_exact_w3 ( X, Y, a, b, c, d, e, f, W )

!*****************************************************************************80
!
!! biharmonic_exact_w3() evaluates exact biharmonic solution W(X,Y) #3.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 August 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real X, Y: the coordinates of the points.
!
!    real a, b, c, d, e, f: parameters.
!
!  Output:
!
!    real W: the solution evaluated at the points (X,Y).
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) d
  real ( kind = rk8 ) e
  real ( kind = rk8 ) f
  real ( kind = rk8 ) R
  real ( kind = rk8 ) W
  real ( kind = rk8 ) X
  real ( kind = rk8 ) Y

  R = sqrt ( ( X - e )**2 + ( Y - f )**2 )

  W =   a     * R**2 * log ( R ) &
      + b     * R**2 &
      + c * log ( R ) &
      + d

  return
end
