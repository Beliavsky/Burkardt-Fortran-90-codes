function c8_abs ( z )

!*****************************************************************************80
!
!! c8_abs() evaluates the absolute value of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the absolute value of a C8 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    real ( kind = rk ) C8_ABS, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_abs
  complex ( kind = ck ) z

  c8_abs = sqrt ( ( real ( z, kind = rk ) )**2 &
                + ( aimag ( z ) )**2 )

  return
end
function c8_acos ( z )

!*****************************************************************************80
!
!! c8_acos() evaluates the inverse cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse cosine of a C8.
!
!    Here we use the relationship:
!
!       C8_ACOS ( Z ) = pi/2 - C8_ASIN ( Z ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ACOS, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8_acos
  complex ( kind = ck ) c8_asin
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  complex ( kind = ck ) z

  c8_acos = 0.5D+00 * r8_pi - c8_asin ( z )

  return
end
function c8_acosh ( z )

!*****************************************************************************80
!
!! c8_acosh() evaluates the inverse hyperbolic cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse hyperbolic cosine of a C8.
!
!    Here we use the relationship:
!
!      C8_ACOSH ( Z ) = i * C8_ACOS ( Z ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ACOSH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_acos
  complex ( kind = ck ) c8_acosh
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_acosh = c8_i * c8_acos ( z )

  return
end
function c8_add ( z1, z2 )

!*****************************************************************************80
!
!! c8_add() adds two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports addition of C8's with the "+" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z1, Z2, the values to add.
!
!  Output:
!
!    complex ( kind = ck ) C8_ADD, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_add
  complex ( kind = ck ) z1
  complex ( kind = ck ) z2

  c8_add = z1 + z2

  return
end
function c8_arg ( x )

!*****************************************************************************80
!
!! c8_arg() returns the argument of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    The value returned by this function is always between 0 and 2*PI.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the complex number.
!
!  Output:
!
!    real ( kind = rk ) C8_ARG, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_arg
  complex ( kind = ck ) x
  real ( kind = rk ) r8_atan

  if ( aimag ( x )           == 0.0D+00 .and. &
       real  ( x, kind = rk ) == 0.0D+00 ) then

    c8_arg = 0.0D+00

  else

    c8_arg = r8_atan ( aimag ( x ), real ( x ) )

  end if

  return
end
function c8_asin ( z )

!*****************************************************************************80
!
!! c8_asin() evaluates the inverse sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse sine of a C8.
!
!    Here we use the relationship:
!
!      C8_ASIN ( Z ) = - i * log ( i * z + sqrt ( 1 - z * z ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ASIN, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_asin
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_asin = - c8_i * log ( c8_i * z + sqrt ( 1.0D+00 - z * z ) )

  return
end
function c8_asinh ( z )

!*****************************************************************************80
!
!! c8_asinh() evaluates the inverse hyperbolic sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse hyperbolic sine of a C8.
!
!    Here we use the relationship:
!
!      C8_ASINH ( Z ) = - i * C8_ASIN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ASINH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_asin
  complex ( kind = ck ) c8_asinh
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_asinh = - c8_i * c8_asin ( c8_i * z )

  return
end
function c8_atan ( z )

!*****************************************************************************80
!
!! c8_atan() evaluates the inverse tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse tangent of a C8.
!
!    Here we use the relationship:
!
!      C8_ATAN ( Z ) = ( i / 2 ) * log ( ( 1 - i * z ) / ( 1 + i * z ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ATAN, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) arg
  complex ( kind = ck ) c8_atan
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  arg = ( 1.0D+00 - c8_i * z ) / ( 1.0D+00 + c8_i * z )

  c8_atan = 0.5D+00 * c8_i * log ( arg )

  return
end
function c8_atanh ( z )

!*****************************************************************************80
!
!! c8_atanh() evaluates the inverse hyperbolic tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the inverse hyperbolic tangent of a C8.
!
!    Here we use the relationship:
!
!      C8_ATANH ( Z ) = - i * C8_ATAN ( i * Z ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_ATANH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_atan
  complex ( kind = ck ) c8_atanh
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_atanh = - c8_i * c8_atan ( c8_i * z )

  return
end
function c8_conj ( z )

!*****************************************************************************80
!
!! c8_conj() evaluates the conjugate of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the conjugate of a C8 with the CONJG function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_CONJ, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  complex ( kind = ck ) c8_conj
  real ( kind = rk ) c8_imag
  real ( kind = rk ) c8_real
  complex ( kind = ck ) z

  a = c8_real ( z )
  b = c8_imag ( z )

  c8_conj = cmplx ( a, -b, kind = ck )

  return
end
function c8_copy ( z )

!*****************************************************************************80
!
!! c8_copy() copies a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the copy of a C8 with the "=" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_COPY, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_copy
  complex ( kind = ck ) z

  c8_copy = z

  return
end
function c8_cos ( z )

!*****************************************************************************80
!
!! c8_cos() evaluates the cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the cosine of a C8 with the COS function.
!
!    We use the relationship:
!
!      C8_COS ( C ) = ( C8_EXP ( i * C ) + C8_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_COS, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_cos
  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_cos = ( c8_exp ( c8_i * z ) + c8_exp ( - c8_i * z ) ) / 2.0D+00

  return
end
function c8_cosh ( z )

!*****************************************************************************80
!
!! c8_cosh() evaluates the hyperbolic cosine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the hyperbolic cosine of a C8.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_COSH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_cosh
  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) z

  c8_cosh = ( c8_exp ( z ) + c8_exp ( - z ) ) / 2.0D+00

  return
end
function c8_cube_root ( x )

!*****************************************************************************80
!
!! c8_cube_root() returns the principal cube root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the cube root of a C8 through the 
!    " ** ( 1.0 / 3.0 )" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_CUBE_ROOT, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg
  real ( kind = rk ) c8_arg
  complex ( kind = ck ) c8_cube_root
  real ( kind = rk ) c8_mag
  real ( kind = rk ) mag
  complex ( kind = ck ) x

  arg = c8_arg ( x )
  mag = c8_mag ( x )

  if ( mag == 0.0D+00 ) then

    c8_cube_root = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  else

    c8_cube_root = mag**( 1.0D+00 / 3.0D+00 ) &
      * cmplx ( cos ( arg / 3.0D+00 ), &
                sin ( arg / 3.0D+00 ), kind = ck )

  end if

  return
end
function c8_div ( z1, z2 )

!*****************************************************************************80
!
!! c8_div() divides two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports division of C8's with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z1, Z2, the arguments.
!
!  Output:
!
!    complex ( kind = ck ) C8_DIV, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  complex ( kind = ck ) c8_div
  real ( kind = rk ) c8_imag
  real ( kind = rk ) c8_real
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ) g
  complex ( kind = ck ) z1
  complex ( kind = ck ) z2

  a = c8_real ( z1 )
  b = c8_imag ( z1 )
  c = c8_real ( z2 )
  d = c8_imag ( z2 )

  e = c * c + d * d

  f = ( a * c + b * d ) / e
  g = ( b * c - a * d ) / e

  c8_div = cmplx ( f, g, kind = ck )

  return
end
function c8_div_r8 ( z1, r )

!*****************************************************************************80
!
!! c8_div_r8() divides a C8 by an R8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!    An R8 is a real ( kind = rk ) value.
!
!    Fortran90 supports division of a C8 by an R8 with the "/" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z1, the value to be divided.
!
!    real ( kind = rk ) R, the divisor.
!
!  Output:
!
!    complex ( kind = ck ) C8_DIV_R8, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  complex ( kind = ck ) c8_div_r8
  real ( kind = rk ) c8_imag
  real ( kind = rk ) c8_real
  real ( kind = rk ) r
  complex ( kind = ck ) z1

  a = c8_real ( z1 ) / r
  b = c8_imag ( z1 ) / r

  c8_div_r8 = cmplx ( a, b, kind = ck )

  return
end
function c8_exp ( z )

!*****************************************************************************80
!
!! c8_exp() evaluates the exponential of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the exponential of a C8 with the EXP function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_EXP, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8_exp
  real ( kind = rk ) c8_imag
  real ( kind = rk ) c8_real
  complex ( kind = ck ) z
  real ( kind = rk ) zi
  real ( kind = rk ) zr

  zr = c8_real ( z )
  zi = c8_imag ( z )

  c8_exp = exp ( zr ) * cmplx ( cos ( zi ), sin ( zi ), kind = ck )

  return
end
subroutine c8_fake_use ( x )

!*****************************************************************************80
!
!! c8_fake_use() pretends to use a variable.
!
!  Discussion:
!
!    Some compilers will issue a warning if a variable is unused.
!    Sometimes there's a good reason to include a variable in a program,
!    but not to use it.  Calling this function with that variable as
!    the argument will shut the compiler up.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the variable to be "used".
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) x

  if ( x /= x ) then
    write ( *, '(a)' ) '  c8_fake_use: variable is NAN.'
  end if

  return
end
function c8_i ( )

!*****************************************************************************80
!
!! c8_i() returns the imaginary unit, i as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    complex ( kind = ck ) C8_I, the value of complex i.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_i

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  return
end
function c8_imag ( z )

!*****************************************************************************80
!
!! c8_imag() evaluates the imaginary part of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the imaginary part of a C8 with the AIMAG function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    real ( kind = rk ) C8_IMAG, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_imag
  complex ( kind = ck ) z

  c8_imag = aimag ( z )

  return
end
function c8_inv ( z )

!*****************************************************************************80
!
!! c8_inv() evaluates the inverse of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the inverse of a C8 with the "1/" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_INV, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8_inv
  complex ( kind = ck ) z
  real ( kind = rk ) z_imag
  real ( kind = rk ) z_norm
  real ( kind = rk ) z_real

  z_real = real ( z, kind = rk )
  z_imag = aimag ( z )

  z_norm = sqrt ( z_real * z_real + z_imag * z_imag )

  c8_inv = cmplx ( z_real, - z_imag, kind = ck ) / z_norm / z_norm

  return
end
function c8_le_l1 ( x, y )

!*****************************************************************************80
!
!! c8_le_l1() := X <= Y for C8 values, and the L1 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, Y, the values to be compared.
!
!  Output:
!
!    logical C8_LE_L1, is TRUE if X <= Y.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c8_le_l1
  complex ( kind = ck ) x
  complex ( kind = ck ) y

  if ( abs ( real ( x, kind = rk ) ) + abs ( aimag ( x ) ) <= &
       abs ( real ( y, kind = rk ) ) + abs ( aimag ( y ) ) ) then
    c8_le_l1 = .true.
  else
    c8_le_l1 = .false.
  end if

  return
end
function c8_le_l2 ( x, y )

!*****************************************************************************80
!
!! c8_le_l2() := X <= Y for C8 values, and the L2 norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag (X) )^2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, Y, the values to be compared.
!
!  Output:
!
!    logical C8_LE_L2, is TRUE if X <= Y.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c8_le_l2
  complex ( kind = ck ) x
  complex ( kind = ck ) y

  if ( ( real ( x, kind = rk ) )**2 + ( aimag ( x ) )**2 <= &
       ( real ( y, kind = rk ) )**2 + ( aimag ( y ) )**2 ) then
    c8_le_l2 = .true.
  else
    c8_le_l2 = .false.
  end if

  return
end
function c8_le_li ( x, y )

!*****************************************************************************80
!
!! c8_le_li() := X <= Y for C8 values, and the L Infinity norm.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    The L Infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, Y, the values to be compared.
!
!  Output:
!
!    logical C8_LE_LI, is TRUE if X <= Y.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  logical c8_le_li
  complex ( kind = ck ) x
  complex ( kind = ck ) y

  if ( max ( abs ( real ( x, kind = rk ) ), abs ( aimag ( x ) ) ) <= &
       max ( abs ( real ( y, kind = rk ) ), abs ( aimag ( y ) ) ) ) then
    c8_le_li = .true.
  else
    c8_le_li = .false.
  end if

  return
end
function c8_log ( z )

!*****************************************************************************80
!
!! c8_log() evaluates the logarithm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the logarithm of a C8 with the LOG function.
!
!    Here we use the relationship:
!
!      C8_LOG ( Z ) = LOG ( MAG ( Z ) ) + i * ARG ( Z )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_LOG, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg
  real ( kind = rk ) c8_arg
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) c8_log
  real ( kind = rk ) c8_mag
  real ( kind = rk ) mag
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  arg = c8_arg ( z )
  mag = c8_mag ( z )

  c8_log = log ( mag ) + c8_i * arg
 
  return
end
function c8_mag ( x )

!*****************************************************************************80
!
!! c8_mag() returns the magnitude of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the magnitude of a C8 with the ABS function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the argument.
!
!  Output:
!
!    real ( kind = rk ) C8_MAG, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_mag
  complex ( kind = ck ) x

  c8_mag = sqrt ( ( real ( x, kind = rk ) )**2 + ( aimag ( x ) )**2 )

  return
end
function c8_mul ( z1, z2 )

!*****************************************************************************80
!
!! c8_mul() multiplies two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports multiplication of C8's with the "*" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z1, Z2, the values to multiply.
!
!  Output:
!
!    complex ( kind = ck ) C8_MUL, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_mul
  complex ( kind = ck ) z1
  complex ( kind = ck ) z2

  c8_mul = z1 * z2

  return
end
function c8_neg ( c1 )

!*****************************************************************************80
!
!! c8_neg() returns the negative of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports negation of a C8 with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) C1, the value to be negated.
!
!  Output:
!
!    complex ( kind = ck ) C8_NEG, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c1
  complex ( kind = ck ) c8_neg

  c8_neg = - c1

  return
end
function c8_nint ( c1 )

!*****************************************************************************80
!
!! c8_nint() returns the nearest complex integer of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) C1, the value to be NINT'ed.
!
!  Output:
!
!    complex ( kind = ck ) C8_NINT, the NINT'ed value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c1
  complex ( kind = ck ) c8_nint
  real ( kind = rk ) r
  real ( kind = rk ) r_min
  real ( kind = rk ) r8_floor
  real ( kind = rk ) x
  real ( kind = rk ) x_min
  real ( kind = rk ) xc
  real ( kind = rk ) y
  real ( kind = rk ) y_min
  real ( kind = rk ) yc

  xc = real ( c1 )
  yc = imag ( c1 )
!
!  Lower left.
!
  x = r8_floor ( real ( c1 ) )
  y = r8_floor ( imag ( c1 ) )
  r = ( x - xc )**2 + ( y - yc )**2
  r_min = r
  x_min = x
  y_min = y
!
!  Lower right.
!
  x = r8_floor ( real ( c1 ) ) + 1.0D+00
  y = r8_floor ( imag ( c1 ) )
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if
!
!  Upper right.
!
  x = r8_floor ( real ( c1 ) ) + 1.0D+00
  y = r8_floor ( imag ( c1 ) ) + 1.0D+00
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if
!
!  Upper left.
!
  x = r8_floor ( real ( c1 ) )
  y = r8_floor ( imag ( c1 ) ) + 1.0D+00
  r = ( x - xc )**2 + ( y - yc )**2
  if ( r < r_min ) then
    r_min = r
    x_min = x
    y_min = y
  end if

  c8_nint = cmplx ( x_min, y_min, kind = ck )

  return
end
function c8_norm_l1 ( x )

!*****************************************************************************80
!
!! c8_norm_l1() evaluates the L1 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Numbers of equal norm lie along diamonds centered at (0,0).
!
!    The L1 norm can be defined here as:
!
!      C8_NORM_L1(X) = abs ( real (X) ) + abs ( imag (X) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the value whose norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8_NORM_L1, the norm of X.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_norm_l1
  complex ( kind = ck ) x

  c8_norm_l1 = abs ( real ( x, kind = rk ) ) + abs ( aimag ( x ) )

  return
end
function c8_norm_l2 ( x )

!*****************************************************************************80
!
!! c8_norm_l2() evaluates the L2 norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Numbers of equal norm lie on circles centered at (0,0).
!
!    The L2 norm can be defined here as:
!
!      C8_NORM_L2(X) = sqrt ( ( real (X) )^2 + ( imag ( X ) )^2 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the value whose norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8_NORM_L2, the 2-norm of X.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_norm_l2
  complex ( kind = ck ) x

  c8_norm_l2 = sqrt ( ( real ( x, kind = rk ) )**2 &
                   + ( aimag ( x ) )**2 )

  return
end
function c8_norm_li ( x )

!*****************************************************************************80
!
!! c8_norm_li() evaluates the L-infinity norm of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Numbers of equal norm lie along squares whose centers are at (0,0).
!
!    The L-infinity norm can be defined here as:
!
!      C8_NORM_LI(X) = max ( abs ( real (X) ), abs ( imag (X) ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the value whose norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8_NORM_LI, the infinity norm of X.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_norm_li
  complex ( kind = ck ) x

  c8_norm_li = max ( abs ( real ( x, kind = rk ) ), abs ( aimag ( x ) ) )

  return
end
function c8_normal_01 ( )

!*****************************************************************************80
!
!! c8_normal_01() returns a unit pseudonormal C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    complex ( kind = ck ) C8_NORMAL_01, a unit pseudornormal value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8_normal_01
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) v1
  real ( kind = rk ) v2
  real ( kind = rk ) x_c
  real ( kind = rk ) x_r

  call random_number ( harvest = v1 )
  call random_number ( harvest = v2 )

  x_r = sqrt ( - 2.0D+00 * log ( v1 ) ) * cos ( 2.0D+00 * r8_pi * v2 )
  x_c = sqrt ( - 2.0D+00 * log ( v1 ) ) * sin ( 2.0D+00 * r8_pi * v2 )

  c8_normal_01 = cmplx ( x_r, x_c, kind = ck )

  return
end
function c8_one ( )

!*****************************************************************************80
!
!! c8_one() returns the value of 1 as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    complex ( kind = ck ) C8_ONE, the value of complex 1.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_one

  c8_one = cmplx ( 1.0D+00, 0.0D+00, kind = ck )

  return
end
subroutine c8_print ( a, title )

!*****************************************************************************80
!
!! c8_print() prints a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) A, the value to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) a
  character ( len = * ) title

  if ( 0 < len_trim ( title ) ) then
    write ( *, '(a,2x,a,g14.6,a,g14.6,a)' ) &
      trim ( title ), '(', real ( a ), ',', imag ( a ), ')'
  else
    write ( *, '(a,g14.6,a,g14.6,a)' ) &
      '(', real ( a ), ',', imag ( a ), ')'
  end if

  return
end
function c8_real ( z )

!*****************************************************************************80
!
!! c8_real() evaluates the real part of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the real part of a C8 with the REAL function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    real ( kind = rk ) C8_REAL, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_real
  complex ( kind = ck ) z

  c8_real = real ( z )

  return
end
function c8_sin ( z )

!*****************************************************************************80
!
!! c8_sin() evaluates the sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the sine of a C8 with the SIN function.
!
!    We use the relationship:
!
!      C8_SIN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_SIN, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) c8_sin
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_sin = - c8_i * ( c8_exp ( c8_i * z ) - c8_exp ( - c8_i * z ) ) / 2.0D+00

  return
end
function c8_sinh ( z )

!*****************************************************************************80
!
!! c8_sinh() evaluates the hyperbolic sine of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the hyperbolic sine of a C8.
!
!    We use the relationship:
!
!      C8_SINH ( C ) = ( C8_EXP ( C ) - C8_EXP ( - C ) ) / 2
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_SINH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) c8_sinh
  complex ( kind = ck ) z

  c8_sinh = ( c8_exp ( z ) - c8_exp ( - z ) ) / 2.0D+00

  return
end
function c8_sqrt ( x )

!*****************************************************************************80
!
!! c8_sqrt() returns the principal square root of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 supports the square root of a C8 with the SQRT function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_SQRT, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) arg
  real ( kind = rk ) c8_arg
  real ( kind = rk ) c8_mag
  complex ( kind = ck ) c8_sqrt
  real ( kind = rk ) mag
  complex ( kind = ck ) x

  arg = c8_arg ( x )
  mag = c8_mag ( x )

  if ( mag == 0.0D+00 ) then

    c8_sqrt = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  else

    c8_sqrt = sqrt ( mag ) &
      * cmplx ( cos ( arg / 2.0D+00 ), &
                sin ( arg / 2.0D+00 ), kind = ck )

  end if

  return
end
function c8_sub ( z1, z2 )

!*****************************************************************************80
!
!! c8_sub() subtracts two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 directly supports C8 subtraction with the "-" operator.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z1, Z2, the values to subtract.
!
!  Output:
!
!    complex ( kind = ck ) C8_SUB, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_sub
  complex ( kind = ck ) z1
  complex ( kind = ck ) z2

  c8_sub = z1 - z2

  return
end
subroutine c8_swap ( x, y )

!*****************************************************************************80
!
!! c8_swap() swaps two C8's.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) X, Y: values to interchange.
!
!  Output:
!
!    complex ( kind = ck ) X, Y: the interchanged values.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) x
  complex ( kind = ck ) y
  complex ( kind = ck ) z

  z = x
  x = y
  y = z

  return
end
function c8_tan ( z )

!*****************************************************************************80
!
!! c8_tan() evaluates the tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the tangent of a C8.
!
!    We use the relationship:
!
!      C8_TAN ( C ) = - i * ( C8_EXP ( i * C ) - C8_EXP ( - i * C ) ) 
!                         / ( C8_EXP ( I * C ) + C8_EXP ( - i * C ) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_TAN, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) c8_i
  complex ( kind = ck ) c8_tan
  complex ( kind = ck ) z

  c8_i = cmplx ( 0.0D+00, 1.0D+00, kind = ck )

  c8_tan =  - c8_i * ( c8_exp ( c8_i * z ) - c8_exp ( - c8_i * z ) ) &
         /           ( c8_exp ( c8_i * z ) + c8_exp ( - c8_i * z ) )

  return
end
function c8_tanh ( z )

!*****************************************************************************80
!
!! c8_tanh() evaluates the hyperbolic tangent of a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    Fortran90 does not support the hyperbolic tangent of a C8.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    complex ( kind = ck ) C8_TANH, the function value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_exp
  complex ( kind = ck ) c8_tanh
  complex ( kind = ck ) z

  c8_tanh = ( c8_exp ( z ) - c8_exp ( - z ) ) &
          / ( c8_exp ( z ) + c8_exp ( - z ) )

  return
end
subroutine c8_to_cartesian ( z, x, y )

!*****************************************************************************80
!
!! c8_to_cartesian() converts a C8 to Cartesian form.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    real ( kind = rk ) X, Y, the Cartesian form.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x
  real ( kind = rk ) y
  complex ( kind = ck ) z

  x = real ( z )
  y = aimag ( z )

  return
end
subroutine c8_to_polar ( z, r, theta )

!*****************************************************************************80
!
!! c8_to_polar() converts a C8 to polar form.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    complex ( kind = ck ) Z, the argument.
!
!  Output:
!
!    real ( kind = rk ) R, THETA, the polar form.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c8_arg
  real ( kind = rk ) c8_mag
  real ( kind = rk ) r
  real ( kind = rk ) theta
  complex ( kind = ck ) z

  r = c8_mag ( z )
  theta = c8_arg ( z )

  return
end
function c8_uniform_01 ( )

!*****************************************************************************80
!
!! c8_uniform_01() returns a unit pseudorandom C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!    The angle should be uniformly distributed between 0 and 2 * PI,
!    the square root of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    complex ( kind = ck ) C8_UNIFORM_01, a pseudorandom complex value.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) c8_uniform_01
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  call random_number ( harvest = r )
  call random_number ( harvest = theta )
  theta = 2.0D+00 * r8_pi * theta

  c8_uniform_01 = r * cmplx ( cos ( theta ), sin ( theta ), kind = ck )

  return
end
function c8_zero ( )

!*****************************************************************************80
!
!! c8_zero() returns the value of 0 as a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    complex ( kind = ck ) C8_ZERO, the value of complex 0.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) c8_zero

  c8_zero = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  return
end
subroutine c8mat_add ( m, n, alpha, a, beta, b, c )

!*****************************************************************************80
!
!! c8mat_add() combines two C8MAT's with complex scalar factors.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    complex ( kind = ck ) ALPHA, the first scale factor.
!
!    complex ( kind = ck ) A(M,N), the first matrix.
!
!    complex ( kind = ck ) BETA, the second scale factor.
!
!    complex ( kind = ck ) B(M,N), the second matrix.
!
!  Output:
!
!    complex ( kind = ck ) C(M,N), the result matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  complex ( kind = ck ) alpha
  complex ( kind = ck ) b(m,n)
  complex ( kind = ck ) beta
  complex ( kind = ck ) c(m,n)

  c(1:m,1:n) = alpha * a(1:m,1:n) + beta * b(1:m,1:n)

  return
end
subroutine c8mat_add_r8 ( m, n, alpha, a, beta, b, c )

!*****************************************************************************80
!
!! c8mat_add_r8() combines two C8MAT's with real scalar factors.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) ALPHA, the first scale factor.
!
!    complex ( kind = ck ) A(M,N), the first matrix.
!
!    real ( kind = rk ) BETA, the second scale factor.
!
!    complex ( kind = ck ) B(M,N), the second matrix.
!
!  Output:
!
!    complex ( kind = ck ) C(M,N), the result matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  real ( kind = rk ) alpha
  complex ( kind = ck ) b(m,n)
  real ( kind = rk ) beta
  complex ( kind = ck ) c(m,n)

  c(1:m,1:n) = alpha * a(1:m,1:n) + beta * b(1:m,1:n)

  return
end
subroutine c8mat_copy ( m, n, a, b )

!*****************************************************************************80
!
!! c8mat_copy() copies a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    complex ( kind = ck ) A(M,N), the matrix.
!
!  Output:
!
!    complex ( kind = ck ) B(M,N), the copied matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  complex ( kind = ck ) b(m,n)

  b(1:m,1:n) = a(1:m,1:n)

  return
end
subroutine c8mat_fss ( n, a, nb, b, info )

!*****************************************************************************80
!
!! c8mat_fss() factors and solves a system with multiple right hand sides.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!    This routine does not save the LU factors of the matrix, and hence cannot
!    be used to efficiently solve multiple linear systems, or even to
!    factor A at one time, and solve a single linear system at a later time.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!    N must be positive.
!
!    complex ( kind = ck ) A(N,N): the coefficient matrix of the linear system.
!
!    integer NB, the number of right hand sides.
!
!    complex ( kind = ck ) B(N,NB): the right hand sides of the linear system.
!
!    integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
!  Output:
!
!    complex ( kind = ck ) A(N,N): the unit upper triangular U factor
!    of an LU factorization of the original coefficient matrix.
!
!    complex ( kind = ck ) B(N,NB): the solutions of the linear systems.
!
!    integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nb

  complex ( kind = ck ) a(n,n)
  complex ( kind = ck ) b(n,nb)
  integer i
  integer info
  integer ipiv
  integer j
  integer jcol
  real ( kind = rk ) piv
  complex ( kind = ck ) row(n)
  complex ( kind = ck ) t(nb)
  complex ( kind = ck ) temp

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a(i,jcol) ) ) then
        piv = abs ( a(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C8MAT_FSS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a(jcol,1:n)
      a(jcol,1:n) = a(ipiv,1:n)
      a(ipiv,1:n) = row(1:n)

      t(1:nb)      = b(jcol,1:nb)
      b(jcol,1:nb) = b(ipiv,1:nb)
      b(ipiv,1:nb) = t(1:nb)

    end if
!
!  Scale the pivot row.
!
    a(jcol,jcol+1:n) = a(jcol,jcol+1:n) / a(jcol,jcol)
    b(jcol,1:nb) = b(jcol,1:nb) / a(jcol,jcol)
    a(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a(i,jcol) /= 0.0D+00 ) then
        temp = - a(i,jcol)
        a(i,jcol) = 0.0D+00
        a(i,jcol+1:n) = a(i,jcol+1:n) + temp * a(jcol,jcol+1:n)
        b(i,1:nb) = b(i,1:nb) + temp * b(jcol,1:nb)
      end if
    end do

  end do
!
!  Back solve.
!
  do j = 1, nb
    do jcol = n, 2, -1
      b(1:jcol-1,j) = b(1:jcol-1,j) - a(1:jcol-1,jcol) * b(jcol,j)
    end do
  end do

  return
end
subroutine c8mat_identity ( n, a )

!*****************************************************************************80
!
!! c8mat_identity() sets a C8MAT to the identity.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    complex ( kind = ck ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n,n)
  integer i

  a(1:n,1:n) = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  do i = 1, n
    a(i,i) = cmplx ( 1.0D+00, 0.0D+00, kind = ck )
  end do

  return
end
subroutine c8mat_indicator ( m, n, a )

!*****************************************************************************80
!
!! c8mat_indicator() returns the C8MAT indicator matrix.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!  Output:
!
!    complex ( kind = ck ) A(M,N), the matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  integer i
  integer j

  do j = 1, n
    do i = 1, m
      a(i,j) = cmplx ( i, j, kind = ck )
    end do
  end do

  return
end
subroutine c8mat_minvm ( n1, n2, a, b, c )

!*****************************************************************************80
!
!! c8mat_minvm() computes inverse(A) * B for C8MAT's.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N1, N2, the order of the matrices.
!
!    complex ( kind = ck ) A(N1,N1), B(N1,N2), the matrices.
!
!  Output:
!
!    complex ( kind = ck ) C(N1,N2), the result, C = inverse(A) * B.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n1
  integer n2

  complex ( kind = ck ) a(n1,n1)
  complex ( kind = ck ) alu(n1,n1)
  complex ( kind = ck ) b(n1,n2)
  complex ( kind = ck ) c(n1,n2)
  integer info

  alu(1:n1,1:n1) = a(1:n1,1:n1)
  c(1:n1,1:n2) = b(1:n1,1:n2)

  call c8mat_fss ( n1, alu, n2, c, info )
 
  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'C8MAT_MINVM - Fatal error!'
    write ( *, '(a)' ) '  The matrix A was numerically singular.'
    stop 1
  end if

  return
end
subroutine c8mat_mm ( n1, n2, n3, a, b, c )

!*****************************************************************************80
!
!! c8mat_mm() multiplies two C8MAT's.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N1, N2, N3, define the orders of
!    the matrices.
!
!    complex ( kind = ck ) A(N1,N2), B(N2,N3), the matrix factors.
!
!  Output:
!
!    complex ( kind = ck ) C(N1,N3), the product matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n1
  integer n2
  integer n3

  complex ( kind = ck ) a(n1,n2)
  complex ( kind = ck ) b(n2,n3)
  complex ( kind = ck ) c(n1,n3)

  c(1:n1,1:n3) = matmul ( a(1:n1,1:n2), b(1:n2,1:n3) )

  return
end
subroutine c8mat_nint ( m, n, a )

!*****************************************************************************80
!
!! c8mat_nint() rounds the entries of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns of A.
!
!    complex ( kind = ck ) A(M,N), the matrix to be NINT'ed.
!
!  Output:
!
!    complex ( kind = ck ) A(M,N), the NINT'ed matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  complex ( kind = ck ) c8_nint
  integer i
  integer j

  do j = 1, n
    do i = 1, m
      a(i,j) = c8_nint ( a(i,j) )
    end do
  end do

  return
end
function c8mat_norm_fro ( m, n, a )

!*****************************************************************************80
!
!! c8mat_norm_fro() returns the Frobenius norm of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!    The Frobenius norm is defined as
!
!      C8MAT_NORM_FRO = sqrt (
!        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J) * A(I,J) )
!
!    The matrix Frobenius norm is not derived from a vector norm, but
!    is compatible with the vector L2 norm, so that:
!
!      c8vec_norm_l2 ( A * x ) <= c8mat_norm_fro ( A ) * c8vec_norm_l2 ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    complex ( kind = ck ) A(M,N), the matrix whose Frobenius
!    norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8MAT_NORM_FRO, the Frobenius norm of A.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  real ( kind = rk ) c8mat_norm_fro
  real ( kind = rk ) value

  value = &
    sqrt & 
    ( &
      sum &
      ( &
        ( &
          abs &
          ( &
            a(1:m,1:n) &
          ) &
        ) ** 2 &
      ) &
    )

  c8mat_norm_fro = value

  return
end
function c8mat_norm_l1 ( m, n, a )

!*****************************************************************************80
!
!! c8mat_norm_l1() returns the matrix L1 norm of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!    The matrix L1 norm is defined as:
!
!      C8MAT_NORM_L1 = max ( 1 <= J <= N )
!        sum ( 1 <= I <= M ) abs ( A(I,J) ).
!
!    The matrix L1 norm is derived from the vector L1 norm, and
!    satisifies:
!
!      c8vec_norm_l1 ( A * x ) <= c8mat_norm_l1 ( A ) * c8vec_norm_l1 ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    complex ( kind = ck ) A(M,N), the matrix whose L1 norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8MAT_NORM_L1, the L1 norm of A.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  real ( kind = rk ) c8mat_norm_l1
  real ( kind = rk ) col_sum
  integer j

  c8mat_norm_l1 = 0.0D+00

  do j = 1, n
    col_sum = sum ( abs ( a(1:m,j) ) )
    c8mat_norm_l1 = max ( c8mat_norm_l1, col_sum )
  end do

  return
end
function c8mat_norm_li ( m, n, a )

!*****************************************************************************80
!
!! c8mat_norm_li() returns the matrix L-oo norm of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is an MxN array of C8's, stored by (I,J) -> [I+J*M].
!
!    The matrix L-oo norm is defined as:
!
!      C8MAT_NORM_LI =  max ( 1 <= I <= M ) sum ( 1 <= J <= N ) abs ( A(I,J) ).
!
!    The matrix L-oo norm is derived from the vector L-oo norm,
!    and satisifies:
!
!      c8vec_norm_li ( A * x ) <= c8mat_norm_li ( A ) * c8vec_norm_li ( x ).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the number of rows in A.
!
!    integer N, the number of columns in A.
!
!    complex ( kind = ck ) A(M,N), the matrix whose L-oo
!    norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8MAT_NORM_LI, the L-oo norm of A.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  real ( kind = rk ) c8mat_norm_li
  integer i
  real ( kind = rk ) row_sum

  c8mat_norm_li = 0.0D+00

  do i = 1, m
    row_sum = sum ( abs ( a(i,1:n) ) )
    c8mat_norm_li = max ( c8mat_norm_li, row_sum )
  end do

  return
end
subroutine c8mat_print ( m, n, a, title )

!*****************************************************************************80
!
!! c8mat_print() prints a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    complex ( kind = ck ) A(M,N), the matrix.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  character ( len = * ) title

  call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! c8mat_print_some() prints some of a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns.
!
!    complex ( kind = ck ) A(M,N), the matrix.
!
!    integer ILO, JLO, IHI, JHI, the first row and
!    column, and the last row and column to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 4
  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  character ( len = 20 ) ctemp(incx)
  integer i
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title
  complex ( kind = ck ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  if ( m <= 0 .or. n <= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  (None)' 
    return
  end if
!
!  Print the columns of the matrix, in strips of INCX.
!
  do j2lo = jlo, min ( jhi, n ), incx

    j2hi = j2lo + incx - 1
    j2hi = min ( j2hi, n )
    j2hi = min ( j2hi, jhi )

    inc = j2hi + 1 - j2lo

    write ( *, '(a)' ) ' '

    do j = j2lo, j2hi
      j2 = j + 1 - j2lo
      write ( ctemp(j2), '(i10,10x)' ) j
    end do

    write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
    write ( *, '(a)' ) '  Row'
    write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
    i2lo = max ( ilo, 1 )
    i2hi = min ( ihi, m )

    do i = i2lo, i2hi
!
!  Print out (up to) INCX entries in row I, that lie in the current strip.
!
      do j2 = 1, inc

        j = j2lo - 1 + j2

        if ( imag ( a(i,j) ) == 0.0D+00 ) then
          write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j), kind = rk )
        else
          write ( ctemp(j2), '(2g10.3)' ) a(i,j)
        end if

      end do

      write ( *, '(i5,a1,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

    end do

  end do

  return
end
subroutine c8mat_scale ( m, n, alpha, a )

!*****************************************************************************80
!
!! c8mat_scale() scales a C8MAT by a complex scalar factor.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    complex ( kind = ck ) ALPHA, the scale factor.
!
!    complex ( kind = ck ) A(M,N), the matrix to be scaled.
!
!  Output:
!
!    complex ( kind = ck ) A(M,N), the scaled matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  complex ( kind = ck ) alpha

  a(1:m,1:n) = alpha * a(1:m,1:n)

  return
end
subroutine c8mat_scale_r8 ( m, n, alpha, a )

!*****************************************************************************80
!
!! c8mat_scale_r8() scales a C8MAT by a real scalar factor.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!    real ( kind = rk ) ALPHA, the scale factor.
!
!    complex ( kind = ck ) A(M,N), the matrix to be scaled.
!
!  Output:
!
!    complex ( kind = ck ) A(M,N): the scaled matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)
  real ( kind = rk ) alpha

  a(1:m,1:n) = alpha * a(1:m,1:n)

  return
end
subroutine c8mat_uniform_01 ( m, n, c )

!*****************************************************************************80
!
!! c8mat_uniform_01() returns a unit pseudorandom C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the number of rows and columns in the matrix.
!
!  Output:
!
!    complex ( kind = ck ) C(M,N), the pseudorandom complex matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  complex ( kind = ck ) c(m,n)
  integer i
  integer j
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do j = 1, n
    do i = 1, m

      call random_number ( harvest = r )
      call random_number ( harvest = theta )
      theta = 2.0D+00 * r8_pi * theta

      c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ), kind = ck )

    end do

  end do

  return
end
subroutine c8mat_zero ( m, n, a )

!*****************************************************************************80
!
!! c8mat_zero() zeroes a C8MAT.
!
!  Discussion:
!
!    A C8MAT is a matrix of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, N, the order of the matrix.
!
!  Output:
!
!    complex ( kind = ck ) A(M,N), the zeroed matrix.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer m
  integer n

  complex ( kind = ck ) a(m,n)

  a(1:m,1:n) = 0.0D+00

  return
end
subroutine c8vec_copy ( n, a, b )

!*****************************************************************************80
!
!! c8vec_copy() copies a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of elements of A.
!
!    complex ( kind = ck ) A(N), the array to be copied.
!
!  Output:
!
!    complex ( kind = ck ) B(N), the copy of the array.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  complex ( kind = ck ) b(n)

  b(1:n) = a(1:n)

  return
end
subroutine c8vec_indicator ( n, a )

!*****************************************************************************80
!
!! c8vec_indicator() sets a C8VEC to the indicator vector.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of elements of A.
!
!  Output:
!
!    complex ( kind = ck ) A(N), the array to be initialized.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i

  do i = 1, n
    a(i) = cmplx ( i, -i, kind = ck )
  end do

  return
end
subroutine c8vec_nint ( n, a )

!*****************************************************************************80
!
!! c8vec_nint() rounds the entries of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the vector.
!
!    complex ( kind = ck ) A(N), the vector to be NINT'ed.
!
!  Output:
!
!    complex ( kind = ck ) A(N): the NINT'ed array.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  complex ( kind = ck ) c8_nint
  integer i

  do i = 1, n
    a(i) = c8_nint ( a(i) )
  end do

  return
end
function c8vec_norm ( n, a )

!*****************************************************************************80
!
!! c8vec_norm() returns the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The vector L2 norm is defined as:
!
!      C8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries.
!
!    complex ( kind = ck ) A(N), the vector.
!
!  Output:
!
!    real ( kind = rk ) C8VEC_NORM, the norm.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  real ( kind = rk ) c8vec_norm

  c8vec_norm = sqrt ( &
               sum ( &
               ( abs ( a(1:n) ) )**2 &
               ) )

  return
end
function c8vec_norm_l1 ( n, a )

!*****************************************************************************80
!
!! c8vec_norm_l1() returns the L1 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The vector L1 norm is defined as:
!
!      C8VEC_NORM_L1 = sum ( 1 <= I <= N ) abs ( A(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries.
!
!    complex ( kind = ck ) A(N), the vector.
!
!  Output:
!
!    real ( kind = rk ) C8VEC_NORM_L1, the norm.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  real ( kind = rk ) c8vec_norm_l1

  c8vec_norm_l1 = sum ( abs ( a(1:n) ) )

  return
end
function c8vec_norm_l2 ( n, a )

!*****************************************************************************80
!
!! c8vec_norm_l2() returns the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The vector L2 norm is defined as:
!
!      C8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries.
!
!    complex ( kind = ck ) A(N), the vector.
!
!  Output:
!
!    real ( kind = rk ) C8VEC_NORM_L2, the norm.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  real ( kind = rk ) c8vec_norm_l2

  c8vec_norm_l2 = sqrt ( &
                  sum ( &
                  ( abs ( a(1:n) ) )**2 &
                  ) )

  return
end
function c8vec_norm_li ( n, a )

!*****************************************************************************80
!
!! c8vec_norm_li() returns the Loo norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The vector Loo norm is defined as:
!
!      C8VEC_NORM_Loo = max ( 1 <= I <= N ) abs ( A(I) )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries.
!
!    complex ( kind = ck ) A(N), the vector.
!
!  Output:
!
!    real ( kind = rk ) C8VEC_NORM_LI, the norm.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  real ( kind = rk ) c8vec_norm_li

  c8vec_norm_li = maxval ( abs ( a(1:n) ) )

  return
end
function c8vec_norm_squared ( n, a )

!*****************************************************************************80
!
!! c8vec_norm_squared() returns the square of the L2 norm of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The square of the vector L2 norm is defined as:
!
!      C8VEC_NORM_SQUARED = sum ( 1 <= I <= N ) conjg ( A(I) ) * A(I).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in A.
!
!    complex ( kind = ck ) A(N), the vector whose L2 norm is desired.
!
!  Output:
!
!    real ( kind = rk ) C8VEC_NORM_SQUARED, the L2 norm of A.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  real ( kind = rk ) c8vec_norm_squared

  c8vec_norm_squared = sum ( ( abs ( a(1:n) ) ) ** 2)

  return
end
subroutine c8vec_normal_01 ( n, x )

!*****************************************************************************80
!
!! c8vec_normal_01() returns a unit pseudonormal C8VEC.
!
!  Discussion:
!
!    A C8VEC is an array of double precision complex values.
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 December 2023
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values desired.
!
!  Output:
!
!    complex ( kind = ck ) X(N), a sample of the standard normal PDF.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) r(n)
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  complex ( kind = ck ) x(n)

  call random_number ( harvest = r )

  x = sqrt ( - 2.0D+00 * log ( r ) ) * &
    cmplx ( cos ( 2.0D+00 * r8_pi * r ), &
            sin ( 2.0D+00 * r8_pi * r ), kind = ck )

  return
end
subroutine c8vec_print ( n, a, title )

!*****************************************************************************80
!
!! c8vec_print() prints a C8VEC, with an optional title.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of components of the vector.
!
!    complex ( kind = ck ) A(N), the vector to be printed.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
  end do

  return
end
subroutine c8vec_print_part ( n, a, max_print, title )

!*****************************************************************************80
!
!! c8vec_print_part() prints "part" of a C8VEC.
!
!  Discussion:
!
!    The user specifies MAX_PRINT, the maximum number of lines to print.
!
!    If N, the size of the vector, is no more than MAX_PRINT, then
!    the entire vector is printed, one entry per line.
!
!    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
!    followed by a line of periods suggesting an omission,
!    and the last entry.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries of the vector.
!
!    complex ( kind = ck ) A(N), the vector to be printed.
!
!    integer MAX_PRINT, the maximum number of lines to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  integer max_print
  character ( len = * ) title

  if ( max_print <= 0 ) then
    return
  end if

  if ( n <= 0 ) then
    return
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  if ( n <= max_print ) then

    do i = 1, n
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do

  else if ( 3 <= max_print ) then

    do i = 1, max_print - 2
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    write ( *, '(a)' ) '  ........  ..............  ..............'
    i = n
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)

  else

    do i = 1, max_print - 1
      write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a(i)
    end do
    i = max_print
    write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6,2x,a)' ) i, ':', a(i), &
      '...more entries...'

  end if

  return
end
subroutine c8vec_print_some ( n, x, i_lo, i_hi, title )

!*****************************************************************************80
!
!! c8vec_print_some() prints some of a C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries of the vector.
!
!    complex ( kind = ck ) X(N), the vector to be printed.
!
!    integer I_LO, I_HI, the first and last entries to print.
!
!    character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  integer i
  integer i_hi
  integer i_lo
  character ( len = * ) title
  complex ( kind = ck ) x(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = max ( 1, i_lo ), min ( n, i_hi )
    write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
  end do

  return
end
subroutine c8vec_sort_a_l1 ( n, x )

!*****************************************************************************80
!
!! c8vec_sort_a_l1() ascending sorts a C8VEC by L1 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L1 norm of A+Bi is abs(A) + abs(B).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the array.
!
!    complex ( kind = ck ) X(N): the array to be sorted.
!
!  Output:
!
!    complex ( kind = ck ) X(N): the sorted array.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  logical c8_le_l1
  integer i
  integer indx
  integer isgn
  integer j
  complex ( kind = ck ) t
  complex ( kind = ck ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      t    = x(i)
      x(i) = x(j)
      x(j) = t
 
    else if ( indx < 0 ) then

      if ( c8_le_l1 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_a_l2 ( n, x )

!*****************************************************************************80
!
!! c8vec_sort_a_l2() ascending sorts a C8VEC by L2 norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L2 norm of A+Bi is sqrt ( A^2 + B^2 ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the array.
!
!    complex ( kind = ck ) X(N), the array to be sorted.
!
!  Output:
!
!    complex ( kind = ck ) X(N): the sorted array.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  logical c8_le_l2
  integer i
  integer indx
  integer isgn
  integer j
  complex ( kind = ck ) t
  complex ( kind = ck ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      t    = x(i)
      x(i) = x(j)
      x(j) = t

    else if ( indx < 0 ) then

      if ( c8_le_l2 ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_sort_a_li ( n, x )

!*****************************************************************************80
!
!! c8vec_sort_a_li() ascending sorts a C8VEC by L-infinity norm.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The L infinity norm of A+Bi is max ( abs ( A ), abs ( B ) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of entries in the array.
!
!    complex ( kind = ck ) X(N): the array to be sorted.
!
!  Output:
!
!    complex ( kind = ck ) X(N): the sorted array.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer n

  logical c8_le_li
  integer i
  integer indx
  integer isgn
  integer j
  complex ( kind = ck ) t
  complex ( kind = ck ) x(n)

  if ( n <= 1 ) then
    return
  end if

  i = 0
  indx = 0
  isgn = 0
  j = 0

  do

    call sort_heap_external ( n, indx, i, j, isgn )

    if ( 0 < indx ) then

      t    = x(i)
      x(i) = x(j)
      x(j) = t

    else if ( indx < 0 ) then

      if ( c8_le_li ( x(i), x(j) ) ) then
        isgn = -1
      else
        isgn = +1
      end if

    else if ( indx == 0 ) then

      exit

    end if

  end do

  return
end
subroutine c8vec_spiral ( n, m, c1, c2, c )

!*****************************************************************************80
!
!! c8vec_spiral() returns N points on a spiral between C1 and C2.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    Let the polar form of C1 be ( R1, T1 ) and the polar form of C2 
!    be ( R2, T2 ) where, if necessary, we increase T2 by 2*PI so that T1 <= T2.
!    
!    Then the polar form of the I-th point C(I) is:
!
!      R(I) = ( ( N - I     ) * R1 
!             + (     I - 1 ) * R2 ) 
!              / ( N    - 1 )
!
!      T(I) = ( ( N - I     ) * T1 
!             + (     I - 1 ) * ( T2 + M * 2 * PI ) ) 
!             / ( N     - 1 )
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of points on the spiral.
!
!    integer M, the number of full circuits the 
!    spiral makes.
!
!    complex ( kind = ck ) C1, C2, the first and last points 
!    on the spiral.
!
!  Output:
!
!    complex ( kind = ck ) C(N), the points.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) c(n)
  complex ( kind = ck ) c1
  complex ( kind = ck ) c2
  real ( kind = rk ) c8_arg
  integer i
  integer m
  real ( kind = rk ) r1
  real ( kind = rk ) r2
  real ( kind = rk ) ri
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) t1
  real ( kind = rk ) t2
  real ( kind = rk ) ti

  r1 = abs ( c1 )
  r2 = abs ( c2 )

  t1 = c8_arg ( c1 )
  t2 = c8_arg ( c2 )

  if ( m == 0 ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0D+00 * r8_pi
    end if

  else if ( 0 < m ) then

    if ( t2 < t1 ) then
      t2 = t2 + 2.0D+00 * r8_pi
    end if

    t2 = t2 + real ( m, kind = rk ) * 2.0D+00 * r8_pi

  else if ( m < 0 ) then

    if ( t1 < t2 ) then
      t2 = t2 - 2.0D+00 * r8_pi
    end if

    t2 = t2 - real ( m, kind = rk ) * 2.0D+00 * r8_pi

  end if

  do i = 1, n

    ri = ( real ( n - i,     kind = rk ) * r1 &
         + real (     i - 1, kind = rk ) * r2 ) &
         / real ( n     - 1, kind = rk )

    ti = ( real ( n - i,     kind = rk ) * t1 &
         + real (     i - 1, kind = rk ) * t2 ) &
         / real ( n     - 1, kind = rk )

    call polar_to_c8 ( ri, ti, c(i) )

  end do

  return
end
subroutine c8vec_uniform_01 ( n, c )

!*****************************************************************************80
!
!! c8vec_uniform_01() returns a unit pseudorandom C8VEC.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    The angles should be uniformly distributed between 0 and 2 * PI,
!    the square roots of the radius uniformly distributed between 0 and 1.
!
!    This results in a uniform distribution of values in the unit circle.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of values to compute.
!
!  Output:
!
!    complex ( kind = ck ) C(N), the pseudorandom complex vector.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) c(n)
  integer i
  real ( kind = rk ) r
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do i = 1, n

    call random_number ( harvest = r )
    call random_number ( harvest = theta )
    theta = 2.0 * r8_pi * theta

    c(i) = r * cmplx ( cos ( theta ), sin ( theta ), kind = ck )

  end do

  return
end
subroutine c8vec_unity ( n, a )

!*****************************************************************************80
!
!! c8vec_unity() returns the N roots of unity.
!
!  Discussion:
!
!    A C8VEC is a vector of C8's.
!
!    X(1:N) = exp ( 2 * PI * (0:N-1) / N )
!
!    X(1:N)^N = ( (1,0), (1,0), ..., (1,0) ).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the number of elements of A.
!
!  Output:
!
!    complex ( kind = ck ) A(N), the N roots of unity.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  complex ( kind = ck ) a(n)
  integer i
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) theta

  do i = 1, n
    theta = r8_pi * real ( 2 * ( i - 1 ), kind = rk ) / real ( n, kind = rk )
    a(i) = cmplx ( cos ( theta ), sin ( theta ), kind = ck )
  end do

  return
end
subroutine cartesian_to_c8 ( x, y, z )

!*****************************************************************************80
!
!! cartesian_to_c8() converts a Cartesian form to a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) X, Y, the Cartesian form.
!
!  Output:
!
!    complex ( kind = ck ) Z, the complex number.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) x
  real ( kind = rk ) y
  complex ( kind = ck ) z

  z = cmplx ( x, y, kind = ck )

  return
end
subroutine polar_to_c8 ( r, theta, z )

!*****************************************************************************80
!
!! polar_to_c8() converts a polar form to a C8.
!
!  Discussion:
!
!    A C8 is a complex ( kind = ck ) value.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) R, THETA, the polar form.
!
!  Output:
!
!    complex ( kind = ck ) Z, the complex number.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) theta
  complex ( kind = ck ) z

  z = r * cmplx ( cos ( theta ), sin ( theta ), kind = ck )

  return
end
subroutine r8poly2_root ( a, b, c, r1, r2 )

!*****************************************************************************80
!
!! r8poly2_root() returns the two roots of a quadratic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X * X + B * X + C = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) A, B, C, the coefficients of the polynomial.
!    A must not be zero.
!
!  Output:
!
!    complex ( kind = ck ) R1, R2, the roots of the polynomial, which
!    might be real and distinct, real and equal, or complex conjugates.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  complex ( kind = ck ) disc
  complex ( kind = ck ) q
  complex ( kind = ck ) r1
  complex ( kind = ck ) r2

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY2_ROOT - Fatal error!'
    write ( *, '(a)' ) '  The coefficient A is zero.'
    stop 1
  end if

  disc = b * b - 4.0D+00 * a * c
  q = -0.5D+00 * ( b + sign ( 1.0D+00, b ) * sqrt ( disc ) )
  r1 = q / a
  r2 = c / q

  return
end
subroutine r8poly3_root ( a, b, c, d, r1, r2, r3 )

!*****************************************************************************80
!
!! r8poly3_root() returns the three roots of a cubic polynomial.
!
!  Discussion:
!
!    The polynomial has the form
!
!      A * X^3 + B * X^2 + C * X + D = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Input:
!
!    real ( kind = rk ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!  Output:
!
!    complex ( kind = ck ) R1, R2, R3, the roots of the polynomial, which
!    will include at least one real root.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  complex ( kind = ck ) i
  complex ( kind = ck ) one
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) q
  real ( kind = rk ) r
  complex ( kind = ck ) r1
  complex ( kind = ck ) r2
  complex ( kind = ck ) r3
  real ( kind = rk ) s1
  real ( kind = rk ) s2
  real ( kind = rk ) temp
  real ( kind = rk ) theta

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY3_ROOT - Fatal error!'
    write ( *, '(a)' ) '  A must not be zero!'
    stop 1
  end if

  one = cmplx ( 1.0d+00, 0.0D+00, kind = ck )
  i = sqrt ( - one )

  q = ( ( b / a ) ** 2 - 3.0D+00 * ( c / a ) ) / 9.0D+00

  r = ( 2.0D+00 * ( b / a ) ** 3 - 9.0D+00 * ( b / a ) * ( c / a ) &
      + 27.0D+00 * ( d / a ) ) / 54.0D+00

  if ( r * r < q * q * q ) then

    theta = acos ( r / sqrt ( q ** 3 ) )
    r1 = -2.0D+00 * sqrt ( q ) * cos (   theta                     / 3.0D+00 )
    r2 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 2.0D+00 * r8_pi ) / 3.0D+00 )
    r3 = -2.0D+00 * sqrt ( q ) * cos ( ( theta + 4.0D+00 * r8_pi ) / 3.0D+00 )

  else if ( q * q * q <= r * r ) then

    temp = -r + sqrt ( r ** 2 - q ** 3 )
    s1 = sign ( 1.0D+00, temp ) * ( abs ( temp ) ) ** ( 1.0D+00 / 3.0D+00 )

    temp = -r - sqrt ( r ** 2 - q ** 3 )
    s2 = sign ( 1.0D+00, temp ) * ( abs ( temp ) ) ** ( 1.0D+00 / 3.0D+00 )

    r1 = s1 + s2
    r2 = -0.5D+00 * ( s1 + s2 ) + i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )
    r3 = -0.5D+00 * ( s1 + s2 ) - i * 0.5D+00 * sqrt ( 3.0D+00 ) * ( s1 - s2 )

  end if

  r1 = r1 - b / ( 3.0D+00 * a )
  r2 = r2 - b / ( 3.0D+00 * a )
  r3 = r3 - b / ( 3.0D+00 * a )

  return
end
subroutine r8poly4_root ( a, b, c, d, e, r1, r2, r3, r4 )

!*****************************************************************************80
!
!! r8poly4_root() returns the four roots of a quartic polynomial.
!
!  Discussion:
!
!    The polynomial has the form:
!
!      A * X^4 + B * X^3 + C * X^2 + D * X + E = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2021
!
!  Input:
!
!    real ( kind = rk ) A, B, C, D, the coefficients of the polynomial.
!    A must not be zero.
!
!  Output:
!
!    complex ( kind = ck ) R1, R2, R3, R4, the roots of the polynomial.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) a3
  real ( kind = rk ) a4
  real ( kind = rk ) b
  real ( kind = rk ) b3
  real ( kind = rk ) b4
  real ( kind = rk ) c
  real ( kind = rk ) c3
  real ( kind = rk ) c4
  real ( kind = rk ) d
  real ( kind = rk ) d3
  real ( kind = rk ) d4
  real ( kind = rk ) e
  complex ( kind = ck ) p
  complex ( kind = ck ) q
  complex ( kind = ck ) r
  complex ( kind = ck ) r1
  complex ( kind = ck ) r2
  complex ( kind = ck ) r3
  complex ( kind = ck ) r4
  complex ( kind = ck ) zero

  zero = cmplx ( 0.0D+00, 0.0D+00, kind = ck )

  if ( a == 0.0D+00 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8POLY4_ROOT - Fatal error!'
    write ( *, '(a)') '  A must not be zero!'
    stop 1
  end if

  a4 = b / a
  b4 = c / a
  c4 = d / a
  d4 = e / a
!
!  Set the coefficients of the resolvent cubic equation.
!
  a3 = 1.0D+00
  b3 = -b4
  c3 = a4 * c4 - 4.0D+00 * d4
  d3 = -a4 * a4 * d4 + 4.0D+00 * b4 * d4 - c4 * c4
!
!  Find the roots of the resolvent cubic.
!
  call r8poly3_root ( a3, b3, c3, d3, r1, r2, r3 )
!
!  Choose one root of the cubic, here R1.
!
!  Set R = sqrt ( 0.25D+00 * A4 ** 2 - B4 + R1 )
!
  r = sqrt ( 0.25D+00 * a4 ** 2 - b4 + r1 )

  if ( r /= zero ) then

    p = sqrt ( 0.75D+00 * a4 ** 2 - r ** 2 - 2.0D+00 * b4 &
        + 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4 ** 3 ) / r )

    q = sqrt ( 0.75D+00 * a4 ** 2 - r ** 2 - 2.0D+00 * b4 &
        - 0.25D+00 * ( 4.0D+00 * a4 * b4 - 8.0D+00 * c4 - a4 ** 3 ) / r )

  else

    p = sqrt ( 0.75D+00 * a4 ** 2 - 2.0D+00 * b4 &
      + 2.0D+00 * sqrt ( r1 ** 2 - 4.0D+00 * d4 ) )

    q = sqrt ( 0.75D+00 * a4 ** 2 - 2.0D+00 * b4 &
      - 2.0D+00 * sqrt ( r1 ** 2 - 4.0D+00 * d4 ) )

  end if
!
!  Set the roots.
!
  r1 = -0.25D+00 * a4 + 0.5D+00 * r + 0.5D+00 * p
  r2 = -0.25D+00 * a4 + 0.5D+00 * r - 0.5D+00 * p
  r3 = -0.25D+00 * a4 - 0.5D+00 * r + 0.5D+00 * q
  r4 = -0.25D+00 * a4 - 0.5D+00 * r - 0.5D+00 * q

  return
end

