program main

!*****************************************************************************80
!
!! toms243_test() tests toms243().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2018
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) fx1
  complex ( kind = ck ) fx2
  integer n_data
  complex ( kind = ck ) x

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOMS243_TEST():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  TOMS243 computes the natural logarithm of a complex value.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '               X                               FX exact'
  write ( *, '(a)' ) '                                               FX computed'
  write ( *, '(a)' ) ''

  n_data = 0

  do

    call c8_log_values ( n_data, x, fx1 )

    if ( n_data == 0 ) then
      exit
    end if

    call toms243 ( x, fx2 )

    write ( *, '(2x,a,f8.4,a,f8.4,a,2x,a,f18.12,a,f18.12,a)' ) &
      '(', real ( x, kind = ck ), ',', aimag ( x ), ' )', &
      '(', real ( fx1, kind = ck ), ',', aimag ( fx1 ), ' )'
    write ( *, '(2x,1x,8x,1x,8x,2x,2x,a,f18.12,a,f18.12,a)' ) &
      '(', real ( fx2, kind = ck ), ',', aimag ( fx2 ), ' )'

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'TOMS243_TEST:'
  write ( *, '(a)' ) '  Normal end of execution:'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop 0
end
subroutine c8_log_values ( n_data, x, fx )

!*****************************************************************************80
!
!! C8_LOG_VALUES: natural logarithm function for complex values.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 January 2018
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    David Collens,
!    Algorithm 243: Logarithm of a Complex Number,
!    Communications of the Association for Computing Machinery,
!    Volume 7, Number 11, November 1964, page 660.
!
!  Parameters:
!
!    Input/output, integer N_DATA.  The user sets N_DATA to 
!    0 before the first call.  On each call, the routine increments N_DATA 
!    by 1, and returns the corresponding data; when there is no more data, the
!    output value of N_DATA will be 0 again.
!
!    Output, complex ( kind = ck ) X, the argument of the function.
!
!    Output, complex ( kind = ck ) FX, the value of the function.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  integer, parameter :: n_max = 12

  complex ( kind = ck ) fx
  complex ( kind = ck ), save, dimension ( n_max ) :: fx_vec = (/ &
  ( 1.039720770839918D+00, - 2.356194490192345D+00 ), &
  ( 0.804718956217050D+00, + 2.677945044588987D+00 ), &
  ( 0.346573590279973D+00, - 2.356194490192345D+00 ), &
  ( 0.000000000000000D+00, + 3.141592653589793D+00 ), &
  ( 0.693147180559945D+00, - 1.570796326794897D+00 ), &
  ( 0.000000000000000D+00, - 1.570796326794897D+00 ), &
  ( 0.000000000000000D+00, + 1.570796326794897D+00 ), &
  ( 0.693147180559945D+00, + 1.570796326794897D+00 ), &
  ( 0.346573590279973D+00, - 0.785398163397448D+00 ), &
  ( 0.000000000000000D+00, + 0.000000000000000D+00 ), &
  ( 1.039720770839918D+00, - 0.785398163397448D+00 ), &
  ( 0.804718956217050D+00, + 0.463647609000806D+00 ) /)
  integer n_data
  complex ( kind = ck ) x
  complex ( kind = ck ), save, dimension ( n_max ) :: x_vec = (/ &
    ( -2.0D+00, - 2.0D+00 ), &
    ( -2.0D+00, + 1.0D+00 ), &
    ( -1.0D+00, - 1.0D+00 ), &
    ( -1.0D+00, + 0.0D+00 ), &
    (  0.0D+00, - 2.0D+00 ), &
    (  0.0D+00, - 1.0D+00 ), &
    (  0.0D+00, + 1.0D+00 ), &
    (  0.0D+00, + 2.0D+00 ), &
    (  1.0D+00, - 1.0D+00 ), &
    (  1.0D+00, + 0.0D+00 ), &
    (  2.0D+00, - 2.0D+00 ), &
    (  2.0D+00, + 1.0D+00 ) /)

  if ( n_data < 0 ) then
    n_data = 0
  end if

  n_data = n_data + 1

  if ( n_max < n_data ) then
    n_data = 0
    x = ( 0.0D+00, 0.0D+00 )
    fx = ( 0.0D+00, 0.0D+00 )
  else
    x = x_vec(n_data)
    fx = fx_vec(n_data)
  end if

  return
end
