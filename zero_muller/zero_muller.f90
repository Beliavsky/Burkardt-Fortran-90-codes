subroutine zero_muller ( func, fatol, itmax, x1, x2, x3, xatol, xrtol, &
  xnew, fxnew )

!*****************************************************************************80
!
!! zero_muller() carries out Muller's method, using complex arithmetic.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 March 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Gisela Engeln-Muellges, Frank Uhlig,
!    Numerical Algorithms with C,
!    Springer, 1996,
!    ISBN: 3-540-60530-4,
!    LC: QA297.E56213.
!
!  Input:
!
!    external FUNC, the name of the routine that evaluates the function.
!    FUNC should have the form:
!      subroutine func ( x, fx )
!      complex ( kind = ck ) fx
!      complex ( kind = ck ) x
!
!    real ( kind = rk ) FATOL, the absolute error tolerance for F(X).
!
!    integer ITMAX, the maximum number of steps allowed.
!
!    complex ( kind = ck ) X1, X2, X3, three distinct points to start the
!    iteration.
!
!    real ( kind = rk ) XATOL, XRTOL, absolute and relative
!    error tolerances for the root.
!
!  Output:
!
!    complex ( kind = ck ) XNEW, the estimated root.
!
!    complex ( kind = ck ) FXNEW, the value of the function at XNEW.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )

  complex ( kind = ck ) a
  complex ( kind = ck ) b
  complex ( kind = ck ) c
  complex ( kind = ck ) c8_temp
  complex ( kind = ck ) discrm
  real ( kind = rk ) fatol
  complex ( kind = ck ) fminus
  complex ( kind = ck ) fplus
  external func
  complex ( kind = ck ) fxmid
  complex ( kind = ck ) fxnew
  complex ( kind = ck ) fxold
  integer iterate
  integer itmax
  real ( kind = rk ) x_ave
  complex ( kind = ck ) x_inc
  complex ( kind = ck ) x1
  complex ( kind = ck ) x2
  complex ( kind = ck ) x3
  real ( kind = rk ) xatol
  complex ( kind = ck ) xlast
  complex ( kind = ck ) xmid
  complex ( kind = ck ) xminus
  complex ( kind = ck ) xnew
  complex ( kind = ck ) xold
  complex ( kind = ck ) xplus
  real ( kind = rk ) xrtol

  xnew = x1
  xmid = x2
  xold = x3

  call func ( xnew, fxnew )
  call func ( xmid, fxmid )
  call func ( xold, fxold )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'zero_muller():'
  write ( *, '(a)' ) '  Muller''s root-finding method (complex root version)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '  Iteration     x_real              x_imag             ||fx||           ||disc||'
  write ( *, '(a)' ) ' '

  iterate = -2
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xold, abs ( fxold )
  iterate = -1
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xmid, abs ( fxmid )
  iterate = 0
  write ( *, '(i6,f20.10,f20.10,f20.10)' ) iterate, xnew, abs ( fxnew )

  if ( abs ( fxnew ) < fatol ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'zero_muller():'
    write ( *, '(a)' ) '  |F(X)| is below the tolerance.'
    return
  end if

  do
!
!  You may need to swap (XMID,FXMID) and (XNEW,FXNEW).
!
    if ( abs ( fxmid ) <= abs ( fxnew ) ) then

      c8_temp = xnew
      xnew = xmid
      xmid = c8_temp

      c8_temp = fxnew
      fxnew = fxmid
      fxmid = c8_temp

    end if

    xlast = xnew
    iterate = iterate + 1

    if ( itmax < iterate ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'zero_muller(): Warning!'
      write ( *, '(a)' ) '  Maximum number of steps taken.'
      exit
    end if

    a =  ( ( xmid - xnew ) * ( fxold - fxnew ) &
         - ( xold - xnew ) * ( fxmid - fxnew ) )

    b = ( ( xold - xnew )**2 * ( fxmid - fxnew ) &
        - ( xmid - xnew )**2 * ( fxold - fxnew ) )

    c = ( ( xold - xnew ) * ( xmid - xnew ) * ( xold - xmid ) * fxnew )

    xold = xmid
    xmid = xnew
!
!  Apply the quadratic formula to get roots XPLUS and XMINUS.
!
    discrm = b**2 - 4.0D+00 * a * c

    if ( a == 0.0D+00 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'zero_muller(): Warning!'
      write ( *, '(a)' ) '  The algorithm has broken down.'
      write ( *, '(a)' ) '  The quadratic coefficient A is zero.'
      exit
    end if

    xplus = xnew + ( ( - b + sqrt ( discrm ) ) / ( 2.0D+00 * a ) )

    call func ( xplus, fplus )

    xminus = xnew + ( ( - b - sqrt ( discrm ) ) / ( 2.0D+00 * a ) )

    call func ( xminus, fminus )
!
!  Choose the root with smallest function value.
!
    if ( abs ( fminus ) < abs ( fplus ) ) then
      xnew = xminus
    else
      xnew = xplus
    end if

    fxold = fxmid
    fxmid = fxnew
    call func ( xnew, fxnew )
    write ( *, '(i6,f20.10,f20.10,f20.10,f20.10)' ) &
      iterate, xnew, abs ( fxnew ), abs ( discrm )
!
!  Check for convergence.
!
    x_ave = abs ( xnew + xmid + xold ) / 3.0D+00
    x_inc = xnew - xmid

    if ( abs ( x_inc ) <= xatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'zero_muller():'
      write ( *, '(a)' ) '  Absolute convergence of the X increment.'
      exit
    end if

    if ( abs ( x_inc ) <= xrtol * x_ave ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'zero_muller():'
      write ( *, '(a)' ) '  Relative convergence of the X increment.'
      exit
    end if

    if ( abs ( fxnew ) <= fatol ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'zero_muller():'
      write ( *, '(a)' ) '  Absolute convergence of |F(X)|.'
      exit
    end if

  end do

  return
end

