program main

!*****************************************************************************80
!
!! roots_rc_test() tests roots_rc().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 4

  real ( kind = rk ) ferr
  real ( kind = rk ) fx(n)
  integer it
  integer, parameter :: it_max = 30
  real ( kind = rk ) q(2*n+2,n+2)
  real ( kind = rk ) x(n)
  real ( kind = rk ) xnew(n)

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'roots_rc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test ROOTS_RC().'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '       FERR          X'
  write ( *, '(a)' ) ' '
!
!  Initialization.
!
  q(1:2*n+2,1:n+2) = 0.0D+00

  xnew(1) = 1.2D+00
  xnew(2:n) = 1.0D+00

  it = 0

  do

    x(1:n) = xnew(1:n)

    fx(1) = 1.0D+00 - x(1)
    fx(2:n) = 10.0D+00 * ( x(2:n) - x(1:n-1)**2 )

    if ( it == 0 ) then
      write ( *, '(2x,14x,5(2x,g14.6))' ) x(1:n)
    else
      write ( *, '(2x,g14.6,5(2x,g14.6))' ) ferr, x(1:n)
    end if

    call roots_rc ( n, x, fx, ferr, xnew, q )

    if ( ferr < 1.0D-07 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Sum of |f(x)| less than tolerance.'
      exit
    end if

    if ( it_max < it ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Too many iterations!'
      exit
    end if

    it = it + 1

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ROOTS_RC_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end

