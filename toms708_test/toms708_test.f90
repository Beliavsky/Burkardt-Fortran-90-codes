program main

!*****************************************************************************80
!
!! toms708_test() tests toms708().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS708_TEST:'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS708.'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call erf_test ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS708_TEST:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests BRATIO.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real b
  integer i
  integer ierr
  real w
  real w1
  real x
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  BRATIO computes the Beta ratio function.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Y       W                         W1 = 1-W'
  write ( *, '(a)' ) ' '

  a = 5.3E+00
  b = 10.1E+00

  do i = 1, 50

    x = real ( i ) / 100.0E+00

    y = 0.5E+00 + ( 0.5E+00 - x )

    call bratio ( a, b, x, y, w, w1, ierr )

    if ( ierr == 0 ) then
      write ( *, '(2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16)' ) x, y, w, w1
    else
      write ( *, '(2x,f6.2,2x,f6.2,2x,a32)' ) &
      x, y, '----FAILURE---- ----FAILURE-----'
    end if

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests BRATIO and BETA_INC_VALUES.
!
!  Modified:
!
!    11 January 2006
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real b
  real fx
  integer ierr
  integer n
  real w
  real w1
  real x
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02:'
  write ( *, '(a)' ) '  BRATIO evaluates the normalized incomplete Beta'
  write ( *, '(a)' ) '    function BETA_INC(A,B,X).'
  write ( *, '(a)' ) '  BETA_INC_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A       B       X       ' // &
    'Exact F                   BRATIO(A,B,X)              DIFF'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call beta_inc_values ( n, a, b, x, fx )

    if ( n == 0 ) then
      exit
    end if

    y = 1.0E+00 - x

    call bratio ( a, b, x, y, w, w1, ierr )

    if ( ierr == 0 ) then

      write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      a, b, x, fx, w, abs ( fx - w )

    else

      write ( *, '(2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16,a14)' ) &
      a, b, x, fx, '---FAILURE----'

    end if

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests BETALN and BETA_LOG_VALUES.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real betaln
  real fxy
  real fxy2
  integer n_data
  real x
  real y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03:'
  write ( *, '(a)' ) '  BETALN evaluates the logarithm of the '
  write ( *, '(a)' ) '    Beta function.'
  write ( *, '(a)' ) '  BETA_LOG_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Y        Exact F                   ' // &
    'BETALN(X,Y)             DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_log_values ( n_data, x, y, fxy )

    if ( n_data == 0 ) then
      exit
    end if

    fxy2 = betaln ( x, y )

    write ( *, '(2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, y, fxy, fxy2, abs ( fxy - fxy2 )

  end do

  return
end
subroutine erf_test ( )

!*****************************************************************************80
!
!! erf_test tests r4_erf and ERF_VALUES.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  integer n_data
  real r4_erf
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'erf_test:'
  write ( *, '(a)' ) '  R4_ERF evaluates the error function.'
  write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '     X        Exact F                   ' // &
    'R4_ERF(X)              DIFF'
  write ( *, '(a)' ) '              (Tabulated)               ' // &
    '(Calculated)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call erf_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = r4_erf ( x )

    write ( *, '(2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests GAMLN and GAMMA_INC_VALUES.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real a
  real eps
  real fx
  real gamln
  integer n_data
  real p
  real q
  real r
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05:'
  write ( *, '(a)' ) '  GRAT1 evaluates the incomplete Gamma function '
  write ( *, '(a)' ) '    for A <= 1.'
  write ( *, '(a)' ) '  GAMMA_INC_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A       X       Exact F                   ' // &
    'GAMMA_INC(A,X)         DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_inc_values ( n_data, a, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    if ( a <= 1.0E+00 ) then

      r = exp ( - x ) * x**a / exp ( gamln ( a ) )

      eps = epsilon ( x )

      call grat1 ( a, x, r, p, q, eps )

      write ( *, '(2x,f6.2,2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
        a, x, fx, p, abs ( p - fx )

    else

      write ( *, '(2x,f6.2,2x,f6.2,2x,g24.16,2x,a)' ) a, x, fx, 'Unavailable'

    end if

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests GAMLN and GAMMA_LOG_VALUES.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  real gamln
  integer n_data
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06:'
  write ( *, '(a)' ) '  GAMLN evaluates the logarithm of the '
  write ( *, '(a)' ) '    Gamma function.'
  write ( *, '(a)' ) '  GAMMA_LOG_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X        Exact F                   ' // &
    'GAMLN(X)               DIFF'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call gamma_log_values ( n_data, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    fx2 = gamln ( x )

    write ( *, '(2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests PSI and PSI_VALUES.
!
!  Modified:
!
!    08 January 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real fx
  real fx2
  real psi
  integer n
  real x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07:'
  write ( *, '(a)' ) '  PSI evaluates the PSI function.'
  write ( *, '(a)' ) '  PSI_VALUES returns some exact values.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X       Exact F                   ' // &
    'PSI(X)                 DIFF'
  write ( *, '(a)' ) ' '

  n = 0

  do

    call psi_values ( n, x, fx )

    if ( n == 0 ) then
      exit
    end if

    if ( x <= 0.0E+00 ) then
      cycle
    end if

    fx2 = psi ( x )

    write ( *, '(2x,f6.2,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
      x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
