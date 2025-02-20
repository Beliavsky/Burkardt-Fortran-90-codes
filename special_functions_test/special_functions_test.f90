program main

!*****************************************************************************80
!
!! special_functions_test() tests special_functions().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2022
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'special_functions_test():'
  write ( *, '(a)' ) '  Fortran90 version'
  write ( *, '(a)' ) '  Test special_functions().'

  call airya_test ( )
  call airyb_test ( )
  call beta_test ( )
  call cchg_test ( )
  call cisia_test ( )
  call cisib_test ( )
  call cjy01_test ( )
  call comelp_test ( )
  call hygfx_test ( )
  call jdzo_test ( )
  call sphj_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'special_functions_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine airya_test ( )

!*****************************************************************************80
!
!! airya_test() tests airya().
!
!  Discussion:
!
!    This program computes Airy functions and their 
!    derivatives using subroutine AIRYA.
!
!  Example:
!
!     x     Ai(x)        Bi(x)        Ai'(x)       Bi'(x)
!    ----------------------------------------------------------------
!     0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
!    10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
!    20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
!    30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
!
!     x     Ai(-x)       Bi(-x)       Ai'(-x)      Bi'(-x)
!    ----------------------------------------------------------------
!     0     .35502805    .61492663     -.25881940      .44828836
!    10     .04024124   -.31467983      .99626504      .11941411
!    20    -.17640613   -.20013931      .89286286     -.79142903
!    30    -.08796819   -.22444694     1.22862060     -.48369473
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    07 October 2013
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  real ( kind = rk8 ) ad
  real ( kind = rk8 ) ai
  real ( kind = rk8 ) bd
  real ( kind = rk8 ) bi
  integer i
  real ( kind = rk8 ) x
  real ( kind = rk8 ), dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 10.0D+00, 20.0D+00, 30.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'airya_test():'
  write ( *, '(a)' ) '  Test airya()'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    x      Ai(x)         Bi(x)' // &
    '         Ai''(x)        Bi''(x)'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    x = x_test(i)
    call airya ( x, ai, bi, ad, bd )
    write ( *, '(1x,f5.1,4g16.8)' ) x, ai, bi, ad, bd
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    x     Ai(-x)        Bi(-x)' // &
    '        Ai''(-x)       Bi''(-x)'
  write ( *, '(a)' ) ' '

  do i = 1, test_num
    x = x_test(i)
    call airya ( -x, ai, bi, ad, bd )
    write ( *, '(1x,f5.1,4g16.8)' ) x, ai, bi, ad, bd
  end do

  return
end
subroutine airyb_test ( )

!*****************************************************************************80
!
!! airyb_test() tests airyb().
!
!  Discussion:
!
!    This program computes Airy functions and their
!    derivatives using subroutine AIRYB
!
!  Modified:
!
!    07 February 2012
!
!       Input:   x  --- Argument of Airy function
!       Output:  AI --- Ai(x)
!                BI --- Bi(x)
!                AD --- Ai'(x)
!                BD --- Bi'(x)
!       Example:
!
!   x       Ai(x)          Bi(x)          Ai'(x)         Bi'(x)
!  ----------------------------------------------------------------
!   0   .35502805D+00  .61492663D+00 -.25881940D+00  .44828836D+00
!  10   .11047533D-09  .45564115D+09 -.35206337D-09  .14292361D+10
!  20   .16916729D-26  .21037650D+26 -.75863916D-26  .93818393D+26
!  30   .32082176D-48  .90572885D+47 -.17598766D-47  .49533045D+48
!
!   x       Ai(-x)         Bi(-x)         Ai'(-x)        Bi'(-x)
!  ----------------------------------------------------------------
!   0       .35502805      .61492663     -.25881940      .44828836
!  10       .04024124     -.31467983      .99626504      .11941411
!  20      -.17640613     -.20013931      .89286286     -.79142903
!  30      -.08796819     -.22444694     1.22862060     -.48369473
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  real ( kind = rk8 ) ad
  real ( kind = rk8 ) ai
  real ( kind = rk8 ) bd
  real ( kind = rk8 ) bi
  integer i
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x_test(test_num)

  save x_test

  data x_test / 0.0D+00, 10.0D+00, 20.0D+00, 30.0D+00 /

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'airyb_test():'
  write ( *, '(a)' ) '  Test airyb()'

  write(*,30)
  write(*,40)

  do i = 1, test_num
    x = x_test(i)
    call AIRYB(X,AI,BI,AD,BD)
    write(*,10)X,AI,BI,AD,BD
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write(*,50)
  write(*,40)

  do i = 1, test_num
    x = x_test(i)
    call AIRYB(-X,AI,BI,AD,BD)
    write(*,20)X,AI,BI,AD,BD
  end do

  return
10  FORMAT(1X,F5.1,4D16.8)
20  FORMAT(1X,F5.1,4D16.8)
30  FORMAT(4X,'x',8X,'Ai(x)',11X,'Bi(x)',11X,'Ai''(x)',10X,'Bi''(x)')
40  FORMAT(2X,'----------------------------------', &
             '-----------------------------------')
50  FORMAT(4X,'x',8X,'Ai(-x)',10X,'Bi(-x)',10X, &
             'Ai''(-x)',9X,'Bi''(-x)')
end
subroutine beta_test ( )

!*****************************************************************************80
!
!! beta_test() tests beta().
!
!  Example:
!
!                 p       q           B(p,q)
!               ---------------------------------
!                1.5     2.0     .2666666667D+00
!                2.5     2.0     .1142857143D+00
!                1.5     3.0     .1523809524D+00
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 3

  real ( kind = rk8 ) bt
  real ( kind = rk8 ) p
  real ( kind = rk8 ), dimension ( test_num ) :: p_test = (/ &
    1.5D+00, 2.5D+00, 1.5D+00 /)
  real ( kind = rk8 ) q
  real ( kind = rk8 ), dimension ( test_num ) :: q_test = (/ &
    2.0D+00, 2.0D+00, 3.0D+00 /)
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BETA_TEST:'
  write ( *, '(a)' ) '  Test beta()'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    p       q           B(p,q)'
  write ( *, '(a)' ) '  ---------------------------------'

  do test = 1, test_num

    p = p_test(test)
    q = q_test(test)

    call beta ( p, q, bt )
    write ( *, '(2x,f5.1,3x,f5.1,g20.10)' ) p, q, bt

  end do

  return
end
subroutine cchg_test ( )

!*****************************************************************************80
!
!! cchg_test() tests cchg().
!
!  Discussion:
!
!    This program computes confluent hypergeometric function M(a,b,z) with 
!    real parameters a, b, and a complex argument z using subroutine CCHG.
!
!  Example:
!
!     a      b        z        Re[M(a,b,z)]   Im[M(a,b,z)]
!    -------------------------------------------------------
!    3.3   4.25    10 + 0i    .61677489D+04    0
!    3.3   4.25    25 + 0i    .95781835D+10    0
!    3.3   4.25     3 -  i    .75828716D+01  -.86815474D+01
!    3.3   4.25    15 +10i   -.58313765D+06  -.48195426D+05
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2021
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Local:
!
!    real ( kind = rk8 ) A, B, parameters.
!
!    complex ( kind = ck ) Z, a parameter.
!
!    complex ( kind = ck ) CHG, the value of the function at (A,B,Z).
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) b
  complex ( kind = ck ) chg
  integer i
  complex ( kind = ck ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cchg_test():'
  write ( *, '(a)' ) '  Test cchg()'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '    a      b        z        M(a,b,z)'
  write ( *, '(a)' ) '  ------------------------------------------------------'


  a = 3.3D+00
  b = 4.25D+00

  do i = 1, 4

    if ( i == 1 ) then
      z = cmplx ( 10.0D+00,  0.0D+00, kind = ck )
    else if ( i == 2 ) then
      z = cmplx ( 25.0D+00,  0.0D+00, kind = ck )
    else if ( i == 3 ) then
      z = cmplx (  3.0D+00, -1.0D+00, kind = ck )
    else if ( i == 4 ) then
      z = cmplx ( 15.0D+00, 10.0D+00, kind = ck )
    end if

    call cchg ( a, b, z, chg )

    write ( *, '(2x,f5.1,2x,f5.1,2x,f5.1,2x,f5.1,2x,g14.6,2x,g14.6)' ) &
      a, b, z, chg

  end do

  return
end
subroutine cisia_test ( )

!*****************************************************************************80
!
!! cisia_test() tests cisia().
!
!  Example:
!
!      x        Ci(x)           Si(x)
!    ------------------------------------
!     0.0    - oo                 0
!     5.0    -.190030D+00      1.549931
!    10.0    -.454563D-01      1.658348
!    20.0     .444201D-01      1.548241
!    30.0    -.330326D-01      1.566757
!    40.0     .190201D-01      1.586985
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 April 2013
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 6

  real ( kind = rk8 ) ci
  real ( kind = rk8 ) si
  integer test
  real ( kind = rk8 ) x
  real ( kind = rk8 ), save, dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CISIA_TEST'
  write ( *, '(a)' ) '  Test cisia()'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  CISIA computes the cosine and sine integrals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   x        ci(x)           si(x)'
  write ( *, '(a)' ) '------------------------------------'

  do test = 1, test_num

    x = x_test(test)

    call cisia ( x, ci, si )

    write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si

  end do

  return
end
subroutine cisib_test ( )

!*****************************************************************************80
!
!! cisib_test() tests cisib().
!
!  Example:
!
!      x        Ci(x)           Si(x)
!    ------------------------------------
!     0.0    - oo                 0
!     5.0    -.190030D+00      1.549931
!    10.0    -.454563D-01      1.658348
!    20.0     .444201D-01      1.548241
!    30.0    -.330326D-01      1.566757
!    40.0     .190201D-01      1.586985
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 July 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: test_num = 6

  real ( kind = rk8 ) ci
  real ( kind = rk8 ) si
  integer test
  real ( kind = rk8 ) x
  real ( kind = rk8 ), save, dimension ( test_num ) :: x_test = (/ &
    0.0D+00, 5.0D+00, 10.0D+00, 20.0D+00, 30.0D+00, 40.0D+00 /)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CISIB_TEST'
  write ( *, '(a)' ) '  Test cisib()'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  CISIB computes the cosine and sine integrals.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   x        ci(x)           si(x)'
  write ( *, '(a)' ) '------------------------------------'

  do test = 1, test_num

    x = x_test(test)

    call cisib ( x, ci, si )

    write ( *, '(1x,f5.1,g16.8,g16.8)' ) x, ci, si

  end do

  return
end
subroutine cjy01_test ( )

!*****************************************************************************80
!
!! cjy01_test() tests cjy01().
!
!  Discussion:
!
!    This program computes Bessel functions J0(z), J1(z), Y0(z), Y1(z), and 
!    their derivatives for a complex argument using subroutine CJY01,
!
!  Example:
!
!    z =  4.0 + i  2.0
!
!     n     Re[Jn(z)]       Im[Jn(z)]       Re[Jn'(z)]      Im[Jn'(z)]
!   --------------------------------------------------------------------
!     0  -.13787022D+01   .39054236D+00   .50735255D+00   .12263041D+01
!     1  -.50735255D+00  -.12263041D+01  -.11546013D+01   .58506793D+00
!
!     n     Re[Yn(z)]       Im[Yn(z)]       Re[Yn'(z)]      Im[Yn'(z)]
!   --------------------------------------------------------------------
!     0  -.38145893D+00  -.13291649D+01  -.12793101D+01   .51220420D+00
!     1   .12793101D+01  -.51220420D+00  -.58610052D+00  -.10987930D+01
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 October 2015
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Local:
!
!    z --- Complex argument
!                CBJ0 --- J0(z)
!                CDJ0 --- J0'(z)
!                CBJ1 --- J1(z)
!                CDJ1 --- J1'(z)
!                CBY0 --- Y0(z)
!                CDY0 --- Y0'(z)
!                CBY1 --- Y1(z)
!                CDY1 --- Y1'(z)
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk8 = kind ( 1.0D+00 )

  complex ( kind = ck ) cbj0
  complex ( kind = ck ) cbj1
  complex ( kind = ck ) cby0
  complex ( kind = ck ) cby1
  complex ( kind = ck ) cdj0
  complex ( kind = ck ) cdj1
  complex ( kind = ck ) cdy0
  complex ( kind = ck ) cdy1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) y
  complex ( kind = ck ) z

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CJY01_TEST'
  write ( *, '(a)' ) '  Test CJY01'
  write ( *, '(a)' ) ' '

  x = 4.0D+00
  y = 2.0D+00

  z = cmplx ( x, y, kind = ck )
  write ( *, '(a,g16.8,a,g16.8)' ) '  Z = ', x, ' + i * ', y
  call cjy01 ( Z, CBJ0, CDJ0, CBJ1, CDJ1, CBY0, CDY0, CBY1, CDY1 )
  write ( *, * )
  write ( *, * ) '  n      Re[Jn(z)]       Im[Jn(z)]       Re[Jn''(z)]      Im[Jn''(z)]'
  write ( *, * ) ' --------------------------------------------------------------------'
  write ( *, '(6x,4d16.8)' ) CBJ0, CDJ0
  write ( *, '(6x,4d16.8)' ) CBJ1, CDJ1
  write ( *, * )
  write ( *, * ) '  n      Re[Yn(z)]       Im[Yn(z)]       Re[Yn''(z)]      Im[Yn''(z)]'
  write ( *, * ) ' --------------------------------------------------------------------'
  write ( *, '(6x,4d16.8)' ) CBY0, CDY0
  write ( *, '(6x,4d16.8)' ) CBY1, CDY1

  return
end
subroutine comelp_test ( )

!*****************************************************************************80
!
!! comelp_test() tests comelp().
!
!  Discussion:
!
!    COMELP computes complete elliptic integrals K(k) and E(k).
!
!  Example:
!
!                  k         K(k)          E(K)
!                ---------------------------------
!                 .00      1.570796      1.570796
!                 .25      1.596242      1.545957
!                 .50      1.685750      1.467462
!                 .75      1.910990      1.318472
!                1.00       Ï            1.000000
!       ===================================================
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    15 January 2021
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Local:
!
!    K  --- Modulus k ( 0 <= k <= 1 )
!    CK --- K(k)
!    CE --- E(k)
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) ce
  real ( kind = rk8 ) ck
  real ( kind = rk8 ) hk
  integer i

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'comelp_test():'
  write ( *, '(a)' ) '  Test comelp().'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  COMELP computes complete elliptic integrals H(K), E(K).'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '    k         K(k)          E(K)'
  write ( *, '(a)' ) '  ---------------------------------'

  do i = 0, 4

    hk = real ( i, kind = rk8 ) / 4.0D+00

    call comelp ( hk, ck, ce )

    if ( hk /= 1.0D+00 ) then
      write ( *, '(2x,f5.2,2f14.6)' ) hk, ck, ce
    else
      write ( *, '(2x,f5.2,3x,a,3x,f14.6)' ) hk, 'Infinity', ce
    end if

  end do

  return
end
subroutine hygfx_test ( )

!*****************************************************************************80
!
!! hygfx_test() tests hygfx().
!
!  Example:
!
!     A    B     C     X     F(A,B,C,X)
!
!   -2.5  3.3   6.7  0.25  0.72356129D+00
!   -0.5  3.3   6.7  0.25  0.93610145D+00    
!    0.5  3.3   6.7  0.25  0.10689695D+01    
!    2.5  3.3   6.7  0.25  0.14051563D+01
!
!   -2.5  3.3   6.7  0.55  0.46961432D+00
!   -0.5  3.3   6.7  0.55  0.85187390D+00
!    0.5  3.3   6.7  0.55  0.11795358D+01
!    2.5  3.3   6.7  0.55  0.23999063D+01
!
!   -2.5  3.3   6.7  0.85  0.29106096D+00
!   -0.5  3.3   6.7  0.85  0.75543187D+00
!    0.5  3.3   6.7  0.85  0.13510497D+00
!    2.5  3.3   6.7  0.85  0.57381566D+01
!
!    3.3  6.7  -5.5  0.25  0.15090670D+05
!    3.3  6.7  -0.5  0.25 -0.21631479D+04
!    3.3  6.7   0.5  0.25  0.26451677D+03
!    3.3  6.7   4.5  0.25  0.41946916D+01
!
!    3.3  6.7  -5.5  0.55  0.10170778D+11
!    3.3  6.7  -0.5  0.55 -0.30854772D+07
!    3.3  6.7   0.5  0.55  0.11967860D+06
!    3.3  6.7   4.5  0.55  0.58092729D+02
!
!    3.3  6.7  -5.5  0.85  0.58682088D+19
!    3.3  6.7  -0.5  0.85 -0.10217370D+13
!    3.3  6.7   0.5  0.85  0.92370648D+10
!    3.3  6.7   4.5  0.85  0.20396914D+05 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 March 2012
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) a
  real ( kind = rk8 ) a_test(4)
  real ( kind = rk8 ) b
  real ( kind = rk8 ) c
  real ( kind = rk8 ) c_test(4)
  real ( kind = rk8 ) hf
  integer i
  integer k
  integer l
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x_test(3)

  save a_test
  save c_test
  save x_test

  data a_test / -2.5D+00, -0.5D+00, 0.5D+00, 2.5D+00 /
  data c_test / -5.5D+00, -0.5D+00, 0.5D+00, 4.5D+00/
  data x_test /  0.25D+00, 0.55D+00, 0.85D+00 /

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'HYGFX_TEST():'
  write ( *, '(a)' ) '  HYGFX() evaluates the hypergeometric function'
  write ( *, '(a)' ) '  2F1(a,b,c;x) for real arguments x.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '     A              B            C            X             2F1(A,B,C,X)'

  do l = 1, 3
    x = x_test(l)
    c = 6.7D+00
    b = 3.3D+00
    write ( *, '(a)' ) ' '
    do i = 1, 4
      a = a_test(i)
      call hygfx ( a, b, c, x, hf )
      write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
    end do
  end do

  do l = 1, 3
    x = x_test(l)
    write ( *, '(a)' ) ' '
    do k = 1, 4
      c = c_test(k)
      b = 6.7D+00
      a = 3.3D+00
      call hygfx ( a, b, c, x, hf )
      write ( *, '(4g14.6,g24.16)' ) a, b, c, x, hf
    end do
  end do

  return
end
subroutine jdzo_test ( )

!*********************************************************************72
!
!! jdzo_test() tests jdzo().
!
!  Discussion:
!
!    This program computes the zeros of Bessel functions
!    Jn(x) and Jn'(x), and arranges them in the order
!    of their values
!
!  Modified:
!
!    08 July 2024
!
!       Input :  NT    --- Number of total zeros ( NT  1200 )
!       Output:  ZO(L) --- Value of the L-th zero of Jn(x) and
!                          Jn'(x)
!                N(L)  --- n, order of Jn(x) or Jn'(x) associated
!                          with the L-th zero
!                M(L)  --- m, serial number of the zeros of Jn(x)
!                          or Jn'(x) associated with the L-th zero
!                          ( L is the serial number of all the
!                            zeros of Jn(x) and Jn'(x) )
!                P(L)  --- TM or TE, a code for designating the
!                          zeros of Jn(x) or Jn'(x)
!                          In the waveguide applications, the zeros
!                          of Jn(x) correspond to TM modes and those
!                          of Jn'(x) correspond to TE modes.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  integer, parameter :: nt = 9

  integer j1
  integer j2
  integer k
  integer k0
  integer ks
  integer m(nt)
  integer n(nt)
  character ( len = 4 ) p(nt)
!
!  It seems that zo() goes slightly out of bounds, and
!  damages the first entry of p.
!
  real ( kind = rk8 ) zo(nt)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'jdzo_test():'
  write ( *, '(a)' ) '  jdzo() evaluates zeros of Bessel functions'
  write ( *, '(a)' ) '  Jn and Jn''.'

  write ( *, '(a,i4)' ) '  Number of zeros = ', nt
  call jdzo ( nt, n, m, p, zo )
  write ( *, '(a)' ) ' '
  ks = nt / 101 + 1

  do k0 = 1, ks

    write(*,*)' Table           Zeros of Bessel', &
                     ' functions Jn(x) and Jn''(x)'
    write(*,*)
    write(*,*)' ----------------------------------', &
                     '----------------------------------'
    do k = 1, 50
      j1=100*(k0-1)+k+1
      j2=j1+50
      if (j1.le.nt+1.and.j2.le.nt+1) then
        write(*,65) j1-1, p(j1),n(j1),m(j1),zo(j1), &
          j2-1,p(j2),n(j2),m(j2),zo(j2)
      else if (j1.le.nt+1.and.j2.gt.nt+1) then
        write(*,65) j1-1, p(j1),n(j1),m(j1),zo(j1)
      end if

    end do
    write(*,*)' ----------------------------------', &
                    '----------------------------------'
  end do

65    FORMAT(1X,I4,3X,A2,I4,2H -,I2,F14.8,3X,1H|,2X,I4, &
              3X,A2,I4,2H -,I2,F14.8)
  return
end
subroutine sphj_test ( )

!*****************************************************************************80
!
!! sphj_test() tests sphj().
!
!  Discussion:
!
!    This routine computes the spherical Bessel functions jn(x) and 
!    jn'(x) using SPHJ.
!
!  Example:
!
!    x = 10.0
!    n          jn(x)               jn(x)
!    --------------------------------------------
!    0    -.5440211109D-01    -.7846694180D-01
!    1     .7846694180D-01    -.7009549945D-01
!    2     .7794219363D-01     .5508428371D-01
!    3    -.3949584498D-01     .9374053162D-01
!    4    -.1055892851D+00     .1329879757D-01
!    5    -.5553451162D-01    -.7226857814D-01
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 January 2016
!
!  Author:
!
!    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
!    This version by John Burkardt.
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!      
  implicit none

  integer, parameter :: rk8 = kind ( 1.0D+00 )

  real ( kind = rk8 ) dj(0:250)
  integer k
  integer n
  integer nm
  integer ns
  real ( kind = rk8 ) sj(0:250)
  real ( kind = rk8 ) x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MSPHJ'
  write ( *, '(a)' ) '  SPHJ evaluates spherical Bessel J functions'

  n = 5
  x = 0.905D+00

  if ( n <= 10 ) then
    ns = 1
  else
    ns = 5
  end if

  call sphj ( n, x, nm, sj, dj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   n      x                   jn(x)               jn''(x)'
  write ( *, '(a)' ) ''

  do k = 0, nm, ns
    write ( *, '(1X,I3,3D20.10)' ) k, x, sj(k), dj(k)
  end do

  n = 5
  x = 10.0D+00

  if ( n <= 10 ) then
    ns = 1
  else
    ns = 5
  end if

  call sphj ( n, x, nm, sj, dj )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   n      x                   jn(x)               jn''(x)'
  write ( *, '(a)' ) ''
  do k = 0, nm, ns
    write ( *, '(1X,I3,3D20.10)' ) k, x, sj(k), dj(k)
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
!    May 31 2001   9:45:54.872 AM
!
!  Modified:
!
!    31 May 2001
!
!  Author:
!
!    John Burkardt
!
  implicit none

  character ( len = 8 ) ampm
  integer d
  character ( len = 8 ) date
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  character ( len = 10 ) time
  integer values(8)
  integer y
  character ( len = 5 ) zone

  call date_and_time ( date, time, zone, values )

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

  write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

