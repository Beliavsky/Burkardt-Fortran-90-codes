program main

!*****************************************************************************80
!
!! toms446_test() tests toms446().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'toms446_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test TOMS446().'

  call cheby_test ( )
  call dfrnt_test ( )
  call echeb_test ( )
  call edcheb_test ( )
  call mltply_test ( )
  call ntgrt_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TOMS446_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine cheby_test ( )

!*****************************************************************************80
!
!! CHEBY_TEST tests CHEBY, which computes Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  integer i
  real ( kind = rk ) x(npl,nf)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'CHEBY_TEST'
  write ( *, '(a)' ) '  CHEBY computes the Chebyshev series for several functions.'

  call cheby ( nf, npl, functn, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine dfrnt_test ( )

!*****************************************************************************80
!
!! DFRNT_TEST tests DFRNT, which computes the Chebyshev series of a derivative.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  integer i
  integer j
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) x2(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'DFRNT_TEST'
  write ( *, '(a)' ) '  DFRNT computes the Chebyshev series for the derivative'
  write ( *, '(a)' ) '  of several functions.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf
    x2(1:npl) = x(1:npl,j)
    call dfrnt ( x2, npl, x2 )
    x(1:npl,j) = x2(1:npl)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Chebyshev series for d/dx of:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine echeb_test ( )

!*****************************************************************************80
!
!! ECHEB_TEST tests ECHEB, which evaluates a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  real ( kind = rk ) fval
  real ( kind = rk ) fxj(nf)
  integer j
  integer k
  integer nx
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) xval

  nx = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'ECHEB_TEST'
  write ( *, '(a)' ) '  ECHEB evaluates a Chebyshev series.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf

    x2(1:npl) = x(1:npl,j)

    write ( *, '(a)' ) ' '
    if ( j == 1 ) then
      write ( *, '(a)' ) '  Sin(x)'
    else if ( j == 2 ) then
      write ( *, '(a)' ) '  Cos(x)'
    else if ( j == 3 ) then
      write ( *, '(a)' ) '  Sin(2x)'
    else if ( j == 4 ) then
      write ( *, '(a)' ) '  Cos(2x)'
    else if ( j == 5 ) then
      write ( *, '(a)' ) '  x^5'
    end if

    write ( *, '(a)' ) ' '

    do k = 1, nx

      xval = 2.0D+00 * real ( k - 1, kind = rk ) / real ( nx - 1, kind = rk ) &
        - 1.0D+00

      call functn ( xval, fxj )

      call echeb ( xval, x2, npl, fval )

      write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

    end do

  end do

  return
end
subroutine edcheb_test ( )

!*****************************************************************************80
!
!! EDCHEB_TEST tests EDCHEB, which evaluates the derivative of a Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  real ( kind = rk ) fval
  real ( kind = rk ) fxj(nf)
  integer j
  integer k
  integer nx
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) xval

  nx = 6

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'EDCHEB_TEST'
  write ( *, '(a)' ) '  EDCHEB evaluates the derivative of a Chebyshev series.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf

    x2(1:npl) = x(1:npl,j)

    write ( *, '(a)' ) ' '
    if ( j == 1 ) then
      write ( *, '(a)' ) '  d/dx Sin(x)'
    else if ( j == 2 ) then
      write ( *, '(a)' ) '  d/dx Cos(x)'
    else if ( j == 3 ) then
      write ( *, '(a)' ) '  d/dx Sin(2x)'
    else if ( j == 4 ) then
      write ( *, '(a)' ) '  d/dx Cos(2x)'
    else if ( j == 5 ) then
      write ( *, '(a)' ) '  d/dx x^5'
    end if

    write ( *, '(a)' ) ' '

    do k = 1, nx

      xval = 2.0D+00 * real ( k - 1, kind = rk ) / real ( nx - 1, kind = rk ) &
        - 1.0D+00

      call functn_d ( xval, fxj )

      call edcheb ( xval, x2, npl, fval )

      write ( *, '(5(2x,f10.4))' ) xval, fxj(j), fval

    end do

  end do

  return
end
subroutine mltply_test ( )

!*****************************************************************************80
!
!! MLTPLY_TEST tests MLTPLY, which multiplies two Chebyshev series.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  integer i
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) x1(npl)
  real ( kind = rk ) x2(npl)
  real ( kind = rk ) x3(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'MLTPLY_TEST'
  write ( *, '(a)' ) '  MLTPLY computes the product of two Chebyshev series.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Multiply series for SIN(X) and COS(X)'
  write ( *, '(a)' ) '  and compare with series for 1/2*SIN(2X).'

  call cheby ( nf, npl, functn, x )

  x1(1:npl) = x(1:npl,1)
  x2(1:npl) = x(1:npl,2)
  x(1:npl,3) = 0.5D+00 * x(1:npl,3)

  call mltply ( x1, x2, npl, x3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)   1/2*Sin(2x)     RESULT'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:3), x3(i)
  end do

  return
end
subroutine ntgrt_test ( )

!*****************************************************************************80
!
!! NTGRT_TEST tests NTGRT, which computes the Chebyshev series of an indefinite integral.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5
  integer, parameter :: npl = 10

  external functn
  integer i
  integer j
  real ( kind = rk ) x(npl,nf)
  real ( kind = rk ) x2(npl)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NTGRT_TEST'
  write ( *, '(a)' ) '  NTGRT computes the Chebyshev series for the indefinite'
  write ( *, '(a)' ) '  integral of several functions.'

  call cheby ( nf, npl, functn, x )

  do j = 1, nf
    x2(1:npl) = x(1:npl,j)
    call ntgrt ( x2, npl, x2 )
    x(1:npl,j) = x2(1:npl)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Series for indefinite integral of:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '        Sin(x)      Cos(x)    Sin(2x)     Cos(2x)       X^5'
  write ( *, '(a)' ) ' '

  do i = 1, npl
    write ( *, '(5(2x,f10.4))' ) x(i,1:nf)
  end do

  return
end
subroutine functn ( x, fxj )

!*****************************************************************************80
!
!! FUNCTN evaluates several functions at X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) FXJ(5), the derivative values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5

  real ( kind = rk ) fxj(nf)
  real ( kind = rk ) x

  fxj(1) = sin ( x )
  fxj(2) = cos ( x )
  fxj(3) = sin ( 2.0D+00 * x )
  fxj(4) = cos ( 2.0D+00 * x )
  fxj(5) = x**5

  return
end
subroutine functn_d ( x, fxj )

!*****************************************************************************80
!
!! FUNCTN_D evaluates the derivatives of several functions at X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 September 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) FXJ(5), the derivative values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nf = 5

  real ( kind = rk ) fxj(nf)
  real ( kind = rk ) x

  fxj(1) =  cos ( x )
  fxj(2) = -sin ( x )
  fxj(3) =  2.0D+00 * cos ( 2.0D+00 * x )
  fxj(4) = -2.0D+00 * sin ( 2.0D+00 * x )
  fxj(5) =  5.0D+00 * x**4

  return
end

