program main

!*****************************************************************************80
!
!! sphere_lebedev_rule_test() tests sphere_lebedev_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_lebedev_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test sphere_lebedev_rule().'

  call test01 ( )
  call test02 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_lebedev_rule_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests AVAILABLE_TABLE, ORDER_TABLE, PRECISION_TABLE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!   12 September 2010
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer available
  integer available_table
  integer order
  integer order_table
  integer precision
  integer precision_table
  integer rule
  integer, parameter :: rule_max = 65

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  List Lebedev rule properties.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Rule Avail Order  Prec'
  write ( *, '(a)' ) ' '
  do rule = 1, rule_max
    available = available_table ( rule )
    order = order_table ( rule )
    precision = precision_table ( rule )
    write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) rule, available, order, precision
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests the SPHERE_LEBEDEV_RULE functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 September 2010
!
!  Author:
!
!    Dmitri Laikov
!
!  Reference:
!
!    Vyacheslav Lebedev, Dmitri Laikov,
!    A quadrature formula for the sphere of the 131st
!    algebraic order of accuracy,
!    Russian Academy of Sciences Doklady Mathematics,
!    Volume 59, Number 3, 1999, pages 477-481.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: nmax = 65
  integer, parameter :: mmax = 5897
! integer, parameter :: mmax = ( ( nmax * 2 + 3 ) * ( nmax * 2 + 3 ) / 3 )

  real ( kind = rk ) alpha
  integer available
  integer available_table
  real ( kind = rk ) beta
  real ( kind = rk ) err
  real ( kind = rk ) err_max
  integer i
  real ( kind = rk ) integral_exact
  real ( kind = rk ) integral_approx
  integer j
  integer k
  integer m
  integer n
  integer order
  integer order_table
  integer precision_table
  real ( kind = rk ) s(0:nmax+1)
  real ( kind = rk ) w(mmax)
  real ( kind = rk ) x(mmax)
  real ( kind = rk ) xn(mmax,0:nmax)
  real ( kind = rk ) y(mmax)
  real ( kind = rk ) yn(mmax,0:nmax)
  real ( kind = rk ) z(mmax)
  real ( kind = rk ) zn(mmax,0:nmax)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Generate each available rule and test for accuracy.'

  do n = 1, nmax

    available = available_table ( n )

    if ( available == 1 ) then

      order = order_table ( n )

      call ld_by_order ( order, x, y, z, w )

      s(0) = 1.0D+00
      do k = 1, n + 1
        s(k) = ( 2 * k - 1 ) * s(k-1)
      end do
!
!  For each abscissa X(M), compute the values 1, X(M)^2, X(M)^4, ..., X(M)^2*N.
!
      do m = 1, order
        xn(m,0) = 1.0D+00
        yn(m,0) = 1.0D+00
        zn(m,0) = 1.0D+00
        do k = 1, n
          xn(m,k) = xn(m,k-1) * x(m) * x(m)
          yn(m,k) = yn(m,k-1) * y(m) * y(m)
          zn(m,k) = zn(m,k-1) * z(m) * z(m)
        end do
      end do

      err_max = 0.0D+00
      do i = 0, n
        do j = 0, n - i
          k = n - i - j
!
!  Apply Lebedev rule to x^2i y^2j z^2k.
!
          integral_approx = 0.0D+00
          do m = 1, order
            integral_approx = integral_approx &
              + w(m) * xn(m,i) * yn(m,j) * zn(m,k)
          end do
!
!  Compute exact value of integral (aside from factor of 4 pi!).
!
          integral_exact = s(i) * s(j) * s(k) / s(1+i+j+k)
!
!  Record the maximum error for this rule.
!
          err = abs ( ( integral_approx - integral_exact ) / integral_exact )
          err_max = max ( err_max, err )

        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,i4,a,g14.6)' ) &
        '  Order = ', order, &
        '  LMAXW = ', precision_table ( n ), &
        '  max error = ', err_max
!
!  Convert (X,Y,Z) to (Theta,Phi) and print the data.
!
      if ( order <= 50 ) then
        do m = 1, order
          call xyz_to_tp ( x(m), y(m), z(m), alpha, beta )
          write ( *, '(g24.15,g24.15,g24.15)' ) alpha, beta, w(m)
        end do
      end if

    end if

  end do

  return
end
