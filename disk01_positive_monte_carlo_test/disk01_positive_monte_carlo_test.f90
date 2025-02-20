program main

!*****************************************************************************80
!
!! disk01_positive_monte_carlo_test() tests disk01_positive_monte_carlo().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 May 2016
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'disk01_positive_monte_carlo_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test disk01_positive_monte_carlo().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'disk01_positive_monte_carlo_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() uses disk01_positive_sample() with an increasing number of points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 January 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00

  real ( kind = rk ) disk01_positive_area
  integer e(2)
  real ( kind = rk ) err
  real ( kind = rk ) exact
  integer i
  integer j
  integer n
  real ( kind = rk ) q
  real ( kind = rk ), allocatable :: value(:)
  real ( kind = rk ), allocatable :: x(:,:)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Use disk01_positive_sample() to estimate integrals'
  write ( *, '(a)' ) '  in the unit quarter disk.'

  do i = 0, 4
    e(1) = i
    do j = 0, 4 - e(1) 
      e(2) = j
      call disk01_positive_monomial_integral ( e, exact )
      write ( *, '(a)' ) ''
      write ( *, '(a,i1,a,i1)' ) '  Estimate integral of X^', e(1), 'Y^', e(2)
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         N        Estimate       Error'
      write ( *, '(a)' ) ''

      n = 1
 
      do while ( n <= 65536 )

        allocate ( x(1:2,1:n) )
        call disk01_positive_sample ( n, x )

        allocate ( value(1:n) )
        call monomial_value ( 2, n, e, x, value )

        q = disk01_positive_area ( ) * sum ( value(1:n) ) / real ( n, kind = rk )

        err = abs ( q - exact )
        write ( *, '(2x,i8,2x,g14.6,2x,e10.2)' ) n, q, err

        deallocate ( value )
        deallocate ( x )

        n = 2 * n

      end do

      write ( *, '(2x,a8,2x,g14.6,2x,e10.2)' ) '  Exact:', exact, 0.0D+00

    end do

  end do

  return
end
