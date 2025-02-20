program main

!*****************************************************************************80
!
!! beta_nc_test() tests beta_nc().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'beta_nc_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test beta_nc().'

  call test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'beta_nc_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests beta_noncentral_cdf().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    03 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) error_max
  real ( kind = rk ) fx
  real ( kind = rk ) fx2
  real ( kind = rk ) lambda
  integer n_data
  real ( kind = rk ) x

  error_max = 1.0D-10

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
  write ( *, '(a)' ) '  incomplete Beta Function against values'
  write ( *, '(a)' ) '  computed by beta_noncentral_cdf().'
  write ( *, '(a)' ) ' '
  write ( *, '(a,a)' ) &
    '      A        B     LAMBDA        X       ', &
    ' CDF                         CDF                    DIFF'
  write ( *, '(a,a)' ) &
    '                                           ', &
    '(tabulated)       (BETA_NC)'
  write ( *, '(a)' ) ' '

  n_data = 0

  do

    call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

    if ( n_data == 0 ) then
      exit
    end if

    call beta_noncentral_cdf ( a, b, lambda, x, error_max, fx2 )

    write ( *, &
    '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) &
    a, b, lambda, x, fx, fx2, abs ( fx - fx2 )

  end do

  return
end
 
