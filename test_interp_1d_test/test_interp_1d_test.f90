program main

!*****************************************************************************80
!
!! test_interp_1d_test() tests test_interp_1d().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 September 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer nd

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_interp_1d_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test test_interp_1d().'
  write ( *, '(a)' ) '  The R8LIB library is needed.'

  call test01 ( )

  nd = 11
  call test02 ( nd )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test_interp_1d_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() simply prints the title of each function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer prob
  integer prob_num
  character ( len = 80 ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01():'
  write ( *, '(a)' ) '  Print the title of each function.'

  call p00_prob_num ( prob_num )
  
  write ( *, '(a)' ) ' '
  write ( *, '(a,i2,a)' ) '  There are ', prob_num, ' functions available:'
  write ( *, '(a)' ) '  Index  Title'
  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_title ( prob, title )

    write ( *, '(2x,i2,2x,a)' ) prob, trim ( title )

  end do

  return
end
subroutine test02 ( nd )

!*****************************************************************************80
!
!! test02() evaluates each test function at ND sample points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 August 2012
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer ND, the number of sample points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer prob
  integer prob_num
  character ( len = 80 ) title
  real ( kind = rk ) x(nd)
  real ( kind = rk ) y(nd)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  p00_f() can sample any test function.'

  call p00_prob_num ( prob_num )

  a = 0.0D+00
  b = 1.0D+00
  call r8vec_linspace ( nd, a, b, x )

  write ( *, '(a)' ) ' '

  do prob = 1, prob_num

    call p00_f ( prob, nd, x, y )
    write ( title, '(a,i2)' ) 'X, Y for problem ', prob
    call r8vec2_print ( nd, x, y, trim ( title ) )

  end do

  return
end
