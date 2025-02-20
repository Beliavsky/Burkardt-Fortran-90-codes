program main

!*****************************************************************************80
!
!! burgers_solution_test() tests burgers_solution().
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
  write ( *, '(a)' ) 'burgers_solution_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test burgers_solution().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'burgers_solution_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests sets up a small test case.
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

  integer, parameter :: vtn = 11
  integer, parameter :: vxn = 11

  character ( len = 80 ) filename
  real ( kind = rk ) nu
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) thi
  real ( kind = rk ) tlo
  real ( kind = rk ) vu(vxn,vtn)
  real ( kind = rk ) vt(vtn)
  real ( kind = rk ) vx(vxn)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo

  nu = 0.01D+00 / r8_pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01()'
  write ( *, '(a)' ) '  burgers_viscous_time_exact1() computes'
  write ( *, '(a)' ) '  exact solution #1 to the Burgers equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
  write ( *, '(a,i4)' ) '  NX = ', vxn
  write ( *, '(a,i4)' ) '  NT = ', vtn

  xlo = -1.0D+00
  xhi = +1.0D+00
  call r8vec_even ( vxn, xlo, xhi, vx )
  call r8vec_print ( vxn, vx, '  X grid points:' )

  tlo = 0.0D+00
  thi = 3.0D+00 / r8_pi
  call r8vec_even ( vtn, tlo, thi, vt )
  call r8vec_print ( vtn, vt, '  T grid points:' )

  call burgers_viscous_time_exact1 ( nu, vxn, vx, vtn, vt, vu )

  call r8mat_print ( vxn, vtn, vu, '  U(X,T) at grid points:' )

  filename = 'burgers_solution_test01.txt'

  call r8mat_write ( filename, vxn, vtn, vu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( filename ) // '".'

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() sets up a finer test case.
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

  integer, parameter :: vtn = 41
  integer, parameter :: vxn = 41

  character ( len = 80 ) filename
  real ( kind = rk ) nu
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) thi
  real ( kind = rk ) tlo
  real ( kind = rk ) vu(vxn,vtn)
  real ( kind = rk ) vt(vtn)
  real ( kind = rk ) vx(vxn)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo

  nu = 0.01D+00 / r8_pi

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02():'
  write ( *, '(a)' ) '  burgers_viscous_time_exact1() computes'
  write ( *, '(a)' ) '  exact solution #1 to the Burgers equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
  write ( *, '(a,i4)' ) '  NX = ', vxn
  write ( *, '(a,i4)' ) '  NT = ', vtn

  xlo = -1.0D+00
  xhi = +1.0D+00
  call r8vec_even ( vxn, xlo, xhi, vx )
  call r8vec_print ( vxn, vx, '  X grid points:' )

  tlo = 0.0D+00
  thi = 3.0D+00 / r8_pi
  call r8vec_even ( vtn, tlo, thi, vt )
  call r8vec_print ( vtn, vt, '  T grid points:' )

  call burgers_viscous_time_exact1 ( nu, vxn, vx, vtn, vt, vu )

  filename = 'burgers_solution_test02.txt'

  call r8mat_write ( filename, vxn, vtn, vu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( filename ) // '".'

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() sets up a small test case.
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

  integer, parameter :: vtn = 11
  integer, parameter :: vxn = 11

  character ( len = 80 ) filename
  real ( kind = rk ) nu
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) thi
  real ( kind = rk ) tlo
  real ( kind = rk ) vu(vxn,vtn)
  real ( kind = rk ) vt(vtn)
  real ( kind = rk ) vx(vxn)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo

  nu = 0.5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  burgers_viscous_time_exact2() computes'
  write ( *, '(a)' ) '  exact solution #2 to the Burgers equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
  write ( *, '(a,i4)' ) '  NX = ', vxn
  write ( *, '(a,i4)' ) '  NT = ', vtn

  xlo = 0.0D+00
  xhi = 2.0D+00 * r8_pi
  call r8vec_even ( vxn, xlo, xhi, vx )
  call r8vec_print ( vxn, vx, '  X grid points:' )

  tlo = 0.0D+00
  thi = 1.0D+00
  call r8vec_even ( vtn, tlo, thi, vt )
  call r8vec_print ( vtn, vt, '  T grid points:' )

  call burgers_viscous_time_exact2 ( nu, vxn, vx, vtn, vt, vu )

  call r8mat_print ( vxn, vtn, vu, '  U(X,T) at grid points:' )

  filename = 'burgers_solution_test03.txt'

  call r8mat_write ( filename, vxn, vtn, vu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( filename ) // '".'

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! test04() sets up a finer test case.
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

  integer, parameter :: vtn = 41
  integer, parameter :: vxn = 41

  character ( len = 80 ) filename
  real ( kind = rk ) nu
  real ( kind = rk ), parameter :: r8_pi = 3.141592653589793D+00
  real ( kind = rk ) thi
  real ( kind = rk ) tlo
  real ( kind = rk ) vu(vxn,vtn)
  real ( kind = rk ) vt(vtn)
  real ( kind = rk ) vx(vxn)
  real ( kind = rk ) xhi
  real ( kind = rk ) xlo

  nu = 0.5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test04()'
  write ( *, '(a)' ) '  burgers_viscous_time_exact2() computes'
  write ( *, '(a)' ) '  exact solution #2 to the Burgers equation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Viscosity NU = ', nu
  write ( *, '(a,i4)' ) '  NX = ', vxn
  write ( *, '(a,i4)' ) '  NT = ', vtn

  xlo = 0.0D+00
  xhi = 2.0D+00 * r8_pi
  call r8vec_even ( vxn, xlo, xhi, vx )
  call r8vec_print ( vxn, vx, '  X grid points:' )

  tlo = 0.0D+00
  thi = 1.0D+00
  call r8vec_even ( vtn, tlo, thi, vt )
  call r8vec_print ( vtn, vt, '  T grid points:' )

  call burgers_viscous_time_exact2 ( nu, vxn, vx, vtn, vt, vu )

  filename = 'burgers_solution_test04.txt'

  call r8mat_write ( filename, vxn, vtn, vu )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Data written to file "' // trim ( filename ) // '".'

  return
end
 
