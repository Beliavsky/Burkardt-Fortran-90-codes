program main

!*****************************************************************************80
!
!! simplex_gm_rule_test() tests simplex_gm_rule().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'simplex_gm_rule_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test SIMPLEX_GM_RULE().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
  call test09 ( )
  call test10 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'SIMPLEX_GM_RULE_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SIMPLEX_UNIT_TO_GENERAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 2
  integer, parameter :: n = 10

  integer j
  real ( kind = rk ) phy(m,n)
  real ( kind = rk ) phy_unit(m,m+1)
  real ( kind = rk ) ref(m,n)
  integer seed
  real ( kind = rk ), dimension(m,m+1) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, &
    2.0D+00, 5.0D+00 /), (/ m, m + 1 /) )
  real ( kind = rk ), dimension(m,m+1) :: t_unit = reshape ( (/ &
    0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00 /), (/ m, m + 1 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_TO_GENERAL'
  write ( *, '(a)' ) '  maps points in the unit simplex to a general simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we consider a simplex in 2D, a triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vertices of the general triangle are:'
  write ( *, '(a)' ) ' '
  do j = 1, m + 1
    write ( *, '(2x,f8.4,2x,f8.4)' ) t(1:m,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI     ETA )   ( X       Y  )'
  write ( *, '(a)' ) ' '

  call simplex_unit_to_general ( m, m+1, t, t_unit, phy_unit )

  do j = 1, m + 1

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) t_unit(1:m,j), phy_unit(1:m,j)

  end do

  call simplex_unit_sample ( m, n, seed, ref )

  call simplex_unit_to_general ( m, n, t, ref, phy )

  do j = 1, n

    write ( *, '(2x,2f8.4,2x,2f8.4)' ) ref(1:m,j), phy(1:m,j)

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SIMPLEX_UNIT_TO_GENERAL.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: n = 10

  integer j
  real ( kind = rk ) phy(m,n)
  real ( kind = rk ) phy_unit(m,m+1)
  real ( kind = rk ) ref(m,n)
  integer seed
  real ( kind = rk ), dimension(m,m+1) :: t = reshape ( (/ &
    1.0D+00, 1.0D+00, 1.0D+00, &
    3.0D+00, 1.0D+00, 1.0D+00, &
    1.0D+00, 4.0D+00, 1.0D+00, &
    1.0D+00, 1.0D+00, 5.0D+00 /), (/ m, m + 1 /) )
  real ( kind = rk ), dimension(m,m+1) :: t_unit = reshape ( (/ &
    0.0D+00, 0.0D+00, 0.0D+00, &
    1.0D+00, 0.0D+00, 0.0D+00, &
    0.0D+00, 1.0D+00, 0.0D+00, &
    0.0D+00, 0.0D+00, 1.0D+00 /), (/ m, m + 1 /) )

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  SIMPLEX_UNIT_TO_GENERAL'
  write ( *, '(a)' ) '  maps points in the unit simplex to a general simplex.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Here we consider a simplex in 3D, a tetrahedron.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The vertices of the general tetrahedron are:'
  write ( *, '(a)' ) ' '
  do j = 1, m + 1
    write ( *, '(2x,f8.4,2x,f8.4,2x,f8.4)' ) t(1:m,j)
  end do

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   (  XSI     ETA     MU )    ( X       Y       Z )'
  write ( *, '(a)' ) ' '

  call simplex_unit_to_general ( m, m+1, t, t_unit, phy_unit )

  do j = 1, m + 1

    write ( *, '(2x,3f8.4,2x,3f8.4)' ) t_unit(1:m,j), phy_unit(1:m,j)

  end do

  call simplex_unit_sample ( m, n, seed, ref )

  call simplex_unit_to_general ( m, n, t, ref, phy )

  do j = 1, n

    write ( *, '(2x,3f8.4,2x,3f8.4)' ) ref(1:m,j), phy(1:m,j)

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests GM_RULE_SIZE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    08 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: test_num = 4

  integer m
  integer, dimension ( test_num ) :: m_test = (/ &
    2, 3, 5, 10 /)
  integer degree
  integer n
  integer rule
  integer test

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  GM_RULE_SIZE returns N, the number of points'
  write ( *, '(a)' ) '  associated with a Grundmann-Moeller quadrature rule'
  write ( *, '(a)' ) '  for the unit simplex of dimension M'
  write ( *, '(a)' ) '  with rule index RULE'
  write ( *, '(a)' ) '  and degree of exactness DEGREE = 2*RULE+1.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   M      RULE    DEGREE N'

  do test = 1, test_num

    m = m_test(test)

    write ( *, '(a)' ) ' '

    do rule = 0, 5

      call gm_rule_size ( rule, m, n )
      degree = 2 * rule + 1

      write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) m, rule, degree, n

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests GM_UNIT_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer point
  integer n
  integer rule
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  GM_UNIT_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional unit simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'

  m = 3
  rule = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use M = ', m
  write ( *, '(a,i8)' ) '  RULE = ', rule
  write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

  call gm_rule_size ( rule, m, n )

  allocate ( w(1:n) )
  allocate ( x(1:m,1:n) )

  call gm_unit_rule_set ( rule, m, n, w, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     POINT        W             X             Y             Z'
  write ( *, '(a)' ) ' '

  do point = 1, n
    write ( *, '(2x,i8,2x,f12.6,2x,f12.6,2x,f12.6,2x,f12.6)' ) &
      point, w(point), x(1:m,point)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests GM_UNIT_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: test_num = 4

  integer m
  integer, dimension ( test_num ) :: m_test = (/ &
    2, 3, 5, 10 /)
  integer n
  integer rule
  integer test
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ) w_sum
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  GM_UNIT_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional unit simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we compute various rules, and simply'
  write ( *, '(a)' ) '  report the number of points, and the sum of weights.'

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   M      RULE    N  WEIGHT SUM'

  do test = 1, test_num

    m = m_test(test)

    write ( *, '(a)' ) ' '

    do rule = 0, 5

      call gm_rule_size ( rule, m, n )

      allocate ( w(1:n) )
      allocate ( x(1:m,1:n) )

      call gm_unit_rule_set ( rule, m, n, w, x )

      w_sum = sum ( w(1:n) )

      write ( *, '(2x,i8,2x,i8,2x,i8,2x,g24.16)' ) &
        m, rule, n, w_sum

      deallocate ( w )
      deallocate ( x )

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests GM_UNIT_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer point
  integer n
  integer rule
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  character ( len = 12 ) w_file
  integer w_unit
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x
  character ( len = 12 ) x_file
  integer x_unit

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  GM_UNIT_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional unit simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, we write a rule to a file.'

  m = 3
  rule = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use M = ', m
  write ( *, '(a,i8)' ) '  RULE = ', rule
  write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

  call gm_rule_size ( rule, m, n )

  allocate ( w(1:n) )
  allocate ( x(1:m,1:n) )

  call gm_unit_rule_set ( rule, m, n, w, x )

  call get_unit ( w_unit )

  write ( w_file, '(a2,i1,a,i1,a7)' ) 'gm', rule, '_', m, 'd_w.txt'

  open ( unit = w_unit, file = w_file, status = 'replace' )

  do point = 1, n
    write ( w_unit, '(f20.16)' ) w(point)
  end do

  close ( unit = w_unit )

  call get_unit ( x_unit )

  write ( x_file, '(a2,i1,a,i1,a7)' ) 'gm', rule, '_', m, 'd_x.txt'

  open ( unit = x_unit, file = x_file, status = 'replace' )

  do point = 1, n
    write ( x_unit, '(3f20.16)' ) x(1:m,point)
  end do

  close ( unit = x_unit )

  write ( *, '(a,i2,a)' ) '  Wrote rule ', rule, ' to "' &
    // trim ( w_file ) // '" and "' // trim ( x_file ) // '".'

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! TEST07 tests GM_UNIT_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 5

  integer degree
  integer, parameter :: degree_max = 4
  integer expon(m)
  integer h
  real ( kind = rk ), allocatable, dimension ( : ) :: mono
  logical more
  integer n
  real ( kind = rk ) quad_error
  integer rule
  integer, parameter :: rule_max = 3
  integer t
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  GM_UNIT_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional unit simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, look at all the monomials up to'
  write ( *, '(a)' ) '  some maximum degree, choose a few low order rules'
  write ( *, '(a)' ) '  and determine the quadrature error for each.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use M = ', m

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      Rule     Order     Quad_Error'
  write ( *, '(a)' ) ' '

  do degree = 0, degree_max

    more = .false.

    do

      call comp_next ( degree, m, expon, more, h, t )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i1,a,i1,a,i1,a,i1,a,i1)' ) &
        '  F(X) = X1^', expon(1), ' * X2^', expon(2), ' * X3^', expon(3), &
              ' * X4^', expon(4), ' * X5^', expon(5)

      write ( *, '(a)' ) ' '

      do rule = 0, rule_max

        call gm_rule_size ( rule, m, n )

        allocate ( mono(1:n) )
        allocate ( w(1:n) )
        allocate ( x(1:m,1:n) )

        call gm_unit_rule_set ( rule, m, n, w, x )

        call simplex_unit_monomial_quadrature ( m, expon, n, &
          x, w, quad_error )

        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) rule, n, quad_error

        deallocate ( mono )
        deallocate ( w )
        deallocate ( x )

      end do

      if ( .not. more ) then
        exit
      end if

    end do

  end do

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! TEST08 tests GM_GENERAL_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 July 2007
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer j
  integer point
  integer n
  integer rule
  real ( kind = rk ) :: t(3,4) = reshape ( (/ &
    1.0, 0.0, 0.0, &
    2.0, 0.0, 0.0, &
    1.0, 2.0, 0.0, &
    1.0, 0.0, 3.0 /), (/ 3, 4 /) )
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  GM_GENERAL_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional general simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'

  m = 3
  rule = 2

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  Here we use M = ', m
  write ( *, '(a,i8)' ) '  RULE = ', rule
  write ( *, '(a,i8)' ) '  DEGREE = ', 2 * rule + 1

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Simplex vertices:'
  write ( *, '(a)' ) ''
  do j = 1, 4
    write ( *, '(3g14.6)' ) t(1:3,j)
  end do
  
  call gm_rule_size ( rule, m, n )

  allocate ( w(1:n) )
  allocate ( x(1:m,1:n) )

  call gm_general_rule_set ( rule, m, n, t, w, x )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '     POINT        W             X             Y             Z'
  write ( *, '(a)' ) ' '

  do point = 1, n
    write ( *, '(2x,i8,2x,f12.6,2x,f12.6,2x,f12.6,2x,f12.6)' ) &
      point, w(point), x(1:m,point)
  end do

  deallocate ( w )
  deallocate ( x )

  return
end
subroutine test09 ( )

!*****************************************************************************80
!
!! TEST09 tests GM_UNIT_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3

  integer e(m)
  integer :: e_test(m,10) = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    1, 1, 0, &
    1, 0, 1, &
    0, 2, 0, &
    0, 1, 1, &
    0, 0, 2 /), (/ m, 10 /) )
  integer k
  integer n
  real ( kind = rk ) result(10)
  integer rule
  real ( kind = rk ), allocatable, dimension ( : ) :: value
  real ( kind = rk ) volume
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  GM_UNIT_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional unit simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, look at all the monomials up to'
  write ( *, '(a)' ) '  some maximum degree, choose a few low order rules'
  write ( *, '(a)' ) '  and determine the quadrature error for each.'

  call simplex_unit_volume ( m, volume )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Simplex volume = ', volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1               X               Y ' // &
    '              Z               X^2              XY             XZ' // &
    '              Y^2             YZ               Z^2'
  write ( *, '(a)' ) ' '

  do rule = 0, 5

    call gm_rule_size ( rule, m, n )
    allocate ( value(1:n) )
    allocate ( w(1:n) )
    allocate ( x(1:m,1:n) )

    call gm_unit_rule_set ( rule, m, n, w, x )

    do k = 1, 10

      e(1:m) = e_test(1:m,k)

      call monomial_value ( m, n, e, x, value )

      result(k) = dot_product ( w(1:n), value(1:n) )

    end do

    write ( *, '(2x,i8,10(2x,g14.6))' ) n, result(1:10)

    deallocate ( value )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
subroutine test10 ( )

!*****************************************************************************80
!
!! TEST10 tests GM_GENERAL_RULE_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 March 2017
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3

  integer e(m)
  integer :: e_test(m,10) = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    1, 1, 0, &
    1, 0, 1, &
    0, 2, 0, &
    0, 1, 1, &
    0, 0, 2 /), (/ m, 10 /) )
  integer j
  integer k
  integer n
  real ( kind = rk ) result(10)
  integer rule
  real ( kind = rk ) :: t(3,4) = reshape ( (/ &
    1.0, 0.0, 0.0, &
    2.0, 0.0, 0.0, &
    1.0, 2.0, 0.0, &
    1.0, 0.0, 3.0 /), (/ 3, 4 /) )
  real ( kind = rk ), allocatable, dimension ( : ) :: value
  real ( kind = rk ) volume
  real ( kind = rk ), allocatable, dimension ( : ) :: w
  real ( kind = rk ), allocatable, dimension ( :, : ) :: x

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  GM_GENERAL_RULE_SET determines the weights and abscissas'
  write ( *, '(a)' ) '  of a Grundmann-Moeller quadrature rule for'
  write ( *, '(a)' ) '  the M dimensional general simplex,'
  write ( *, '(a)' ) '  using a rule of index RULE,'
  write ( *, '(a)' ) '  which will have degree of exactness 2*RULE+1.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  In this test, look at all the monomials up to'
  write ( *, '(a)' ) '  some maximum degree, choose a few low order rules'
  write ( *, '(a)' ) '  and determine the quadrature error for each.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Simplex vertices:'
  write ( *, '(a)' ) ''
  do j = 1, 4
    write ( *, '(3g14.6)' ) t(1:3,j)
  end do

  call simplex_general_volume ( m, t, volume )
  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Simplex volume = ', volume

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '         N        1               X               Y ' // &
    '              Z               X^2              XY             XZ' // &
    '              Y^2             YZ               Z^2'
  write ( *, '(a)' ) ' '

  do rule = 0, 5

    call gm_rule_size ( rule, m, n )
    allocate ( value(1:n) )
    allocate ( w(1:n) )
    allocate ( x(1:m,1:n) )

    call gm_general_rule_set ( rule, m, n, t, w, x )

    do k = 1, 10

      e(1:m) = e_test(1:m,k)

      call monomial_value ( m, n, e, x, value )

      result(k) = dot_product ( w(1:n), value(1:n) )

    end do

    write ( *, '(2x,i8,10(2x,g14.6))' ) n, result(1:10)

    deallocate ( value )
    deallocate ( w )
    deallocate ( x )

  end do

  return
end
