program main

!*****************************************************************************80
!
!! sphere_triangle_quad_test() tests sphere_triangle_quad().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 November 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_triangle_quad_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test sphere_triangle_quad().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'sphere_triangle_quad_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test01 ( )

!*****************************************************************************80
!
!! TEST01 tests SPHERE01_TRIANGLE_QUAD_01, 02, 03.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  real ( kind = rk ) result_01
  real ( kind = rk ) result_02
  real ( kind = rk ) result_03
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) &
     '  Approximate the integral of a function on a random spherical triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QUAD_01 uses centroids of spherical triangles.'
  write ( *, '(a)' ) '  QUAD_02 uses vertices of spherical triangles.'
  write ( *, '(a)' ) '  QUAD_03 uses midsides of spherical triangles.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QUAD_01      QUAD_02      QUAD_03'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    call sphere01_triangle_quad_01 ( v1, v2, v3, polyterm_value_3d, result_01 )

    call sphere01_triangle_quad_02 ( v1, v2, v3, polyterm_value_3d, result_02 )

    call sphere01_triangle_quad_03 ( v1, v2, v3, polyterm_value_3d, result_03 )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      result_01, result_02, result_03

  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! TEST02 tests SPHERE01_TRIANGLE_QUAD_00.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  integer, parameter :: n_mc1 = 1000
  integer, parameter :: n_mc2 = 10000
  integer, parameter :: n_mc3 = 100000
  real ( kind = rk ) result_mc1
  real ( kind = rk ) result_mc2
  real ( kind = rk ) result_mc3
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Approximate an integral on a random spherical triangle.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QUAD_MC1 uses a Monte Carlo method with', n_mc1, ' points.'
  write ( *, '(a)' ) '  QUAD_MC2 uses a Monte Carlo method with', n_mc2, ' points.'
  write ( *, '(a)' ) '  QUAD_MC3 uses a Monte Carlo method with', n_mc3, ' points.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  QUAD_MC1      QUAD_MC2      QUAD_MC3'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    call sphere01_triangle_quad_00 ( n_mc1, v1, v2, v3, polyterm_value_3d, &
      seed, result_mc1 )
    call sphere01_triangle_quad_00 ( n_mc2, v1, v2, v3, polyterm_value_3d, &
      seed, result_mc2 )
    call sphere01_triangle_quad_00 ( n_mc3, v1, v2, v3, polyterm_value_3d, &
      seed, result_mc3 )

    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) &
      result_mc1, result_mc2, result_mc3

  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! TEST03 tests SPHERE01_TRIANGLE_QUAD_ICOS1C.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) best
  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ) error
  integer factor
  integer factor_log
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  integer node_num
  real ( kind = rk ) result
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SPHERE01_TRIANGLE_QUAD_ICOS1C approximates the'
  write ( *, '(a)' ) '  integral of a function over a spherical triangle on'
  write ( *, '(a)' ) '  the surface of the unit sphere using a centroid rule.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  We do not have an exact result, so we compare each'
  write ( *, '(a)' ) '  estimate to the final one.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FACTOR   N   RESULT'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 2**11
    call sphere01_triangle_quad_icos1c ( v1, v2, v3, factor, &
      polyterm_value_3d, node_num, best )

    factor = 1
    do factor_log = 0, 9

      call sphere01_triangle_quad_icos1c ( v1, v2, v3, factor, &
        polyterm_value_3d, node_num, result )

      error = abs ( result - best )

      write ( *, '(2x,i4,2x,i8,2x,g16.8,2x,g10.2)' ) &
        factor, node_num, result, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! TEST04 tests SPHERE01_TRIANGLE_QUAD_ICOS1M.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) best
  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ) error
  integer factor
  integer factor_log
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  integer node_num
  real ( kind = rk ) result
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  SPHERE01_TRIANGLE_QUAD_ICOS1M approximate the'
  write ( *, '(a)' ) '  integral of a function over a spherical triangle'
  write ( *, '(a)' ) '  on the surface of the unit sphere using a'
  write ( *, '(a)' ) '  midpoint rule.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  We do not have an exact result, so we compare each'
  write ( *, '(a)' ) '  estimate to the final one.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  FACTOR   N   RESULT'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 2**11
    call sphere01_triangle_quad_icos1m ( v1, v2, v3, factor, &
      polyterm_value_3d, node_num, best )

    factor = 1
    do factor_log = 0, 9

      call sphere01_triangle_quad_icos1m ( v1, v2, v3, factor, &
        polyterm_value_3d, node_num, result )

      error = abs ( result - best )

      write ( *, '(2x,i4,2x,i8,2x,g16.8,2x,g10.2)' ) &
        factor, node_num, result, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! TEST05 tests SPHERE01_TRIANGLE_QUAD_ICOS1V.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) best
  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ) error
  integer factor
  integer factor_log
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  integer node_num
  real ( kind = rk ) result
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  SPHERE01_TRIANGLE_QUAD_ICOS1V approximates the'
  write ( *, '(a)' ) '  integral of a function over a spherical triangle'
  write ( *, '(a)' ) '  on the surface of the unit sphere using a vertex rule.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  We do not have an exact result, so we compare each'
  write ( *, '(a)' ) '  estimate to the final one.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACTOR   N   RESULT'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 2**11
    call sphere01_triangle_quad_icos1v ( v1, v2, v3, factor, &
      polyterm_value_3d, node_num, best )

    factor = 1
    do factor_log = 0, 9

      call sphere01_triangle_quad_icos1v ( v1, v2, v3, factor, &
        polyterm_value_3d, node_num, result )

      error = abs ( result - best )

      write ( *, '(2x,i4,2x,i8,2x,g16.8,2x,g10.2)' ) &
        factor, node_num, result, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! TEST06 tests SPHERE01_TRIANGLE_QUAD_ICOS2V.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    23 April 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) best
  integer e(3)
  integer, dimension ( 3, 17 ) :: e_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 1, 0, &
    0, 0, 1, &
    2, 0, 0, &
    0, 2, 2, &
    2, 2, 2, &
    0, 2, 4, &
    0, 0, 6, &
    1, 2, 4, &
    2, 4, 2, &
    6, 2, 0, &
    0, 0, 8, &
    6, 0, 4, &
    4, 6, 2, &
    2, 4, 8, &
   16, 0, 0 /), (/ 3, 17 /) )
  real ( kind = rk ) error
  integer factor
  integer factor_log
  real ( kind = rk ), external :: polyterm_value_3d
  integer j
  integer node_num
  real ( kind = rk ) result
  integer seed
  real ( kind = rk ) v1(3)
  real ( kind = rk ) v2(3)
  real ( kind = rk ) v3(3)

  seed = 123456789

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  SPHERE01_TRIANGLE_QUAD_ICOS2V approximates the'
  write ( *, '(a)' ) '  integral of a function over a spherical triangle'
  write ( *, '(a)' ) '  on the surface of the unit sphere using a vertex rule.'
  write ( *, '(a)' ) ' ' 
  write ( *, '(a)' ) '  We do not have an exact result, so we compare each'
  write ( *, '(a)' ) '  estimate to the final one.'
!
!  Choose three points at random to define a spherical triangle.
!
  call sphere01_sample ( 1, seed, v1 )
  call sphere01_sample ( 1, seed, v2 )
  call sphere01_sample ( 1, seed, v3 )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Vertices of random spherical triangle:'
  write ( *, '(a)' ) ' '
  call r8vec_transpose_print ( 3, v1, '  V1:' )
  call r8vec_transpose_print ( 3, v2, '  V2:' )
  call r8vec_transpose_print ( 3, v3, '  V3:' )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FACTOR   N   RESULT'

  do j = 1, 17

    e(1:3) = e_test(1:3,j)

    call polyterm_exponent ( 'SET', e )

    call polyterm_exponent ( 'PRINT', e )

    factor = 2**11
    call sphere01_triangle_quad_icos2v ( v1, v2, v3, factor, &
      polyterm_value_3d, node_num, best )

    factor = 1
    do factor_log = 0, 9

      call sphere01_triangle_quad_icos2v ( v1, v2, v3, factor, &
        polyterm_value_3d, node_num, result )

      error = abs ( result - best )

      write ( *, '(2x,i4,2x,i8,2x,g16.8,2x,g10.2)' ) &
        factor, node_num, result, error

      factor = factor * 2

    end do

  end do

  return
end
subroutine polyterm_exponent ( action, e )

!*****************************************************************************80
!
!! POLYTERM_EXPONENT gets or sets the exponents for the polynomial term.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = 3 ) ACTION.
!    'GET' asks the routine to return the current values in E.
!    'SET' asks the routine to set the current values to E.
!
!    Input/output, integer E(3), storage used to set or get values.
!
  implicit none

  character ( len = * )  action
  integer e(3)
  integer, save, dimension ( 3 ) :: e_save = (/ 0, 0, 0 /)
  character ( len = 80 ) text
  character ( len = 80 ) text2

  if ( action(1:1) == 'G' ) then

    e(1:3) = e_save(1:3)

  else if ( action(1:1) == 'P' ) then

    write ( *, * ) ' '

    if ( all ( e_save(1:3) == 0 ) ) then

      text = 'P(X,Y,Z) = 1'

    else

      text = 'P(X,Y,Z) = '

      if ( e_save(1) == 0 ) then

      else if ( e_save(1) == 1 ) then

        call s_cat ( text, ' X', text )

      else

        call s_cat ( text, ' X^', text )

        write ( text2, '(i2)' ) e_save(1)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if

      if ( e_save(2) == 0 ) then

      else if ( e_save(2) == 1 ) then

        call s_cat ( text, ' Y', text )

      else

        call s_cat ( text, ' Y^', text )

        write ( text2, '(i2)' ) e_save(2)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
       
      if ( e_save(3) == 0 ) then

      else if ( e_save(3) == 1 ) then

        call s_cat ( text, ' Z', text )

      else

        call s_cat ( text, ' Z^', text )

        write ( text2, '(i2)' ) e_save(3)
        text2 = adjustl ( text2 )
        call s_cat ( text, text2, text )

      end if
 
    end if

    write ( *, '(a)' ) trim ( text )
    
  else if ( action(1:1) == 'S' ) then

    e_save(1:3) = e(1:3)

  end if

  return
end
function polyterm_value_3d ( x )

!*****************************************************************************80
!
!! POLYTERM_VALUE_3D evaluates a single polynomial term in 3D.
!
!  Discussion:
!
!    The polynomial term has the form:
!
!      X(1)^E(1) * X(2)^E(2) * X(3)^E(3)
!
!    The exponents E(1:3) are set by calling POLYTERM_EXPONENT_SET.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    08 June 2002
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X(3), the point where the polynomial term 
!    is to be evaluated.
!
!    Output, real ( kind = rk ) POLYTERM_VALUE_3D, the value of the 
!    polynomial term.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer e(3)
  real ( kind = rk ) factor
  integer i
  real ( kind = rk ) polyterm_value_3d
  real ( kind = rk ) value
  real ( kind = rk ) x(3)

  call polyterm_exponent ( 'GET', e )

  value = 1.0D+00

  do i = 1, 3

    if ( e(i) == 0 ) then
      factor = 1.0D+00
    else if ( e(i) == 1 ) then
      factor = x(i)
    else if ( x(i) == 0.0D+00 ) then
      factor = 0.0D+00
    else
      factor = x(i)**e(i)
    end if

    value = value * factor

  end do

  polyterm_value_3d = value
  
  return
end
