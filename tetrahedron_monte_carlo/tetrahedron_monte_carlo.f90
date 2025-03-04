function r8mat_det_4d ( a )

!*****************************************************************************80
!
!! R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
!
!  Discussion:
!
!    An R8MAT is a two dimensional matrix of double precision real values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    01 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) A(4,4), the matrix whose determinant is desired.
!
!    Output, real ( kind = rk ) R8MAT_DET_4D, the determinant of the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d

  r8mat_det_4d = &
         a(1,1) * ( &
             a(2,2) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           + a(2,4) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) ) &
       - a(1,2) * ( &
             a(2,1) * ( a(3,3) * a(4,4) - a(3,4) * a(4,3) ) &
           - a(2,3) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) ) &
       + a(1,3) * ( &
             a(2,1) * ( a(3,2) * a(4,4) - a(3,4) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,4) - a(3,4) * a(4,1) ) &
           + a(2,4) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) ) &
       - a(1,4) * ( &
             a(2,1) * ( a(3,2) * a(4,3) - a(3,3) * a(4,2) ) &
           - a(2,2) * ( a(3,1) * a(4,3) - a(3,3) * a(4,1) ) &
           + a(2,3) * ( a(3,1) * a(4,2) - a(3,2) * a(4,1) ) )

  return
end
subroutine r8mat_transpose_print ( m, n, a, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = * ) title

  call r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title )

  return
end
subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

!*****************************************************************************80
!
!! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
!
!  Discussion:
!
!    An R8MAT is an array of R8 values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    14 June 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, N, the number of rows and columns.
!
!    Input, real ( kind = rk ) A(M,N), an M by N matrix to be printed.
!
!    Input, integer ILO, JLO, the first row and column to print.
!
!    Input, integer IHI, JHI, the last row and column to print.
!
!    Input, character ( len = * ) TITLE, an optional title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: incx = 5
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  character ( len = 14 ) ctemp(incx)
  integer i
  integer i2
  integer i2hi
  integer i2lo
  integer ihi
  integer ilo
  integer inc
  integer j
  integer j2hi
  integer j2lo
  integer jhi
  integer jlo
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )

  do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

    i2hi = i2lo + incx - 1
    i2hi = min ( i2hi, m )
    i2hi = min ( i2hi, ihi )

    inc = i2hi + 1 - i2lo

    write ( *, '(a)' ) ' '

    do i = i2lo, i2hi
      i2 = i + 1 - i2lo
      write ( ctemp(i2), '(i8,6x)' ) i
    end do

    write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
    write ( *, '(a)' ) '  Col'
    write ( *, '(a)' ) ' '

    j2lo = max ( jlo, 1 )
    j2hi = min ( jhi, n )

    do j = j2lo, j2hi

      do i2 = 1, inc
        i = i2lo - 1 + i2
        write ( ctemp(i2), '(g14.6)' ) a(i,j)
      end do

      write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

    end do

  end do

  return
end
subroutine r8vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    For now, the input quantity SEED is an integer variable.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input/output, integer SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = rk ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  integer, parameter :: i4_huge = 2147483647
  integer k
  integer seed
  real ( kind = rk ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = rk ) * 4.656612875D-10

  end do

  return
end
subroutine reference_to_physical_tet4 ( t, n, ref, phy )

!*****************************************************************************80
!
!! REFERENCE_TO_PHYSICAL_TET4 maps TET4 reference points to physical points.
!
!  Discussion:
!
!    Given the vertices of an order 4 physical tetrahedron and a point 
!    (R,S,T) in the reference tetrahedron, the routine computes the value 
!    of the corresponding point (X,Y,Z) in the physical tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,4), the coordinates of the vertices.  
!    The vertices are assumed to be the images of (1,0,0), (0,1,0),
!    (0,0,1) and (0,0,0) respectively.
!
!    Input, integer N, the number of points to transform.
!
!    Input, real ( kind = rk ) REF(3,N), points in the reference tetrahedron.
!
!    Output, real ( kind = rk ) PHY(3,N), corresponding points in the
!    physical tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  integer i
  real ( kind = rk ) phy(3,n)
  real ( kind = rk ) ref(3,n)
  real ( kind = rk ) t(3,4)

  do i = 1, 3
    phy(i,1:n) =                                                    &
        t(i,1) *             ref(1,1:n)                             &
      + t(i,2) *                          ref(2,1:n)                &
      + t(i,3) *                                       ref(3,1:n)   &
      + t(i,4) * ( 1.0D+00 - ref(1,1:n) - ref(2,1:n) - ref(3,1:n) ) 
  end do

  return
end
subroutine tetrahedron_integrand_01 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_01 evaluates 1 integrand function.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input, real ( kind = rk ) P(3,P_NUM), the evaluation points.
!
!    Input, integer F_NUM, the number of integrands.
!
!    Output, real ( kind = rk ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  real ( kind = rk ) p(3,p_num)

  fp(1,1:p_num) = 1.0D+00

  return
end
subroutine tetrahedron_integrand_02 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_02 evaluates 3 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input, real ( kind = rk ) P(3,P_NUM), the evaluation points.
!
!    Input, integer F_NUM, the number of integrands.
!
!    Output, real ( kind = rk ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  real ( kind = rk ) p(3,p_num)

  fp(1,1:p_num) = p(1,1:p_num)
  fp(2,1:p_num) = p(2,1:p_num)
  fp(3,1:p_num) = p(3,1:p_num)

  return
end
subroutine tetrahedron_integrand_03 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_03 evaluates 6 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input, real ( kind = rk ) P(3,P_NUM), the evaluation points.
!
!    Input, integer F_NUM, the number of integrands.
!
!    Output, real ( kind = rk ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  real ( kind = rk ) p(3,p_num)

  fp(1,1:p_num) = p(1,1:p_num) * p(1,1:p_num)
  fp(2,1:p_num) = p(1,1:p_num) * p(2,1:p_num)
  fp(3,1:p_num) = p(1,1:p_num) * p(3,1:p_num)
  fp(4,1:p_num) = p(2,1:p_num) * p(2,1:p_num)
  fp(5,1:p_num) = p(2,1:p_num) * p(3,1:p_num)
  fp(6,1:p_num) = p(3,1:p_num) * p(3,1:p_num)

  return
end
subroutine tetrahedron_integrand_04 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_04 evaluates 10 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input, real ( kind = rk ) P(3,P_NUM), the evaluation points.
!
!    Input, integer F_NUM, the number of integrands.
!
!    Output, real ( kind = rk ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  real ( kind = rk ) p(3,p_num)

  fp( 1,1:p_num) = p(1,1:p_num)**3
  fp( 2,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)
  fp( 3,1:p_num) = p(1,1:p_num)**2                   * p(3,1:p_num)
  fp( 4,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)**2
  fp( 5,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)    * p(3,1:p_num)
  fp( 6,1:p_num) = p(1,1:p_num)                      * p(3,1:p_num)**2
  fp( 7,1:p_num) =                   p(2,1:p_num)**3
  fp( 8,1:p_num) =                   p(2,1:p_num)**2 * p(3,1:p_num)
  fp( 9,1:p_num) =                   p(2,1:p_num)    * p(3,1:p_num)**2
  fp(10,1:p_num) =                                     p(3,1:p_num)**3

  return
end
subroutine tetrahedron_integrand_05 ( p_num, p, f_num, fp )

!*****************************************************************************80
!
!! TETRAHEDRON_INTEGRAND_05 evaluates 15 integrand functions.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input, real ( kind = rk ) P(3,P_NUM), the evaluation points.
!
!    Input, integer F_NUM, the number of integrands.
!
!    Output, real ( kind = rk ) FP(F_NUM,P_NUM), the integrand values.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  real ( kind = rk ) p(3,p_num)

  fp( 1,1:p_num) = p(1,1:p_num)**4
  fp( 2,1:p_num) = p(1,1:p_num)**3 * p(2,1:p_num)
  fp( 3,1:p_num) = p(1,1:p_num)**3                   * p(3,1:p_num)
  fp( 4,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)**2
  fp( 5,1:p_num) = p(1,1:p_num)**2 * p(2,1:p_num)    * p(3,1:p_num)
  fp( 6,1:p_num) = p(1,1:p_num)**2                   * p(3,1:p_num)**2
  fp( 7,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)**3
  fp( 8,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)**2 * p(3,1:p_num)
  fp( 9,1:p_num) = p(1,1:p_num)    * p(2,1:p_num)    * p(3,1:p_num)**2
  fp(10,1:p_num) = p(1,1:p_num)                      * p(3,1:p_num)**3
  fp(11,1:p_num) =                   p(2,1:p_num)**4
  fp(12,1:p_num) =                   p(2,1:p_num)**3 * p(3,1:p_num)
  fp(13,1:p_num) =                   p(2,1:p_num)**2 * p(3,1:p_num)**2
  fp(14,1:p_num) =                   p(2,1:p_num)    * p(3,1:p_num)**3
  fp(15,1:p_num) =                                     p(3,1:p_num)**4
  
  return
end
subroutine tetrahedron_monte_carlo ( t, p_num, f_num, tetrahedron_unit_sample, &
  tetrahedron_integrand, seed, result )

!*****************************************************************************80
!
!! TETRAHEDRON_MONTE_CARLO applies the Monte Carlo rule to integrate a function.
!
!  Discussion:
!
!    The function f(x,y,z) is to be integrated over a tetrahedron.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) T(3,4), the vertices.
!
!    Input, integer P_NUM, the number of sample points.
!
!    Input, integer F_NUM, the number of functions to integrate.
!
!    Input, external TETRAHEDRON_UNIT_SAMPLE, the sampling routine.
!
!    Input, external TETRAHEDRON_INTEGRAND, the integrand routine.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) RESULT(F_NUM), the approximate integrals.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer f_num
  integer p_num

  real ( kind = rk ) fp(f_num,p_num)
  integer i
  real ( kind = rk ) p(3,p_num)
  real ( kind = rk ) p2(3,p_num)
  real ( kind = rk ) result(f_num)
  integer seed
  real ( kind = rk ) t(3,4)
  external tetrahedron_sample
  external tetrahedron_integrand
  real ( kind = rk ) volume

  call tetrahedron_volume ( t, volume )

  call tetrahedron_unit_sample ( p_num, seed, p )

  call reference_to_physical_tet4 ( t, p_num, p, p2 )

  call tetrahedron_integrand ( p_num, p2, f_num, fp )

  do i = 1, f_num
    result(i) = volume * sum ( fp(i,1:p_num) ) / real ( p_num, kind = rk )
  end do

  return
end
subroutine tetrahedron_unit_sample_01 ( p_num, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_UNIT_SAMPLE_01 selects points from the unit tetrahedron.
!
!  Discussion:
!
!    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
!
!    Any point in the unit tetrahedron CAN be chosen by this algorithm.
!
!    However, the points that are chosen tend to be clustered near
!    the centroid.
!
!    This routine is supplied as an example of "bad" sampling.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) P(3,P_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer p_num

  real ( kind = rk ) e(4)
  real ( kind = rk ) e_sum
  integer j
  real ( kind = rk ) p(3,p_num)
  integer seed

  do j = 1, p_num

    call r8vec_uniform_01 ( 4, seed, e )

    e_sum = sum ( e(1:4) )

    e(1:4) = e(1:4) / e_sum
!
!  We may take the values E(1:3) as being the barycentric
!  coordinates of the point.
!
    p(1:3,j) = e(1:3)

  end do

  return
end
subroutine tetrahedron_unit_sample_02 ( p_num, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_UNIT_SAMPLE_02 selects points from the unit tetrahedron.
!
!  Discussion:
!
!    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
!
!    The sampling is uniform.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Claudio Rocchini, Paolo Cignoni,
!    Generating Random Points in a Tetrahedron,
!    Journal of Graphics Tools,
!    Volume 5, Number 5, 2000, pages 9-12.
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) P(3,P_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer p_num

  real ( kind = rk ) c(3)
  integer j
  integer seed
  real ( kind = rk ) t
  real ( kind = rk ) p(3,p_num)

  do j = 1, p_num

    call r8vec_uniform_01 ( 3, seed, c )

    if ( 1.0D+00 < c(1) + c(2) ) then
      c(1) = 1.0D+00 - c(1)
      c(2) = 1.0D+00 - c(2)
    end if

    if ( 1.0D+00 < c(2) + c(3) ) then
      t = c(3)
      c(3) = 1.0D+00 - c(1) - c(2)
      c(2) = 1.0D+00 - t
    else if ( 1.0D+00 < c(1) + c(2) + c(3) ) then
      t = c(3)
      c(3) = c(1) + c(2) + c(3) - 1.0D+00
      c(1) = 1.0D+00 - c(2) - t
    end if

    p(1:3,j) = c(1:3)

  end do

  return
end
subroutine tetrahedron_unit_sample_03 ( p_num, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_UNIT_SAMPLE_03 selects points from the unit tetrahedron.
!
!  Discussion:
!
!    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Greg Turk,
!    Generating Random Points in a Triangle,
!    in Graphics Gems,
!    edited by Andrew Glassner,
!    AP Professional, 1990, pages 24-28.
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) P(3,P_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer p_num

  real ( kind = rk ) a
  real ( kind = rk ) b
  real ( kind = rk ) c
  real ( kind = rk ) d
  real ( kind = rk ) e
  real ( kind = rk ) f
  real ( kind = rk ) g
  integer j
  real ( kind = rk ) p(3,p_num)
  real ( kind = rk ) r(3)
  integer seed

  do j = 1, p_num

    call r8vec_uniform_01 ( 3, seed, r )

    e = r(1)**(1.0D+00/3.0D+00)
    f = sqrt ( r(2) )
    g = r(3)

    a =   1.0D+00 - e
    b = ( 1.0D+00 - f )       * e
    c = ( 1.0D+00 - g ) * f   * e
    d =             g   * f   * e

    p(1,j) = a
    p(2,j) = b
    p(3,j) = c

  end do

  return
end
subroutine tetrahedron_unit_sample_04 ( p_num, seed, p )

!*****************************************************************************80
!
!! TETRAHEDRON_UNIT_SAMPLE_04 selects points from the unit tetrahedron.
!
!  Discussion:
!
!    The unit tetrahedron has vertices (1,0,0), (0,1,0), (0,0,1), (0,0,0).
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    16 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Reuven Rubinstein,
!    Monte Carlo Optimization, Simulation, and Sensitivity 
!    of Queueing Networks,
!    Krieger, 1992,
!    ISBN: 0894647644,
!    LC: QA298.R79.
!
!  Parameters:
!
!    Input, integer P_NUM, the number of points.
!
!    Input/output, integer SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = rk ) P(3,P_NUM), the points.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer p_num

  real ( kind = rk ) e(4)
  integer j
  real ( kind = rk ) p(3,p_num)
  integer seed
!
!  The construction begins by sampling DIM_NUM+1 points from the
!  exponential distribution with parameter 1.
!
  do j = 1, p_num

    call r8vec_uniform_01 ( 4, seed, e )

    e(1:4) = - log ( e(1:4) )

    p(1:3,j) = e(1:3) / sum ( e(1:4) )

  end do

  return
end
subroutine tetrahedron_volume ( tet_xyz, volume )

!*****************************************************************************80
!
!! TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    30 December 2004
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) TET_XYZ(3,4), the coordinates of the vertices.
!
!    Output, real ( kind = rk ) VOLUME, the volume of the tetrahedron.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: dim_num = 3

  real ( kind = rk ) a(4,4)
  real ( kind = rk ) r8mat_det_4d
  real ( kind = rk ) tet_xyz(dim_num,4)
  real ( kind = rk ) volume

  a(1:dim_num,1:4) = tet_xyz( 1:dim_num,1:4)
  a(4,1:4) = 1.0D+00

  volume = abs ( r8mat_det_4d ( a ) ) / 6.0D+00

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  character ( len = 8 ) ampm
  integer d
  integer h
  integer m
  integer mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer n
  integer s
  integer values(8)
  integer y

  call date_and_time ( values = values )

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

  write ( *, '(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
