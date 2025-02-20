program main

!*****************************************************************************80
!
!! MAIN is the main program for TEST_NLS_TEST.
!
!  Discussion:
!
!    TEST_NLS_TEST tests the TEST_NLS library.
!
!    This sample program demonstrates how the problems can
!    be used to test an algorithm or package that solves nonlinear
!    least squares problems.
!
!    In this case, we use the problems to test the performance of one of
!    the MINPACK routines, LMDER1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 December 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: abstol = 0.00001D+00
  real ( kind = rk ), allocatable, dimension ( : ) :: f
  external fcn
  real ( kind = rk ), allocatable, dimension (:,:) :: fjac
  real ( kind = rk ) fnrm
  real ( kind = rk ), allocatable, dimension ( : ) :: g
  integer info
  integer known
  integer m
  integer, parameter :: maxit = 30
  integer mprob
  integer n
  integer nprob
  integer nprob2
  character ( len = 80 ) title
  real ( kind = rk ), allocatable, dimension ( : ) :: x
  real ( kind = rk ) xnrm

  common /comprb/ nprob2

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NLS_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test the TEST_NLS library.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Show how the sample problems can be used.'
  write ( *, '(a)' ) '  In this example, we use the sample problems with the'
  write ( *, '(a)' ) '  MINPACK routine LMDER1.'
!
!  Get the number of problems.
!
  call p00_mprob ( mprob )
!
!  Solve each problem
!
  do nprob = 1, mprob

    nprob2 = nprob
!
!  Get the size of the problem.
!
    call p00_mn ( nprob, m, n )

    if ( m == 0 .or. n == 0 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'There are no more problems!'
      exit
    end if

    if ( m < 0 ) then
      m = abs ( m )
    end if

    if ( n < 0 ) then
      n = abs ( n )
    end if
!
!  Get the title of the problem.
!
    call p00_title ( nprob, title )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) trim ( title )
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) '  Number of equations M = ', m
    write ( *, '(a,i6)' ) '  Number of variables N = ', n
    write ( *, '(a)' ) ' '
!
!  Allocate space.
!
    allocate ( f(1:m) )
    allocate ( fjac(1:m,1:n) )
    allocate ( g(1:n) )
    allocate ( x(1:n) )
!
!  Get, and describe, starting point
!
    call p00_start ( nprob, n, x )

    if ( n <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting point X:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) x(1:n)
    else
      xnrm = sqrt ( sum ( x(1:n)**2 ) )
      write ( *, '(a)') ' '
      write ( *, '(a,g14.6)' ) '  ||X|| = ', xnrm
    end if

    call p00_f ( nprob, m, n, x, f )

    if ( m <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X):'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) f(1:m)
    else
      fnrm = sqrt ( sum ( f(1:m)**2 ) )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ||F(X)|| = ', fnrm
    end if

    if ( n <= 10 ) then
      call p00_g ( nprob, m, n, x, g )
      call r8vec_print ( n, g, '  The least squares gradient:' )
    end if
!
!  Call MINPACK routine LMDER1 to minimize F.
!
    call lmder1 ( fcn, m, n, x, f, fjac, m, abstol, info )
!
!  Report results.
!
    write ( *, '(a)' ) ' '
    write ( *, '(a,i6)' ) 'LMDER1 return flag INFO = ', info
    if ( info == 0 ) then
      write ( *, '(a)' ) 'Improper input parameters.'
    else if ( info == 1 ) then
      write ( *, '(a)' ) 'Relative error in the sum of squares is at most TOL.'
    else if ( info == 2 ) then
      write ( *, '(a)' ) 'Relative error in X is at most TOL.'
    else if ( info == 3 ) then
      write ( *, '(a)' ) 'Relative error in X and sum of squares at most TOL.'
    else if ( info == 4 ) then
      write ( *, '(a)' ) 'FVEC is orthogonal to the columns of the jacobian.'
    else if ( info == 5 ) then
      write ( *, '(a)' ) 'Too many calls to FCN with IFLAG = 1.'
    else if ( info == 6 ) then
      write ( *, '(a)' ) 'TOL is too small, can''t improve sum of squares.'
    else if ( info == 7 ) then
      write ( *, '(a)' ) 'TOL is too small, can''t improve X further.'
    end if
!
!  Describe final point
!
    if ( n <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Final point X:'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) x(1:n)
    else
      xnrm = sqrt ( sum ( x(1:n)**2 ) )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ||X|| = ',xnrm
    end if

    if ( m <= 10 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X):'
      write ( *, '(a)' ) ' '
      write ( *, '(5g14.6)' ) f(1:m)
    else
      fnrm = sqrt ( sum ( f(1:m)**2 ) )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  ||F(X)|| = ', fnrm
    end if

    if ( n <= 10 ) then
      call p00_g ( nprob, m, n, x, g )
      call r8vec_print ( n, g, '  The least squares gradient:' )
    end if

    if ( m <= 5 ) then
      call p00_j ( nprob, m, n, x, fjac )
      call r8mat_print ( m, n, fjac, '  Jacobian matrix:' )
    end if
!
!  Report true solution, if known.
!
    call p00_sol ( nprob, m, n, known, x )

    if ( known == 1 ) then

      if ( n <= 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Solution X:'
        write ( *, '(a)' ) ' '
        write ( *, '(5g14.6)' ) x(1:n)
      else
        xnrm = sqrt ( sum ( x(1:n)**2 ) )
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  ||X|| = ',xnrm
      end if

      call p00_f ( nprob, m, n, x, f )

      if ( m <= 10 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  F(X):'
        write ( *, '(a)' ) ' '
        write ( *, '(5g14.6)' ) f(1:m)
      else
        fnrm = sqrt ( sum ( f(1:m)**2 ) )
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  ||F(X)|| = ', fnrm
      end if

      if ( n <= 10 ) then
        call p00_g ( nprob, m, n, x, g )
        call r8vec_print ( n, g, '  The least squares gradient:' )
      end if

      if ( m <= 5 ) then
        call p00_j ( nprob, m, n, x, fjac )
        call r8mat_print ( m, n, fjac, '  Jacobian matrix:' )
      end if

    end if

    deallocate ( f )
    deallocate ( fjac )
    deallocate ( g )
    deallocate ( x )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_NLS_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine fcn ( m, n, x, f, fjac, ldfjac, iflag )

!*****************************************************************************80
!
!! FCN is a user-supplied routine to evaluate the function and jacobian.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    16 December 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer M, the number of equations.
!
!    Input, integer N, the number of variables.
!
!    Input, real ( kind = rk ) X(N), the point where evaluation is desired.
!
!    Output, real ( kind = rk ) F(M), the value of the function, if requested.
!
!    Output, real ( kind = rk ) FJAC(LDFJAC,N), the value of the M by N
!    jacobian matrix, if requested.
!
!    Input, integer LDFJAC, the leading dimension of FJAC.
!
!    Input, integer IFLAG, is 1 if F is to be evaluated; otherwise,
!    FJAC is to be evaluated.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ldfjac
  integer m
  integer n

  real ( kind = rk ) f(m)
  real ( kind = rk ) fjac(ldfjac,n)
  integer iflag
  integer nprob
  real ( kind = rk ) x(n)
!
!  This common block is just a trick to allow the main program to
!  specify a particular problem.
!
  common /comprb/ nprob

  if ( iflag == 1 ) then
    call p00_f ( nprob, m, n, x, f )
  else
    call p00_j ( nprob, m, n, x, fjac )
  end if

  return
end
