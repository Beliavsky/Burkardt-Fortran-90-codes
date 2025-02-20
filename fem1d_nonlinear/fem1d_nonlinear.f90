program main

!*****************************************************************************80
!
!! MAIN is the main program for FEM1D_NONLINEAR.
!
!  Discussion:
!
!    FEM1D_NONLINLEAR solves a nonlinear one dimensional boundary value problem.
!
!    The differential equation has the form:
!
!      -d/dx (p(x) du/dx) + q(x)*u +u*u' =  f(x)
!
!    The finite-element method uses piecewise linear basis functions.
!
!    Here U is an unknown scalar function of X defined on the
!    interval [XL,XR], and P, Q and F are given functions of X.
!
!    The values of U or U' at XL and XR are also specified.
!
!    Sample problem #1:
!
!    u(x)  = x,
!    p(x)  = 1,
!    q(x)  = 0,
!    f(x)  = x,
!    u(0)  = 0,
!    u'(1) = 1.
!    The code should solve this problem exactly.
!
!    Sample problem #2:
!
!    u(x)  = 2*(1-cos(0.5*pi*x))/pi,
!    p(x)  = 1,
!    q(x)  = 0,
!    f(x)  = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0)  = 0,
!    u'(1) = 1.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    real ( kind = rk ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    real ( kind = rk ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    real ( kind = rk ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    real ( kind = rk ) F(NU).
!    ASSEMBLE stores into F the right hand side of the linear
!    equations.
!    SOLVE replaces those values of F by the solution of the
!    linear equations.
!
!    real ( kind = rk ) FOLD(NU).
!    FOLD contains the value of F from the previous iteration,
!    and is used in ASSEMBLE to add correction terms to the
!    matrix and right hand side.
!
!    real ( kind = rk ) H(N), the length of the subintervals.
!
!    integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    integer IMAX.
!    The number of Newton iterations to carry out.
!
!    integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out when compared to the exact solution.
!
!    integer NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    integer N.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    real ( kind = rk ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    real ( kind = rk ) XQUAD(N)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
!    real ( kind = rk ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: n = 10
  integer, parameter :: nl = 2

  real ( kind = rk ) adiag(n+1)
  real ( kind = rk ) aleft(n+1)
  real ( kind = rk ) arite(n+1)
  real ( kind = rk ) f(n+1)
  real ( kind = rk ) h(n)
  integer i
  integer ibc
  integer imax
  integer indx(0:n)
  integer node(nl,n)
  integer nprint
  integer nquad
  integer nu
  integer problem
  real ( kind = rk ) u(n+1)
  real ( kind = rk ) ul
  real ( kind = rk ) ur
  real ( kind = rk ) xl
  real ( kind = rk ) xn(0:n)
  real ( kind = rk ) xquad(n)
  real ( kind = rk ) xr

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_NONLINEAR'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Solve a nonlinear boundary value problem:'
  write ( *, '(a)' ) '    -d/dx (p(x) du/dx) + q(x)*u + u*u'' =  f(x)'
  write ( *, '(a)' ) '  on an interval [xl,xr], with the values of'
  write ( *, '(a)' ) '  u or u'' specified at xl and xr.'
!
!  Do two problems
!
  do problem = 1, 2
!
!  Initialize variables that define the problem.
!
    call init ( ibc, imax, nprint, nquad, problem, ul, ur, xl, xr )
!
!  Compute the quantities that describe the geometry of the problem.
!
    call geometry ( h, ibc, indx, nl, node, n, nu, xl, xn, xquad, xr )

    u(1:nu) = 0.0
!
!  Begin the Newton iteration.
!
    do i = 1, imax
!
!  Assemble the linear system.
!
      if ( i <= 3 ) then
        call assemble_picard ( adiag, aleft, arite, f, h, indx, n, nl, &
          node, nquad, nu, problem, u, ul, ur, xn, xquad )
      else
        call assemble_newton ( adiag, aleft, arite, f, h, indx, n, nl, &
          node, nquad, nu, problem, u, ul, ur, xn, xquad )
      end if
!
!  Print out the linear system, just once.
!
      if ( i == 1 ) then
        call prsys ( adiag, aleft, arite, f, nu )
      end if
!
!  Solve the linear system.
!
      call solve ( adiag, aleft, arite, f, nu, u )
!
!  Print the current solution.
!
      call output ( u, ibc, indx, n, nu, ul, ur, xn )

    end do
!
!  Compare the solution to the exact solution.
!
    call compare ( u, indx, n, nl, node, nprint, nu, problem, ul, &
      ur, xl, xn, xr )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'FEM1D_NONLINEAR:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine assemble_newton ( adiag, aleft, arite, f, h, indx, n, nl, &
  node, nquad, nu, problem, u, ul, ur, xn, xquad )

!*****************************************************************************80
!
!! ASSEMBLE_NEWTON assembles the Newton linear system.
!
!  Discussion:
!
!    The linear system being solved here is for the Newton correction
!    to an approximate solution of a nonlinear system.
!
!    Thus, we suppose that we have a nonlinear function F(X),
!    and an approximate solution X0.  If we can suppose there is an
!    exact solution X* that is "nearby", and in fact close enough
!    that Taylor's theorem gives us a useful estimate, then we
!    may write:
!
!      F(X*) = F(X0) + F'(X0) * ( X* - X0 ) + Order ( X* - X0 )^2
!
!    and by rearranging, we get the Newton system (which is only
!    approximately correct):
!
!      F'(X0) * ( X* - X0 ) = - F(X0)
!
!    We solve this system and add the solution to X0 to get a
!    new approximate solution that, we hope, is much improved.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    29 April 2007
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Output, real ( kind = rk ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Output, real ( kind = rk ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Output, real ( kind = rk ) F(NU), the right hand side of the linear
!    equations.
!
!    Input, real ( kind = rk ) H(N), the length of the subintervals.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer N.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Input, real ( kind = rk ) U(NU), the solution value
!    from the previous iteration,
!
!    Input, real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Input, real ( kind = rk ) XQUAD(N)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nu

  real ( kind = rk ) aij
  real ( kind = rk ) adiag(nu)
  real ( kind = rk ) aleft(nu)
  real ( kind = rk ) arite(nu)
  real ( kind = rk ) f(nu)
  real ( kind = rk ) ff
  real ( kind = rk ) h(n)
  real ( kind = rk ) he
  integer ie
  integer ig
  integer il
  integer indx(0:n)
  integer iu
  integer iul
  integer iur
  integer iq
  integer jg
  integer jl
  integer jr
  integer ju
  integer nl
  integer node(nl,n)
  integer nquad
  real ( kind = rk ) phii
  real ( kind = rk ) phiix
  real ( kind = rk ) phij
  real ( kind = rk ) phijx
  real ( kind = rk ) pp
  integer problem
  real ( kind = rk ) qq
  real ( kind = rk ) total
  real ( kind = rk ) u(nu)
  real ( kind = rk ) ul
  real ( kind = rk ) uold
  real ( kind = rk ) uoldx
  real ( kind = rk ) ur
  real ( kind = rk ) x
  real ( kind = rk ) xleft
  real ( kind = rk ) xn(0:n)
  real ( kind = rk ) xqe
  real ( kind = rk ) xquad(n)
  real ( kind = rk ) xrite

  f(1:nu) = 0.0D+00
  adiag(1:nu) = 0.0D+00
  aleft(1:nu) = 0.0D+00
  arite(1:nu) = 0.0D+00
!
!  For element IE...
!
  do ie = 1, n

    he = h(ie)
    xleft = xn(node(1,ie))
    xrite = xn(node(2,ie))
!
!  For quadrature point IQ...
!
    do iq = 1, nquad

      xqe = xquad(ie)
!
!  Compute value of U for previous solution.
!
      total = 0.0D+00

      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( iu <= 0 ) then

          if ( il == 1 ) then
            total = total + ul
          else
            total = total + ur
          end if

        else
          total = total + u(iu)
        end if

      end do

      uold = total / real ( nl, kind = rk )
!
!  Compute value of U' for previous solution.
!
      jl = node(1,ie)
      jr = node(2,ie)
      iul = indx(jl)
      iur = indx(jr)

      if ( iul <= 0 ) then
        uoldx = ( u(iur) - ul ) / he
      else if ( iur <= 0 ) then
        uoldx = ( ur - u(iul) ) / he
      else
        uoldx = ( u(iur) - u(iul) ) / he
      end if
!
!  For basis function IL...
!
      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( 0 < iu ) then

          call phi ( il, xqe, phii, phiix, xleft, xrite )

          f(iu) = f(iu) + he * phii * ( ff ( xqe, problem ) + uold * uoldx )
!
!  Handle boundary conditions that prescribe the value of U'.
!
          if ( ig == 0 ) then

            x = 0.0D+00
            f(iu) = f(iu) - pp ( problem ) * ul

          else if ( ig == n ) then

            x = 1.0D+00
            f(iu) = f(iu) + pp ( problem ) * ur

          end if
!
!  For basis function JL...
!
          do jl = 1, nl

            jg = node(jl,ie)
            ju = indx(jg)

            call phi ( jl, xqe, phij, phijx, xleft, xrite )

            aij = he * ( pp ( problem ) * phiix * phijx &
                  + qq ( problem ) * phii * phij &
                  + uold * phii * phijx &
                  + uoldx * phij * phii )

            if ( ju <= 0 ) then

              if ( jg == 0 ) then
                f(iu) = f(iu) - aij * ul
              else if ( jg == n ) then
                f(iu) = f(iu) - aij * ur
              end if

            else if ( iu == ju ) then
              adiag(iu) = adiag(iu) + aij
            else if ( ju < iu ) then
              aleft(iu) = aleft(iu) + aij
            else
              arite(iu) = arite(iu) + aij
            end if

          end do

        end if

      end do

    end do

  end do

  return
end
subroutine assemble_picard ( adiag, aleft, arite, f, h, indx, n, nl, &
  node, nquad, nu, problem, u, ul, ur, xn, xquad )

!*****************************************************************************80
!
!! ASSEMBLE_PICARD assembles the Picard linear system.
!
!  Discussion:
!
!    The equation we are trying to solve has the form:
!
!      -d/dx ( p(x) du/dx ) + q(x) * u + u * u' = f(x)
!
!    For the Picard iteration, we need to modify the nonlinear term u * u'
!    so that it is linear in the unknown u, and any other factors of u are
!    lagged.  One way to do this gives us the following equation:
!
!      -d/dx ( p(x) du/dx ) + q(x) * u + u * uold' = f(x)
!
!    where uold is the previous iterate.
!
!    Now we can formulate this system as a (linear) finite element problem
!
!      A * u = rhs
!
!    to be solved for the new approximate solution u.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) ADIAG(NU).
!    ADIAG(I) is the "diagonal" coefficient of the I-th
!    equation in the linear system.  That is, ADIAG(I) is
!    the coefficient of the I-th unknown in the I-th equation.
!
!    Output, real ( kind = rk ) ALEFT(NU).
!    ALEFT(I) is the "left hand" coefficient of the I-th
!    equation in the linear system.  That is, ALEFT(I) is the
!    coefficient of the (I-1)-th unknown in the I-th equation.
!    There is no value in ALEFT(1), since the first equation
!    does not refer to a "0-th" unknown.
!
!    Output, real ( kind = rk ) ARITE(NU).
!    ARITE(I) is the "right hand" coefficient of the I-th
!    equation in the linear system.  ARITE(I) is the coefficient
!    of the (I+1)-th unknown in the I-th equation.  There is
!    no value in ARITE(NU) because the NU-th equation does not
!    refer to an "NU+1"-th unknown.
!
!    Output, real ( kind = rk ) F(NU), the right hand side of the linear
!    equations.
!
!    Input, real ( kind = rk ) H(N), the length of the subintervals.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer N.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Input, real ( kind = rk ) U(NU), the solution value
!    from the previous iteration,
!
!    Input, real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Input, real ( kind = rk ) XQUAD(N)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nu

  real ( kind = rk ) aij
  real ( kind = rk ) adiag(nu)
  real ( kind = rk ) aleft(nu)
  real ( kind = rk ) arite(nu)
  real ( kind = rk ) f(nu)
  real ( kind = rk ) ff
  real ( kind = rk ) h(n)
  real ( kind = rk ) he
  integer ie
  integer ig
  integer il
  integer indx(0:n)
  integer iu
  integer iul
  integer iur
  integer iq
  integer jg
  integer jl
  integer jr
  integer ju
  integer nl
  integer node(nl,n)
  integer nquad
  real ( kind = rk ) phii
  real ( kind = rk ) phiix
  real ( kind = rk ) phij
  real ( kind = rk ) phijx
  real ( kind = rk ) pp
  integer problem
  real ( kind = rk ) qq
  real ( kind = rk ) total
  real ( kind = rk ) u(nu)
  real ( kind = rk ) ul
  real ( kind = rk ) uold
  real ( kind = rk ) uoldx
  real ( kind = rk ) ur
  real ( kind = rk ) x
  real ( kind = rk ) xleft
  real ( kind = rk ) xn(0:n)
  real ( kind = rk ) xqe
  real ( kind = rk ) xquad(n)
  real ( kind = rk ) xrite

  f(1:nu) = 0.0D+00
  adiag(1:nu) = 0.0D+00
  aleft(1:nu) = 0.0D+00
  arite(1:nu) = 0.0D+00
!
!  For element IE...
!
  do ie = 1, n

    he = h(ie)
    xleft = xn(node(1,ie))
    xrite = xn(node(2,ie))
!
!  For quadrature point IQ...
!
    do iq = 1, nquad

      xqe = xquad(ie)
!
!  Compute value of U for previous solution.
!
      total = 0.0D+00

      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( iu <= 0 ) then

          if ( il == 1 ) then
            total = total + ul
          else
            total = total + ur
          end if

        else
          total = total + u(iu)
        end if

      end do

      uold = total / real ( nl, kind = rk )
!
!  Compute value of U' for previous solution.
!
      jl = node(1,ie)
      jr = node(2,ie)
      iul = indx(jl)
      iur = indx(jr)

      if ( iul <= 0 ) then
        uoldx = ( u(iur) - ul ) / he
      else if ( iur <= 0 ) then
        uoldx = ( ur - u(iul) ) / he
      else
        uoldx = ( u(iur) - u(iul) ) / he
      end if
!
!  For basis function IL...
!
      do il = 1, nl

        ig = node(il,ie)
        iu = indx(ig)

        if ( 0 < iu ) then

          call phi ( il, xqe, phii, phiix, xleft, xrite )

          f(iu) = f(iu) + he * phii * ( ff ( xqe, problem ) )
!
!  Handle boundary conditions that prescribe the value of U'.
!
          if ( ig == 0 ) then

            x = 0.0
            f(iu) = f(iu) - pp ( problem ) * ul

          else if ( ig == n ) then

            x = 1.0D+00
            f(iu) = f(iu) + pp ( problem ) * ur

          end if
!
!  For basis function JL...
!
          do jl = 1, nl

            jg = node(jl,ie)
            ju = indx(jg)

            call phi ( jl, xqe, phij, phijx, xleft, xrite )

            aij = he * ( pp ( problem ) * phiix * phijx &
                  + qq ( problem ) * phii * phij &
                  + uold * phii * phijx )

            if ( ju <= 0 ) then

              if ( jg == 0 ) then
                f(iu) = f(iu) - aij * ul
              else if ( jg == n ) then
                f(iu) = f(iu) - aij * ur
              end if

            else if ( iu == ju ) then
              adiag(iu) = adiag(iu) + aij
            else if ( ju < iu ) then
              aleft(iu) = aleft(iu) + aij
            else
              arite(iu) = arite(iu) + aij
            end if

          end do

        end if

      end do

    end do

  end do

  return
end
subroutine compare ( u, indx, n, nl, node, nprint, nu, problem, ul, ur, &
  xl, xn, xr )

!*****************************************************************************80
!
!! COMPARE compares the computed and exact solutions.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) U(NU), the solution of the linear equations.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer N.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Input, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out when compared to the exact solution.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Input, real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = rk ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    Input, real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Input, real ( kind = rk ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n
  integer nl
  integer nu

  real ( kind = rk ) u(nu)
  integer i
  integer ig
  integer indx(0:n)
  integer iu
  integer j
  integer k
  integer node(nl,n)
  integer nprint
  real ( kind = rk ) phii
  real ( kind = rk ) phiix
  integer problem
  real ( kind = rk ) u_exact
  real ( kind = rk ) u_value
  real ( kind = rk ) ul
  real ( kind = rk ) ur
  real ( kind = rk ) ux
  real ( kind = rk ) x
  real ( kind = rk ) xl
  real ( kind = rk ) xleft
  real ( kind = rk ) xn(0:n)
  real ( kind = rk ) xr
  real ( kind = rk ) xrite

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Compare computed and exact solutions:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      X      Computed U      Exact U'
  write ( *, '(a)' ) ' '

  do i = 1, nprint

    x = ( real ( nprint - i,     kind = rk ) * xl   &
        + real (          i - 1, kind = rk ) * xr ) &
        / real ( nprint     - 1, kind = rk )

    ux = u_exact ( x, problem )

    do j = 1, n

      xleft = xn(j-1)
      xrite = xn(j)
!
!  Search for the interval that X lies in.
!
      if ( xleft <= x .and. x <= xrite ) then

        u_value = 0.0D+00

        do k = 1, nl

          ig = node(k,j)
          iu = indx(ig)
          call phi ( k, x, phii, phiix, xleft, xrite )

          if ( iu <= 0 ) then
            if ( j == 1 .and. k == 1 ) then
              u_value = u_value + ul * phii
            else if ( j == n .and. k == nl ) then
              u_value = u_value + ur * phii
            end if
          else
            u_value = u_value + u(iu) * phii
          end if

        end do

        exit

      end if

    end do

    write(*,'(3g14.6)') x, u_value, ux

  end do

  return
end
function ff ( x, problem )

!*****************************************************************************80
!
!! FF returns the right hand side of the differential equation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Output, real ( kind = rk ) FF, the value of F(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ff
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  integer problem
  real ( kind = rk ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    ff = x
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    ff = - 0.5D+00 * pi * cos ( 0.5D+00 * pi * x ) &
      + 2.0D+00 * sin ( 0.5D+00 * pi * x ) &
      * ( 1.0D+00 - cos ( 0.5D+00 * pi * x ) ) / pi

  end if

  return
end
subroutine geometry ( h, ibc, indx, nl, node, n, nu, xl, xn, xquad, xr )

!*****************************************************************************80
!
!! GEOMETRY sets up the geometry for the interval [XL,XR].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = rk ) H(N), the length of the subintervals.
!
!    Input, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Output, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    Input, integer NL.
!    The number of basis functions used in a single
!    subinterval.  (NL-1) is the degree of the polynomials
!    used.  For this code, NL is fixed at 2, meaning that
!    piecewise linear functions are used as the basis.
!
!    Output, integer NODE(NL,N).
!    For each subinterval I:
!    NODE(1,I) is the number of the left node, and
!    NODE(2,I) is the number of the right node.
!
!    Input, integer N.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Output, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = rk ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    Output, real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
!    Output, real ( kind = rk ) XQUAD(N)
!    XQUAD(I) is the location of the single quadrature point
!    in interval I.
!
!    Input, real ( kind = rk ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nl
  integer n

  real ( kind = rk ) h(n)
  integer i
  integer ibc
  integer indx(0:n)
  integer node(nl,n)
  integer nu
  real ( kind = rk ) xl
  real ( kind = rk ) xn(0:n)
  real ( kind = rk ) xquad(n)
  real ( kind = rk ) xr
!
!  Set the value of XN, the locations of the nodes.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node      Location'
  write ( *, '(a)' ) ' '
  do i = 0, n
    xn(i) =  ( real ( n - i, kind = rk ) * xl   &
             + real (     i, kind = rk ) * xr ) &
             / real ( n,     kind = rk )
    write(*,'(i6,g14.6)') i, xn(i)
  end do
!
!  Set the lengths of each subinterval.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint    Length'
  write ( *, '(a)' ) ' '
  do i = 1, n
    h(i) = xn(i) - xn(i-1)
    write(*,'(i6,g14.6)') i, h(i)
  end do
!
!  Set the quadrature points, each of which is the midpoint of its subinterval.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint    Quadrature point'
  write ( *, '(a)' ) ' '
  do i = 1, n
    xquad(i) = 0.5D+00 * ( xn(i-1) + xn(i) )
    write(*,'(i6,g14.6)') i, xquad(i)
  end do
!
!  Set the value of NODE, which records, for each interval,
!  the node numbers at the left and right.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Subint  Left Node  Right Node'
  write ( *, '(a)' ) ' '
  do i = 1, n
    node(1,i) = i - 1
    node(2,i) = i
    write(*,'(3i6)') i, node(1,i), node(2,i)
  end do
!
!  Starting with node 0, see if an unknown is associated with
!  the node.  If so, give it an index.
!
  nu = 0
!
!  Handle first node.
!
  i = 0
  if ( ibc == 1 .or. ibc == 3 ) then
    indx(i) = -1
  else
    nu = nu + 1
    indx(i) = nu
  end if
!
!  Handle nodes 1 through nsub-1
!
  do i = 1, n - 1
    nu = nu + 1
    indx(i) = nu
  end do
!
!  Handle the last node.
!
  i = n
  if ( ibc == 2 .or. ibc == 3 ) then
    indx(i) = -1
  else
    nu = nu + 1
    indx(i) = nu
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Node  Unknown'
  write ( *, '(a)' ) ' '
  do i = 0, n
    write(*,'(2i6)') i, indx(i)
  end do

  return
end
subroutine init ( ibc, imax, nprint, nquad, problem, ul, ur, xl, xr )

!*****************************************************************************80
!
!! INIT initializes variables that define the problem.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Output, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Output, integer IMAX.
!    The number of Newton iterations to carry out.
!
!    Output, integer NPRINT.
!    The number of points at which the computed solution
!    should be printed out when compared to the exact solution.
!
!    Output, integer NQUAD.
!    The number of quadrature points used in a subinterval.
!    This code uses NQUAD = 1.
!
!    Output, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Output, real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Output, real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Output, real ( kind = rk ) XL.
!    XL is the left endpoint of the interval over which the
!    differential equation is being solved.
!
!    Output, real ( kind = rk ) XR.
!    XR is the right endpoint of the interval over which the
!    differential equation is being solved.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ibc
  integer imax
  integer nprint
  integer nquad
  integer problem
  real ( kind = rk ) ul
  real ( kind = rk ) ur
  real ( kind = rk ) xl
  real ( kind = rk ) xr

  ibc = 1
  imax = 10
  nprint = 9
  nquad = 1
  ul = 0.0D+00
  ur = 1.0D+00
  xl = 0.0D+00
  xr = 1.0D+00
!
!  Print out the values that have been set.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The equation is to be solved for'
  write ( *, '(a,g14.6)' ) 'X greater than XL  =  ', xl
  write ( *, '(a,g14.6)' ) ' and less than XR = ', xr
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'The boundary conditions are:'
  write ( *, '(a)' ) ' '

  if ( ibc == 1 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XL, U=', ul
  else
    write ( *, '(a,g14.6)' ) '  At X = XL, U''=', ul
  end if

  if ( ibc == 2 .or. ibc == 3 ) then
    write ( *, '(a,g14.6)' ) '  At X = XR, U=', ur
  else
    write ( *, '(a,g14.6)' ) '  At X = XR, U''=', ur
  end if

  if ( problem == 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is test problem #1:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  P(X) = 1, Q(X) = 0, F(X) = X.'
    write ( *, '(a)' ) '  Boundary conditions: U(0) = 0, U''(1) = 1.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The exact solution is U(X) = X'
  else if ( problem == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  This is test problem #2:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  P(X) = 1, Q(X) = 0, '
    write ( *, '(a)' ) '  F(X) = -0.5*pi*cos(0.5*pi*X)'
    write ( *, '(a)' ) '        + 2*sin(0.5*pi*X)*(1-cos(0.5*pi*X)/pi.'
    write ( *, '(a)' ) '  Boundary conditions: U(0) = 0, U''(1) = 1.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The exact solution is U(X) = 2*(1-cos(pi*x/2))/pi'
  end if

  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) 'Number of quadrature points per element is ', nquad
  write ( *, '(a,i8)' ) 'Number of iterations is ', imax

  return
end
subroutine output ( u, ibc, indx, nsub, nu, ul, ur, xn )

!*****************************************************************************80
!
!! OUTPUT prints out the computed solution at the nodes.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) U(NU), the solution of the linear equations.
!
!    Input, integer IBC.
!    IBC declares what the boundary conditions are.
!    1, at the left endpoint, U has the value UL,
!       at the right endpoint, U' has the value UR.
!    2, at the left endpoint, U' has the value UL,
!       at the right endpoint, U has the value UR.
!    3, at the left endpoint, U has the value UL,
!       and at the right endpoint, U has the value UR.
!    4, at the left endpoint, U' has the value UL,
!       at the right endpoint U' has the value UR.
!
!    Input, integer INDX(0:N).
!    For a node I, INDX(I) is the index of the unknown
!    associated with node I.
!    If INDX(I) is equal to -1, then no unknown is associated
!    with the node, because a boundary condition fixing the
!    value of U has been applied at the node instead.
!    Unknowns are numbered beginning with 1.
!    If IBC is 2 or 4, then there is an unknown value of U
!    at node 0, which will be unknown number 1.  Otherwise,
!    unknown number 1 will be associated with node 1.
!    If IBC is 1 or 4, then there is an unknown value of U
!    at node N, which will be unknown N or N+1,
!    depending on whether there was an unknown at node 0.
!
!    integer NSUB.
!    The number of subintervals into which the interval
!    [XL,XR] is broken.
!
!    Input, integer NU.
!    NU is the number of unknowns in the linear system.
!    Depending on the value of IBC, there will be N-1,
!    N, or N+1 unknown values, which are the coefficients
!    of basis functions.
!
!    Input, real ( kind = rk ) UL.
!    If IBC is 1 or 3, UL is the value that U is required
!    to have at X = XL.
!    If IBC is 2 or 4, UL is the value that U' is required
!    to have at X = XL.
!
!    Input, real ( kind = rk ) UR.
!    If IBC is 2 or 3, UR is the value that U is required
!    to have at X = XR.
!    If IBC is 1 or 4, UR is the value that U' is required
!    to have at X = XR.
!
!    Input, real ( kind = rk ) XN(0:N).
!    XN(I) is the location of the I-th node.  XN(0) is XL,
!    and XN(N) is XR.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nsub
  integer nu

  integer i
  integer ibc
  integer indx(0:nsub)
  real ( kind = rk ) u(nu)
  real ( kind = rk ) u_value
  real ( kind = rk ) ul
  real ( kind = rk ) ur
  real ( kind = rk ) xn(0:nsub)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Computed solution:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Node    X(I)        U(X(I))'
  write ( *, '(a)' ) ' '

  do i = 0, nsub

    if ( i == 0 ) then

      if ( ibc == 1 .or. ibc == 3 ) then
        u_value = ul
      else
        u_value = u(indx(i))
      end if

    else if ( i == nsub ) then

      if ( ibc == 2 .or. ibc == 3 ) then
        u_value = ur
      else
        u_value = u(indx(i))
      end if

    else

      u_value = u(indx(i))

    end if

    write(*,'(i4,2g14.6)')i, xn(i), u_value

  end do

  return
end
subroutine phi ( il, x, phii, phiix, xleft, xrite )

!*****************************************************************************80
!
!! PHI evaluates a linear basis function and its derivative.
!
!  Discussion:
!
!    In any interval, there are just two basis functions.  The first
!    basis function is a line which is 1 at the left endpoint
!    and 0 at the right.  The second basis function is 0 at
!    the left endpoint and 1 at the right.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer IL, the index of the basis function.
!    1, the function which is 1 at XLEFT and 0 at XRITE.
!    2, the function which is 0 at XLEFT and 1 at XRITE.
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Output, real ( kind = rk ) PHII, PHIIX, the value of the
!    basis function and its derivative at X.
!
!    Input, real ( kind = rk ) XLEFT, XRITE, the left and right
!    endpoints of the interval.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer il
  real ( kind = rk ) phii
  real ( kind = rk ) phiix
  real ( kind = rk ) x
  real ( kind = rk ) xleft
  real ( kind = rk ) xrite

  if ( xleft <= x .and. x <= xrite ) then

    if ( il == 1 ) then
      phii = ( xrite - x ) / ( xrite - xleft )
      phiix = -1.0D+00 / ( xrite - xleft )
    else
      phii = ( x - xleft ) / ( xrite - xleft )
      phiix = 1.0D+00 / ( xrite - xleft )
    end if
!
!  If X is outside of the interval, then the basis function
!  is always zero.
!
  else

    phii = 0.0D+00
    phiix = 0.0D+00

  end if

  return
end
function pp ( problem )

!*****************************************************************************80
!
!! PP returns the value of the coefficient function P(X).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, 
!      u = x, 
!      p=1, 
!      q=0,
!      f=x,
!      u(0)=0,
!      u'(1)=1.
!
!    * 2, 
!      u = 2*(1-cos(0.5*pi*x))/pi, 
!      p=1, 
!      q=0,
!      f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!      u(0) = 0,
!      u'(1)=1.
!
!    Output, real ( kind = rk ) PP, the value of P(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) pp
  integer problem
!
!  Test problem 1
!
  if ( problem == 1 ) then

    pp = 1.0D+00
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    pp = 1.0D+00

  end if

  return
end
subroutine prsys ( adiag, aleft, arite, f, nu )

!*****************************************************************************80
!
!! PRSYS prints out the tridiagonal linear system to be solved.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) ADIAG(NU), ALEFT(NU), ARITE(NU),
!    the diagonal, left and right entries of the equations.
!
!    Input, real ( kind = rk ) F(NU), the right hand side of the linear
!    system to be solved.
!
!    Input, integer NU.
!    NU is the number of equations to be solved.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nu

  real ( kind = rk ) adiag(nu)
  real ( kind = rk ) aleft(nu)
  real ( kind = rk ) arite(nu)
  real ( kind = rk ) f(nu)
  integer i

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Printout of tridiagonal linear system:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Equation  ALEFT  ADIAG  ARITE  RHS'
  write ( *, '(a)' ) ' '
  do i = 1, nu
    if ( i == 1 ) then
      write(*,'(i3,14x,3g14.6)') i, adiag(i), arite(i), f(i)
    else if ( i < nu ) then
      write(*,'(i3,4g14.6)') i, aleft(i), adiag(i), arite(i), f(i)
    else
      write(*,'(i3,2g14.6,14x,g14.6)') i, aleft(i), adiag(i), f(i)
    end if
  end do

  return
end
function qq ( problem )

!*****************************************************************************80
!
!! QQ returns the value of the coefficient function Q(X).
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Output, real ( kind = rk ) QQ, the value of Q(X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer problem
  real ( kind = rk ) qq
!
!  Test problem 1
!
  if ( problem == 1 ) then

    qq = 0.0D+00
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    qq = 0.0D+00

  end if

  return
end
subroutine solve ( adiag, aleft, arite, f, nu, u )

!*****************************************************************************80
!
!! SOLVE solves a tridiagonal matrix system of the form A*x = b.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) ADIAG(NU), ALEFT(NU), ARITE(NU).
!    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
!    left and right entries of the equations.
!    On output, ADIAG and ARITE have been changed in order
!    to compute the solution.
!    Note that for the first equation, there is no ALEFT
!    coefficient, and for the last, there is no ARITE.
!    So there is no need to store a value in ALEFT(1), nor
!    in ARITE(NU).
!
!    Input, real ( kind = rk ) F(NU), the right hand side of the linear
!    system to be solved.
!
!    Input, integer NU.
!    NU is the number of equations to be solved.
!
!    Output, real ( kind = rk ) U(NU), the solution of the linear system.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nu

  real ( kind = rk ) adiag(nu)
  real ( kind = rk ) aleft(nu)
  real ( kind = rk ) arite(nu-1)
  real ( kind = rk ) f(nu)
  integer i
  real ( kind = rk ) u(nu)

  arite(1) = arite(1) / adiag(1)
  do i = 2, nu-1
    adiag(i) = adiag(i) - aleft(i) * arite(i-1)
    arite(i) = arite(i) / adiag(i)
  end do
  adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)

  u(1:nu) = f(1:nu)

  u(1) = u(1) / adiag(1)
  do i = 2, nu
    u(i) = ( u(i) - aleft(i) * u(i-1) ) / adiag(i)
  end do

  do i = nu-1, 1, -1
    u(i) = u(i) - arite(i) * u(i+1)
  end do

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

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
function u_exact ( x, problem )

!*****************************************************************************80
!
!! U_EXACT returns the value of the exact solution at a point X.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    01 November 2006
!
!  Author:
!
!    FORTRAN90 version by John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = rk ) X, the evaluation point.
!
!    Input, integer PROBLEM, indicates which problem to be solved.
!    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
!    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
!    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
!    u(0) = 0, u'(1)=1.
!
!    Output, real ( kind = rk ) U_EXACT, the value of the exact solution at X.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  integer problem
  real ( kind = rk ) u_exact
  real ( kind = rk ) x
!
!  Test problem 1
!
  if ( problem == 1 ) then

    u_exact = x
!
!  Test problem 2
!
  else if ( problem == 2 ) then

    u_exact = 2.0D+00 * ( 1.0D+00 - cos ( 0.5D+00 * pi * x ) ) / pi

  end if

  return
end
