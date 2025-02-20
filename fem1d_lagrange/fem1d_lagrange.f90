subroutine fem1d_lagrange_stiffness ( x_num, x, q_num, f, a, m, b )

!*****************************************************************************80
!
!! FEM1D_LAGRANGE_STIFFNESS evaluates the Lagrange polynomial stiffness matrix.
!
!  Discussion:
!
!    The finite element method is to be applied over a given interval that
!    has been meshed with X_NUM points X.
!
!    The finite element basis functions are to be the X_NUM Lagrange
!    basis polynomials L(i)(X), such that
!      L(i)(X(j)) = delta(i,j).
!
!    The following items are computed:
!    * A, the stiffness matrix, with A(I,J) = integral L'(i)(x) L'(j)(x)
!    * M, the mass matrix, with M(I,J) = integral L(i)(x) L(j)(x)
!    * B, the load matrix, with B(I) = integral L(i)(x) F(x)
!
!    The integrals are approximated by quadrature.
!
!    Boundary conditions are not handled by this routine.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer X_NUM, the number of nodes.
!
!    Input, real ( kind = rk ) X(X_NUM), the coordinates of the nodes.
!
!    Input, integer Q_NUM, the number of quadrature points to use.
!
!    Input, external, real ( kind = rk ) F(X), the right hand side function.
!
!    Output, real ( kind = rk ) A(X_NUM,X_NUM), the stiffness matrix.
!
!    Output, real ( kind = rk ) M(X_NUM,X_NUM), the mass matrix.
!
!    Output, real ( kind = rk ) B(X_NUM), the right hand side vector.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer x_num

  real ( kind = rk ) a(x_num,x_num)
  real ( kind = rk ) b(x_num)
  real ( kind = rk ), external :: f
  real ( kind = rk ), allocatable :: l(:,:)
  real ( kind = rk ) li
  real ( kind = rk ) lj
  real ( kind = rk ), allocatable :: lp(:,:)
  real ( kind = rk ) lpi
  real ( kind = rk ) lpj
  real ( kind = rk ) m(x_num,x_num)
  integer q_i
  integer q_num
  real ( kind = rk ), allocatable :: q_w(:)
  real ( kind = rk ), allocatable :: q_x(:)
  real ( kind = rk ) x(x_num)
  integer x_i
  integer x_j
!
!  Get the quadrature rule for [-1,+1].
!
  allocate ( q_x(1:q_num) )
  allocate ( q_w(1:q_num) )

  call legendre_set ( q_num, q_x, q_w )
!
!  Adjust the quadrature rule to the interval [ x(1), x(x_num) }
!
  q_x(1:q_num) =  ( ( 1.0D+00 - q_x(1:q_num) ) * x(1) &
                  + ( 1.0D+00 + q_x(1:q_num) ) * x(x_num) ) &
                  /   2.0D+00

  q_w(1:q_num) = q_w(1:q_num) * ( x(x_num) - x(1) ) / 2.0D+00
!
!  Evaluate all the Lagrange basis polynomials at all the quadrature points.
!
  allocate ( l(1:q_num,1:x_num) )

  call lagrange_value ( x_num, x, q_num, q_x, l )
!
!  Evaluate all the Lagrange basis polynomial derivatives at all 
!  the quadrature points.
!
  allocate ( lp(1:q_num,1:x_num) )

  call lagrange_derivative ( x_num, x, q_num, q_x, lp )
!
!  Assemble the matrix and right hand side.
!
  a(1:x_num,1:x_num) = 0.0D+00
  m(1:x_num,1:x_num) = 0.0D+00
  b(1:x_num) = 0.0D+00

  do x_i = 1, x_num
    do q_i = 1, q_num
      li = l(q_i,x_i)
      lpi = lp(q_i,x_i)
      do x_j = 1, x_num
        lj = l(q_i,x_j)
        lpj = lp(q_i,x_j)
        a(x_i,x_j) = a(x_i,x_j) + q_w(q_i) * lpi * lpj
        m(x_i,x_j) = m(x_i,x_j) + q_w(q_i) * li * lj
      end do
      b(x_i) = b(x_i) + q_w(q_i) * li * f ( q_x(q_i) )
    end do
  end do

  deallocate ( l )
  deallocate ( lp )
  deallocate ( q_w )
  deallocate ( q_x )

  return
end
subroutine lagrange_derivative ( nd, xd, ni, xi, lpi ) 

!*****************************************************************************80
!
!! LAGRANGE_DERIVATIVE evaluates the Lagrange basis derivative.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 November 2014
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the interpolation nodes.
!
!    Input, integer NI, the number of evaluation points.
!
!    Input, real ( kind = rk ) XI(NI), the evaluation points.
!
!    Output, real ( kind = rk ) LPI(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function derivative.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  integer i
  integer j
  integer j1
  integer j2
  real ( kind = rk ) lpi(ni,nd)
  real ( kind = rk ) p
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  
  lpi(1:ni,1:nd) = 0.0D+00

  do i = 1, ni
    do j = 1, nd

      do j1 = 1, nd

        if ( j1 /= j ) then
          p = 1.0D+00
          do j2 = 1, nd
            if ( j2 == j1 ) then
              p = p / ( xd(j) - xd(j2) )
            else if ( j2 /= j ) then
              p = p * ( xi(i) - xd(j2) ) / ( xd(j) - xd(j2) )
            end if
          end do
          lpi(i,j) = lpi(i,j) + p
        end if

      end do

    end do
  end do

  return
end
subroutine lagrange_value ( nd, xd, ni, xi, li ) 

!*****************************************************************************80
!
!! LAGRANGE_VALUE evaluates the Lagrange basis polynomials.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    09 October 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ND, the number of data points.
!
!    Input, real ( kind = rk ) XD(ND), the interpolation nodes.
!
!    Input, integer NI, the number of evaluation points.
!
!    Input, real ( kind = rk ) XI(NI), the evaluation points.
!
!    Output, real ( kind = rk ) LI(NI,ND), the value, at the I-th point XI, 
!    of the Jth basis function.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer nd
  integer ni

  integer i
  integer j
  real ( kind = rk ) li(ni,nd)
  real ( kind = rk ) xd(nd)
  real ( kind = rk ) xi(ni)
  
  do i = 1, ni
    do j = 1, nd
      li(i,j) = product ( ( xi(i) - xd(1:j-1)  ) / ( xd(j) - xd(1:j-1)  ) ) &
              * product ( ( xi(i) - xd(j+1:nd) ) / ( xd(j) - xd(j+1:nd) ) )
    end do
  end do

  return
end
subroutine legendre_set ( n, x, w )

!*****************************************************************************80
!
!! LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
!
!  Discussion:
!
!    The integral:
!
!      integral ( -1 <= x <= 1 ) f(x) dx
!
!    The quadrature rule:
!
!      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
!    The quadrature rule is exact for polynomials through degree 2*N-1.
!
!    The abscissas are the zeroes of the Legendre polynomial P(N)(X).
!
!    Mathematica can compute the abscissas and weights of a Gauss-Legendre
!    rule of order N for the interval [A,B] with P digits of precision
!    by the commands:
!
!    Needs["NumericalDifferentialEquationAnalysis`"]
!    GaussianQuadratureWeights [ n, a, b, p ]
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    20 April 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    Dover, 2006,
!    ISBN: 0486445798,
!    LC: QA311.K713.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Stephen Wolfram,
!    The Mathematica Book,
!    Fourth Edition,
!    Cambridge University Press, 1999,
!    ISBN: 0-521-64314-7,
!    LC: QA76.95.W65.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
!  Parameters:
!
!    Input, integer N, the order.
!    N must be between 1 and 16.
!
!    Output, real ( kind = rk ) X(N), the abscissas.
!
!    Output, real ( kind = rk ) W(N), the weights.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) w(n)
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = 0.000000000000000000000000000000D+00

    w(1) = 2.000000000000000000000000000000D+00

  else if ( n == 2 ) then

    x(1) = -0.577350269189625764509148780502D+00
    x(2) = 0.577350269189625764509148780502D+00

    w(1) = 1.000000000000000000000000000000D+00
    w(2) = 1.000000000000000000000000000000D+00

  else if ( n == 3 ) then

    x(1) = -0.774596669241483377035853079956D+00
    x(2) = 0.000000000000000000000000000000D+00
    x(3) = 0.774596669241483377035853079956D+00

    w(1) = 0.555555555555555555555555555556D+00
    w(2) = 0.888888888888888888888888888889D+00
    w(3) = 0.555555555555555555555555555556D+00

  else if ( n == 4 ) then

    x(1) = -0.861136311594052575223946488893D+00
    x(2) = -0.339981043584856264802665759103D+00
    x(3) = 0.339981043584856264802665759103D+00
    x(4) = 0.861136311594052575223946488893D+00

    w(1) = 0.347854845137453857373063949222D+00
    w(2) = 0.652145154862546142626936050778D+00
    w(3) = 0.652145154862546142626936050778D+00
    w(4) = 0.347854845137453857373063949222D+00

  else if ( n == 5 ) then

    x(1) = -0.906179845938663992797626878299D+00
    x(2) = -0.538469310105683091036314420700D+00
    x(3) = 0.000000000000000000000000000000D+00
    x(4) = 0.538469310105683091036314420700D+00
    x(5) = 0.906179845938663992797626878299D+00

    w(1) = 0.236926885056189087514264040720D+00
    w(2) = 0.478628670499366468041291514836D+00
    w(3) = 0.568888888888888888888888888889D+00
    w(4) = 0.478628670499366468041291514836D+00
    w(5) = 0.236926885056189087514264040720D+00

  else if ( n == 6 ) then

    x(1) = -0.932469514203152027812301554494D+00
    x(2) = -0.661209386466264513661399595020D+00
    x(3) = -0.238619186083196908630501721681D+00
    x(4) = 0.238619186083196908630501721681D+00
    x(5) = 0.661209386466264513661399595020D+00
    x(6) = 0.932469514203152027812301554494D+00

    w(1) = 0.171324492379170345040296142173D+00
    w(2) = 0.360761573048138607569833513838D+00
    w(3) = 0.467913934572691047389870343990D+00
    w(4) = 0.467913934572691047389870343990D+00
    w(5) = 0.360761573048138607569833513838D+00
    w(6) = 0.171324492379170345040296142173D+00

  else if ( n == 7 ) then

    x(1) = -0.949107912342758524526189684048D+00
    x(2) = -0.741531185599394439863864773281D+00
    x(3) = -0.405845151377397166906606412077D+00
    x(4) = 0.000000000000000000000000000000D+00
    x(5) = 0.405845151377397166906606412077D+00
    x(6) = 0.741531185599394439863864773281D+00
    x(7) = 0.949107912342758524526189684048D+00

    w(1) = 0.129484966168869693270611432679D+00
    w(2) = 0.279705391489276667901467771424D+00
    w(3) = 0.381830050505118944950369775489D+00
    w(4) = 0.417959183673469387755102040816D+00
    w(5) = 0.381830050505118944950369775489D+00
    w(6) = 0.279705391489276667901467771424D+00
    w(7) = 0.129484966168869693270611432679D+00

  else if ( n == 8 ) then

    x(1) = -0.960289856497536231683560868569D+00
    x(2) = -0.796666477413626739591553936476D+00
    x(3) = -0.525532409916328985817739049189D+00
    x(4) = -0.183434642495649804939476142360D+00
    x(5) = 0.183434642495649804939476142360D+00
    x(6) = 0.525532409916328985817739049189D+00
    x(7) = 0.796666477413626739591553936476D+00
    x(8) = 0.960289856497536231683560868569D+00

    w(1) = 0.101228536290376259152531354310D+00
    w(2) = 0.222381034453374470544355994426D+00
    w(3) = 0.313706645877887287337962201987D+00
    w(4) = 0.362683783378361982965150449277D+00
    w(5) = 0.362683783378361982965150449277D+00
    w(6) = 0.313706645877887287337962201987D+00
    w(7) = 0.222381034453374470544355994426D+00
    w(8) = 0.101228536290376259152531354310D+00

  else if ( n == 9 ) then

    x(1) = -0.968160239507626089835576203D+00
    x(2) = -0.836031107326635794299429788D+00
    x(3) = -0.613371432700590397308702039D+00
    x(4) = -0.324253423403808929038538015D+00
    x(5) = 0.000000000000000000000000000D+00
    x(6) = 0.324253423403808929038538015D+00
    x(7) = 0.613371432700590397308702039D+00
    x(8) = 0.836031107326635794299429788D+00
    x(9) = 0.968160239507626089835576203D+00

    w(1) = 0.081274388361574411971892158111D+00
    w(2) = 0.18064816069485740405847203124D+00
    w(3) = 0.26061069640293546231874286942D+00
    w(4) = 0.31234707704000284006863040658D+00
    w(5) = 0.33023935500125976316452506929D+00
    w(6) = 0.31234707704000284006863040658D+00
    w(7) = 0.26061069640293546231874286942D+00
    w(8) = 0.18064816069485740405847203124D+00
    w(9) = 0.081274388361574411971892158111D+00

  else if ( n == 10 ) then

    x(1) = -0.973906528517171720077964012D+00
    x(2) = -0.865063366688984510732096688D+00
    x(3) = -0.679409568299024406234327365D+00
    x(4) = -0.433395394129247190799265943D+00
    x(5) = -0.148874338981631210884826001D+00
    x(6) = 0.148874338981631210884826001D+00
    x(7) = 0.433395394129247190799265943D+00
    x(8) = 0.679409568299024406234327365D+00
    x(9) = 0.865063366688984510732096688D+00
    x(10) = 0.973906528517171720077964012D+00

    w(1) = 0.066671344308688137593568809893D+00
    w(2) = 0.14945134915058059314577633966D+00
    w(3) = 0.21908636251598204399553493423D+00
    w(4) = 0.26926671930999635509122692157D+00
    w(5) = 0.29552422471475287017389299465D+00
    w(6) = 0.29552422471475287017389299465D+00
    w(7) = 0.26926671930999635509122692157D+00
    w(8) = 0.21908636251598204399553493423D+00
    w(9) = 0.14945134915058059314577633966D+00
    w(10) = 0.066671344308688137593568809893D+00

  else if ( n == 11 ) then

    x(1) = -0.978228658146056992803938001D+00
    x(2) = -0.887062599768095299075157769D+00
    x(3) = -0.730152005574049324093416252D+00
    x(4) = -0.519096129206811815925725669D+00
    x(5) = -0.269543155952344972331531985D+00
    x(6) = 0.000000000000000000000000000D+00
    x(7) = 0.269543155952344972331531985D+00
    x(8) = 0.519096129206811815925725669D+00
    x(9) = 0.730152005574049324093416252D+00
    x(10) = 0.887062599768095299075157769D+00
    x(11) = 0.978228658146056992803938001D+00

    w(1) = 0.055668567116173666482753720443D+00
    w(2) = 0.12558036946490462463469429922D+00
    w(3) = 0.18629021092773425142609764143D+00
    w(4) = 0.23319376459199047991852370484D+00
    w(5) = 0.26280454451024666218068886989D+00
    w(6) = 0.27292508677790063071448352834D+00
    w(7) = 0.26280454451024666218068886989D+00
    w(8) = 0.23319376459199047991852370484D+00
    w(9) = 0.18629021092773425142609764143D+00
    w(10) = 0.12558036946490462463469429922D+00
    w(11) = 0.055668567116173666482753720443D+00

  else if ( n == 12 ) then

    x(1) = -0.981560634246719250690549090D+00
    x(2) = -0.904117256370474856678465866D+00
    x(3) = -0.769902674194304687036893833D+00
    x(4) = -0.587317954286617447296702419D+00
    x(5) = -0.367831498998180193752691537D+00
    x(6) = -0.125233408511468915472441369D+00
    x(7) = 0.125233408511468915472441369D+00
    x(8) = 0.367831498998180193752691537D+00
    x(9) = 0.587317954286617447296702419D+00
    x(10) = 0.769902674194304687036893833D+00
    x(11) = 0.904117256370474856678465866D+00
    x(12) = 0.981560634246719250690549090D+00

    w(1) = 0.047175336386511827194615961485D+00
    w(2) = 0.10693932599531843096025471819D+00
    w(3) = 0.16007832854334622633465252954D+00
    w(4) = 0.20316742672306592174906445581D+00
    w(5) = 0.23349253653835480876084989892D+00
    w(6) = 0.24914704581340278500056243604D+00
    w(7) = 0.24914704581340278500056243604D+00
    w(8) = 0.23349253653835480876084989892D+00
    w(9) = 0.20316742672306592174906445581D+00
    w(10) = 0.16007832854334622633465252954D+00
    w(11) = 0.10693932599531843096025471819D+00
    w(12) = 0.047175336386511827194615961485D+00

  else if ( n == 13 ) then

    x(1) = -0.984183054718588149472829449D+00
    x(2) = -0.917598399222977965206547837D+00
    x(3) = -0.801578090733309912794206490D+00
    x(4) = -0.642349339440340220643984607D+00
    x(5) = -0.448492751036446852877912852D+00
    x(6) = -0.230458315955134794065528121D+00
    x(7) = 0.000000000000000000000000000D+00
    x(8) = 0.230458315955134794065528121D+00
    x(9) = 0.448492751036446852877912852D+00
    x(10) = 0.642349339440340220643984607D+00
    x(11) = 0.80157809073330991279420649D+00
    x(12) = 0.91759839922297796520654784D+00
    x(13) = 0.98418305471858814947282945D+00

    w(1) = 0.040484004765315879520021592201D+00
    w(2) = 0.092121499837728447914421775954D+00
    w(3) = 0.13887351021978723846360177687D+00
    w(4) = 0.17814598076194573828004669200D+00
    w(5) = 0.20781604753688850231252321931D+00
    w(6) = 0.22628318026289723841209018604D+00
    w(7) = 0.23255155323087391019458951527D+00
    w(8) = 0.22628318026289723841209018604D+00
    w(9) = 0.20781604753688850231252321931D+00
    w(10) = 0.17814598076194573828004669200D+00
    w(11) = 0.13887351021978723846360177687D+00
    w(12) = 0.092121499837728447914421775954D+00
    w(13) = 0.040484004765315879520021592201D+00

  else if ( n == 14 ) then

    x(1) = -0.986283808696812338841597267D+00
    x(2) = -0.928434883663573517336391139D+00
    x(3) = -0.827201315069764993189794743D+00
    x(4) = -0.687292904811685470148019803D+00
    x(5) = -0.515248636358154091965290719D+00
    x(6) = -0.319112368927889760435671824D+00
    x(7) = -0.108054948707343662066244650D+00
    x(8) = 0.108054948707343662066244650D+00
    x(9) = 0.31911236892788976043567182D+00
    x(10) = 0.51524863635815409196529072D+00
    x(11) = 0.68729290481168547014801980D+00
    x(12) = 0.82720131506976499318979474D+00
    x(13) = 0.92843488366357351733639114D+00
    x(14) = 0.98628380869681233884159727D+00

    w(1) = 0.035119460331751863031832876138D+00
    w(2) = 0.08015808715976020980563327706D+00
    w(3) = 0.12151857068790318468941480907D+00
    w(4) = 0.15720316715819353456960193862D+00
    w(5) = 0.18553839747793781374171659013D+00
    w(6) = 0.20519846372129560396592406566D+00
    w(7) = 0.21526385346315779019587644332D+00
    w(8) = 0.21526385346315779019587644332D+00
    w(9) = 0.20519846372129560396592406566D+00
    w(10) = 0.18553839747793781374171659013D+00
    w(11) = 0.15720316715819353456960193862D+00
    w(12) = 0.12151857068790318468941480907D+00
    w(13) = 0.08015808715976020980563327706D+00
    w(14) = 0.035119460331751863031832876138D+00

  else if ( n == 15 ) then

    x(1) = -0.987992518020485428489565719D+00
    x(2) = -0.937273392400705904307758948D+00
    x(3) = -0.848206583410427216200648321D+00
    x(4) = -0.724417731360170047416186055D+00
    x(5) = -0.570972172608538847537226737D+00
    x(6) = -0.394151347077563369897207371D+00
    x(7) = -0.201194093997434522300628303D+00
    x(8) = 0.00000000000000000000000000D+00
    x(9) = 0.20119409399743452230062830D+00
    x(10) = 0.39415134707756336989720737D+00
    x(11) = 0.57097217260853884753722674D+00
    x(12) = 0.72441773136017004741618605D+00
    x(13) = 0.84820658341042721620064832D+00
    x(14) = 0.93727339240070590430775895D+00
    x(15) = 0.98799251802048542848956572D+00

    w(1) = 0.030753241996117268354628393577D+00
    w(2) = 0.070366047488108124709267416451D+00
    w(3) = 0.107159220467171935011869546686D+00
    w(4) = 0.13957067792615431444780479451D+00
    w(5) = 0.16626920581699393355320086048D+00
    w(6) = 0.18616100001556221102680056187D+00
    w(7) = 0.19843148532711157645611832644D+00
    w(8) = 0.20257824192556127288062019997D+00
    w(9) = 0.19843148532711157645611832644D+00
    w(10) = 0.18616100001556221102680056187D+00
    w(11) = 0.16626920581699393355320086048D+00
    w(12) = 0.13957067792615431444780479451D+00
    w(13) = 0.107159220467171935011869546686D+00
    w(14) = 0.070366047488108124709267416451D+00
    w(15) = 0.030753241996117268354628393577D+00

  else if ( n == 16 ) then

    x(1) = -0.989400934991649932596154173D+00
    x(2) = -0.944575023073232576077988416D+00
    x(3) = -0.865631202387831743880467898D+00
    x(4) = -0.755404408355003033895101195D+00
    x(5) = -0.617876244402643748446671764D+00
    x(6) = -0.458016777657227386342419443D+00
    x(7) = -0.281603550779258913230460501D+00
    x(8) = -0.09501250983763744018531934D+00
    x(9) = 0.09501250983763744018531934D+00
    x(10) = 0.28160355077925891323046050D+00
    x(11) = 0.45801677765722738634241944D+00
    x(12) = 0.61787624440264374844667176D+00
    x(13) = 0.75540440835500303389510119D+00
    x(14) = 0.86563120238783174388046790D+00
    x(15) = 0.94457502307323257607798842D+00
    x(16) = 0.98940093499164993259615417D+00

    w(1) = 0.027152459411754094851780572456D+00
    w(2) = 0.062253523938647892862843836994D+00
    w(3) = 0.09515851168249278480992510760D+00
    w(4) = 0.12462897125553387205247628219D+00
    w(5) = 0.14959598881657673208150173055D+00
    w(6) = 0.16915651939500253818931207903D+00
    w(7) = 0.18260341504492358886676366797D+00
    w(8) = 0.18945061045506849628539672321D+00
    w(9) = 0.18945061045506849628539672321D+00
    w(10) = 0.18260341504492358886676366797D+00
    w(11) = 0.16915651939500253818931207903D+00
    w(12) = 0.14959598881657673208150173055D+00
    w(13) = 0.12462897125553387205247628219D+00
    w(14) = 0.09515851168249278480992510760D+00
    w(15) = 0.062253523938647892862843836994D+00
    w(16) = 0.027152459411754094851780572456D+00

  else

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'LEGENDRE_SET - Fatal error!'
    write ( *, '(a,i8)' ) '  Illegal value of N = ', n
    write ( *, '(a)' ) &
      '  Legal values are 1 through 16.'
    stop 1

  end if

  return
end
subroutine r8mat_fs ( n, a, b, info )

!*****************************************************************************80
!
!! R8MAT_FS factors and solves a system with one right hand side.
!
!  Discussion:
!
!    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
!
!    This routine differs from R8MAT_FSS in two ways:
!    * only one right hand side is allowed;
!    * the input matrix A is not modified.
!
!    This routine uses partial pivoting, but no pivot vector is required.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    21 January 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the order of the matrix.
!    N must be positive.
!
!    Input/output, real ( kind = rk ) A(N,N).
!    On input, A is the coefficient matrix of the linear system.
!    On output, A is in unit upper triangular form, and
!    represents the U factor of an LU factorization of the
!    original coefficient matrix.
!
!    Input/output, real ( kind = rk ) B(N).
!    On input, the right hand side of the linear system.
!    On output, the solution of the linear systems.
!
!    Output, integer INFO, singularity flag.
!    0, no singularity detected.
!    nonzero, the factorization failed on the INFO-th step.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  real ( kind = rk ) a2(n,n)
  real ( kind = rk ) b(n)
  integer i
  integer info
  integer ipiv
  integer jcol
  real ( kind = rk ) piv
  real ( kind = rk ) row(n)
  real ( kind = rk ) t
  real ( kind = rk ) temp

  a2(1:n,1:n) = a(1:n,1:n)

  info = 0

  do jcol = 1, n
!
!  Find the maximum element in column I.
!
    piv = abs ( a2(jcol,jcol) )
    ipiv = jcol
    do i = jcol + 1, n
      if ( piv < abs ( a2(i,jcol) ) ) then
        piv = abs ( a2(i,jcol) )
        ipiv = i
      end if
    end do

    if ( piv == 0.0D+00 ) then
      info = jcol
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8MAT_FS - Fatal error!'
      write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      stop 1
    end if
!
!  Switch rows JCOL and IPIV, and B.
!
    if ( jcol /= ipiv ) then

      row(1:n) = a2(jcol,1:n)
      a2(jcol,1:n) = a2(ipiv,1:n)
      a2(ipiv,1:n) = row(1:n)

      t       = b(jcol)
      b(jcol) = b(ipiv)
      b(ipiv) = t

    end if
!
!  Scale the pivot row.
!
    a2(jcol,jcol+1:n) = a2(jcol,jcol+1:n) / a2(jcol,jcol)
    b(jcol) = b(jcol) / a2(jcol,jcol)
    a2(jcol,jcol) = 1.0D+00
!
!  Use the pivot row to eliminate lower entries in that column.
!
    do i = jcol + 1, n
      if ( a2(i,jcol) /= 0.0D+00 ) then
        temp = - a2(i,jcol)
        a2(i,jcol) = 0.0D+00
        a2(i,jcol+1:n) = a2(i,jcol+1:n) + temp * a2(jcol,jcol+1:n)
        b(i) = b(i) + temp * b(jcol)
      end if
    end do

  end do
!
!  Back solve.
!
  do jcol = n, 2, -1
    b(1:jcol-1) = b(1:jcol-1) - a2(1:jcol-1,jcol) * b(jcol)
  end do

  return
end
subroutine r8vec_linspace ( n, a, b, x )

!*****************************************************************************80
!
!! R8VEC_LINSPACE creates a vector of linearly spaced values.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
!
!    In other words, the interval is divided into N-1 even subintervals,
!    and the endpoints of intervals are used as the points.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    14 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real ( kind = rk ) A, B, the first and last entries.
!
!    Output, real ( kind = rk ) X(N), a vector of linearly spaced data.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a
  real ( kind = rk ) b
  integer i
  real ( kind = rk ) x(n)

  if ( n == 1 ) then

    x(1) = ( a + b ) / 2.0D+00

  else

    do i = 1, n
      x(i) = ( real ( n - i,     kind = rk ) * a   &
             + real (     i - 1, kind = rk ) * b ) &
             / real ( n     - 1, kind = rk )
    end do

  end if

  return
end
subroutine r8vec_print ( n, a, title )

!*****************************************************************************80
!
!! R8VEC_PRINT prints an R8VEC.
!
!  Discussion:
!
!    An R8VEC is a vector of R8's.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    22 August 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer N, the number of components of the vector.
!
!    Input, real ( kind = rk ) A(N), the vector to be printed.
!
!    Input, character ( len = * ) TITLE, a title.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n)
  integer i
  character ( len = * ) title

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '

  do i = 1, n
    write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
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
