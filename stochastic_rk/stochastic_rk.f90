subroutine rk1_ti_step ( x, t, h, q, fi, gi, xstar )

!*****************************************************************************80
!
!! rk1_ti_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is first-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FI, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GI, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ), external :: fi
  real ( kind = rk8 ), external :: gi
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) xstar

  a21 =   1.0D+00

  q1 = 1.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  xstar = x1 + a21 * k1

  return
end
subroutine rk2_ti_step ( x, t, h, q, fi, gi, xstar )

!*****************************************************************************80
!
!! rk2_ti_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is second-order, and suitable for time-invariant
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FI, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GI, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ) a31
  real ( kind = rk8 ) a32
  real ( kind = rk8 ), external :: fi
  real ( kind = rk8 ), external :: gi
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) k2
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) q2
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) w2
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) xstar

  a21 =   1.0D+00
  a31 =   0.5D+00
  a32 =   0.5D+00

  q1 = 2.0D+00
  q2 = 2.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  xstar = x1 + a31 * k1 + a32 * k2

  return
end
subroutine rk3_ti_step ( x, t, h, q, fi, gi, xstar )

!*****************************************************************************80
!
!! rk3_ti_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is third-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FI, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GI, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ) a31
  real ( kind = rk8 ) a32
  real ( kind = rk8 ) a41
  real ( kind = rk8 ) a42
  real ( kind = rk8 ) a43
  real ( kind = rk8 ), external :: fi
  real ( kind = rk8 ), external :: gi
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) k2
  real ( kind = rk8 ) k3
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) q2
  real ( kind = rk8 ) q3
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) t3
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) w2
  real ( kind = rk8 ) w3
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) x3
  real ( kind = rk8 ) xstar

  a21 =   1.52880952525675D+00
  a31 =   0.0D+00
  a32 =   0.51578733443615D+00
  a41 =   0.53289582961739D+00
  a42 =   0.25574324768195D+00
  a43 =   0.21136092270067D+00

  q1 = 1.87653936176981D+00
  q2 = 3.91017166264989D+00
  q3 = 4.73124353935667D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( ) * sqrt ( q3 * q / h )
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

  xstar = x1 + a41 * k1 + a42 * k2 + a43 * k3

  return
end
subroutine rk4_ti_step ( x, t, h, q, fi, gi, xstar )

!*****************************************************************************80
!
!! rk4_ti_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
!    systems in which F and G do not depend explicitly on time.
!
!    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FI, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GI, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ) a31
  real ( kind = rk8 ) a32
  real ( kind = rk8 ) a41
  real ( kind = rk8 ) a42
  real ( kind = rk8 ) a43
  real ( kind = rk8 ) a51
  real ( kind = rk8 ) a52
  real ( kind = rk8 ) a53
  real ( kind = rk8 ) a54
  real ( kind = rk8 ), external :: fi
  real ( kind = rk8 ), external :: gi
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) k2
  real ( kind = rk8 ) k3
  real ( kind = rk8 ) k4
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) q2
  real ( kind = rk8 ) q3
  real ( kind = rk8 ) q4
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) t3
  real ( kind = rk8 ) t4
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) w2
  real ( kind = rk8 ) w3
  real ( kind = rk8 ) w4
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) x3
  real ( kind = rk8 ) x4
  real ( kind = rk8 ) xstar

  a21 =   2.71644396264860D+00
  a31 = - 6.95653259006152D+00
  a32 =   0.78313689457981D+00
  a41 =   0.0D+00
  a42 =   0.48257353309214D+00
  a43 =   0.26171080165848D+00
  a51 =   0.47012396888046D+00
  a52 =   0.36597075368373D+00
  a53 =   0.08906615686702D+00
  a54 =   0.07483912056879D+00

  q1 =   2.12709852335625D+00
  q2 =   2.73245878238737D+00
  q3 =  11.22760917474960D+00
  q4 =  13.36199560336697D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( ) * sqrt ( q2 * q / h )
  k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( ) * sqrt ( q3 * q / h )
  k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

  t4 = t1 + a41 * h  + a42 * h + a43 * h
  x4 = x1 + a41 * k1 + a42 * k2
  w4 = r8_normal_01 ( ) * sqrt ( q4 * q / h )
  k4 = h * fi ( x4 ) + h * gi ( x4 ) * w4

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

  return
end
subroutine rk1_tv_step ( x, t, h, q, fv, gv, xstar )

!*****************************************************************************80
!
!! rk1_tv_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is first-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FV, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GV, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ), external :: fv
  real ( kind = rk8 ), external :: gv
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) xstar

  a21 =   1.0D+00

  q1 = 1.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  xstar = x1 + a21 * k1

  return
end
subroutine rk2_tv_step ( x, t, h, q, fv, gv, xstar )

!*****************************************************************************80
!
!! rk2_tv_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is second-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FV, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GV, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ) a31
  real ( kind = rk8 ) a32
  real ( kind = rk8 ), external :: fv
  real ( kind = rk8 ), external :: gv
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) k2
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) q2
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) w2
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) xstar

  a21 =   1.0D+00
  a31 =   0.5D+00
  a32 =   0.5D+00

  q1 = 2.0D+00
  q2 = 2.0D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( ) * sqrt ( q2 * q / h )
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

  xstar = x1 + a31 * k1 + a32 * k2

  return
end
subroutine rk4_tv_step ( x, t, h, q, fv, gv, xstar )

!*****************************************************************************80
!
!! rk4_tv_step() takes one step of a stochastic Runge Kutta scheme.
!
!  Discussion:
!
!    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
!    systems.
!
!    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 June 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Jeremy Kasdin,
!    Runge-Kutta algorithm for the numerical integration of
!    stochastic differential equations,
!    Journal of Guidance, Control, and Dynamics,
!    Volume 18, Number 1, January-February 1995, pages 114-120.
!
!    Jeremy Kasdin,
!    Discrete Simulation of Colored Noise and Stochastic Processes
!    and 1/f^a Power Law Noise Generation,
!    Proceedings of the IEEE,
!    Volume 83, Number 5, 1995, pages 802-827.
!
!  Input:
!
!    real ( kind = rk8 ) X, the value at the current time.
!
!    real ( kind = rk8 ) T, the current time.
!
!    real ( kind = rk8 ) H, the time step.
!
!    real ( kind = rk8 ) Q, the spectral density of the input white noise.
!
!    external real ( kind = rk8 ) FV, the name of the deterministic
!    right hand side function.
!
!    external real ( kind = rk8 ) GV, the name of the stochastic
!    right hand side function.
!
!  Output:
!
!    real ( kind = rk8 ) XSTAR, the value at time T+H.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ) a21
  real ( kind = rk8 ) a31
  real ( kind = rk8 ) a32
  real ( kind = rk8 ) a41
  real ( kind = rk8 ) a42
  real ( kind = rk8 ) a43
  real ( kind = rk8 ) a51
  real ( kind = rk8 ) a52
  real ( kind = rk8 ) a53
  real ( kind = rk8 ) a54
  real ( kind = rk8 ), external :: fv
  real ( kind = rk8 ), external :: gv
  real ( kind = rk8 ) h
  real ( kind = rk8 ) k1
  real ( kind = rk8 ) k2
  real ( kind = rk8 ) k3
  real ( kind = rk8 ) k4
  real ( kind = rk8 ) q
  real ( kind = rk8 ) q1
  real ( kind = rk8 ) q2
  real ( kind = rk8 ) q3
  real ( kind = rk8 ) q4
  real ( kind = rk8 ) r8_normal_01
  real ( kind = rk8 ) t
  real ( kind = rk8 ) t1
  real ( kind = rk8 ) t2
  real ( kind = rk8 ) t3
  real ( kind = rk8 ) t4
  real ( kind = rk8 ) w1
  real ( kind = rk8 ) w2
  real ( kind = rk8 ) w3
  real ( kind = rk8 ) w4
  real ( kind = rk8 ) x
  real ( kind = rk8 ) x1
  real ( kind = rk8 ) x2
  real ( kind = rk8 ) x3
  real ( kind = rk8 ) x4
  real ( kind = rk8 ) xstar

  a21 =   0.66667754298442D+00
  a31 =   0.63493935027993D+00
  a32 =   0.00342761715422D+00
  a41 = - 2.32428921184321D+00
  a42 =   2.69723745129487D+00
  a43 =   0.29093673271592D+00
  a51 =   0.25001351164789D+00
  a52 =   0.67428574806272D+00
  a53 = - 0.00831795169360D+00
  a54 =   0.08401868181222D+00

  q1 = 3.99956364361748D+00
  q2 = 1.64524970733585D+00
  q3 = 1.59330355118722D+00
  q4 = 0.26330006501868D+00

  t1 = t
  x1 = x
  w1 = r8_normal_01 ( ) * sqrt ( q1 * q / h )
  k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

  t2 = t1 + a21 * h
  x2 = x1 + a21 * k1
  w2 = r8_normal_01 ( ) * sqrt ( q2 * q / h )
  k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

  t3 = t1 + a31 * h  + a32 * h
  x3 = x1 + a31 * k1 + a32 * k2
  w3 = r8_normal_01 ( ) * sqrt ( q3 * q / h )
  k3 = h * fv ( t3, x3 ) + h * gv ( t3, x3 ) * w3

  t4 = t1 + a41 * h  + a42 * h  + a43 * h
  x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3
  w4 = r8_normal_01 ( ) * sqrt ( q4 * q / h )
  k4 = h * fv ( t4, x4 ) + h * gv ( t4, x4 ) * w4

  xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

  return
end
function r8_normal_01 ( )

!*****************************************************************************80
!
!! r8_normal_01() returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!    Because this routine uses the Box Muller method, it requires pairs
!    of uniform random values to generate a pair of normal random values.
!    This means that on every other call, the code can use the second
!    value that it calculated.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Output:
!
!    real ( kind = rk8 ) R8_NORMAL_01, a normally distributed random value.
!
  implicit none

  integer, parameter :: rk8 = kind ( 1.0d+00 )

  real ( kind = rk8 ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk8 ) r1
  real ( kind = rk8 ) r2
  real ( kind = rk8 ) r8_normal_01
  integer, parameter :: two = 2
  integer, save :: used = 0
  real ( kind = rk8 ) v1
  real ( kind = rk8 ), save :: v2 = 0.0D+00
!
!  If USED is even, generate two uniforms, create two normals,
!  return the first normal.
!
  if ( mod ( used, two ) == 0 ) then

    call random_number ( harvest = r1 )
    call random_number ( harvest = r2 )

    v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
    v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

    r8_normal_01 = v1
!
!  If USED is odd, return the second normal.
!
  else

    r8_normal_01 = v2

  end if

  used = used + 1

  return
end

