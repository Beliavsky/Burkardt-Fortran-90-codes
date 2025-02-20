subroutine f_exact ( z, f, fz, fzz )

!*****************************************************************************80
!
!! f_exact() computes a function related to an exact solution of the KPP Fisher PDE.
!
!  Discussion:
!
!    ut = uxx + u * ( 1 - u )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Ablowitz, Anthony Zeppetella,
!    Explicit solutions of Fisher's equation for a special wave speed,
!    Bulletin of Mathematical Biology,
!    Volume 41, pages 835-840, 1979.
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) Z: the evaluation point.
!
!  Output:
!
!    real ( kind = rk ) F, FZ, FZZ: the function and derivatives at Z.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: c_in
      real ( kind = rk ), optional :: c_out
      real ( kind = rk ), optional :: k_in
      real ( kind = rk ), optional :: k_out
    end subroutine
  end interface

  real ( kind = rk ) a
  real ( kind = rk ) f
  real ( kind = rk ) fz
  real ( kind = rk ) fzz
  real ( kind = rk ) k
  real ( kind = rk ) z

  call fisher_parameters ( a_out = a, k_out = k )

  f =     1.0 / ( 1.0 + a * exp ( k * z ) )**2
  fz =  - 2.0 / ( 1.0 + a * exp ( k * z ) )**3 * a    * k    * exp (       k * z )
  fzz = + 6.0 / ( 1.0 + a * exp ( k * z ) )**4 * a**2 * k**2 * exp ( 2.0 * k * z ) &
        - 2.0 / ( 1.0 + a * exp ( k * z ) )**3 * a    * k**2 * exp (       k * z ) 

  return
end
subroutine fisher_exact ( t, x, u, ut, ux, uxx )

!*****************************************************************************80
!
!! fisher_exact() computes an exact solution of the KPP Fisher equation.
!
!  Discussion:
!
!    ut = uxx + u * ( 1 - u )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Ablowitz, Anthony Zeppetella,
!    Explicit solutions of Fisher's equation for a special wave speed,
!    Bulletin of Mathematical Biology,
!    Volume 41, pages 835-840, 1979.
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) T, X: the time and position.
!
!  Output:
!
!    real ( kind = rk ) U, UT, UX, UXX: the solution and derivatives at (T,X).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: c_in
      real ( kind = rk ), optional :: c_out
      real ( kind = rk ), optional :: k_in
      real ( kind = rk ), optional :: k_out
    end subroutine
  end interface

  real ( kind = rk ) a
  real ( kind = rk ) c
  real ( kind = rk ) k
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) x
  real ( kind = rk ) z

  call fisher_parameters ( a_out = a, c_out = c, k_out = k )

  z = x - c * t

  u =     1.0     / ( 1.0 + a * exp ( k * z ) )**2
  ut =    2.0 * c / ( 1.0 + a * exp ( k * z ) )**3 * a    * k    * exp (       k * z )
  ux =  - 2.0     / ( 1.0 + a * exp ( k * z ) )**3 * a    * k    * exp (       k * z )
  uxx = + 6.0     / ( 1.0 + a * exp ( k * z ) )**4 * a**2 * k**2 * exp ( 2.0 * k * z ) &
        - 2.0     / ( 1.0 + a * exp ( k * z ) )**3 * a    * k**2 * exp (       k * z ) 

  return
end
subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )

!*****************************************************************************80
!
!! fisher_parameters() returns parameters for the Fisher ODE.
!
!  Discussion:
!
!    If input values are specified, this resets the default parameters.
!    Otherwise, the output will be the current defaults.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    real ( kind = rk ) a_in: a new value for A.
!
!    real ( kind = rk ) c_in: a new value for C.
!
!    real ( kind = rk ) k_in: a new value for K.
!
!  Output:
!
!    real ( kind = rk ) a_out: the default or new value for A.
!
!    real ( kind = rk ) c_out: the default or new value for C.
!
!    real ( kind = rk ) k_out: the default or new value for K.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), save     :: a_default = 2.0D+00
  real ( kind = rk ), optional :: a_in
  real ( kind = rk ), optional :: a_out
  real ( kind = rk ), save     :: c_default = 5.0D+00 / sqrt ( 6.0D+00 )
  real ( kind = rk ), optional :: c_in
  real ( kind = rk ), optional :: c_out
  real ( kind = rk ), save     :: k_default = 1.0D+00 / sqrt ( 6.0D+00 )
  real ( kind = rk ), optional :: k_in
  real ( kind = rk ), optional :: k_out
!
!  Update defaults if input was supplied.
!
  if ( present ( a_in ) ) then
    a_default = a_in
  end if

  if ( present ( c_in ) ) then
    c_default = c_in
  end if

  if ( present ( k_in ) ) then
    k_default = k_in
  end if
!
!  Return values.
!
  if ( present ( a_out ) ) then
    a_out = a_default
  end if

  if ( present ( c_out ) ) then
    c_out = c_default
  end if

  if ( present ( k_out ) ) then
    k_out = k_default
  end if

  return
end
subroutine fisher_residual ( t, x, r )

!*****************************************************************************80
!
!! fisher_residual() computes the residual of the KPP Fisher equation.
!
!  Discussion:
!
!    ut = uxx + u * ( 1 - u )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
 !  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Ablowitz, Anthony Zeppetella,
!    Explicit solutions of Fisher's equation for a special wave speed,
!    Bulletin of Mathematical Biology,
!    Volume 41, pages 835-840, 1979.
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) T, X: the time and position.
!
!  Output:
!
!    real ( kind = rk ) R: the residual at that time and position.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) r
  real ( kind = rk ) t
  real ( kind = rk ) u
  real ( kind = rk ) ut
  real ( kind = rk ) ux
  real ( kind = rk ) uxx
  real ( kind = rk ) x

  call fisher_exact ( t, x, u, ut, ux, uxx )
  r = ut - uxx - u * ( 1.0 - u )

  return
end
subroutine f_residual ( f, fz, fzz, r )

!*****************************************************************************80
!
!! f_residual() computes a residual related to the KPP Fisher equation.
!
!  Discussion:
!
!    The KPP Fisher equation is
!
!      ut = uxx + u * ( 1 - u )
!
!    For a particular value of c, and function f(z), the KPP Fisher 
!    equation admits a traveling wave solution of the form
!
!      u(x,t) = f(x-ct) = f(z)
!
!    where the function f(z) satisfies the ODE
!
!      f'' + cf' + f - f^2 = 0
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 May 2024
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Mark Ablowitz, Anthony Zeppetella,
!    Explicit solutions of Fisher's equation for a special wave speed,
!    Bulletin of Mathematical Biology,
!    Volume 41, pages 835-840, 1979.
!
!    Daniel Arrigo,
!    Analytical Techniques for Solving Nonlinear Partial Differential Equations,
!    Morgan and Clayfoot, 2019,
!    ISBN: 978 168 173 5351.
!
!  Input:
!
!    real ( kind = rk ) F, FZ, FZZ: the associated function and derivatives.
!
!  Output:
!
!    real ( kind = rk ) R: the residual of the f(z) solution.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  interface
    subroutine fisher_parameters ( a_in, c_in, k_in, a_out, c_out, k_out )
      integer, parameter :: rk = kind ( 1.0D+00 )
      real ( kind = rk ), optional :: a_in
      real ( kind = rk ), optional :: a_out
      real ( kind = rk ), optional :: c_in
      real ( kind = rk ), optional :: c_out
      real ( kind = rk ), optional :: k_in
      real ( kind = rk ), optional :: k_out
    end subroutine
  end interface

  real ( kind = rk ) c
  real ( kind = rk ) f
  real ( kind = rk ) fz
  real ( kind = rk ) fzz
  real ( kind = rk ) r

  call fisher_parameters ( c_out = c )

  r = fzz + c * fz + f - f**2

  return
end
