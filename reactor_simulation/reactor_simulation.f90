program main

!*****************************************************************************80
!
!! reactor_simulation() simulations the reactor shielding problem.
!
!  Discussion:
!
!    This is a Monte Carlo simulation, using
!    uniform random numbers, which investigates the
!    effectiveness of a shield intended to absorb the
!    neutrons emitted from a nuclear reactor.
!   
!    The reactor is modeled as a point source,
!    located at (0,0,0).
!   
!    A particle emitted from the reactor has a random
!    initial direction, and an energy selected from
!    [Emin,Emax] with a 1/Sqrt(E) distribution.
!   
!    The shield is modeled as a wall of thickness THICK,
!    extending from 0 to THICK in the X direction, and
!    extending forever in the Y and Z directions.
!   
!    Based on the particle energy, a distance D is computed
!    which measures how far the particle could travel through
!    the shield before colliding.
!   
!    Based on the particle direction, the position is updated
!    by D units.
!   
!    If the particle is now to the left of the shield, it is
!    counted as being REFLECTED.
!   
!    If the particle is to the right of the shield, it is 
!    counted as being ABSORBED.
!   
!    If the particle is inside the shield, it has COLLIDED.
!    A particle that collides is either absorbed (end of story)
!    or SCATTERED with a new random direction and a new (lower)
!    energy.
!   
!    Every particle is followed from origin to its final fate,
!    which is reflection, transmission, or absorption.
!    At the end, a summary is printed, giving the number of
!    particles with each fate, and the average energy of each
!    group of particles.
!   
!    Increasing NTOT, the number of particles used, will improve the
!    expected reliability of the results.
!   
!    Increasing THICK, the thickness of the shield, should 
!    result in more absorptions and reflections.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Local Parameters:
!
!    Local, real ( kind = rk ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Local, real ( kind = rk ) D, the distance that the particle can
!    travel through the slab, given its current energy.
!
!    Local, real ( kind = rk ) E, the energy of the particle.
!
!    Local, real ( kind = rk ) EA, energy absorbed by the slab.
!
!    Local, real ( kind = rk ) ER, energy reflected by the slab.
!
!    Local, real ( kind = rk ) ET, energy transmitted through the slab.
!
!    Local, real ( kind = rk ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Local, integer NA, number of particles absorbed by the slab.
!
!    Local, integer NPART, the index of the current particle.
!
!    Local, integer NR, number of particles reflected by the slab.
!
!    Local, integer NT, number of particles transmitted 
!    by the slab.
!
!    Local, integer NTOT, the total number of particles to be
!    emitted from the neutron source.
!
!    Local, real ( kind = rk ) SA, standard deviation of absorbed energy.
!
!    Local, real ( kind = rk ) SR, standard deviation of reflected energy.
!
!    Local, real ( kind = rk ) ST, standard deviation of transmitted energy.
!
!    Local, real ( kind = rk ) THICK, the thickness of the slab that is
!    intended to absorb most of the particles.
!
!    Local, real ( kind = rk ) X, Y, Z, the current position of the particle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical absorb
  real ( kind = rk ) azm
  real ( kind = rk ) d
  real ( kind = rk ) dist2c
  real ( kind = rk ) e
  real ( kind = rk ) ea
  real ( kind = rk ) er
  real ( kind = rk ) et
  real ( kind = rk ) mu
  integer na
  integer npart
  integer nr
  integer nt
  integer, parameter :: ntot = 100000
  real ( kind = rk ) sa
  real ( kind = rk ) sr
  real ( kind = rk ) st
  integer test
  integer, parameter :: test_num = 5
  real ( kind = rk ), parameter :: thick = 2.0D+00
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'reactor_simulation():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  The reactor shielding simulation.'
  write ( *, '(a)' ) ' '
  write ( *, '(a,g14.6)' ) '  Shield thickness is THICK = ', thick
  write ( *, '(a,i8)' ) '  Number of simulated particles is NTOT = ', ntot
  write ( *, '(a,i8)' ) '  Number of tests TEST_NUM = ', test_num

  do test = 1, test_num

    write ( *, '(a)' ) ' '
    write ( *, '(a,i2)' ) '  Test # ', test
!
!  Initialize.
!
    ea = 0.0D+00
    er = 0.0D+00
    et = 0.0D+00
    na = 0
    nr = 0
    nt = 0
    sa = 0.0D+00
    sr = 0.0D+00
    st = 0.0D+00
!
!  Loop over the particles.
!
    do npart = 1, ntot
!
!  Generate a new particle.
!
      call source ( e, mu, azm, x, y, z )

      do
!
!  Compute the distance that the particle can travel through the slab,
!  based on its current energy.
!
        d = dist2c ( e )
!
!  Update the particle's position by D units.
!
        call update ( mu, azm, d, x, y, z )
!
!  The particle was reflected by the shield, and this path is complete.
!
        if ( x < 0.0D+00 ) then

          nr = nr + 1
          er = er + e
          sr = sr + e * e
          exit
!
!  The particle was transmitted through the shield, and this path is complete.
!
        else if ( thick < x ) then

          nt = nt + 1
          et = et + e
          st = st + e * e
          exit
!
!  The particle collided with the shield, and was absorbed.  This path is done.
!
        else if ( absorb ( ) ) then

          na = na + 1
          ea = ea + e
          sa = sa + e * e
          exit
!
!  The particle collided with the shield and was scattered.
!  Find the scattering angle and energy, and continue along the new path.
!
        else

          call scatter ( e, mu, azm )

        end if

      end do

    end do
!
!  Print the results of the simulation.
!
    call output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

  end do
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'REACTOR_SIMULATION:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
function absorb ( )

!*****************************************************************************80
!
!! ABSORB determines if a colliding particle is absorbed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, logical ABSORB, is TRUE if the particle is absorbed.
!
!  Local:
!
!    real ( kind = rk ) PA, the probability of absorption.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  logical absorb
  real ( kind = rk ), parameter :: pa = 0.1D+00
  real ( kind = rk ) u

  call random_number ( harvest = u )

  if ( u <= pa ) then
    absorb = .true.
  else
    absorb = .false.
  end if

  return
end
function cross ( e )

!*****************************************************************************80
!
!! CROSS returns the "cross section" of a particle based on its energy.
!
!  Discussion:
!
!    The particle's cross section is a measure of its likelihood to collide
!    with the material of the slab.  This quantity typically depends on both
!    the particle's energy and the kind of medium through which it is traveling.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = rk ) E, the energy of the particle.
!
!    Output, real ( kind = rk ) CROSS, the cross section.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cross
  real ( kind = rk ) e
  real ( kind = rk ) s
  real ( kind = rk ) y

  s = abs ( sin ( 100.0D+00 * ( exp ( e ) - 1.0D+00 ) ) &
    + sin ( 18.81D+00 * ( exp ( e ) - 1.0D+00 ) ) )

  y = max ( 0.02D+00, s )

  cross = 10.0D+00 * exp ( -0.1D+00 / y )

  return
end
function dist2c ( e )

!*****************************************************************************80
!
!! DIST2C returns the distance to collision.
!
!  Discussion:
!
!    Assuming the particle has a given energy, and assuming it is currently
!    somewhere inside the shield, it is possible to determine a typical distance
!    which the particle can travel before it collides with the material of
!    the shield.
!
!    The computation of the collision distance is made by estimating a
!    "cross section" (as though having more energy made the particle "bigger"
!    and hence more likely to collide) and then randomly selecting a distance
!    that is logarithmically distributed.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = rk ) E, the energy of the particle.
!
!    Output, real ( kind = rk ) DIST2C, the distance the particle can travel
!    through the slab before colliding.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cross
  real ( kind = rk ) dist2c
  real ( kind = rk ) e
  real ( kind = rk ) u

  call random_number ( harvest = u )

  dist2c = - log ( u ) / cross ( e )

  return
end
function energy ( )

!*****************************************************************************80
!
!! ENERGY assigns an energy to an emitted particle.
!
!  Discussion:
!
!    The energy E is in the range [EMIN,EMAX], with distribution
!    const/sqrt(energy).
!
!    An inverse function approach is used to compute this.
!
!    The energies are measured in MeV.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, real ( kind = rk ) ENERGY, a randomly chosen energy that is
!    distributed as described above.
!
!  Local:
!
!    real ( kind = rk ) EMIN, EMAX, the minimum and maximum
!    energies.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) c
  real ( kind = rk ), parameter :: emax = 2.5D+00
  real ( kind = rk ), parameter :: emin = 1.0D-03
  real ( kind = rk ) energy
  real ( kind = rk ) u

  call random_number ( harvest = u )

  c = 1.0D+00 / ( 2.0D+00 * ( sqrt ( emax ) - sqrt ( emin ) ) )

  energy = ( u / ( 2.0D+00 * c ) + sqrt ( emin ) )
  energy = energy * energy

  return
end
subroutine output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

!*****************************************************************************80
!
!! OUTPUT prints the results of the reactor shielding simulation.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, integer NA, number of particles absorbed by the slab.
!
!    Input, real ( kind = rk ) EA, energy absorbed by the slab.
!
!    Input, real ( kind = rk ) SA, the sum of the squares of the 
!    absorbed energies.
!
!    Input, integer NR, number of particles reflected by the slab.
!
!    Input, real ( kind = rk ) ER, energy reflected by the slab.
!
!    Input, real ( kind = rk ) SR, the sum of the squares of the 
!    reflected energies.
!
!    Input, integer NT, number of particles transmitted
!    by the slab.
!
!    Input, real ( kind = rk ) ET, energy transmitted through the slab.
!
!    Input, real ( kind = rk ) ST, the sum of the squares of the 
!    transmitted energies.
!
!    Input, integer NTOT, the total number of particles.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ea
  real ( kind = rk ) ea_ave
  real ( kind = rk ) er
  real ( kind = rk ) er_ave
  real ( kind = rk ) et
  real ( kind = rk ) et_ave
  real ( kind = rk ) etot
  integer na
  integer nr
  integer nt
  integer ntot
  real ( kind = rk ) pa
  real ( kind = rk ) pr
  real ( kind = rk ) pt
  real ( kind = rk ) ptot
  real ( kind = rk ) sa
  real ( kind = rk ) sr
  real ( kind = rk ) st

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  The Reactor Shielding Problem:'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) &
    '                           Total                   Average'
  write ( *, '(a)' ) '                    #      Energy      ' &
    // 'Percent     Energy         StDev'
  write ( *, '(a)' ) ' '

  etot = ea + er + et

  if ( 0 < na ) then
    ea_ave = ea / real ( na, kind = rk )
    sa = sqrt ( sa / real ( na, kind = rk ) - ea_ave * ea_ave )
  else
    ea_ave = 0.0D+00
  end if

  pa = real ( na * 100, kind = rk ) / real ( ntot, kind = rk )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' ) &
    'Absorbed   ', na, ea, pa, ea_ave, sa

  if ( 0 < nr ) then
    er_ave = er / real ( nr, kind = rk )
    sr = sqrt ( sr / real ( nr, kind = rk ) - er_ave * er_ave )
  else
    er_ave = 0.0D+00
  end if

  pr = real ( nr * 100, kind = rk ) / real ( ntot, kind = rk )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Reflected  ', nr, er, pr, er_ave, sr

  if ( 0 < nt ) then
    et_ave = et / real ( nt, kind = rk )
    st = sqrt ( st / real ( nt, kind = rk ) - et_ave * et_ave )
  else
    et_ave = 0.0D+00
  end if

  pt = real ( nt * 100, kind = rk ) / real ( ntot, kind = rk )

  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Transmitted', nt, et, pt, et_ave, st

  ptot = 100.0D+00

  write ( *, '(a)' ) ' '
  write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  &
    'Total      ', ntot, etot, ptot

  return
end
subroutine scatter ( e, mu, azm )

!*****************************************************************************80
!
!! SCATTER returns the new direction and energy of a particle that is scattered.
!
!  Discussion:
!
!    The scattering direction is chosen uniformly on the sphere.
!
!    The energy of the scattered particle is chosen uniformly in
!    [ 0.3*E, E ].
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input/output, real ( kind = rk ) E.  On input, the particle energy
!    before collision.  On output, the particle energy after collision
!    and scattering.
!
!    Output, real ( kind = rk ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Output, real ( kind = rk ) AZM, the azimuthal angle of the particle's
!    direction.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) azm
  real ( kind = rk ) e
  real ( kind = rk ) mu
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) u

  call random_number ( harvest = u )
  mu = - 1.0D+00 + 2.0D+00 * u

  call random_number ( harvest = u )
  azm = u * 2.0D+00 * pi

  call random_number ( harvest = u )
  e = ( u * 0.7D+00 + 0.3D+00 ) * e

  return
end
subroutine source ( e, mu, azm, x, y, z )

!*****************************************************************************80
!
!! SOURCE generates a new particle from the neutron source.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Output, real ( kind = rk ) E, the initial energy of the particle.
!
!    Output, real ( kind = rk ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Output, real ( kind = rk ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Output, real ( kind = rk ) X, Y, Z, the initial coordinates of the particle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) azm
  real ( kind = rk ) e
  real ( kind = rk ) energy
  real ( kind = rk ) mu
  real ( kind = rk ), parameter :: pi = 3.141592653589793D+00
  real ( kind = rk ) u
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  call random_number ( harvest = u )
  mu = u

  call random_number ( harvest = u )
  azm = u * 2.0D+00 * pi

  x = 0.0D+00
  y = 0.0D+00
  z = 0.0D+00

  e = energy ( )

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
  implicit none

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
subroutine update ( mu, azm, d, x, y, z )

!*****************************************************************************80
!
!! UPDATE determines the position of the particle after it has traveled D units.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    02 April 2008
!
!  Author:
!
!    Original FORTRAN77 version by Kahaner, Moler, Nash.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Kahaner, Cleve Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1989,
!    ISBN: 0-13-627258-4,
!    LC: TA345.K34.
!
!  Parameters:
!
!    Input, real ( kind = rk ) MU, the cosine of the angle between the
!    particle's direction and the X axis.
!
!    Input, real ( kind = rk ) AZM, the azimuthal angle of the particle's
!    direction.
!
!    Input, real ( kind = rk ) D, the distance the particle traveled.
!
!    Input/output, real ( kind = rk ) X, Y, Z.  On input, the previous
!    coordinates of the particle.  On output, the updated coordinates of the
!    particle.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) azm
  real ( kind = rk ) d
  real ( kind = rk ) mu
  real ( kind = rk ) s
  real ( kind = rk ) x
  real ( kind = rk ) y
  real ( kind = rk ) z

  s = sqrt ( 1.0D+00 - mu * mu )

  x = x + d * mu
  y = y + d * s * cos ( azm )
  z = z + d * s * sin ( azm )

  return
end
