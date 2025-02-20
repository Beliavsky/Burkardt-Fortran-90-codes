program main

!*****************************************************************************80
!
!! duel_simulation() simulates a duel.
!
!  Discussion:
!
!    Player 1 fires at player 2, and hits with a probability of P(1).
!    If Player 2 misses, then Player 2 fires at Player 1, hitting with
!    a probability of P(2).
!
!    The duel continues with alternating shots until only one player 
!    survives.
!
!    The simulation is intended to estimate the probabilities that a
!    player will survive, and the number of turns required.
!
!    The exact probability that player 1 will survive is
!
!      P(1) / ( P(1) + P(2) - P(1) * P(2) )
!
!    Player 2's chance is
!
!     P(2) * ( 1 - P(1) ) / ( P(1) + P(2) - P(1) * P(2) )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    17 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Nahin,
!    Duelling Idiots and Other Probability Puzzlers,
!    Princeton University Press, 2000,
!    ISBN13: 978-0691009797,
!    LC: QA273.N29.
!
!    Martin Shubik,
!    "Does the Fittest Necessarily Survive?",
!    in Readings in Game Theory and Political Behavior,
!    edited by Martin Shubik,
!    Doubleday, 1954,
!    LC: H61.S53.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_ACCURACY, B_ACCURACY, the probabilities that 
!    players A and B will hit their opponent in a single shot.
!
!    Input, integer DUEL_NUM, the number of duels to run.
!
!    Output, real ( kind = rk ) A_PROB, B_PROB, the estimated probablities that 
!    players A and B will survive.
!
!    Output, real ( kind = rk ) TURN_AVERAGE, the average number of turns 
!    required to complete the duel.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_accuracy
  real ( kind = rk ) a_prob
  integer a_wins
  real ( kind = rk ) b_accuracy
  real ( kind = rk ) b_prob
  integer b_wins
  integer duel
  integer duel_num
  integer winner

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'DUEL_SIMULATION:'
  write ( *, '(a)' ) '  FORTRAN90 version'

  write ( *, '(a)' ) 'Enter number of duels to run:'
  read ( *, * ) duel_num

  write ( *, '(a)' ) 'Enter player A''s accuracy between 0.0 and 1.0:'
  read ( *, * ) a_accuracy

  write ( *, '(a)' ) 'Enter player B''s accuracy between 0.0 and 1.0:'
  read ( *, * ) b_accuracy

  a_wins = 0
  b_wins = 0

  do duel = 1, duel_num

    call duel_result ( a_accuracy, b_accuracy, winner )

    if ( winner == 1 ) then
      a_wins = a_wins + 1
    else
      b_wins = b_wins + 1
    end if

  end do

  a_prob = real ( a_wins, kind = rk ) / real ( duel_num, kind = rk )
  b_prob = real ( b_wins, kind = rk ) / real ( duel_num, kind = rk )

  write ( *, '(a)' ) ''
  write ( *, '(a,g14.6)' ) '  Player A wins with probability ', a_prob
  write ( *, '(a,g14.6)' ) '  Player B wins with probability ', b_prob
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'duel_simulation()gg:'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  stop
end
subroutine duel_result ( a_accuracy, b_accuracy, winner )

!*****************************************************************************80
!
!! DUEL_RESULT returns the outcome of a single duel.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 September 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Martin Shubik,
!    “Does the Fittest Necessarily Survive?”,
!    in Readings in Game Theory and Political Behavior,
!    edited by Martin Shubik,
!    Doubleday, 1954,
!    LC: H61.S53.
!
!  Parameters:
!
!    Input, real ( kind = rk ) A_ACCURACY, B_ACCURACY, the probabilities that 
!    player A and B will hit their opponent in a single shot.
!
!    Output, integer WINNER, the survivor of the duel.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a_accuracy
  real ( kind = rk ) b_accuracy
  real ( kind = rk ) r
  integer winner

  do

    call random_number ( harvest = r )

    if ( r <= a_accuracy ) then
      winner = 1
      exit
    end if

    call random_number ( harvest = r )

    if ( r <= b_accuracy ) then
      winner = 2
      exit
    end if

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
