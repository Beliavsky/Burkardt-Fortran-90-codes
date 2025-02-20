program main

!*****************************************************************************80
!
!! apportionment_test() tests apportionment().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'apportionment_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test apportionment().'

  call test01 ( )
  call test02 ( )
  call test03 ( )
  call test04 ( )
  call test05 ( )
  call test06 ( )
  call test07 ( )
  call test08 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'apportionment_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end
subroutine test01 ( )

!*****************************************************************************80
!
!! test01() tests state_num_year() and rep_num_year().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer rep_num
  integer state_num
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test01()'
  write ( *, '(a)' ) '  state_num_year() returns the number of states in'
  write ( *, '(a)' ) '  the union at the end of a given year.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  rep_num_year() returns the number of reps in'
  write ( *, '(a)' ) '  the House of Representatives (based only on the'
  write ( *, '(a)' ) '  decennial census.)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Year  States  reps'
  write ( *, '(a)' ) ' '

  do year = 1790, 2010, 10
    call state_num_year ( year, state_num )
    call rep_num_year ( year, rep_num )
    write ( *, '(2x,i4,6x,i2,4x,i3)' ) year, state_num, rep_num
  end do

  return
end
subroutine test02 ( )

!*****************************************************************************80
!
!! test02() tests state_statehood().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer d
  integer m
  character ( len = 9 ) month
  integer state
  character ( len = 20 ) state_name
  integer y

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test02()'
  write ( *, '(a)' ) '  state_statehood() returns the statehood date of a state.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '   #  Name                     Statehood Date'
  write ( *, '(a)' ) ' '

  do state = 1, 51
    call state_statehood ( state, y, m, d )
    call i4_to_month_name ( m, month )
    write ( *, '(2x,i2,2x,a20,2x,i2,2x,a9,2x,i4)' ) &
      state, state_name(state), d, month, y
  end do

  return
end
subroutine test03 ( )

!*****************************************************************************80
!
!! test03() gets the historic representation values.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: state_num = 51

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  integer rep_num
  integer state
  character ( len = 2 ) state_id
  integer state_pop(state_num)
  integer state_rep(state_num)
  character ( len = 12 ) string
  integer us_pop
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test03():'
  write ( *, '(a)' ) '  Get the historic representation values.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  What were the state representations based on the last decennial census?
!
  call state_rep_year ( year, state_rep )

  call rep_num_year ( year, rep_num )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%'
  write ( *, '(a)' ) ' '

  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4)' ) &
      state_id (state), string, state_rep(state), f1, f2

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------'

  call i4_to_s_right_comma ( us_pop, string )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4)' ) &
    'US', string, rep_num, 100.0, 100.0

  return
end
subroutine test04 ( )

!*****************************************************************************80
!
!! test04() tests apportion_hamilton().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  integer indx(51)
  integer rep_num
  integer state
  character ( len = 2 ) state_id
  integer state_num
  integer state_pop(51)
  integer, allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer us_pop
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test04()'
  write ( *, '(a)' ) '  apportion_hamilton() uses Hamilton''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_hamilton ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )
    f3 = real ( state_pop(state), kind = rk ) / real ( state_rep(state), kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = rk ) / real ( rep_num, kind = rk )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test05 ( )

!*****************************************************************************80
!
!! test05() tests apportion_jefferson().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  integer indx(51)
  integer rep_num
  integer state
  character ( len = 2 ) state_id
  integer state_num
  integer state_pop(51)
  integer, allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer us_pop
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test05()'
  write ( *, '(a)' ) '  apportion_jefferson() uses Jefferson''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_jefferson ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )
    f3 = real ( state_pop(state), kind = rk ) / real ( state_rep(state), kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = rk ) / real ( rep_num, kind = rk )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test06 ( )

!*****************************************************************************80
!
!! test06() tests apportion_adams().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  integer indx(51)
  integer rep_num
  integer state
  character ( len = 2 ) state_id
  integer state_num
  integer state_pop(51)
  integer, allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer us_pop
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test06()'
  write ( *, '(a)' ) '  apportion_adams() uses Adams''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_adams ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )
    f3 = real ( state_pop(state), kind = rk ) / real ( state_rep(state), kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = rk ) / real ( rep_num, kind = rk )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test07 ( )

!*****************************************************************************80
!
!! test07() tests apportion_webster().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  integer indx(51)
  integer rep_num
  integer state
  character ( len = 2 ) state_id
  integer state_num
  integer state_pop(51)
  integer, allocatable, dimension ( : ) :: state_rep
  character ( len = 12 ) string
  integer us_pop
  integer year

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test07():'
  write ( *, '(a)' ) '  apportion_webster() uses Webster''s method to'
  write ( *, '(a)' ) '  apportion representatives.'
!
!  Pick a year.
!
  year = 1960
!
!  What were the state populations in the last decennial census?
!
  call state_pop_year ( year, state_pop )
!
!  Make an index vector.
!
  call i4vec_nonzero_first ( 51, state_pop, state_num, indx )
!
!  "Squeeze" the population vector.
!
  state_pop(1:state_num) = state_pop(indx(1:state_num))

  allocate ( state_rep(1:state_num) )

  call rep_num_year ( year, rep_num )

  call apportion_webster ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )
    f3 = real ( state_pop(state), kind = rk ) / real ( state_rep(state), kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id ( indx(state) ), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = rk ) / real ( rep_num, kind = rk )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  deallocate ( state_rep )

  return
end
subroutine test08 ( )

!*****************************************************************************80
!
!! test08() tests apportion_huntington_hill().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    26 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: state_num = 15

  real ( kind = rk ) f1
  real ( kind = rk ) f2
  real ( kind = rk ) f3
  integer, parameter :: rep_num = 105
  integer state
  character ( len = 2 ), dimension ( state_num ) :: state_id = (/ &
    'CT', 'DE', 'GA', 'KY', 'MD', &
    'MA', 'NH', 'NJ', 'NY', 'NC', &
    'PA', 'RI', 'SC', 'VT', 'VA' /) 
  integer, dimension ( state_num ) :: state_pop = (/ &
    236841,  55540,  70835,  68705, 278514, &
    475327, 141822, 179570, 331589, 353523, &
    432879,  68446, 206236,  85533, 630560 /)
  integer state_rep(state_num)
  character ( len = 12 ) string
  integer us_pop
  integer, parameter :: year = 1790

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'test08()'
  write ( *, '(a)' ) '  apportion_huntingon_hill() uses the Huntington-Hill'
  write ( *, '(a)' ) '  apportionment method.'

  call apportion_huntington_hill ( state_num, state_pop, rep_num, state_rep )

  write ( *, '(a)' ) ' '
  write ( *, '(a,i4)' ) '  Year: ', year
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  ID    Population   Rep     Pop%      Rep%      Pop/Rep'
  write ( *, '(a)' ) ' '
  us_pop = sum ( state_pop(1:state_num) )

  do state = 1, state_num

    f1 = real ( 100.0 * state_pop(state), kind = rk ) / real ( us_pop, kind = rk )
    f2 = real ( 100.0 * state_rep(state), kind = rk ) / real ( rep_num, kind = rk )
    f3 = real ( state_pop(state), kind = rk ) / real ( state_rep(state), kind = rk )

    call i4_to_s_right_comma ( state_pop(state), string )

    write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
      state_id(state), string, state_rep(state), f1, f2, f3

  end do

  write ( *, '(a)' ) '  --  ------------  ---  --------  --------  ------------'

  call i4_to_s_right_comma ( us_pop, string )

  f3 = real ( us_pop, kind = rk ) / real ( rep_num, kind = rk )

  write ( *, '(2x,a2,2x,a12,2x,i3,2x,f8.4,2x,f8.4,2x,f12.0)' ) &
    'US', string, rep_num, 100.0, 100.0, f3

  return
end
 
