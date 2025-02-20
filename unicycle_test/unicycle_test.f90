program main

!*****************************************************************************80
!
!! unicycle_test() tests unicycle().
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'unicycle_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  Test unicycle().'

  call perm_is_unicycle_test ( )
  call unicycle_enum_test ( )
  call unicycle_index_test ( )
  call unicycle_index_to_sequence_test ( )
  call unicycle_inverse_test ( )
  call unicycle_next_test ( )
  call unicycle_random_test ( )
  call unicycle_rank_test ( )
  call unicycle_unrank_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'unicycle_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine perm_is_unicycle_test ( )

!*****************************************************************************80
!
!! PERM_IS_UNICYCLE_TEST tests PERM_IS_UNICYCLE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5
  integer, parameter :: test_num = 10

  integer p(n)
  logical perm_is_unicycle
  integer test
  integer u(n)
  logical value

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'PERM_IS_UNICYCLE_TEST'
  write ( *, '(a)' ) '  PERM_IS_UNICYCLE determines whether a permutation'
  write ( *, '(a)' ) '  is a unicyle'

  do test = 1, test_num

    call perm_random ( n, p )

    value = perm_is_unicycle ( n, p )

    if ( value ) then

      call perm_print ( n, p, '  This permutation is a unicycle' )
      call unicycle_index_to_sequence ( n, p, u )
      call unicycle_print ( n, u, '  The permutation in sequence form' )

    else

      call perm_print ( n, p, '  This permutation is NOT a unicycle' )

    end if

  end do

  return
end
subroutine unicycle_enum_test ( )

!*****************************************************************************80
!
!! UNICYCLE_ENUM_TEST tests UNICYCLE_ENUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n_max = 10

  integer n
  integer num

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_ENUM_TEST'
  write ( *, '(a)' ) '  UNICYCLE_ENUM enumerates the unicycles of N objects.'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  N    Number'
  write ( *, '(a)' ) ' '

  do n = 0, n_max

    call unicycle_enum ( n, num )
    write ( *, '(2x,i3,2x,i8)' ) n, num

  end do

  return
end
subroutine unicycle_index_test ( )

!*****************************************************************************80
!
!! UNICYCLE_INDEX_TEST tests UNICYCLE_INDEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 6

  integer u(n)
  integer u_index(n)
  integer u2(n)
  integer test
  integer :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_INDEX_TEST'
  write ( *, '(a)' ) '  UNICYCLE_INDEX converts a unicycle to index form.'

  do test = 1, test_num 

    call unicycle_random ( n, u )

    call unicycle_print ( n, u, '  The unicycle:' )

    call unicycle_index ( n, u, u_index )
    
    call unicycle_index_print ( n, u_index, '  The index form:' )

    call unicycle_index_to_sequence ( n, u_index, u2 )

    call unicycle_print ( n, u2, '  The unicycle recovered:' )

  end do

  return
end
subroutine unicycle_index_to_sequence_test ( )

!*****************************************************************************80
!
!! UNICYCLE_INDEX_TO_SEQUENCE_TEST tests UNICYCLE_INDEX_TO_SEQUENCE.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 6

  integer u(n)
  integer u_index(n)
  integer u2(n)
  integer test
  integer :: test_num = 5

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_INDEX_TO_SEQUENCE_TEST'
  write ( *, '(a)' ) '  UNICYCLE_INDEX_TO_SEQUENCE converts an index to unicycle form.'

  do test = 1, test_num 

    call unicycle_random ( n, u )

    call unicycle_print ( n, u, '  The unicycle:' )

    call unicycle_index ( n, u, u_index )
    
    call unicycle_index_print ( n, u_index, '  The index form:' )

    call unicycle_index_to_sequence ( n, u_index, u2 )

    call unicycle_print ( n, u2, '  The unicycle recovered:' )

  end do

  return
end
subroutine unicycle_inverse_test ( )

!*****************************************************************************80
!
!! UNICYCLE_INVERSE_TEST tests UNICYCLE_INVERSE;
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 January 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 7

  integer, dimension ( n ) :: u = (/ 1, 7, 6, 2, 4, 3, 5 /)
  integer u_inverse(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_INVERSE_TEST'
  write ( *, '(a)' ) '  UNICYCLE_INVERSE inverts a unicycle;'

  call unicycle_print ( n, u, '  The original unicycle:' )
 
  call unicycle_inverse ( n, u, u_inverse )
 
  call unicycle_print ( n, u_inverse, '  The inverse unicycle:' )
 
  return
end
subroutine unicycle_next_test ( )

!*****************************************************************************80
!
!! UNICYCLE_NEXT_TEST tests UNICYCLE_NEXT.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    17 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  integer rank
  integer u(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_NEXT_TEST'
  write ( *, '(a)' ) '  UNICYCLE_NEXT generates unicycles in lex order.'
  write ( *, '(a)' ) ' '
  rank = -1
 
  do

    call unicycle_next ( n, u, rank )

    if ( rank == - 1 ) then
      exit
    end if

    write ( *, '(2x,i3,a1,2x,10i2)' ) rank, ':', u(1:n)

  end do
 
  return
end
subroutine unicycle_random_test ( )

!*****************************************************************************80
!
!! UNICYCLE_RANDOM_TEST tests UNICYCLE_RANDOM;
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  integer i
  integer u(n)

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'unicycle_random_test():'
  write ( *, '(a)' ) '  unicycle_random() produces a random unicyle;'
  write ( *, '(a,i8)' ) '  For this test, N = ', n
  write ( *, '(a)' ) ' '

  do i = 1, 5
    call unicycle_random ( n, u )
    call unicycle_print ( n, u, ' ' )
  end do
 
  return
end
subroutine unicycle_rank_test ( )

!*****************************************************************************80
!
!! UNICYCLE_RANK_TEST tests UNICYCLE_RANK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  integer, save, dimension ( n ) :: u = (/ 1, 5, 2, 3, 4 /)
  integer rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_RANK_TEST'
  write ( *, '(a)' ) '  UNICYCLE_RANK ranks a unicycle.'

  call unicycle_print ( n, u, '  The unicycle:' )
 
  call unicycle_rank ( n, u, rank )
 
  write ( *, '(a)' ) ' '
  write ( *, '(a,i8)' ) '  The rank is:', rank
 
  return
end
subroutine unicycle_unrank_test ( )

!*****************************************************************************80
!
!! UNICYCLE_UNRANK_TEST tests UNICYCLE_UNRANK.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    13 June 2012
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 5

  integer u(n)
  integer rank

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'UNICYCLE_UNRANK_TEST'
  write ( *, '(a)' ) '  UNICYCLE_UNRANK, given a rank, computes the'
  write ( *, '(a)' ) '  corresponding unicycle.'
  write ( *, '(a)' ) ' '
  rank = 6
  write ( *, '(a,i8)' ) '  The requested rank is ', rank
 
  call unicycle_unrank ( n, rank, u )
 
  call unicycle_print ( n, u, '  The unicycle:' )
 
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

