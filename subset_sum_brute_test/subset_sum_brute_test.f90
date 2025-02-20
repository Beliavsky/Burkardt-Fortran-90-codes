  program main

!*****************************************************************************80
!
!! subset_sum_brute_test() tests subset_sum_brute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_brute_test():'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test subset_sum_brute().'

  call subset_sum_brute_test01 ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_brute_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine subset_sum_brute_test01 ( )

!*****************************************************************************80
!
!! subset_sum_brute_test01() tests subset_sum_brute().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    25 October 2022
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: n = 21

  integer choice(n)
  integer i
  integer target
  integer w_sum
  integer, save, dimension ( n ) :: weight = (/ &
      518533, 1037066, 2074132, 1648264, 796528, &
     1593056,  686112, 1372224,  244448, 488896, &
      977792, 1955584, 1411168,  322336, 644672, &
     1289344,   78688,  157376,  314752, 629504, &
     1259008 /)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'subset_sum_brute_test01():'
  write ( *, '(a)' ) '  subset_sum_brute() looks for a selection '
  write ( *, '(a)' ) '  from a set of weights that adds up to a '
  write ( *, '(a)' ) '  given target.'
!
!  Define the problem data.
!
  target = 2463098
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Target value:'
  write ( *, '(2x,i8)' ) target
 
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '   I  W(I)'
  write ( *, '(a)' ) ''
  do i = 1, n
    write ( *, '(2x,i2,2x,i8)' ) i, weight(i)
  end do

  call subset_sum_brute ( n, weight, target, choice )

  if ( choice(1) == -1 ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '  No solution was found.'
  else
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) '   I*     W*'
    write ( *, '(a)' ) ''
    w_sum = 0
    do i = 1, n
      if ( choice(i) == 1 ) then
        w_sum = w_sum + weight(i)
        write ( *, '(2x,i2,2x,i8)' ) i, weight(i)
      end if
    end do
    write ( *, '(a)' ) ''
    write ( *, '(a,i12)' ) '  Sum:    ', w_sum
   write ( *, '(a,i12)' ) '  Target: ', target
  end if

  return
end

