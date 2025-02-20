subroutine subset_sum_backtrack ( s, n, v, more, u, t )
 
!*****************************************************************************80
!
!! subset_sum_backtrack() seeks, one at a time, subsets of V that sum to S.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    06 November 2022
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer S, the desired sum.
!
!    Input, integer N, the number of values.
!
!    Input, integer V(N), the values.
!    These must be nonnegative, and sorted in ascending order.  
!    Duplicate values are allowed.
!
!    Input/output, logical MORE, should be set to FALSE before the first call.
!    Thereafter, on output, MORE is TRUE if a solution is being returned,
!    and FALSE if there are no more solutions.
!
!    Input/output, integer U(N), should be set to 0 before the 
!    first call.  On output with MORE TRUE, U indexes the selected entries
!    of V that form a solution.
!
!    Input/output, integer T, should be set to 0 before the first 
!    call.  On output, if MORE is true, T is the highest index of the selected 
!    values, although this is of little interest to the user.
!
  implicit none

  integer n

  integer i
  logical more
  integer s
  integer su
  integer t
  integer told
  integer u(n)
  integer v(n)

  if ( .not. more ) then
  
    t = 0;
    u(1:n) = 0
    
  else
  
    more = .false.
    u(t) = 0

    told = t
    t = -1
    do i = told - 1, 1, -1
      if ( u(i) == 1 ) then
        t = i
        exit
      end if
    end do
    
    if ( t < 1 ) then
      return
    end if

    u(t) = 0
    t = t + 1
    u(t) = 1
      
  end if
    
  do

    su = dot_product ( u, v )
  
    if ( su < s .and. t < n ) then

      t = t + 1;
      u(t) = 1;        

    else if ( su == s ) then

      more = .true.
      return

    else

      u(t) = 0

      told = t
      t = -1
      do i = told - 1, 1, -1
        if ( u(i) == 1 ) then
          t = i
          exit
        end if
      end do
      
      if ( t < 1 ) then
        exit
      end if

      u(t) = 0
      t = t + 1
      u(t) = 1
      
    end if
        
  end do

  return
end

