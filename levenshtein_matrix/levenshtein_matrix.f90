subroutine levenshtein_matrix ( m, s, n, t, d )

!*****************************************************************************80
!
!! levenshtein_matrix() computes the Levenshtein distance matrix between strings.
!
!  Discussion:
!
!    Let S and T be source and target strings.  Consider the task of
!    converting S to T in the minimal number of steps, involving
!    * Insertion: insert a new character
!    * Deletion: delete a character
!    * Substitution: swap one character for another.
!    The Levenshtein distance from S to T is the minimal number of such
!    steps required to transform S into T.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 March 2018
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer M, the length of string S.
!
!    character ( len = * ) S, the first string.
!
!    integer N, the length of string T.
!
!    character ( len = * ) T, the second string.
!
!  Output:
!
!    integer D(0:M,0:N), the Levenshtein matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer m
  integer n

  integer d(0:m,0:n)
  integer i
  integer j
  character ( len = * ) s
  integer substitution_cost
  character ( len = * ) t

  d(0,0) = 0
!
!  Source prefixes can be transformed into empty string by
!  dropping all characters,
!
  do i = 1, m
    d(i,0) = i
  end do
!
!  Target prefixes can be reached from empty source prefix
!  by inserting every character.
!
  do j = 1, n
    d(0,j) = j
  end do

  do j = 1, n
    do i = 1, m
      if ( s(i:i) == t(j:j) ) then
        substitution_cost = 0
      else
        substitution_cost = 1
      end if
      d(i,j) = min ( d(i-1,j) + 1, &
               min ( d(i,j-1) + 1, &
                     d(i-1,j-1) + substitution_cost ) )
    end do
  end do
 
  return
end

