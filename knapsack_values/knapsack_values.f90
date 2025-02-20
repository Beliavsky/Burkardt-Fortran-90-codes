subroutine knapsack_values ( n_data, n, v, w, s, k )

!*****************************************************************************80
!
!! knapsack_values() returns samples of the knapsack problem.
!
!  Discussion:
!
!    For each value of n_data, starting with 0, the user makes two calls.
!    The first returns the appropriate value of n, so that the user 
!    can allocate arrays v, w, and s of the appropriate size.
!    The second call returns v, w, s and k.
!
!  Licensing:
!
!    This code is distributed under the MIT license. 
!
!  Modified:
!
!    19 November 2024
!
!  Author:
!
!    John Burkardt
!
!  Input on first call:
!
!    integer n_data: the user sets n_data to 0 before the first call.  
!
!  Output on first call:
!
!    integer n: the size of the arrays for this problem.
!    But if n is returned as 0, then there are no more data examples.
!
!  Input on second call:
!
!    integer n_data: the same value as on the previous call.  
!
!    integer n: the number of items
!
!  Output on second call:
!
!    integer n_data: on each call, the routine increments n_data by 1, and
!    returns the corresponding data when there is no more data, the
!    output value of n_data will be 0 again.
!
!    integer v(n): the value of each item
!
!    integer w(n): the weight of each item
!
!    integer s(n): 1 if the item is used in the optimal solution, 0 otherwise.
!
!    integer k: the weight capacity of the knapsack.
!
  implicit none

  integer n

  integer k
  integer n_data
  integer, parameter :: n_max = 10
  integer s(n)
  integer v(n)
  integer w(n)

  if ( n_data == n_max ) then
    n_data = 0
    n = 0
    k = 0
    return
  end if
!
!  If n == 0, this is the first call for this problem.
!  Return n, so user can allocate array sizes.
!
  if ( n == 0 ) then

    if ( n_data == 0 ) then
      n = 5
    else if ( n_data == 1 ) then
      n = 6
    else if ( n_data == 2 ) then
      n = 6
    else if ( n_data == 3 ) then
      n = 7
    else if ( n_data == 4 ) then
      n = 7
    else if ( n_data == 5 ) then
      n = 8
    else if ( n_data == 6 ) then
      n = 10
    else if ( n_data == 7 ) then
      n = 10
    else if ( n_data == 8 ) then
      n = 15
    else if ( n_data == 9 ) then
      n = 24
    end if

    return
  end if
!
!  On second call, set arrays and increment n_data.
!
  if ( n_data == 0 ) then

    v = (/ 24,  13,  23,  15,  16 /)
    w = (/ 12,   7,  11,   8,   9 /)
    s = (/  0,   1,   1,   1,   0 /)
    k = 26

  else if ( n_data == 1 ) then

    v = (/ 50,  50,  64,  46,  50,   5 /)
    w = (/ 56,  59,  80,  64,  75,  17 /)
    s = (/  1,   1,   0,   0,   1,   0 /)
    k = 190

  else if ( n_data == 2 ) then

    v = (/ 175, 90, 20, 50, 10, 200 /)
    w = (/ 10,   9,  4,  2,  1,  20 /)
    s = (/  1,   1,  0,  0,  1,   0 /)
    k = 20

  else if ( n_data == 3 ) then
 
    v = (/ 70,  20,  39,  37,   7,   5,  10 /)
    w = (/ 31,  10,  20,  19,   4,   3,   6 /)
    s = (/  1,   0,   0,   1,   0,   0,   0 /)
    k = 50

  else if ( n_data == 4 ) then
 
    v = (/442, 525, 511, 593, 546, 564, 617 /)
    w = (/ 41,  50,  49,  59,  55,  57,  60 /)
    s = (/  1,   0,   0,   1,   0,   0,   1 /)
    k = 170

  else if ( n_data == 5 ) then

    v = (/350, 400, 450,  20,  70,   8,   5,   5 /)
    w = (/ 25,  35,  45,   5,  25,   3,   2,   2 /)
    s = (/  1,   0,   1,   1,   1,   0,   1,   1 /)
    k = 104

  else if ( n_data == 6 ) then

    v = (/505, 352, 458, 220, 354, 414, 498, 545, 473, 543 /)
    w = (/ 23,  26,  20,  18,  32,  27,  29,  26,  30,  27 /)
    s = (/  1,   0,   0,   1,   0,   0,   0,   1,   0,   0 /)
    k = 67

  else if ( n_data == 7 ) then

    v = (/ 92,  57,  49,  68,  60,  43,  67,  84,  87,  72 /)
    w = (/ 23,  31,  29,  44,  53,  38,  63,  85,  89,  82 /)
    s = (/  1,   1,   1,   1,   0,   1,   0,   0,   0,   0 /)
    k = 165

  else if ( n_data == 8 ) then

    v = (/ 135, 139, 149, 150, 156, 163, 173, 184, 192, 201, &
           210, 214, 221, 229, 240 /)
    w = (/  70,  73,  77,  80,  82,  87,  90,  94,  98, 106, &
           110, 113, 115, 118, 120 /)
    s = (/   1,   0,   1,   0,   1,   0,   1,   1,   1,   0, &
             0,   0,   0,   1,   1 /)
    k = 750

  else if ( n_data == 9 ) then

    v = (/ 825594, 1677009, 1676628, 1523970,  943972, &
            97426,   69666, 1296457, 1679693, 1902996, &
          1844992, 1049289, 1252836, 1319836,  953277, &
          2067538,  675367,  853655, 1826027,   65731, &
           901489,  577243,  466257,  369261 /)
    w = (/ 382745,  799601,  909247,  729069,  467902, &
            44328,   34610,  698150,  823460,  903959, &
           853665,  551830,  610856,  670702,  488960, &
           951111,  323046,  446298,  931161,   31385, &
           496951,  264724,  224916,  169684 /)
    s = (/      1,       1,       0,       1,       1, &
                1,       0,       0,       0,       1, &
                1,       0,       1,       0,       0, &
                1,       0,       0,       0,       0, &
                0,       1,       1,       1 /)
    k = 6404180

  end if

  n_data = n_data + 1

  return
end
