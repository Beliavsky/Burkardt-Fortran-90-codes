program main

!*****************************************************************************80
!
!! MONOMIAL_TEST tests the MONOMIAL library.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONOMIAL_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the MONOMIAL library.'

  call mono_between_enum_test ( )
  call mono_between_next_grevlex_test ( )
  call mono_between_next_grlex_test ( )
  call mono_between_random_test ( )

  call mono_next_grevlex_test ( )
  call mono_next_grlex_test ( )
  call mono_print_test ( )
  call mono_rank_grlex_test ( )

  call mono_total_enum_test ( )
  call mono_total_next_grevlex_test ( )
  call mono_total_next_grlex_test ( )
  call mono_total_random_test ( )

  call mono_unrank_grlex_test ( )

  call mono_upto_enum_test ( )
  call mono_upto_next_grevlex_test ( )
  call mono_upto_next_grlex_test ( )
  call mono_upto_random_test ( )

  call mono_value_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONOMIAL_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ''
  call timestamp ( )

  return
end
subroutine mono_between_enum_test ( )

!*****************************************************************************80
!
!! MONO_BETWEEN_ENUM_TEST tests MONO_BETWEEN_ENUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer m
  integer mono_between_enum
  integer n1
  integer n2
  integer v

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_BETWEEN_ENUM_TEST'
  write ( *, '(a)' ) '  MONO_BETWEEN_ENUM can enumerate the number of monomials'
  write ( *, '(a)' ) '  in M variables, of total degree between N1 and N2.'

  m = 3
  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Using spatial dimension M = ', m
  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '   N2:'
  do n2 = 0, 8
    write ( *, '(2x,i4)', advance = 'no' ) n2
  end do
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '  N1 +------------------------------------------------------'
  do n1 = 0, 8
    write ( *, '(a,i2,a)', advance = 'no' ) '  ', n1, ' |'
    do n2 = 0, 8
      v = mono_between_enum ( m, n1, n2 )
      write ( *, '(2x,i4)', advance = 'no' ) v
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine mono_between_next_grevlex_test ( )

!*****************************************************************************80
!
!! MONO_BETWEEN_NEXT_GREVLEX_TEST tests MONO_BETWEEN_NEXT_GREVLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n1
  integer n2
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_BETWEEN_NEXT_GREVLEX_TEST'
  write ( *, '(a)' ) '  MONO_BETWEEN_NEXT_GREVLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree N between N1 and N2,'
  write ( *, '(a)' ) '  one at a time, in graded reverse lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,N1).'
  write ( *, '(a)' ) '  The process ends with (N2,0,...,0,0)'

  n1 = 2
  n2 = 3

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M =  ', m
  write ( *, '(a,i2)' ) '      N1 = ', n1
  write ( *, '(a,i2)' ) '      N2 = ', n2
  write ( *, '(a)' ) ''

  x = (/ 0, 0, n1 /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n2 ) then
      exit
    end if
 
    call mono_between_next_grevlex ( m, n1, n2, x )
    i = i + 1

  end do

  return
end
subroutine mono_between_next_grlex_test ( )

!*****************************************************************************80
!
!! MONO_BETWEEN_NEXT_GRELEX_TEST tests MONO_BETWEEN_NEXT_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n1
  integer n2
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_BETWEEN_NEXT_GRLEX_TEST'
  write ( *, '(a)' ) '  MONO_BETWEEN_NEXT_GRLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree N between N1 and N2,'
  write ( *, '(a)' ) '  one at a time, in graded lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,N1).'
  write ( *, '(a)' ) '  The process ends with (N2,0,...,0,0)'

  n1 = 2
  n2 = 3

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M =  ', m
  write ( *, '(a,i2)' ) '      N1 = ', n1
  write ( *, '(a,i2)' ) '      N2 = ', n2
  write ( *, '(a)' ) ''

  x = (/ 0, 0, n1 /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n2 ) then
      exit
    end if
 
    call mono_between_next_grlex ( m, n1, n2, x )
    i = i + 1

  end do

  return
end
subroutine mono_between_random_test ( )

!*****************************************************************************80
!
!! MONO_BETWEEN_RANDOM_TEST tests MONO_BETWEEN_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer , parameter :: m = 3

  integer n1
  integer n2
  integer rank
  integer test
  integer test_num
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_BETWEEN_RANDOM_TEST'
  write ( *, '(a)' ) '  MONO_BETWEEN_RANDOM selects at random a monomial'
  write ( *, '(a)' ) '  in M dimensions of total degree between N1 and N2.'

  n1 = 2
  n2 = 3

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Let M =  ', m
  write ( *, '(a,i3)' ) '      N1 = ', n1
  write ( *, '(a,i3)' ) '      N2 = ', n2
  write ( *, '(a)' ) ''

  test_num = 5

  do test = 1, test_num
    call mono_between_random ( m, n1, n2, rank, x )
    write ( *, '(2x,i3,4x,3i2)' ) rank, x(1:m)
  end do

  return
end

subroutine mono_next_grevlex_test ( )

!*****************************************************************************80
!
!! MONO_NEXT_GREVLEX_TEST tests MONO_NEXT_GREVLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4

  integer d
  integer k
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_NEXT_GREVLEX_TEST'
  write ( *, '(a)' ) '  MONO_NEXT_GREVLEX returns the next monomial'
  write ( *, '(a)' ) '  in graded reverse lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M =  ', m

  k = 0
  x(1:m) = 0

  do
    d = sum ( x(1:m) )
    write ( *, '(2x,i2,2x,i2,2x,a,2x,4i2)' ) k, d, '|', x(1:m)
    if ( x(1) == 3 ) then
      exit
    end if
    k = k + 1
    call mono_next_grevlex ( m, x )
  end do

  return
end
subroutine mono_next_grlex_test ( )

!*****************************************************************************80
!
!! MONO_NEXT_GRLEX_TEST tests MONO_NEXT_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    05 February 2015
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 4

  integer d
  integer k
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_NEXT_GRLEX_TEST'
  write ( *, '(a)' ) '  MONO_NEXT_GRLEX returns the next monomial'
  write ( *, '(a)' ) '  in graded lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M =  ', m

  k = 0
  x(1:m) = 0

  do
    d = sum ( x(1:m) )
    write ( *, '(2x,i2,2x,i2,2x,a,2x,4i2)' ) k, d, '|', x(1:m)
    if ( x(1) == 3 ) then
      exit
    end if
    k = k + 1
    call mono_next_grlex ( m, x )
  end do

  return
end
subroutine mono_print_test ( )

!*****************************************************************************80
!
!! MONO_PRINT_TEST tests MONO_PRINT.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 November 2014
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, dimension ( 1 ) :: f1 = (/ 5 /)
  integer, dimension ( 1 ) :: f2 = (/ - 5 /)
  integer, dimension ( 4 ) :: f3 = (/ 2, 1, 0, 3 /)
  integer, dimension ( 3 ) :: f4 = (/ 17, -3, 199 /)
  integer m

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_PRINT_TEST'
  write ( *, '(a)' ) '  MONO_PRINT can print out a monomial.'
  write ( *, '(a)' ) ''

  m = 1
  call mono_print ( m, f1, '  Monomial [5]:' )

  m = 1
  call mono_print ( m, f2, '  Monomial [5]:' )

  m = 4
  call mono_print ( m, f3, '  Monomial [2,1,0,3]:' )

  m = 3
  call mono_print ( m, f4, '  Monomial [17,-3,199]:' )

  return
end
subroutine mono_rank_grlex_test ( )

!*****************************************************************************80
!
!! MONO_RANK_GRLEX_TEST tests MONO_RANK_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3
  integer, parameter :: test_num = 8

  integer i
  integer n
  integer rank
  integer test
  integer x(m)
  integer, dimension ( m, test_num ) :: x_test = reshape ( (/ &
    0, 0, 0, &
    1, 0, 0, &
    0, 0, 1, &
    0, 2, 0, &
    1, 0, 2, &
    0, 3, 1, &
    3, 2, 1, &
    5, 2, 1 /), (/ m, test_num /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_RANK_GRLEX_TEST'
  write ( *, '(a)' ) '  MONO_RANK_GRLEX returns the rank of a monomial in '
  write ( *, '(a)' ) '  the sequence of all monomials in M dimensions.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Print a monomial sequence with ranks assigned.'

  n = 4

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M = ', m
  write ( *, '(a,i2)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x = (/ 0, 0, 0 /)
  i = 1

  do

    write ( *, '(2x,i3,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_upto_next_grlex ( m, n, x )
    i = i + 1

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  Now, given a monomial, retrieve its rank in the sequence:'
  write ( *, '(a)' ) ''

  do test = 1, test_num
    x(1:m) = x_test(1:m,test)
    call mono_rank_grlex ( m, x, rank )
    write ( *, '(2x,i3,4x,3i2)' ) rank, x(1:m)
  end do

  return
end
subroutine mono_total_enum_test ( )

!*****************************************************************************80
!
!! MONO_TOTAL_ENUM_TEST tests MONO_TOTAL_ENUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer m
  integer mono_total_enum
  integer n
  integer v

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_TOTAL_ENUM_TEST'
  write ( *, '(a)' ) '  MONO_TOTAL_ENUM can enumerate the number of monomials'
  write ( *, '(a)' ) '  in M variables, of total degree N.'

  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '    N:'
  do n = 0, 8
    write ( *, '(2x,i4)', advance = 'no' ) n
  end do
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '   M +------------------------------------------------------'
  do m = 1, 8
    write ( *, '(2x,i2,a)', advance = 'no' ) m, ' |'
    do n = 0, 8
      v = mono_total_enum ( m, n )
      write ( *, '(2x,i4)', advance = 'no' ) v
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine mono_total_next_grevlex_test ( )

!*****************************************************************************80
!
!! MONO_TOTAL_NEXT_GREVLEX_TEST tests MONO_TOTAL_NEXT_GREVLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_TOTAL_NEXT_GREVLEX_TEST'
  write ( *, '(a)' ) '  MONO_TOTAL_NEXT_GREVLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree N,'
  write ( *, '(a)' ) '  one at a time, in graded reverse lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,N).'
  write ( *, '(a)' ) '  The process ends with (N,0,...,0,0)'

  n = 3

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M = ', m
  write ( *, '(a,i2)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x = (/ 0, 0, n /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_total_next_grevlex ( m, n, x )
    i = i + 1

  end do

  return
end
subroutine mono_total_next_grlex_test ( )

!*****************************************************************************80
!
!! MONO_TOTAL_NEXT_GRLEX_TEST tests MONO_TOTAL_NEXT_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_TOTAL_NEXT_GRLEX_TEST'
  write ( *, '(a)' ) '  MONO_TOTAL_NEXT_GRLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree N,'
  write ( *, '(a)' ) '  one at a time, in graded lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,N).'
  write ( *, '(a)' ) '  The process ends with (N,0,...,0,0)'

  n = 3

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M = ', m
  write ( *, '(a,i2)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x = (/ 0, 0, n /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_total_next_grlex ( m, n, x )
    i = i + 1

  end do

  return
end
subroutine mono_total_random_test ( )

!*****************************************************************************80
!
!! MONO_TOTAL_RANDOM_TEST tests MONO_TOTAL_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer , parameter :: m = 3

  integer n
  integer rank
  integer test
  integer test_num
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_TOTAL_RANDOM_TEST'
  write ( *, '(a)' ) '  MONO_TOTAL_RANDOM selects at random a monomial'
  write ( *, '(a)' ) '  in M dimensions of total degree N.'

  n = 4

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Let M = ', m
  write ( *, '(a,i3)' ) '      N = ', n
  write ( *, '(a)' ) ''

  test_num = 5

  do test = 1, test_num
    call mono_total_random ( m, n, rank, x )
    write ( *, '(2x,i3,4x,3i2)' ) rank, x(1:m)
  end do

  return
end
subroutine mono_unrank_grlex_test ( )

!*****************************************************************************80
!
!! MONO_RANK_GRLEX_TEST tests MONO_UNRANK_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer , parameter :: m = 3

  integer i
  integer i4_uniform_ab
  integer mono_upto_enum
  integer n
  integer rank
  integer rank_max
  integer test
  integer test_num
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_UNRANK_GRLEXT_TEST'
  write ( *, '(a)' ) '  MONO_UNRANK_GRLEX is given a rank, and returns'
  write ( *, '(a)' ) '  the corresponding monomial in the sequence of '
  write ( *, '(a)' ) '  all monomials in M dimensions in grlex order.'

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  For reference, print a monomial sequence with ranks.'

  n = 4
  rank_max = mono_upto_enum ( m, n )

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Let M = ', m
  write ( *, '(a,i3)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x(1:m) = 0

  i = 1

  do

    write ( *, '(2x,i3,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_upto_next_grlex ( m, n, x )
    i = i + 1

  end do

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Now choose random ranks between 1 and ', rank_max
  write ( *, '(a)' ) ''

  test_num = 5

  do test = 1, test_num

    rank = i4_uniform_ab ( 1, rank_max )    
    call mono_unrank_grlex ( m, rank, x )
    write ( *, '(2x,i3,4x,3i2)' ) rank, x(1:m)

  end do

  return
end
subroutine mono_upto_enum_test ( )

!*****************************************************************************80
!
!! MONO_UPTO_ENUM_TEST tests MONO_UPTO_ENUM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    18 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer m
  integer mono_upto_enum
  integer n
  integer v

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_UPTO_ENUM_TEST'
  write ( *, '(a)' ) '  MONO_UPTO_ENUM can enumerate the number of monomials'
  write ( *, '(a)' ) '  in M variables, of total degree 0 up to N.'

  write ( *, '(a)' ) ''
  write ( *, '(a)', advance = 'no' ) '    N:'
  do n = 0, 8
    write ( *, '(2x,i4)', advance = 'no' ) n
  end do
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) &
    '   M +------------------------------------------------------'
  do m = 1, 8
    write ( *, '(2x,i2,a)', advance = 'no' ) m, ' |'
    do n = 0, 8
      v = mono_upto_enum ( m, n )
      write ( *, '(1x,i5)', advance = 'no' ) v
    end do
    write ( *, '(a)' ) ''
  end do

  return
end
subroutine mono_upto_next_grevlex_test ( )

!*****************************************************************************80
!
!! MONO_UPTO_NEXT_GREVLEX_TEST tests MONO_UPTO_NEXT_GREVLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_UPTO_NEXT_GREVLEX_TEST'
  write ( *, '(a)' ) '  MONO_UPTO_NEXT_GREVLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree up to N,'
  write ( *, '(a)' ) '  one at a time, in graded reverse lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,0).'
  write ( *, '(a)' ) '  The process ends with (N,0,...,0,0)'

  n = 4

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M = ', m
  write ( *, '(a,i2)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x = (/ 0, 0, 0 /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_upto_next_grevlex ( m, n, x )
    i = i + 1

  end do

  return
end
subroutine mono_upto_next_grlex_test ( )

!*****************************************************************************80
!
!! MONO_UPTO_NEXT_GRLEX_TEST tests MONO_UPTO_NEXT_GRLEX.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    04 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: m = 3

  integer i
  integer n
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_UPTO_NEXT_GRLEX_TEST'
  write ( *, '(a)' ) '  MONO_UPTO_NEXT_GRLEX can list the monomials'
  write ( *, '(a)' ) '  in M variables, of total degree up to N,'
  write ( *, '(a)' ) '  one at a time, in graded lexicographic order.'
  write ( *, '(a)' ) ''
  write ( *, '(a)' ) '  We start the process with (0,0,...,0,0).'
  write ( *, '(a)' ) '  The process ends with (N,0,...,0,0)'

  n = 4

  write ( *, '(a)' ) ''
  write ( *, '(a,i2)' ) '  Let M = ', m
  write ( *, '(a,i2)' ) '      N = ', n
  write ( *, '(a)' ) ''

  x = (/ 0, 0, 0 /)
  i = 1

  do

    write ( *, '(2x,i2,4x,3i2)' ) i, x(1:m)

    if ( x(1) == n ) then
      exit
    end if

    call mono_upto_next_grlex ( m, n, x )
    i = i + 1

  end do

  return
end
subroutine mono_upto_random_test ( )

!*****************************************************************************80
!
!! MONO_UPTO_RANDOM_TEST tests MONO_UPTO_RANDOM.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    21 November 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer , parameter :: m = 3

  integer n
  integer rank
  integer test
  integer test_num
  integer x(m)

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_UPTO_RANDOM_TEST'
  write ( *, '(a)' ) '  MONO_UPTO_RANDOM selects at random a monomial'
  write ( *, '(a)' ) '  in M dimensions of total degree no greater than N.'

  n = 4

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Let M = ', m
  write ( *, '(a,i3)' ) '      N = ', n
  write ( *, '(a)' ) ''

  test_num = 5

  do test = 1, test_num
    call mono_upto_random ( m, n, rank, x )
    write ( *, '(2x,i3,4x,3i2)' ) rank, x(1:m)
  end do

  return
end
subroutine mono_value_test ( )

!*****************************************************************************80
!
!! MONO_VALUE_TEST tests MONO_VALUE.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    10 December 2013
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: m = 3
  integer, parameter :: nx = 2

  integer f(m)
  integer j
  integer n
  integer rank
  integer test
  integer test_num
  real ( kind = rk ) v(nx)
  real ( kind = rk ), dimension ( m, nx ) :: x = reshape ( (/ &
     1.0D+00, 2.0D+00, 3.0D+00, &
    -2.0D+00, 4.0D+00, 1.0D+00 /), (/ m, nx /) )

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'MONO_VALUE_TEST'
  write ( *, '(a)' ) '  MONO_VALUE evaluates a monomial.'

  n = 6

  write ( *, '(a)' ) ''
  write ( *, '(a,i3)' ) '  Let M = ', m
  write ( *, '(a,i3)' ) '      N = ', n

  test_num = 5

  do test = 1, test_num

    call mono_upto_random ( m, n, rank, f )
    write ( *, '(a)' ) ''
    call mono_print ( m, f, '  M(X) = ' )
    call mono_value ( m, nx, f, x, v )
    do j = 1, nx
      write ( *, '(a,f4.0,a,f4.0,a,f4.0,a,g14.6)' ) &
        '  M(', x(1,j), ',', x(2,j), ',', x(3,j), ') = ', v(j)
    end do
    
  end do

  return
end

