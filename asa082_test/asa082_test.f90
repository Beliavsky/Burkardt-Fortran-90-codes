program main

!*****************************************************************************80
!
!! asa082_test() tests asa082().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa082_test():'
  write ( *, '(a)' ) '  FORTRAN90 version'
  write ( *, '(a)' ) '  asa082() computes the determinant of '
  write ( *, '(a)' ) '  an orthogonal matrix.'

  call detq_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'asa082_test():'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine detq_test ( )

!*****************************************************************************80
!
!! detq_test() tests detq().
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ), allocatable :: a(:,:)
  real ( kind = rk ) d1
  real ( kind = rk ) d2
  integer ifault
  integer n

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'detq_test():'
  write ( *, '(a)' ) '  detq() finds the determinant of an orthogonal matrix.'

  do n = 5, 10
    allocate ( a(1:n,1:n) )
    call helmert_matrix ( n, a )
    write ( *, '(a)' ) ''
    write ( *, '(a,i4)' ) '  Helmert matrix of order ', n
    if ( .false. ) then
      call r8mat_print ( n, n, a, '  Helmert matrix:' )
    end if
    call helmert_determinant ( n, d1 )
    write ( *, '(a,g14.6)' ) '  determinant =      ', d1
    call detq ( a, n, d2, ifault )
    if ( ifault == 1 ) then
      write ( *, '(a)' ) '  DETQ failed for this case.'
    else
      write ( *, '(a,g14.6)' ) '  DETQ determinant = ', d2
    end if
    deallocate ( a )
  end do

  return
end
subroutine helmert_matrix ( n, a )

!*****************************************************************************80
!
!! helmert_matrix() returns the Helmert matrix.
!
!  Discussion:
!
!    A matrix is a (standard) Helmert matrix if it is orthogonal,
!    and the elements which are above the diagonal and below the
!    first row are zero.
!
!    If the elements of the first row are all strictly positive,
!    then the matrix is a strictly Helmertian matrix.
!
!    It is possible to require in addition that all elements below
!    the diagonal be strictly positive, but in the reference, this
!    condition is discarded as being cumbersome and not useful.
!
!    A Helmert matrix can be regarded as a change of basis matrix
!    between a pair of orthonormal coordinate bases.  The first row
!    gives the coordinates of the first new basis vector in the old
!    basis.  Each later row describes combinations of (an increasingly
!    extensive set of) old basis vectors that constitute the new
!    basis vectors.
!
!    Helmert matrices have important applications in statistics.
!
!  Formula:
!
!    If I = 1 then
!      A(I,J) = 1 / sqrt ( N )
!    else if J < I then
!      A(I,J) = 1 / sqrt ( I * ( I - 1 ) )
!    else if J = I then
!      A(I,J) = - sqrt (I-1) / sqrt ( I )
!    else
!      A(I,J) = 0
!
!  Example:
!
!    N = 5
!
!    0.4472    0.4472    0.4472    0.4472    0.4472
!    0.7071   -0.7071         0         0         0
!    0.4082    0.4082   -0.8165         0         0
!    0.2887    0.2887    0.2887   -0.8660         0
!    0.2236    0.2236    0.2236    0.2236   -0.8944
!
!  Properties:
!
!    A is generally not symmetric: A' /= A.
!
!    A is orthogonal: A' * A = A * A' = I.
!
!    Because A is orthogonal, it is normal: A' * A = A * A'.
!
!    A is not symmetric: A' /= A.
!
!    det ( A ) = (-1)^(N+1)
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    HO Lancaster,
!    The Helmert Matrices,
!    American Mathematical Monthly,
!    Volume 72, 1965, pages 4-12.
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) A(N,N), the matrix.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer n

  real ( kind = rk ) a(n,n)
  integer i
  integer j
!
!  A begins with the first row, diagonal, and lower triangle set to 1.
!
  do i = 1, n
    do j = 1, n

      if ( i == 1 ) then
        a(i,j) = 1.0D+00 / sqrt ( real ( n, kind = rk ) )
      else if ( j < i ) then
        a(i,j) = 1.0D+00 / sqrt ( real ( i * ( i - 1 ), kind = rk ) )
      else if ( i == j ) then
        a(i,j) = - sqrt ( real ( i - 1, kind = rk ) ) &
                 / sqrt ( real ( i,     kind = rk ) )
      else
        a(i,j) = 0.0D+00
      end if

    end do
  end do

  return
end
subroutine helmert_determinant ( n, determ )

!*****************************************************************************80
!
!! helmert_determinant() returns the determinant of the Helmert matrix.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    27 August 2021
!
!  Author:
!
!    John Burkardt
!
!  Input:
!
!    integer N, the order of the matrix.
!
!  Output:
!
!    real ( kind = rk ) DETERM, the determinant.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) determ
  integer n

  if ( mod ( n, 2 ) == 0 ) then
    determ = -1.0D+00
  else
    determ = +1.0D+00
  end if

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! timestamp() prints the current YMDHMS date as a time stamp.
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
!    27 August 2021
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

