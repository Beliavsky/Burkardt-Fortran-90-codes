subroutine swap ( varval, class, clsize, in, ik, iv, critvl, ntrans, ifault )

!*****************************************************************************80
!
!! swap() interchanges objects between different classes to improve a criterion.
!
!  Discussion:
!
!    This routine is given a classification of objects, including the
!    number of objects in each class, and the current value of some criterion
!    which is desired to be minimized.
!
!    The routine calculates the change in criterion for all possible swaps,
!    that is, operations in which two objects in different classes exchange
!    places. Each swap that would result in a lowering of the criterion is
!    executed, and the related quantities are updated.
!
!    When no more advantageous swaps can be found, the routine returns.
!
!    The routine relies on a user-supplied routine, CRSWAP, to report the
!    expected change in the criterion for a given swap, and to carry
!    out that transfer if requested.
!
!    The variables CLASS and CRITVL have been added to the argument list
!    of CRSWAP.
!
!    Also, the order of the two classes "L" and "M" was interchanged in
!    the call to CRSWAP.  The original order was counterintuitive.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Banfield, Bassill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Input:
!
!    real ( kind = rk ) VARVAL(IN,IV), the data values.  There are
!    IN objects, each having spatial dimension IV.
!
!    integer CLASS(IN), the classification of each object.
!
!    integer CLSIZE(IK), the number of objects in each class.
!
!    integer IN, the number of objects.
!
!    integer IK, the number of classes.
!
!    integer IV, the number of spatial dimensions,
!    or variates, of the objects.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!  Output:
!
!    integer CLASS(IN), the updated classification of each object.
!
!    integer CLSIZE(IK), the updated number of objects in each class.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!    integer NTRANS, the number of transfers executed.
!
!    integer IFAULT, error indicator.
!    0, no error detected.
!    1, the number of classes was less than 2.
!    2, the number of objects was less than the number of classes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ik
  integer in
  integer iv

  integer class(in)
  integer clsize(ik)
  real ( kind = rk ) critvl
  real ( kind = rk ), parameter :: eps = 1.0D-38
  integer i
  integer icount
  integer ifault
  real ( kind = rk ) inc
  integer iswitch
  integer it
  integer itop
  integer j
  integer k
  integer l
  integer m
  integer ntrans
  real ( kind = rk ) varval(in,iv)

  if ( ik <= 1 ) then
    ifault = 1
    return
  end if

  if ( in <= ik ) then
    ifault = 2
    return
  end if

  ifault = 0
  icount = 0
  ntrans = 0
  itop = ( in * ( in - 1 ) ) / 2

  i = 1

  do

    i = i + 1

    if ( itop <= icount ) then
      exit
    end if

    if ( in < i ) then
      i = 1
      cycle
    end if

    l = class(i)
    k = l
    it = i - 1
!
!  Test the swap of object I from class M to L,
!  and object J from class L to M.
!
    do j = 1, it

      icount = icount + 1
      m = class(j)

      if ( l /= j ) then

        if ( clsize(l) /= 1 .or. clsize(m) /= 1 ) then

          iswitch = 1
          call crswap ( varval, class, clsize, in, ik, iv, critvl, &
            i, j, l, m, iswitch, inc )

          if ( inc < - eps ) then

            critvl = critvl + inc
            icount = 0

            iswitch = 2
            call crswap ( varval, class, clsize, in, ik, iv, critvl, &
              i, j, l, m, iswitch, inc )

            ntrans = ntrans + 1
            class(i) = m
            class(j) = l
            l = m

          end if

        end if
      end if

    end do

  end do

  return
end
subroutine trnsfr ( varval, class, clsize, in, ik, iv, critvl, ntrans, ifault )

!*****************************************************************************80
!
!! trnsfr() transfers objects between classes to improve a criterion.
!
!  Discussion:
!
!    This routine is given a classification of objects, including the
!    number of objects in each class, and the current value of some criterion
!    which is desired to be minimized.
!
!    The routine calculates the change in criterion for all possible transfers
!    of any object from its current class to a different class.  Each transfer
!    that would result in a lowering of the criterion is executed, and the
!    related quantities are updated.
!
!    When no more advantageous transfers can be found, the routine returns.
!
!    The routine relies on a user-supplied routine, CRTRAN, to report the
!    expected change in the criterion for a given transfer, and to carry
!    out that transfer if requested.
!
!    The variables CLASS and CRITVL have been added to the argument list
!    of CRTRAN.
!
!    Also, the order of the two classes "L" and "M" was interchanged in
!    the call to CRTRAN.  The original order was counterintuitive.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    28 August 2021
!
!  Author:
!
!    Original FORTRAN77 version by Banfield, Bassill.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Colin Banfield, LC Bassill,
!    Algorithm AS 113:
!    A transfer for non-hierarchichal classification,
!    Applied Statistics,
!    Volume 26, Number 2, 1977, pages 206-210.
!
!  Input:
!
!    real ( kind = rk ) VARVAL(IN,IV), the data values.  There are IN
!    objects, each having spatial dimension IV.
!
!    integer CLASS(IN), the classification of each object.
!
!    integer CLSIZE(IK), the number of objects in each class.
!
!    integer IN, the number of objects.
!
!    integer IK, the number of classes.
!
!    integer IV, the number of spatial dimensions, or
!    variates, of the objects.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!  Output:
!
!    integer CLASS(IN), the updated classification of each object.
!
!    integer CLSIZE(IK), the updated number of objects in each class.
!
!    real ( kind = rk ) CRITVL, the current value of the criterion.
!
!    integer NTRANS, the number of transfers executed.
!
!    integer IFAULT, error indicator.
!    0, no error detected.
!    1, the number of classes was less than 2.
!    2, the number of objects was less than the number of classes.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer ik
  integer in
  integer iv

  integer class(in)
  integer clsize(ik)
  real ( kind = rk ) critvl
  real ( kind = rk ), parameter :: eps = 1.0D-38
  integer i
  integer icount
  integer ifault
  real ( kind = rk ) inc
  real ( kind = rk ) inco
  integer iswitch
  integer l
  integer lo
  integer m
  integer ntrans
  real ( kind = rk ) varval(in,iv)

  if ( ik <= 1 ) then
    ifault = 1
    return
  end if

  if ( in <= ik ) then
    ifault = 2
    return
  end if

  ifault = 0
  ntrans = 0
  i = 0
  icount = 0

  do

    i = i + 1

    if ( in <= icount ) then
      exit
    end if

    if ( in < i ) then
      i = 0
      icount = 0
      cycle
    end if

    m = class(i)
    if ( clsize(m) <= 1 ) then
      icount = icount + 1
      cycle
    end if

    inco = - eps
    lo = m
!
!  Test the transfer of object I from class M to class L.
!
    do l = 1, ik

      if ( l /= m ) then

        iswitch = 1
        call crtran ( varval, class, clsize, in, ik, iv, critvl, &
          i, m, l, iswitch, inc )
!
!  Remember the values of L and INC.
!
        if ( inc < inco ) then
          lo = l
          inco = inc
        end if

      end if

    end do

    icount = icount + 1
!
!  Execute the transfer of object I from class M to class LO.
!
    if ( lo /= m ) then

      l = lo
      critvl = critvl + inco
      icount = 0

      iswitch = 2
      call crtran ( varval, class, clsize, in, ik, iv, critvl, &
        i, m, l, iswitch, inc )

      ntrans = ntrans + 1
      class(i) = l
      clsize(l) = clsize(l) + 1
      clsize(m) = clsize(m) - 1

    end if

  end do

  return
end
 
