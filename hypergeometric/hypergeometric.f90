function pfq ( a, b, ip, iq, z, lnpfq, ix, nsigfig )

!*****************************************************************************80
!
!! pfq() evaluates generalized hypergeometric functions for complex arguments.
!
!  Discussion:
!
!    A numerical evaluator for the generalized hypergeometric
!    function for complex arguments with large magnitudes
!    using a direct summation of the Gauss series.
!
!    pFq is defined by (borrowed from Maple):
!      pFq = sum(z^k / k! * product(pochhammer(n[i], k), i=1..p) /
!            product(pochhammer(d[j], k), j=1..q), k=0..infinity )
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
!  Reference:
!
!    Warren Perger, Atul Bhalla, Mark Nardin,
!    ACPAPFQ, a numerical evaluator for the generalized hypergeometric series,
!    Computational Physics Communications,
!    Volume 77, page 249, 1993.
!
!  Input:
!
!    complex ( kind = ck ) a(ip): numerator parameters
!
!    complex ( kind = ck ) b(iq): denominator parameters
!
!    integer ip: number of numerator parameters
!
!    integer iq: number of denominator parameters
!
!    complex ( kind = ck ) z: the argument
!
!    integer lnpfq: set to 1 if desired result is the natural
!                    log of pfq (default is 0)
!
!    integer ix: maximum number of terms in a,b (see below)
!
!    integer nsigfig: number of desired significant figures (default=10)
!
!  Output:
!
!    complex ( kind = ck ) pfq: the value of pFq(z).
!
!  Examples:
!
!    a=[1+i,1]; b=[2-i,3,3]; z=1.5;
!               >> genHyper(a,b,z)
!               ans =
!                          1.02992154295955 +     0.106416425916656i
!               or with more precision,
!               >> genHyper(a,b,z,0,0,15)
!               ans =
!                          1.02992154295896 +     0.106416425915575i
!               using the log option,
!               >> genHyper(a,b,z,1,0,15)
!               ans =
!                        0.0347923403326305 +     0.102959427435454i
!               >> exp(ans)
!               ans =
!                          1.02992154295896 +     0.106416425915575i
!
!
!     ****************************************************************
!     *                                                              *
!     *    SOLUTION TO THE GENERALIZED HYPERGEOMETRIC FUNCTION       *
!     *                                                              *
!     *                           by                                 *
!     *                                                              *
!     *                      W. F. PERGER,                           *
!     *                                                              *
!     *              MARK NARDIN  and ATUL BHALLA                    *
!     *                                                              *
!     *                                                              *
!     *            Electrical Engineering Department                 *
!     *            Michigan Technological University                 *
!     *                  1400 Townsend Drive                         *
!     *                Houghton, MI  49931-1295   USA                *
!     *                     Copyright 1993                           *
!     *                                                              *
!     *               e-mail address: wfp@mtu.edu                    *
!     *                                                              *
!     *  Description : A numerical evaluator for the generalized     *
!     *    hypergeometric function for complex arguments with large  *
!     *    magnitudes using a direct summation of the Gauss series.  *
!     *    The method used allows an accuracy of up to thirteen      *
!     *    decimal places through the use of large integer arrays    *
!     *    and a single final division.                              *
!     *    (original subroutines for the confluent hypergeometric    *
!     *    written by Mark Nardin, 1989; modifications made to cal-  *
!     *    culate the generalized hypergeometric function were       *
!     *    written by W.F. Perger and A. Bhalla, June, 1990)         *
!     *                                                              *
!     *  The evaluation of the pFq series is accomplished by a func- *
!     *  ion call to PFQ, which is a real ( kind = rk ) complex func-  *
!     *  tion.  The required input is:                               *
!     *  1. real ( kind = rk ) complex arrays A and B.  These are the  *
!     *     arrays containing the parameters in the numerator and de-*
!     *     nominator, respectively.                                 *
!     *  2. Integers IP and IQ.  These integers indicate the number  *
!     *     of numerator and denominator terms, respectively (these  *
!     *     are p and q in the pFq function).                        *
!     *  3. real ( kind = rk ) complex argument Z.                     *
!     *  4. Integer LNPFQ.  This integer should be set to '1' if the *
!     *     result from PFQ is to be returned as the natural logaritm*
!     *     of the series, or '0' if not.  The user can generally set*
!     *     LNPFQ = '0' and change it if required.                   *
!     *  5. Integer IX.  This integer should be set to '0' if the    *
!     *     user desires the program PFQ to estimate the number of   *
!     *     array terms (in A and B) to be used, or an integer       *
!     *     greater than zero specifying the number of integer pos-  *
!     *     itions to be used.  This input parameter is escpecially  *
!     *     useful as a means to check the results of a given run.   *
!     *     Specificially, if the user obtains a result for a given  *
!     *     set of parameters, then changes IX and re-runs the eval- *
!     *     uator, and if the number of array positions was insuffi- *
!     *     cient, then the two results will likely differ.  The rec-*
!     *     commended would be to generally set IX = '0' and then set*
!     *     it to 100 or so for a second run.  Note that the LENGTH  *
!     *     parameter currently sets the upper limit on IX to 777,   *
!     *     but that can easily be changed (it is a single PARAMETER *
!     *     statement) and the program recompiled.                   *
!     *  6. Integer NSIGFIG.  This integer specifies the requested   *
!     *     number of significant figures in the final result.  If   *
!     *     the user attempts to request more than the number of bits*
!     *     in the mantissa allows, the program will abort with an   *
!     *     appropriate error message.  The recommended value is 10. *
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) a(ip)
  complex ( kind = ck ) a1(2)
  real ( kind = rk ) argi
  real ( kind = rk ) argr
  complex ( kind = ck ) b(iq)
  complex ( kind = ck ) b1(1)
  complex ( kind = ck ) cgamma
  real ( kind = rk ) diff
  real ( kind = rk ) dnum
  complex ( kind = ck ) gam1
  complex ( kind = ck ) gam2
  complex ( kind = ck ) gam3
  complex ( kind = ck ) gam4
  complex ( kind = ck ) gam5
  complex ( kind = ck ) gam6
  complex ( kind = ck ) gam7
  complex ( kind = ck ) hyper
  complex ( kind = ck ) hyper1
  complex ( kind = ck ) hyper2
  integer i
  integer ip
  integer iq
  integer ix
  integer lnpfq
  integer nsigfig
  complex ( kind = ck ) pfq
  real ( kind = rk ) precis
  complex ( kind = ck ) z
  complex ( kind = ck ) z1

  if ((lnpfq/=0) .and. (lnpfq/=1)) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'pfq(): Fatal error!'
    write ( *, '(a)' ) '  error in input arguments: lnpfq /= 0 or 1'
    stop 1
  endif

  if ((ip>iq) .and. (abs(z)> 1.0D+00 )) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'pfq(): Fatal error!'
    write ( *, '(a)' ) '  iq < ip and 1.0 < |z|.'
    write ( *, '(a)' ) '  The series does not converge.'
    stop 1
  endif

  if (ip==2 .and. iq==1 .and. abs(z)>0.9) then

    if (lnpfq/=1) then
!
!  check to see if the gamma function arguments are o.k.; if not,
!  then the series will have to be used.
!      precis - machine precision
!
      precis= 1.0D+00 
      precis=precis/ 2.0d+00
      dnum=precis+ 1.0D+00 
      do while (dnum> 1.0D+00 )
        precis=precis/ 2.0d+00
        dnum=precis+ 1.0D+00 
      end do
      precis= 2.0d+00 *precis

      do i=1 , 6

        if (i==1) then
          argi=aimag(b(1))
          argr=real(b(1), kind = rk )
        elseif (i==2) then
          argi=aimag(b(1)-a(1)-a(2))
          argr=real(b(1)-a(1)-a(2), kind = rk )
        elseif (i==3) then
          argi=aimag(b(1)-a(1))
          argr=real(b(1)-a(1), kind = rk )
        elseif (i==4) then
          argi=aimag(a(1)+a(2)-b(1))
          argr=real(a(1)+a(2)-b(1), kind = rk )
        elseif (i==5) then
          argi=aimag(a(1))
          argr=real(a(1), kind = rk )
        elseif (i==6) then
          argi=aimag(a(2))
          argr=real(a(2), kind = rk )
        endif
!
!  cases where the argument is real
!
        if (argi==0.0) then
!
!  cases where the argument is real and negative
!
!
!  use the series expansion if the argument is too near a pole
!
          if (argr<=0.0) then
            diff = abs ( real ( nint(argr), kind = rk )-argr )
            if (diff<= 2.0d+00 *precis) then
              pfq=hyper(a,b,ip,iq,z,lnpfq,ix,nsigfig)
              return
            endif
          endif

        endif

      enddo

      gam1=cgamma(b(1),lnpfq)
      gam2=cgamma(b(1)-a(1)-a(2),lnpfq)
      gam3=cgamma(b(1)-a(1),lnpfq)
      gam4=cgamma(b(1)-a(2),lnpfq)
      gam5=cgamma(a(1)+a(2)-b(1),lnpfq)
      gam6=cgamma(a(1),lnpfq)
      gam7=cgamma(a(2),lnpfq)
      a1(1)=a(1)
      a1(2)=a(2)
      b1(1)=a(1)+a(2)-b(1)+ 1.0D+00 
      z1= 1.0D+00 -z
      hyper1=hyper(a1,b1,ip,iq,z1,lnpfq,ix,nsigfig)
      a1(1)=b(1)-a(1)
      a1(2)=b(1)-a(2)
      b1(1)=b(1)-a(1)-a(2)+ 1.0D+00 
      hyper2=hyper(a1,b1,ip,iq,z1,lnpfq,ix,nsigfig)
      pfq=gam1*gam2*hyper1/(gam3*gam4)+( 1.0D+00 -z)**(b(1)-a(1)-a(2)) &
        *gam1*gam5*hyper2/(gam6*gam7)
      return
    endif

  endif

  pfq = hyper ( a, b, ip, iq, z, lnpfq, ix, nsigfig )

  return
end

subroutine aradd ( a, b, c, z, l, rmax )

!*****************************************************************************80
!
!! aradd() computes the sum of two arrays.
!
!  Discussion:
!
!    Accepts two arrays of numbers and returns the sum of the array.  
!    Each array is holding the value of one number in the series.  
!    The parameter L is the size of the array representing the number 
!    and RMAX is the actual number of digits needed to give the numbers
!    the desired accuracy.   
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(-1:*)
  real ( kind = rk ) b(-1:*)
  real ( kind = rk ) c(-1:*)
  integer ediff
  integer goon190
  integer goon300
  integer i
  integer j
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) z(-1:*)

  do i=0 , l+1
    z(i)=0.0
  enddo
  ediff = int ( anint(a(l+1)-b(l+1)) )
  if ( abs(a(1))< 0.5d+00 .or. ediff<=-l ) then
    do i=-1 , l+1
      c(i)=b(i)
    enddo
    if (c(1)< 0.5d+00 ) then
      c(-1)= 1.0d+00 
      c(l+1)=0.0
    endif
    return
  else
    if (abs(b(1))< 0.5d+00 .or. ediff>=l) then
      do i=-1 , l+1
        c(i)=a(i)
      enddo
      if (c(1)< 0.5d+00 ) then
        c(-1)= 1.0d+00 
        c(l+1)=0.0
      endif
      return
    else
      z(-1)=a(-1)
      goon300=1
      goon190=1
      if (abs(a(-1)-b(-1))>= 0.5d+00 ) then
        goon300=0
        if (ediff>0) then
          z(l+1)=a(l+1)
        elseif (ediff<0) then
          z(l+1)=b(l+1)
          z(-1)=b(-1)
          goon190=0
        else
          do i=1 , l
            if (a(i)>b(i)) then
              z(l+1)=a(l+1)
              exit
            endif
            if (a(i)<b(i)) then
              z(l+1)=b(l+1)
              z(-1)=b(-1)
              goon190=0
            endif
          enddo
        end if

      elseif (ediff>0) then
        z(l+1)=a(l+1)
        do i=l , 1+ediff , -1
          z(i)=a(i)+b(i-ediff)+z(i)
          if (z(i)>=rmax) then
            z(i)=z(i)-rmax
            z(i-1)= 1.0d+00 
          endif
        enddo
        do i=ediff , 1 , -1
          z(i)=a(i)+z(i)
          if (z(i)>=rmax) then
            z(i)=z(i)-rmax
            z(i-1)= 1.0d+00 
          endif
        enddo
        if (z(0)> 0.5d+00 ) then
          do i=l , 1 , -1
            z(i)=z(i-1)
          enddo
          z(l+1)=z(l+1)+1
          z(0)=0.0
        endif
      elseif (ediff<0) then
        z(l+1)=b(l+1)
        do i=l , 1-ediff , -1
          z(i)=a(i+ediff)+b(i)+z(i)
          if (z(i)>=rmax) then
            z(i)=z(i)-rmax
            z(i-1)= 1.0d+00 
          endif
        enddo
        do i=0-ediff , 1 , -1
          z(i)=b(i)+z(i)
          if (z(i)>=rmax) then
            z(i)=z(i)-rmax
            z(i-1)= 1.0d+00 
          endif
        enddo
        if (z(0)> 0.5d+00 ) then
          do i=l , 1 , -1
            z(i)=z(i-1)
          enddo
          z(l+1)=z(l+1)+ 1.0d+00 
          z(0)=0.0
        endif
      else
        z(l+1)=a(l+1)
        do i=l , 1 , -1
          z(i)=a(i)+b(i)+z(i)
          if (z(i)>=rmax) then
            z(i)=z(i)-rmax
            z(i-1)= 1.0d+00 
          endif
        enddo
        if (z(0)> 0.5d+00 ) then
          do i=l , 1 , -1
            z(i)=z(i-1)
          enddo
          z(l+1)=z(l+1)+ 1.0d+00 
          z(0)=0.0
        endif
      endif
      if (goon300==1) then
        i=i+1
        do while (z(i)< 0.5d+00 .and. i<l+1)
          i=i+1
        end do
        if (i==l+1) then
          z(-1)= 1.0d+00 
          z(l+1)=0.0
          do i=-1 , l+1
            c(i)=z(i)
          enddo
          if (c(1)< 0.5d+00 ) then
            c(-1)= 1.0d+00 
            c(l+1)=0.0
          endif
          return
        endif
        do j=1 , l+1-i
          z(j)=z(j+i-1)
        enddo
        do j=l+2-i , l
          z(j)=0.0
        enddo
        z(l+1)=z(l+1)-i+1
        do i=-1 , l+1
          c(i)=z(i)
        enddo
        if (c(1)< 0.5d+00 ) then
          c(-1)= 1.0d+00 
          c(l+1)=0.0
        endif
        return
      end if

      if (goon190==1) then
        if (ediff>0) then
          do i=l , 1+ediff , -1
            z(i)=a(i)-b(i-ediff)+z(i)
            if (z(i)<0.0) then
              z(i)=z(i)+rmax
              z(i-1)=- 1.0d+00 
            endif
          enddo
          do i=ediff , 1 , -1
            z(i)=a(i)+z(i)
            if (z(i)<0.0) then
              z(i)=z(i)+rmax
              z(i-1)=- 1.0d+00 
            endif
          enddo
        else
          do i=l , 1 , -1
            z(i)=a(i)-b(i)+z(i)
            if (z(i)<0.0) then
              z(i)=z(i)+rmax
              z(i-1)=- 1.0d+00 
            endif
          enddo
        endif
        if (z(1)> 0.5d+00 ) then
          do i=-1 , l+1
            c(i)=z(i)
          enddo
          if (c(1)< 0.5d+00 ) then
            c(-1)= 1.0d+00 
            c(l+1)=0.0
          endif
          return
        endif
        i=1
        i=i+1
        do while (z(i)< 0.5d+00 .and. i<l+1)
          i=i+1
        end do
        if (i==l+1) then
          z(-1)= 1.0d+00 
          z(l+1)=0.0
          do i=-1 , l+1
            c(i)=z(i)
          enddo
          if (c(1)< 0.5d+00 ) then
            c(-1)= 1.0d+00 
            c(l+1)=0.0
          endif
          return
        endif
        do j=1 , l+1-i
          z(j)=z(j+i-1)
        enddo
        do j=l+2-i , l
          z(j)=0.0
        enddo
        z(l+1)=z(l+1)-i+1
        do i=-1 , l+1
          c(i)=z(i)
        enddo
        if (c(1)< 0.5d+00 ) then
          c(-1)= 1.0d+00 
          c(l+1)=0.0
        endif
        return
      end if
    endif

    if (ediff<0) then
      do i=l , 1-ediff , -1
        z(i)=b(i)-a(i+ediff)+z(i)
        if (z(i)<0.0) then
          z(i)=z(i)+rmax
          z(i-1)=- 1.0d+00 
        endif
      enddo
      do i=0-ediff , 1 , -1
        z(i)=b(i)+z(i)
        if (z(i)<0.0) then
          z(i)=z(i)+rmax
          z(i-1)=- 1.0d+00 
        endif
      enddo
    else
      do i=l , 1 , -1
        z(i)=b(i)-a(i)+z(i)
        if (z(i)<0.0) then
          z(i)=z(i)+rmax
          z(i-1)=- 1.0d+00 
        endif
      enddo
    endif
  endif

  if (z(1)> 0.5d+00 ) then
    do i=-1 , l+1
      c(i)=z(i)
    enddo
    if (c(1)< 0.5d+00 ) then
      c(-1)= 1.0d+00 
      c(l+1)=0.0
    endif
    return
  endif

  i=1
  i=i+1
  do while (z(i)< 0.5d+00 .and. i<l+1)
    i=i+1
  end do

  if (i==l+1) then
    z(-1)= 1.0d+00 
    z(l+1)=0.0
    do i=-1 , l+1
      c(i)=z(i)
    enddo
    if (c(1)< 0.5d+00 ) then
      c(-1)= 1.0d+00 
      c(l+1)=0.0
    endif
    return
  endif

  do j=1 , l+1-i
    z(j)=z(j+i-1)
  enddo
  do j=l+2-i , l
    z(j)=0.0
  enddo
  z(l+1)=z(l+1)-i+1
  do i=-1 , l+1
    c(i)=z(i)
  enddo

  if (c(1)< 0.5d+00 ) then
    c(-1)= 1.0d+00 
    c(l+1)=0.0
  endif

  return
end

subroutine armult ( a, b, c, z, l, rmax )

!*****************************************************************************80
!
!! armult() returns the product of two arrays.
!
!  Discussion:
!
!    Accepts two arrays and returns the product.
!    L and RMAX are the size of the arrays and the number of
!    digits needed to represent the numbers with the required
!    accuracy.     
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(-1:*)
  real ( kind = rk ) b
  real ( kind = rk ) b2
  real ( kind = rk ) c(-1:*)
  real ( kind = rk ) carry
  integer i
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) z(-1:*)

  z(-1) = sign ( 1.0D+00, b ) * a(-1)
  b2=abs(b)
  z(l+1)=a(l+1)
  do i=0 , l
    z(i)=0.0
  enddo

  if ( b2 <= 1.0D-10 .or. a(1) <= 1.0D-10 ) then
    z(-1) = 1.0D+00
    z(l+1)=0.0
  else
    do i=l , 1 , -1
      z(i)=a(i)*b2+z(i)
      if (z(i)>=rmax) then
        carry=aint(z(i)/rmax)
        z(i)=z(i)-carry*rmax
        z(i-1)=carry
      endif
    enddo
    if (z(0)>= 0.5d+00 ) then
      do i=l , 1 , -1
        z(i)=z(i-1)
      enddo
      z(l+1)=z(l+1)+ 1.0D+00
      if (z(1)>=rmax) then
        do i=l , 1 , -1
          z(i)=z(i-1)
        enddo
        carry=aint(z(1)/rmax)
        z(2)=z(2)-carry*rmax
        z(1)=carry
        z(l+1)=z(l+1)+ 1.0D+00
      endif
      z(0)=0.0
    endif
  end if

  do i=-1 , l+1
    c(i)=z(i)
  enddo

  if (c(1)< 0.5d+00 ) then
    c(-1)= 1.0D+00
    c(l+1)=0.0
  endif

  return
end

subroutine arsub ( a, b, c, wk1, wk2, l, rmax )

!*****************************************************************************80
!
!! arsub() subtracts one array from another.
!
!  Discussion:
!
!    Accepts two arrays and subtracts each element
!    in the second array from the element in the first array
!    and returns the solution.  The parameters L and RMAX are
!    the size of the array and the number of digits needed for
!    the accuracy, respectively.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(-1:*)
  real ( kind = rk ) b(-1:*)
  real ( kind = rk ) c(-1:*)
  integer i
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) wk1(-1:*)
  real ( kind = rk ) wk2(-1:*)
 
  do i=-1 , l+1
    wk2(i)=b(i)
  enddo
  wk2(-1)=- wk2(-1)
  call aradd(a,wk2,c,wk1,l,rmax)

  return
end

subroutine arydiv ( ar, ai, br, bi, c, l, lnpfq, rmax, ibit )

!*****************************************************************************80
!
!! arydiv() divides arrays.
!
!  Discussion:
!
!    Returns the real ( kind = rk ) complex number resulting from the 
!    division of four arrays, representing two complex numbers.  The 
!    number returned will be in one of two different forms:  either 
!    standard scientific or as the log (base 10) of the number.  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ae(2,2)
  real ( kind = rk ) ai(-1:*)
  real ( kind = rk ) ar(-1:*)
  real ( kind = rk ) be(2,2)
  real ( kind = rk ) bi(-1:*)
  real ( kind = rk ) br(-1:*)
  complex ( kind = ck ) c
  complex ( kind = ck ) cdum
  real ( kind = rk ) ce(2,2)
  real ( kind = rk ) dnum
  real ( kind = rk ) dum1
  real ( kind = rk ) dum2
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  integer ibit
  integer ii10
  integer ir10
  integer itnmax
  integer l
  integer lnpfq
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) n3
  real ( kind = rk ) phi
  integer rexp
  real ( kind = rk ) ri10
  real ( kind = rk ) rmax
  real ( kind = rk ) rr10
  real ( kind = rk ) tenmax
  real ( kind = rk ) x
  real ( kind = rk ) x1
  real ( kind = rk ) x2

  rexp=int(ibit/2)
  x=rexp*(ar(l+1)-2)
  rr10=x*log10( 2.0d+00 )/log10 ( 10.0d+00 )
  ir10=int(rr10)
  rr10=rr10-ir10
  x=rexp*(ai(l+1)-2)
  ri10=x*log10( 2.0d+00 )/log10 ( 10.0d+00 )
  ii10=int(ri10)
  ri10=ri10-ii10
  dum1=sign(ar(1)*rmax*rmax+ar(2)*rmax+ar(3),ar(-1))
  dum2=sign(ai(1)*rmax*rmax+ai(2)*rmax+ai(3),ai(-1))
  dum1=dum1*10**rr10
  dum2=dum2*10**ri10
  cdum=cmplx(dum1,dum2,kind = ck)
  call conv12(cdum,ae)
  ae(1,2)=ae(1,2)+ir10
  ae(2,2)=ae(2,2)+ii10
  x=rexp*(br(l+1)-2)
  rr10=x*log10( 2.0d+00 )/log10 ( 10.0d+00 )
  ir10=int(rr10)
  rr10=rr10-ir10
  x=rexp*(bi(l+1)-2)
  ri10=x*log10( 2.0d+00 )/log10 ( 10.0d+00 )
  ii10=int(ri10)
  ri10=ri10-ii10
  dum1=sign(br(1)*rmax*rmax+br(2)*rmax+br(3),br(-1))
  dum2=sign(bi(1)*rmax*rmax+bi(2)*rmax+bi(3),bi(-1))
  dum1=dum1*10**rr10
  dum2=dum2*10**ri10
  cdum=cmplx(dum1,dum2,kind = ck)
  call conv12(cdum,be)
  be(1,2)=be(1,2)+ir10
  be(2,2)=be(2,2)+ii10
  call ecpdiv(ae,be,ce)

  if (lnpfq==0) then
    call conv21(ce,c)
  else
    call emult(ce(1,1),ce(1,2),ce(1,1),ce(1,2),n1,e1)
    call emult(ce(2,1),ce(2,2),ce(2,1),ce(2,2),n2,e2)
    call eadd(n1,e1,n2,e2,n3,e3)
    n1=ce(1,1)
    e1=ce(1,2)-ce(2,2)
    x2=ce(2,1)
!
!  tenmax - maximum size of exponent of 10
!  the following code can be used to determine tenmax, but it
!  will likely generate an ieee floating point underflow error
!  on a sun workstation.  replace tenmax with the value appropriate
!  for your machine.
!
    tenmax=320
    itnmax=1
    dnum=0.1d0
    itnmax=itnmax+1
    dnum=dnum*0.1d0
    do while (dnum>0.0) 
      itnmax=itnmax+1
      dnum=dnum*0.1d0
    end do

    itnmax=itnmax-1
    tenmax= real (itnmax, kind = rk )

    if (e1>tenmax) then
      x1=tenmax
    elseif (e1<-tenmax) then
      x1=0.0
    else
      x1=n1*( 10.0d+00 **e1)
    endif

    if (x2/=0.0) then
      phi=atan2(x2,x1)
    else
      phi=0.0
    endif

    c = cmplx ( 0.5d+00 *(log(n3)+e3*log(10.0d+00)), phi, kind = ck)

  endif

  return
end

function bits ( )

!*****************************************************************************80
!
!! bits() determines number of significant figures needed.
!
!  Discussion:
!
!    Determines the number of significant figures  
!    of machine precision to arrive at the size of the array
!    the numbers must be stored in to get the accuracy of the
!    solution.   
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
! 
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) bit
  real ( kind = rk ) bit2
  real ( kind = rk ) bits
  integer count
 
  bit=1.0
  count=0
  count=count+1
  bit2=bit*2.0
  bit=bit2+1.0

  do while ((bit-bit2)/=0.0)
    count=count+1
    bit2=bit*2.0
    bit=bit2+1.0
  end do

  bits = count - 3

  return
end

function cgamma ( arg, lnpfq )

!*****************************************************************************80
!
!! cgamma() calculates the complex gamma function.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
!  Reference:
!
!    Farid Parpia, Charlotte Froese Fischer, I. P. Grant,
!    GRASP92: a package for large-scale relativistic atomic 
!    structure calculations,
!    Computer Physics Communications,
!    Volume 175, 2006, pages 745-747.
!
!  Local:
!
!    real ( kind = rk ) expmax: maximum size of exponent of e.
!
!    real ( kind = rk ) precis: machine precision.
!
!    real ( kind = rk ) tenmax: maximum exponent of 10.
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) arg
  real ( kind = rk ) argi
  real ( kind = rk ) argr
  real ( kind = rk ) argui
  real ( kind = rk ) argui2
  real ( kind = rk ) argum
  real ( kind = rk ) argur
  real ( kind = rk ) argur2
  complex ( kind = ck ) cgamma
  real ( kind = rk ) clngi
  real ( kind = rk ) clngr
  real ( kind = rk ) diff
  real ( kind = rk ) dnum
  real ( kind = rk ) expmax
  real ( kind = rk ) fac
  real ( kind = rk ) facneg
  real ( kind = rk ), save, dimension ( 7 ) :: fd = (/ &
    6.0d00, 30.0d00, 42.0d00, 30.0d00, 66.0d00, 2730.0d00, 6.0d00 /)
  real ( kind = rk ), save, dimension ( 7 ) :: fn = (/ &
    1.0d00, -1.0d00, 1.0d00, -1.0d00, 5.0d00, -691.0d00, 7.0d00 /)
  logical first
  real ( kind = rk ) hlntpi
  integer i
  integer itnmax
  integer lnpfq
  logical negarg
  real ( kind = rk ) obasq
  real ( kind = rk ) obasqi
  real ( kind = rk ) obasqr
  real ( kind = rk ) ovlfac
  real ( kind = rk ) ovlfi
  real ( kind = rk ) ovlfr
  real ( kind = rk ) pi
  real ( kind = rk ) precis
  real ( kind = rk ) resi
  real ( kind = rk ) resr
  real ( kind = rk ) tenmax
  real ( kind = rk ), parameter :: tenth = 0.1D+00
  real ( kind = rk ) termi
  real ( kind = rk ) termr
  real ( kind = rk ) twoi
  real ( kind = rk ) zfaci
  real ( kind = rk ) zfacr

  argr = real (arg, kind = rk )
  argi = aimag(arg)
!
!  set the machine-dependent parameters:
!
  first=.true.

  if (first) then

    pi=4.0d0*atan( 1.0d+00 )

    itnmax=1
    dnum=tenth
    itnmax=itnmax+1
    dnum=dnum*tenth
    do while (dnum>0.0) 
      itnmax=itnmax+1
      dnum=dnum*tenth
    end do
    itnmax=itnmax-1
    tenmax=real ( itnmax, kind = rk )
    dnum=tenth**itnmax
    expmax=-log(dnum)

    precis= 1.0d+00 
    precis = precis / 2.0d+00
    dnum=precis+ 1.0d+00 
    do while (dnum> 1.0d+00 ) 
      precis=precis / 2.0d+00
      dnum=precis+ 1.0d+00 
    end do
    precis = 2.0d+00 * precis

    hlntpi= 0.5d+00 *log(2.0d+00*pi)

    do i=1 , 7
      fn(i)=fn(i)/fd(i)
      twoi = 2.0d+00 * real ( i, kind = rk )
      fn(i)=fn(i)/(twoi*(twoi- 1.0d+00 ))
    enddo

    first=.false.

  endif
!
!  cases where the argument is real
!
  if (argi==0.0) then
!
!  cases where the argument is real and negative
!
    if (argr<=0.0) then
!
!  stop with an error message if the argument is too near a pole
!
      diff=abs( real (nint(argr), kind = rk )-argr)

      if (diff<= 2.0d+00 *precis) then

        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'cgamma(): fatal error!'
        write ( *, '(a)' ) '  Argument too close to a pole.'
        stop 1

      else
!
!  otherwise use the reflection formula (abramowitz and stegun 6.1.17)
!  to ensure that the argument is suitable for stirling's formula
!
        argum=pi/(-argr*sin(pi*argr))
        if (argum<0.0) then
          argum=-argum
          clngi=pi
        else
          clngi=0.0
        endif

        facneg=log(argum)
        argur=-argr
        negarg=.true.

      endif
!
!  cases where the argument is real and positive
!
    else

     clngi=0.0
     argur=argr
     negarg=.false.

    endif
!
!  use abramowitz and stegun formula 6.1.15 to ensure that
!  the argument in stirling's formula is greater than 10
!
    ovlfac= 1.0d+00 
    do while (argur< 10.0d+00 )
      ovlfac=ovlfac*argur
      argur=argur+ 1.0d+00 
    end do
!
!  Use stirling's formula to compute log (gamma (argum))
!
    clngr=(argur- 0.5d+00 )*log(argur)-argur+hlntpi
    fac=argur
    obasq= 1.0d+00 /(argur*argur)
    do i=1 , 7
      fac=fac*obasq
      clngr=clngr+fn(i)*fac
    enddo
!
!  include the contributions from the recurrence and reflection formulas.
!
    clngr=clngr-log(ovlfac)
    if (negarg) then
      clngr=facneg-clngr
    endif

  else
!
!  cases where the argument is complex
!
    argur=argr
    argui=argi
    argui2=argui*argui
!
!  use the recurrence formula (abramowitz and stegun 6.1.15)
!  to ensure that the magnitude of the argument in stirling's
!  formula is greater than 10
!
    ovlfr= 1.0d+00 
    ovlfi=0.0
    argum=sqrt(argur*argur+argui2)
    do while (argum< 10.0d+00 )
      termr=ovlfr*argur-ovlfi*argui
      termi=ovlfr*argui+ovlfi*argur
      ovlfr=termr
      ovlfi=termi
      argur=argur+ 1.0d+00 
      argum=sqrt(argur*argur+argui2)
    end do
!
!  use stirling's formula to compute log (gamma (argum))
!
    argur2=argur*argur
    termr= 0.5d+00 *log(argur2+argui2)
    termi=atan2(argui,argur)
    clngr=(argur- 0.5d+00 )*termr-argui*termi-argur+hlntpi
    clngi=(argur- 0.5d+00 )*termi+argui*termr-argui
    fac=(argur2+argui2)**(-2)
    obasqr=(argur2-argui2)*fac
    obasqi=- 2.0d+00 *argur*argui*fac
    zfacr=argur
    zfaci=argui
    do i=1 , 7
      termr=zfacr*obasqr-zfaci*obasqi
      termi=zfacr*obasqi+zfaci*obasqr
      fac=fn(i)
      clngr=clngr+termr*fac
      clngi=clngi+termi*fac
      zfacr=termr
      zfaci=termi
    enddo
!
!  add the relevant pieces from the recurrence formula
!
    clngr=clngr- 0.5d+00 *log(ovlfr*ovlfr+ovlfi*ovlfi)
    clngi=clngi-atan2(ovlfi,ovlfr)

   endif

   if (lnpfq==1) then
     cgamma = cmplx ( clngr, clngi, kind = ck )
     return
   endif
!
!  exponentiate the complex log gamma function to get
!  the complex gamma function
!
  if ((clngr<=expmax) .and. (clngr>=-expmax)) then
    fac=exp(clngr)
  else
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'cgamma(): fatal error!'
    write ( *, '(a)' ) '  exponential argument is out of range.'
    stop 1
  endif

  resr=fac * cos(clngi)
  resi=fac * sin(clngi)
  cgamma = cmplx ( resr, resi, kind = ck )

  return
end

subroutine cmpadd ( ar, ai, br, bi, cr, ci, wk1, l, rmax )

!*****************************************************************************80
!
!! cmpadd() computes the complex sum of two arrays.
!
!  Discussion:
!
!    Takes two arrays representing one real and
!    one imaginary part, and adds two arrays representing
!    another complex number and returns two array holding the
!    complex sum.
!
!      (CR,CI) = (AR+BR, AI+BI) 
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ai(-1:*)
  real ( kind = rk ) ar(-1:*) 
  real ( kind = rk ) bi(-1:*)
  real ( kind = rk ) br(-1:*)
  real ( kind = rk ) ci(-1:*)
  real ( kind = rk ) cr(-1:*)
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) wk1(-1:*)

  call aradd ( ar, br, cr, wk1, l, rmax )
  call aradd ( ai, bi, ci, wk1, l, rmax )

  return
end

subroutine cmpmul ( ar, ai, br, bi, cr, ci, wk1, wk2, cr2, d1, d2, wk6, l, &
  rmax )

!*****************************************************************************80
!
!! cmpmul() multiplies two complex arrays.
!
!  Discussion:
!
!    Takes two arrays representing one real and one imaginary part, 
!    and multiplies it with two arrays representing another complex 
!    number and returns the complex product.   
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ai(-1:*)
  real ( kind = rk ) ar(-1:*)
  real ( kind = rk ) bi
  real ( kind = rk ) br
  real ( kind = rk ) ci(-1:*)
  real ( kind = rk ) cr(-1:*)
  real ( kind = rk ) cr2(-1:*)
  real ( kind = rk ) d1(-1:*)
  real ( kind = rk ) d2(-1:*)
  integer i
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) wk1(-1:*)
  real ( kind = rk ) wk2(-1:*)
  real ( kind = rk ) wk6(-1:*)

  call armult(ar,br,d1,wk6,l,rmax)
  call armult(ai,bi,d2,wk6,l,rmax)
  call arsub(d1,d2,cr2,wk1,wk2,l,rmax)
  call armult(ar,bi,d1,wk6,l,rmax)
  call armult(ai,br,d2,wk6,l,rmax)
  call aradd(d1,d2,ci,wk1,l,rmax)
  do i = -1, l + 1
    cr(i)=cr2(i)
  enddo

  return
end

subroutine cmpsub ( ar, ai, br, bi, cr, ci, wk1, wk2, l, rmax )

!*****************************************************************************80
!
!! cmpsub() takes the complex difference of two arrays.
!
!  Discussion:
!
!    Takes two arrays representing one real and
!    one imaginary part, and subtracts two arrays representing
!    another complex number and returns two array holding the
!    complex sum.
!
!    (CR,CI) = (AR+BR, AI+BI)   
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) ai(-1:*)
  real ( kind = rk ) ar(-1:*)
  real ( kind = rk ) bi(-1:*)
  real ( kind = rk ) br(-1:*)
  real ( kind = rk ) ci(-1:*)
  real ( kind = rk ) cr(-1:*)
  integer l
  real ( kind = rk ) rmax
  real ( kind = rk ) wk1(-1:*)
  real ( kind = rk ) wk2(-1:*)

  call arsub(ar,br,cr,wk1,wk2,l,rmax)
  call arsub(ai,bi,ci,wk1,wk2,l,rmax)

  return
end

subroutine conv12 ( cn, cae )

!*****************************************************************************80
!
!! conv12() converts a number from complex notation to a 2x2 real array.    
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cae(2,2)
  complex ( kind = ck ) cn

  cae(1,1)=real(cn, kind = rk )
  cae(1,2)=0.0
  do
    if (abs(cae(1,1))< 10.0d+00 ) then
      do
        if ((abs(cae(1,1))>= 1.0d+00 ) .or. (cae(1,1)==0.0)) then
          cae(2,1) = aimag(cn)
          cae(2,2) = 0.0
          do
            if (abs(cae(2,1))< 10.0d+00 ) then
              do while ((abs(cae(2,1))< 1.0d+00 ) .and. (cae(2,1)/=0.0)) 
                cae(2,1)=cae(2,1)* 10.0d+00
                cae(2,2)=cae(2,2)- 1.0d+00
              end do
              exit
            else
              cae(2,1)=cae(2,1)/ 10.0d+00
              cae(2,2)=cae(2,2)+ 1.0d+00
            endif
          end do
          exit
        else
          cae(1,1)=cae(1,1)* 10.0d+00
          cae(1,2)=cae(1,2)- 1.0d+00
        endif
      end do
      exit
    else
      cae(1,1)=cae(1,1)/ 10.0d+00
      cae(1,2)=cae(1,2)+ 1.0d+00
    endif
  end do

  return
end

subroutine conv21 ( cae, cn )

!*****************************************************************************80
!
!! conv21() converts a number in a 2x2 real array to a complex number.  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
!  Local:
!
!    real ( kind = rk ) tenmax: maximum exponent of 10
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) cae(2,2)
  complex ( kind = ck ) cn
  real ( kind = rk ) dnum
  integer itnmax
  real ( kind = rk ) tenmax

  itnmax=1
  dnum=0.1d0
  itnmax=itnmax+1
  dnum=dnum*0.1d0
  do while  (dnum>0.0) 
    itnmax=itnmax+1
    dnum=dnum*0.1d0
  end do
  itnmax=itnmax-2
  tenmax = real (itnmax, kind = rk )

  if (cae(1,2)>tenmax .or. cae(2,2)>tenmax) then

    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'conv21(): Fatal error!'
    write ( *, '(a)' ) '  The exponent required for summation is larger'
    write ( *, '(a)' ) '  than the machine can handle.'
    write ( *, '(a)' ) '  Try repeating the calculation with'
    write ( *, '(a)' ) '  lnpfq = 1.'
    stop 1
  elseif (cae(2,2)<-tenmax) then
    cn = cmplx ( cae(1,1)*( 10.0d+00 **cae(1,2)), 0.0, kind = ck )
  else
    cn = cmplx ( cae(1,1)*( 10.0d+00 **cae(1,2)), &
      cae(2,1)*(10.0d+00**cae(2,2)), kind = ck )
  endif

 return
end

subroutine eadd ( n1, e1, n2, e2, nf, ef )

!*****************************************************************************80
!
!! eadd() returns the sum of two numbers as a base and an exponent.    
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) ef
  real ( kind = rk ) ediff
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) nf

  ediff=e1-e2

  if (ediff>36.0d0) then

    nf=n1
    ef=e1

  elseif (ediff<-36.0d0) then

    nf=n2
    ef=e2

  else

    nf=n1*( 10.0d+00 **ediff)+n2
    ef=e2
    do
      if (abs(nf) < 10.0d+00 ) then
        do while ((abs(nf)< 1.0d+00 ) .and. (nf/=0.0)) 
          nf=nf* 10.0d+00
          ef=ef- 1.0d+00
        end do
        exit
      else
        nf=nf/ 10.0d+00
        ef=ef+ 1.0d+00
      endif
    end do

  endif

  return
end

subroutine ecpdiv ( a, b, c )

!*****************************************************************************80
!
!! ecpdiv() divides two numbers represented as 2x2 arrays.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) b(2,2)
  real ( kind = rk ) b2(2,2)
  real ( kind = rk ) c(2,2)
  real ( kind = rk ) c2(2,2)
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) e3
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) n3
 
  b2(1,1)=b(1,1)
  b2(1,2)=b(1,2)
  b2(2,1)=- b(2,1)
  b2(2,2)=b(2,2)

  call ecpmul(a,b2,c2)
  call emult(b(1,1),b(1,2),b(1,1),b(1,2),n1,e1)
  call emult(b(2,1),b(2,2),b(2,1),b(2,2),n2,e2)
  call eadd(n1,e1,n2,e2,n3,e3)
  call ediv(c2(1,1),c2(1,2),n3,e3,c(1,1),c(1,2))
  call ediv(c2(2,1),c2(2,2),n3,e3,c(2,1),c(2,2))

  return
end

subroutine ecpmul ( a, b, c )

!*****************************************************************************80
!
!! ecpmul() multiplies numbers stored as 2x2 arrays.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) a(2,2)
  real ( kind = rk ) b(2,2)
  real ( kind = rk ) c(2,2)
  real ( kind = rk ) c2(2,2)
  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) n1
  real ( kind = rk ) n2

  call emult(a(1,1),a(1,2),b(1,1),b(1,2),n1,e1)
  call emult(a(2,1),a(2,2),b(2,1),b(2,2),n2,e2)
  call esub(n1,e1,n2,e2,c2(1,1),c2(1,2))
  call emult(a(1,1),a(1,2),b(2,1),b(2,2),n1,e1)
  call emult(a(2,1),a(2,2),b(1,1),b(1,2),n2,e2)
  call eadd(n1,e1,n2,e2,c(2,1),c(2,2))
  c(1,1)=c2(1,1)
  c(1,2)=c2(1,2)

  return
end

subroutine ediv ( n1, e1, n2, e2, nf, ef )

!*****************************************************************************80
!
!! ediv() computes the ratio of two exponential numbers.
!
!  Discussion:
!
!    The solution is returned as a base and an exponent.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) ef
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) nf

  nf = n1 / n2
  ef = e1 - e2

  if ((abs(nf)< 1.0D+00 ) .and. (nf/= 0.0d+00 )) then
    nf=nf * 10.0d+00
    ef=ef- 1.0D+00
  endif

  return
end

subroutine emult ( n1, e1, n2, e2, nf, ef )

!*****************************************************************************80
!
!! emult() multiplies using the base and exponent format.
!
!  Discussion:
!
!    Takes one base and exponent and multiplies it
!    by another numbers base and exponent to give the product
!    in the form of base and exponent.    
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) ef
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) nf

  nf=n1*n2
  ef=e1+e2
  if (abs(nf)>= 10.0d+00 ) then
    nf=nf / 10.0d+00
    ef=ef+ 1.0d+00
  endif

  return
end

subroutine esub ( n1, e1, n2, e2, nf, ef )

!*****************************************************************************80
!
!! esub() subtracts numbers in base and exponent format.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  real ( kind = rk ) e1
  real ( kind = rk ) e2
  real ( kind = rk ) ef
  real ( kind = rk ) n1
  real ( kind = rk ) n2
  real ( kind = rk ) nf
 
  call eadd ( n1, e1, - n2, e2, nf, ef )

  return
end

function factorial_log ( z )

!*****************************************************************************80
!
!! factorial_log() computes the log of the factorial.
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none
 
  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) factorial_log
  real ( kind = rk ) pi
  complex ( kind = ck ) z

  if (((real ( z, kind = rk )== 1.0d+00 ) .and. (aimag(z)==0.0d+00)) .or. &
    (abs(z)==0.0d+00)) then
    factorial_log = cmplx ( 0.0d+00, 0.0d+00, kind = ck )
    return
  endif

  pi = 4.0d+00 * atan( 1.0d+00 )

  factorial_log = 0.5d+00 * log(2.0d+00*pi) &
    + (z+ 0.5d+00 )*log(z)-z+( 1.0d+00 /(12.0d0*z)) &
    *( 1.0d+00 -( 1.0d+00 /(30.d0*z*z))*( 1.0d+00 -(2.0d+00/(7.0d0*z*z))))

  return
end

function hyper ( a, b, ip, iq, z, lnpfq, ix, nsigfig )

!*****************************************************************************80
!
!! hyper() sums the Gauss series.   
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  integer, parameter :: length = 777

  complex ( kind = ck ) a(ip)
  real ( kind = rk ) accy
  real ( kind = rk ) ai(10)
  real ( kind = rk ) ai2(10)
  real ( kind = rk ) ar(10)
  real ( kind = rk ) ar2(10)
  complex ( kind = ck ) b(iq)
  real ( kind = rk ) bar1(-1:length)
  real ( kind = rk ) bar2(-1:length)
  real ( kind = rk ) bits
  complex ( kind = ck ) cdum1
  complex ( kind = ck ) cdum2
  real ( kind = rk ) ci(10)
  real ( kind = rk ) ci2(10)
  real ( kind = rk ) cnt
  real ( kind = rk ) cr(10)
  real ( kind = rk ) cr2(10)
  real ( kind = rk ) creal
  real ( kind = rk ) denomi(-1:length)
  real ( kind = rk ) denomr(-1:length)
  real ( kind = rk ) dum1
  real ( kind = rk ) dum2
  real ( kind = rk ) expon
  complex ( kind = ck ) factorial_log
  complex ( kind = ck ) final
  real ( kind = rk ) foo1(-1:length)
  real ( kind = rk ) foo2(-1:length)
  integer goon1
  complex ( kind = ck ) hyper
  integer i
  integer i1
  integer ibit
  integer icount
  integer ii10
  integer ir10
  integer ip
  integer ipremax
  integer iq
  integer ix
  integer ixcnt
  integer l
  integer lmax
  integer lnpfq
  real ( kind = rk ) log2
  real ( kind = rk ) mx1
  real ( kind = rk ) mx2
  integer nmach
  integer nsigfig
  real ( kind = rk ) numi(-1:length)
  real ( kind = rk ) numr(-1:length)
  complex ( kind = ck ) oldtemp
  real ( kind = rk ) qi1(-1:length)
  real ( kind = rk ) qi2(-1:length)
  real ( kind = rk ) qr1(-1:length)
  real ( kind = rk ) qr2(-1:length)
  integer rexp
  real ( kind = rk ) ri10
  real ( kind = rk ) rmax
  real ( kind = rk ) rr10
  real ( kind = rk ) sigfig
  real ( kind = rk ) sumi(-1:length)
  real ( kind = rk ) sumr(-1:length)
  complex ( kind = ck ) temp
  complex ( kind = ck ) temp1
  real ( kind = rk ) wk1(-1:length)
  real ( kind = rk ) wk2(-1:length)
  real ( kind = rk ) wk3(-1:length)
  real ( kind = rk ) wk4(-1:length)
  real ( kind = rk ) wk5(-1:length)
  real ( kind = rk ) wk6(-1:length)
  real ( kind = rk ) x
  real ( kind = rk ) xi
  real ( kind = rk ) xi2
  real ( kind = rk ) xl
  real ( kind = rk ) xr
  real ( kind = rk ) xr2
  complex ( kind = ck ) z

  log2=log10 ( 2.0d+00 )
  ibit=int(bits())
  rmax = 2.0d+00 ** (int(ibit/2))
  sigfig = 2.0d+00 ** (int(ibit/4))

  do i1=1 , ip
    ar2(i1)=real(a(i1), kind = rk )*sigfig
    ar(i1)=aint(ar2(i1))
    ar2(i1)=anint((ar2(i1)-ar(i1))*rmax)
    ai2(i1)=aimag(a(i1))*sigfig
    ai(i1)=aint(ai2(i1))
    ai2(i1)=anint((ai2(i1)-ai(i1))*rmax)
  enddo

  do i1=1 , iq
    cr2(i1)=real(b(i1), kind = rk )*sigfig
    cr(i1)=aint(cr2(i1))
    cr2(i1)=anint((cr2(i1)-cr(i1))*rmax)
    ci2(i1)=aimag(b(i1))*sigfig
    ci(i1)=aint(ci2(i1))
    ci2(i1)=anint((ci2(i1)-ci(i1))*rmax)
  enddo

  xr2=real(z, kind = rk )*sigfig
  xr=aint(xr2)
  xr2=anint((xr2-xr)*rmax)
  xi2=aimag(z)*sigfig
  xi=aint(xi2)
  xi2=anint((xi2-xi)*rmax)
!
!  warn the user that the input value was so close to zero that it
!  was set equal to zero.
!
  do i1 = 1, ip

    if ( ( real(a(i1), kind = rk )/=0.0) .and. &
          (ar(i1)==0.0) .and. (ar2(i1)==0.0)) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'hyper(): warning.'
      write ( *, '(a,i2,a)' ) '  Real part of a(', i1, ') set to zero.'
    endif

    if ((aimag(a(i1))/=0.0) .and. (ai(i1)==0.0) .and. (ai2(i1)==0.0)) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'hyper(): warning.'
      write ( *, '(a,i2,a)' ) '  Imaginary part of a(', i1, ') set to zero'
    endif
  enddo

  do i1=1 , iq

    if ((real(b(i1), kind = rk )/=0.0) .and. (cr(i1)==0.0) &
      .and. (cr2(i1)==0.0)) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'hyper(): warning.'
      write ( *, '(a,i2,a)' ) '  Real part of b(', i1, ') set to zero.'
    endif

    if ((aimag(b(i1))/=0.0) .and. (ci(i1)==0.0) .and. (ci2(i1)==0.0)) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'hyper(): warning.'
      write ( *, '(a,i2,a)' ) '  Imaginary part of b(', i1, ') set to zero.'
    endif
  enddo

  if ((real(z, kind = rk ) /= 0.0) .and. (xr==0.0) .and. (xr2==0.0)) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'hyper(): warning.'
    write ( *, '(a)' ) '  Real part of z set to zero.'
    z = cmplx ( 0.0, aimag(z), kind = ck )
  endif

  if ((aimag(z)/=0.0) .and. (xi==0.0) .and. (xi2==0.0)) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'hyper(): warning.'
    write ( *, '(a)' ) '  Imginary part of z set to zero.'
    z = cmplx ( real(z, kind = rk ), 0.0, kind = ck )
  endif
!
!  screening of numerator arguments for negative integers or zero.
!  icount will force the series to terminate correctly.
!
  nmach=int(log10(2.0d+00**int(bits())))
  icount=-1
  do i1 = 1, ip

    if ((ar2(i1)==0.0) .and. (ar(i1)==0.0) .and. (ai2(i1)==0.0) .and.  &
      (ai(i1)==0.0)) then
      hyper = cmplx ( 1.0d+00, 0.0, kind = ck )
      return
    endif

    if ((ai(i1)==0.0) .and. (ai2(i1)==0.0) .and. (real(a(i1))<0.0)) then
      if (abs(real(a(i1), kind = rk ) &
        - real ( nint ( real ( a(i1), kind = rk ) ), kind = rk ) ) &
        < 10.0d+00 **(-nmach)) then
        if (icount/=-1) then
          icount=min ( icount, -nint ( real ( a(i1), kind = rk ) ) )
        else
          icount=-nint ( real ( a(i1), kind = rk ) )
        endif
      endif
    endif
  enddo
!
!  screening of denominator arguments for zeroes or negative integers
!
  do i1=1 , iq

    if ((cr(i1)==0.0) .and. (cr2(i1)==0.0) .and. (ci(i1)==0.0) .and.   &
      (ci2(i1)==0.0)) then
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'hyper(): Fatal error.'
      write ( *, '(a,i2,a)' ) '  Argument b(', i1, ') was zero.'
      stop 1
    endif

    if ((ci(i1)==0.0) .and. (ci2(i1)==0.0) .and. &
      ( real ( b(i1), kind = rk ) < 0.0) ) then
      if ( ( abs ( real ( b(i1), kind = rk ) &
        - real ( nint ( real ( b(i1), kind = rk ) ), kind = rk ) ) &
        < 10.0d+00 **(-nmach)) .and. &
        (icount>=-nint(real ( b(i1), kind = rk ) ) .or. icount==-1)) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'hyper(): Fatal error.'
        write ( *, '(a,i2,a)' ) '  Argument b(', i1, ') was a negative integer.'
        stop 1
      endif
    endif
  enddo

  nmach = int(log10(2.0d+00**ibit))
  nsigfig = min(nsigfig,int(log10(2.0d+00**ibit)))
  accy = 10.0d+00 **(-nsigfig)
  l = ipremax(a,b,ip,iq,z)
  if ( l /= 1 ) then
!
!  estimate the exponent of the maximum term in the pfq series
!     .
    expon=0.0
    xl = real(l, kind = rk )
    do i=1 , ip
      expon=expon + real (factorial_log(a(i)+xl- 1.0d+00 ), kind = rk ) &
        - real ( factorial_log ( a(i)- 1.0d+00 ), kind = rk )
    enddo

    do i=1 , iq
      expon = expon - real ( factorial_log(b(i)+xl- 1.0d+00 ), kind = rk ) &
        + real ( factorial_log(b(i)- 1.0d+00 ), kind = rk )
    enddo

    expon=expon + xl * real ( log(z), kind = rk ) &
      - real ( factorial_log ( cmplx ( xl, 0.0d+00, kind = ck ) ), kind = rk )
    lmax=int(log10(exp( 1.0d+00 ))*expon)
    l=lmax
!
!  estimate the exponent of where the pfq series will terminate.
!
    temp1= cmplx ( 1.0d+00, 0.0, kind = ck )
    creal= 1.0d+00 
    do i1=1 , ip
      temp1 = temp1 * cmplx ( ar(i1), ai(i1), kind = ck ) / sigfig
    enddo
    do i1=1 , iq
      temp1=temp1 / ( cmplx ( cr(i1), ci(i1), kind = ck ) / sigfig)
      creal=creal*cr(i1)
    enddo
    temp1 = temp1 * cmplx ( xr, xi, kind = ck )
!
!  triple it to make sure.
!
    l=3*l
!
!  divide the number of significant figures necessary by the number of
!  digits available per array position.
!
    l=int((2*l+nsigfig)/nmach)+2
  endif
!
!  make sure there are at least 5 array positions used.
!
  l=max(l,5)
  l=max(l,ix)

  if ( l > length ) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'hyper(): Fatal error!'
    write ( *, '(a)' ) '  Value of L exceeds maximum.'
    stop 1
  endif
 
  if (nsigfig>nmach) then
    write ( *, '(a)' ) ''
    write ( *, '(a)' ) 'hyper(): Warning.'
    write ( *, '(a)' ) '  The number of significant figures requested'
    write ( *, '(a)' ) '  is more than the machine precision.'
  endif

  sumr(-1)= 1.0d+00 
  sumi(-1)= 1.0d+00 
  numr(-1)= 1.0d+00 
  numi(-1)= 1.0d+00 
  denomr(-1)= 1.0d+00 
  denomi(-1)= 1.0d+00 
  do i=0 , l+1
    sumr(i)=0.0
    sumi(i)=0.0
    numr(i)=0.0
    numi(i)=0.0
    denomr(i)=0.0
    denomi(i)=0.0
  enddo

  sumr(1)= 1.0d+00 
  numr(1)= 1.0d+00 
  denomr(1)= 1.0d+00 
  cnt=sigfig
  temp= cmplx(0.0,0.0, kind = ck )
  oldtemp=temp
  ixcnt=0
  rexp=int(ibit/2)
  x=rexp*(sumr(l+1)-2)
  rr10=x*log2
  ir10=int(rr10)
  rr10=rr10-ir10
  x=rexp*(sumi(l+1)-2)
  ri10=x*log2
  ii10=int(ri10)
  ri10=ri10-ii10
  dum1=sign(sumr(1)*rmax*rmax+sumr(2)*rmax+sumr(3),sumr(-1))
  dum2=sign(sumi(1)*rmax*rmax+sumi(2)*rmax+sumi(3),sumi(-1))
  dum1=dum1*10**rr10
  dum2=dum2*10**ri10
  cdum1= cmplx(dum1,dum2, kind = ck )
  x=rexp*(denomr(l+1)-2)
  rr10=x*log2
  ir10=int(rr10)
  rr10=rr10-ir10
  x=rexp*(denomi(l+1)-2)
  ri10=x*log2
  ii10=int(ri10)
  ri10=ri10-ii10
  dum1=sign(denomr(1)*rmax*rmax+denomr(2)*rmax+denomr(3),denomr(-1))
  dum2=sign(denomi(1)*rmax*rmax+denomi(2)*rmax+denomi(3),denomi(-1))
  dum1=dum1*10**rr10
  dum2=dum2*10**ri10
  cdum2= cmplx(dum1,dum2, kind = ck )
  temp=cdum1/cdum2
 
  goon1=1
  do while (goon1==1)
 
    goon1=0
    if (ip<0) then

      if (sumr(1)< 0.5d+00 ) then
        mx1=sumi(l+1)
      elseif (sumi(1)< 0.5d+00 ) then
        mx1=sumr(l+1)
      else
        mx1=max(sumr(l+1),sumi(l+1))
      endif

      if (numr(1)< 0.5d+00 ) then
        mx2=numi(l+1)
      elseif (numi(1)< 0.5d+00 ) then
        mx2=numr(l+1)
      else
        mx2=max(numr(l+1),numi(l+1))
      endif

      if (mx1-mx2>2.0) then

        if (creal>=0.0) then
          if (abs(temp1/cnt)<= 1.0d+00 ) then
            call arydiv(sumr,sumi,denomr,denomi,final,l,lnpfq,rmax,ibit)
            hyper=final
            return
          endif
        endif
      endif

    else

      call arydiv(sumr,sumi,denomr,denomi,temp,l,lnpfq,rmax,ibit)
!
!  estimate the exponent of the maximum term in the pfq series.
!
      expon=0.0
      xl= real ( ixcnt, kind = rk )
      do i = 1, ip
        expon=expon + real (factorial_log(a(i)+xl- 1.0d+00 ), kind = rk ) &
          -real ( factorial_log(a(i)- 1.0d+00 ), kind = rk)
      enddo

      do i = 1, iq
        expon=expon - real (factorial_log(b(i)+xl- 1.0d+00 ), kind = rk ) &
          + real ( factorial_log(b(i)- 1.0d+00 ), kind = rk )
      enddo

      expon=expon+xl * real ( log(z), kind = rk ) &
        - real ( factorial_log ( cmplx(xl,0.0, kind = ck )), kind = rk )
      lmax=int(log10(exp( 1.0d+00 ))*expon)

      if (abs(oldtemp-temp)<abs(temp*accy)) then
        call arydiv(sumr,sumi,denomr,denomi,final,l,lnpfq,rmax,ibit)
        hyper=final
        return
      endif

      oldtemp = temp

    endif

    if ( ixcnt /= icount ) then
      ixcnt = ixcnt + 1
      do i1 = 1, iq
!
!  take the current sum and multiply by the denominator of the next
!  term, for both the most significant half (cr,ci) and the least
!  significant half (cr2,ci2).
!
        call cmpmul(sumr,sumi,cr(i1),ci(i1),qr1,qi1,wk1,wk2,wk3,wk4,wk5,   &
          wk6,l,rmax)
        call cmpmul(sumr,sumi,cr2(i1),ci2(i1),qr2,qi2,wk1,wk2,wk3,wk4,wk5, &
          wk6,l,rmax)
        qr2(l+1)=qr2(l+1)-1
        qi2(l+1)=qi2(l+1)-1
!
!  store this temporarily in the sum arrays.
!
        call cmpadd(qr1,qi1,qr2,qi2,sumr,sumi,wk1,l,rmax)

      enddo
!
!  multiply by the factorial term.
!
      foo1=sumr
      foo2=sumr
      call armult(foo1,cnt,foo2,wk6,l,rmax)
      sumr=foo2
      foo1=sumi
      foo2=sumi
      call armult(foo1,cnt,foo2,wk6,l,rmax)
      sumi=foo2
!
!  multiply by the scaling factor, sigfig, to keep the scale correct.
!
      do i1=1 , ip-iq
        foo1=sumr
        foo2=sumr
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        sumr=foo2
        foo1=sumi
        foo2=sumi
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        sumi=foo2
      enddo

      do i1=1 , iq
!
!  update the denominator.
!
        call cmpmul(denomr,denomi,cr(i1),ci(i1),qr1,qi1,wk1,wk2,wk3,wk4,   &
          wk5,wk6,l,rmax)
        call cmpmul(denomr,denomi,cr2(i1),ci2(i1),qr2,qi2,wk1,wk2,wk3,wk4, &
          wk5,wk6,l,rmax)
        qr2(l+1)=qr2(l+1)-1
        qi2(l+1)=qi2(l+1)-1
        call cmpadd(qr1,qi1,qr2,qi2,denomr,denomi,wk1,l,rmax)
      enddo
!
!  multiply by the factorial term.
!
      foo1=denomr
      foo2=denomr
      call armult(foo1,cnt,foo2,wk6,l,rmax)
      denomr=foo2
      foo1=denomi
      foo2=denomi
      call armult(foo1,cnt,foo2,wk6,l,rmax)
      denomi=foo2
!
!  multiply by the scaling factor, sigfig, to keep the scale correct.
!
      do i1=1 , ip-iq
        foo1=denomr
        foo2=denomr
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        denomr=foo2
        foo1=denomi
        foo2=denomi
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        denomi=foo2
      enddo
!
!  form the next numerator term by multiplying the current
!  numerator term (an array) with the a argument (a scalar).
!
      do i1=1 , ip
        call cmpmul(numr,numi,ar(i1),ai(i1),qr1,qi1,wk1,wk2,wk3,wk4,wk5,   &
          wk6,l,rmax)
        call cmpmul(numr,numi,ar2(i1),ai2(i1),qr2,qi2,wk1,wk2,wk3,wk4,wk5, &
          wk6,l,rmax)
        qr2(l+1)=qr2(l+1)-1
        qi2(l+1)=qi2(l+1)-1
        call cmpadd(qr1,qi1,qr2,qi2,numr,numi,wk1,l,rmax)
      enddo
!
!  finish the new numerator term by multiplying by the z argument.
!
      call cmpmul(numr,numi,xr,xi,qr1,qi1,wk1,wk2,wk3,wk4,wk5,wk6,l,rmax)
      call cmpmul(numr,numi,xr2,xi2,qr2,qi2,wk1,wk2,wk3,wk4,wk5,wk6,l,rmax)
      qr2(l+1)=qr2(l+1)-1
      qi2(l+1)=qi2(l+1)-1
      call cmpadd(qr1,qi1,qr2,qi2,numr,numi,wk1,l,rmax)
!
!  multiply by the scaling factor, sigfig, to keep the scale correct.
!
      do i1=1 , iq-ip
        foo1=numr
        foo2=numr
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        numr=foo2
        foo1=numi
        foo2=numi
        call armult(foo1,sigfig,foo2,wk6,l,rmax)
        numi=foo2
      enddo
!
!  Add the new numerator term with the current running
!  sum of the numerator and store the new running sum in sumr, sumi.
!
      foo1=sumr
      foo2=sumr
      bar1=sumi
      bar2=sumi

      call cmpadd ( foo1, bar1, numr, numi, foo2, bar2, wk1, l, rmax )
      sumi=bar2
      sumr=foo2
!
!  because sigfig represents "one" on the new scale, add sigfig
!  to the current count and, consequently, to the ip arguments
!  in the numerator and the iq arguments in the denominator.
!
      cnt=cnt+sigfig
      do i1=1 , ip
        ar(i1)=ar(i1)+sigfig
      enddo
      do i1=1 , iq
        cr(i1)=cr(i1)+sigfig
      enddo
      goon1=1
    endif
  end do

  call arydiv(sumr,sumi,denomr,denomi,final,l,lnpfq,rmax,ibit)
  hyper = final

  return
end

function ipremax ( a, b, ip, iq, z )

!*****************************************************************************80
!
!! ipremax() predicts the maximum term in the pFq series.  
!
!  Licensing:
!
!    This code is distributed under the MIT license.
!
!  Modified:
!
!    12 July 2022
!
!  Author:
!
!    Warren Perger, Atul Bhalla, Mark Nardin
!    Modifications by John Burkardt
!
  implicit none

  integer, parameter :: ck = kind ( ( 1.0D+00, 1.0D+00 ) )
  integer, parameter :: rk = kind ( 1.0D+00 )

  complex ( kind = ck ) a(ip)
  complex ( kind = ck ) b(iq)
  real ( kind = rk ) expon
  complex ( kind = ck ) factorial_log
  integer i
  integer ip
  integer ipremax
  integer iq
  integer j
  real ( kind = rk ) xl
  real ( kind = rk ) xmax
  real ( kind = rk ) xterm
  complex ( kind = ck ) z
!
!  estimate the exponent of the maximum term in the pfq series.
!
  xterm=0
  do j = 1, 100000 
    expon = 0.0d+00
    xl=real (j, kind = rk )
    do i=1 , ip
      expon = expon + real ( factorial_log(a(i)+xl- 1.0d+00 ), kind = rk ) &
        - real (factorial_log(a(i)- 1.0d+00 ), kind = rk )
    enddo
    do i = 1, iq
      expon=expon - real ( factorial_log(b(i)+xl- 1.0d+00 ), kind = rk ) &
        + real (factorial_log(b(i)- 1.0d+00 ), kind = rk )
    enddo
    expon=expon + xl * real ( log(z), kind = rk ) &
      - real ( factorial_log(cmplx(xl,0.0d+00,kind = ck) ), kind = rk )
    xmax=log10(exp( 1.0d+00 ))*expon
    if ((xmax<xterm) .and. (j>2)) then
      ipremax=j
      return
    endif
    xterm=max(xmax,xterm)
  enddo

  write ( *, '(a)' ) ''
  write ( *, '(a)' ) 'ipremax(): Fatal error!'
  write ( *, '(a)' ) '  Did not find maximum exponent.'
  stop 1

end

