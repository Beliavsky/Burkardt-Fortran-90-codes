subroutine kmns ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, live, &
  iter, wss, ifault )

!*****************************************************************************80
!
!! kmns() carries out the K-means algorithm.
!
!  Discussion:
!
!    This routine attempts to divide M points in N-dimensional space into
!    K clusters so that the within cluster sum of squares is minimized.
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
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Input:
!
!    real ( kind = rk ) A(M,N), the points.
!
!    integer M, the number of points.
!
!    integer N, the number of spatial dimensions.
!
!    real ( kind = rk ) C(K,N), the initial cluster centers.
!
!    integer K, the number of clusters.
!
!    integer ITER, the maximum number of iterations allowed.
!
!  Output:
!
!    real ( kind = rk ) C(K,N), the updated cluster centers.
!
!    integer IC1(M), the cluster to which each point
!    is assigned.
!
!    integer NC(K), the number of points in each cluster.
!
!    real ( kind = rk ) WSS(K), the within-cluster sum of squares
!    of each cluster.
!
!    integer IFAULT, error indicator.
!    0, no error was detected.
!    1, at least one cluster is empty after the initial assignment.  A better
!       set of initial cluster centers is needed.
!    2, the allowed maximum number off iterations was exceeded.
!    3, K is less than or equal to 1, or greater than or equal to M.
!
!  Workspace:
!
!    integer IC2(M), used to store the cluster which
!    each point is most likely to be transferred to at each step.
!
!    real ( kind = rk ) AN1(K).
!
!    real ( kind = rk ) AN2(K).
!
!    integer NCP(K).
!
!    real ( kind = rk ) D(M).
!
!    integer ITRAN(K).
!
!    integer LIVE(K).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) aa
  real ( kind = rk ) an1(k)
  real ( kind = rk ) an2(k)
  real ( kind = rk ) c(k,n)
  real ( kind = rk ) d(m)
  real ( kind = rk ) da
  real ( kind = rk ) db
  real ( kind = rk ) dc
  real ( kind = rk ) dt(2)
  integer i
  integer ic1(m)
  integer ic2(m)
  integer ifault
  integer ii
  integer ij
  integer il
  integer indx
  integer iter
  integer itran(k)
  integer j
  integer l
  integer live(k)
  integer nc(k)
  integer ncp(k)
  real ( kind = rk ) temp
  real ( kind = rk ) wss(k)

  ifault = 0

  if ( k <= 1 .or. m <= k ) then
    ifault = 3
    return
  end if
!
!  For each point I, find its two closest centers, IC1(I) and
!  IC2(I).  Assign the point to IC1(I).
!
  do i = 1, m

    ic1(i) = 1
    ic2(i) = 2

    do il = 1, 2
      dt(il) = 0.0D+00
      do j = 1, n
        da = a(i,j) - c(il,j)
        dt(il) = dt(il) + da * da
      end do
    end do

    if ( dt(2) < dt(1) ) then
      ic1(i) = 2
      ic2(i) = 1
      temp = dt(1)
      dt(1) = dt(2)
      dt(2) = temp
    end if

    do l = 3, k

      db = 0.0D+00
      do j = 1, n
        dc = a(i,j) - c(l,j)
        db = db + dc * dc
      end do

      if ( db < dt(2) ) then

        if ( dt(1) <= db ) then
          dt(2) = db
          ic2(i) = l
        else
          dt(2) = dt(1)
          ic2(i) = ic1(i)
          dt(1) = db
          ic1(i) = l
        end if

      end if

    end do

  end do
!
!  Update cluster centers to be the average of points contained within them.
!
  do l = 1, k
    nc(l) = 0
    do j = 1, n
      c(l,j) = 0.0D+00
    end do
  end do

  do i = 1, m
    l = ic1(i)
    nc(l) = nc(l) + 1
    c(l,1:n) = c(l,1:n) + a(i,1:n)
  end do
!
!  Check to see if there is any empty cluster at this stage.
!
  ifault = 1

  do l = 1, k

    if ( nc(l) == 0 ) then
      ifault = 1
      return
    end if

  end do

  ifault = 0

  do l = 1, k

    aa = real ( nc(l), kind = rk )
    c(l,1:n) = c(l,1:n) / aa
!
!  Initialize AN1, AN2, ITRAN and NCP.
!
!  AN1(L) = NC(L) / (NC(L) - 1)
!  AN2(L) = NC(L) / (NC(L) + 1)
!  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
!           = 0 otherwise
!
!  In the optimal-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated.
!
!  In the quick-transfer stage, NCP(L) stores the step at which
!  cluster L is last updated plus M.
!
    an2(l) = aa / ( aa + 1.0D+00 )

    if ( 1.0D+00 < aa ) then
      an1(l) = aa / ( aa - 1.0D+00 )
    else
      an1(l) = huge ( an1(l) )
    end if

    itran(l) = 1
    ncp(l) = -1

  end do

  indx = 0
  ifault = 2

  do ij = 1, iter
!
!  In this stage, there is only one pass through the data.   Each
!  point is re-allocated, if necessary, to the cluster that will
!  induce the maximum reduction in within-cluster sum of squares.
!
    call optra ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
      live, indx )
!
!  Stop if no transfer took place in the last M optimal transfer steps.
!
    if ( indx == m ) then
      ifault = 0
      exit
    end if
!
!  Each point is tested in turn to see if it should be re-allocated
!  to the cluster to which it is most likely to be transferred,
!  IC2(I), from its present cluster, IC1(I).   Loop through the
!  data until no further change is to take place.
!
    call qtran ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, &
      itran, indx )
!
!  If there are only two clusters, there is no need to re-enter the
!  optimal transfer stage.
!
    if ( k == 2 ) then
      ifault = 0
      exit
    end if
!
!  NCP has to be set to 0 before entering OPTRA.
!
    do l = 1, k
      ncp(l) = 0
    end do

  end do
!
!  If the maximum number of iterations was taken without convergence,
!  IFAULT is 2 now.  This may indicate unforeseen looping.
!
  if ( ifault == 2 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'KMNS - Warning!'
    write ( *, '(a)' ) '  Maximum number of iterations reached'
    write ( *, '(a)' ) '  without convergence.'
  end if
!
!  Compute the within-cluster sum of squares for each cluster.
!
  wss(1:k) = 0.0D+00
  c(1:k,1:n) = 0.0D+00

  do i = 1, m
    ii = ic1(i)
    c(ii,1:n) = c(ii,1:n) + a(i,1:n)
  end do

  do j = 1, n
    do l = 1, k
      c(l,j) = c(l,j) / real ( nc(l), kind = rk )
    end do
    do i = 1, m
      ii = ic1(i)
      da = a(i,j) - c(ii,j)
      wss(ii) = wss(ii) + da * da
    end do
  end do

  return
end
subroutine optra ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
  live, indx )

!*****************************************************************************80
!
!! optra() carries out the optimal transfer stage.
!
!  Discussion:
!
!    This is the optimal transfer stage.
!
!    Each point is re-allocated, if necessary, to the cluster that
!    will induce a maximum reduction in the within-cluster sum of
!    squares.
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
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Input:
!
!    real ( kind = rk ) A(M,N), the points.
!
!    integer M, the number of points.
!
!    integer N, the number of spatial dimensions.
!
!    real ( kind = rk ) C(K,N), the cluster centers.
!
!    integer K, the number of clusters.
!
!    integer IC1(M), the cluster to which each point is assigned.
!
!    integer IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    integer NC(K), the number of points in each cluster.
!
!    real ( kind = rk ) AN1(K).
!
!    real ( kind = rk ) AN2(K).
!
!    integer NCP(K).
!
!    real ( kind = rk ) D(M).
!
!    integer ITRAN(K).
!
!    integer LIVE(K).
!
!    integer INDX, the number of steps since a transfer took place.
!
!  Output:
!
!    real ( kind = rk ) C(K,N), the cluster centers.
!
!    integer IC1(M), the cluster to which each point is assigned.
!
!    integer IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    integer NC(K), the number of points in each cluster.
!
!    real ( kind = rk ) AN1(K).
!
!    real ( kind = rk ) AN2(K).
!
!    integer NCP(K).
!
!    real ( kind = rk ) D(M).
!
!    integer ITRAN(K).
!
!    integer LIVE(K).
!
!    integer INDX, the number of steps since a transfer took place.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) al1
  real ( kind = rk ) al2
  real ( kind = rk ) alt
  real ( kind = rk ) alw
  real ( kind = rk ) an1(k)
  real ( kind = rk ) an2(k)
  real ( kind = rk ) c(k,n)
  real ( kind = rk ) d(m)
  real ( kind = rk ) da
  real ( kind = rk ) db
  real ( kind = rk ) dc
  real ( kind = rk ) dd
  real ( kind = rk ) de
  real ( kind = rk ) df
  integer i
  integer ic1(m)
  integer ic2(m)
  integer indx
  integer itran(k)
  integer j
  integer l
  integer l1
  integer l2
  integer live(k)
  integer ll
  integer nc(k)
  integer ncp(k)
  real ( kind = rk ) r2
  real ( kind = rk ) rr
!
!  If cluster L is updated in the last quick-transfer stage, it
!  belongs to the live set throughout this stage.   Otherwise, at
!  each step, it is not in the live set if it has not been updated
!  in the last M optimal transfer steps.
!
  do l = 1, k
    if ( itran(l) == 1) then
      live(l) = m + 1
    end if
  end do

  do i = 1, m

    indx = indx + 1
    l1 = ic1(i)
    l2 = ic2(i)
    ll = l2
!
!  If point I is the only member of cluster L1, no transfer.
!
    if ( 1 < nc(l1)  ) then
!
!  If L1 has not yet been updated in this stage, no need to
!  re-compute D(I).
!
      if ( ncp(l1) /= 0 ) then
        de = 0.0D+00
        do j = 1, n
          df = a(i,j) - c(l1,j)
          de = de + df * df
        end do
        d(i) = de * an1(l1)
      end if
!
!  Find the cluster with minimum R2.
!
     da = 0.0D+00
      do j = 1, n
        db = a(i,j) - c(l2,j)
        da = da + db * db
      end do
      r2 = da * an2(l2)

      do l = 1, k
!
!  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
!  true, we only need to consider clusters that are in the live set
!  for possible transfer of point I.   Otherwise, we need to consider
!  all possible clusters.
!
        if ( ( i < live(l1) .or. i < live(l2) ) .and. &
               l /= l1 .and. l /= ll ) then

          rr = r2 / an2(l)

          dc = 0.0D+00
          do j = 1, n
            dd = a(i,j) - c(l,j)
            dc = dc + dd * dd
          end do

          if ( dc < rr ) then
            r2 = dc * an2(l)
            l2 = l
          end if

        end if

      end do
!
!  If no transfer is necessary, L2 is the new IC2(I).
!
      if ( d(i) <= r2 ) then

        ic2(i) = l2
!
!  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
!  L2, and update IC1(I) and IC2(I).
!
      else

        indx = 0
        live(l1) = m + i
        live(l2) = m + i
        ncp(l1) = i
        ncp(l2) = i
        al1 = real ( nc(l1), kind = rk )
        alw = al1 - 1.0D+00
        al2 = real ( nc(l2), kind = rk )
        alt = al2 + 1.0D+00
        do j = 1, n
          c(l1,j) = ( c(l1,j) * al1 - a(i,j) ) / alw
          c(l2,j) = ( c(l2,j) * al2 + a(i,j) ) / alt
        end do
        nc(l1) = nc(l1) - 1
        nc(l2) = nc(l2) + 1
        an2(l1) = alw / al1
        if ( 1.0D+00 < alw ) then
          an1(l1) = alw / ( alw - 1.0D+00 )
        else
          an1(l1) = huge ( an1(l1) )
        end if
        an1(l2) = alt / al2
        an2(l2) = alt / ( alt + 1.0D+00 )
        ic1(i) = l2
        ic2(i) = l1

      end if

    end if

    if ( indx == m ) then
      return
    end if

  end do
!
!  ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
!  decreased by M before re-entering OPTRA.
!
  do l = 1, k
    itran(l) = 0
    live(l) = live(l) - m
  end do

  return
end
subroutine qtran ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d, itran, &
  indx )

!*****************************************************************************80
!
!! qtran() carries out the quick transfer stage.
!
!  Discussion:
!
!    This is the quick transfer stage.
!
!    IC1(I) is the cluster which point I belongs to.
!    IC2(I) is the cluster which point I is most likely to be
!    transferred to.
!
!    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
!    reduce within-cluster sum of squares.  The cluster centers are
!    updated after each step.
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
!    Original FORTRAN77 version by John Hartigan, Manchek Wong.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    John Hartigan, Manchek Wong,
!    Algorithm AS 136:
!    A K-Means Clustering Algorithm,
!    Applied Statistics,
!    Volume 28, Number 1, 1979, pages 100-108.
!
!  Input:
!
!    real ( kind = rk ) A(M,N), the points.
!
!    integer M, the number of points.
!
!    integer N, the number of spatial dimensions.
!
!    real ( kind = rk ) C(K,N), the cluster centers.
!
!    integer K, the number of clusters.
!
!    integer IC1(M), the cluster to which each
!    point is assigned.
!
!    integer IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    integer NC(K), the number of points in each cluster.
!
!    real ( kind = rk ) AN1(K).
!
!    real ( kind = rk ) AN2(K).
!
!    integer NCP(K).
!
!    real ( kind = rk ) D(M).
!
!    integer ITRAN(K).
!
!    integer INDX, counts the number of steps since the last transfer.
!
!  Output:
!
!    real ( kind = rk ) C(K,N), the cluster centers.
!
!    integer IC1(M), the cluster to which each
!    point is assigned.
!
!    integer IC2(M), used to store the cluster
!    which each point is most likely to be transferred to at each step.
!
!    integer NC(K), the number of points in each cluster.
!
!    real ( kind = rk ) AN1(K).
!
!    real ( kind = rk ) AN2(K).
!
!    integer NCP(K).
!
!    real ( kind = rk ) D(M).
!
!    integer ITRAN(K).
!
!    integer INDX, counts the number of steps since the last transfer.
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer k
  integer m
  integer n

  real ( kind = rk ) a(m,n)
  real ( kind = rk ) al1
  real ( kind = rk ) al2
  real ( kind = rk ) alt
  real ( kind = rk ) alw
  real ( kind = rk ) an1(k)
  real ( kind = rk ) an2(k)
  real ( kind = rk ) c(k,n)
  real ( kind = rk ) d(m)
  real ( kind = rk ) da
  real ( kind = rk ) db
  real ( kind = rk ) dd
  real ( kind = rk ) de
  integer i
  integer ic1(m)
  integer ic2(m)
  integer icoun
  integer indx
  integer istep
  integer itran(k)
  integer j
  integer l1
  integer l2
  integer nc(k)
  integer ncp(k)
  real ( kind = rk ) r2
!
!  In the optimal transfer stage, NCP(L) indicates the step at which
!  cluster L is last updated.   In the quick transfer stage, NCP(L)
!  is equal to the step at which cluster L is last updated plus M.
!
  icoun = 0
  istep = 0

  do

    do i = 1, m

      icoun = icoun + 1
      istep = istep + 1
      l1 = ic1(i)
      l2 = ic2(i)
!
!  If point I is the only member of cluster L1, no transfer.
!
      if ( 1 < nc(l1) ) then
!
!  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
!  cluster L1.   Note that if cluster L1 is last updated exactly M
!  steps ago, we still need to compute the distance from point I to
!  cluster L1.
!
        if ( istep <= ncp(l1) ) then

          da = 0.0D+00
          do j = 1, n
            db = a(i,j) - c(l1,j)
            da = da + db * db
          end do

          d(i) = da * an1(l1)

        end if
!
!  If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of
!  point I at this step.
!
        if ( istep < ncp(l1) .or. istep < ncp(l2) ) then

          r2 = d(i) / an2(l2)

          dd = 0.0D+00
          do j = 1, n
            de = a(i,j) - c(l2,j)
            dd = dd + de * de
          end do
!
!  Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters
!  L1 and L2.   Also update IC1(I) and IC2(I).   Note that if any
!  updating occurs in this stage, INDX is set back to 0.
!
          if ( dd < r2 ) then

            icoun = 0
            indx = 0
            itran(l1) = 1
            itran(l2) = 1
            ncp(l1) = istep + m
            ncp(l2) = istep + m
            al1 = real ( nc(l1), kind = rk )
            alw = al1 - 1.0D+00
            al2 = real ( nc(l2), kind = rk )
            alt = al2 + 1.0D+00
            do j = 1, n
              c(l1,j) = ( c(l1,j) * al1 - a(i,j) ) / alw
              c(l2,j) = ( c(l2,j) * al2 + a(i,j) ) / alt
            end do
            nc(l1) = nc(l1) - 1
            nc(l2) = nc(l2) + 1
            an2(l1) = alw / al1
            if ( 1.0D+00 < alw ) then
              an1(l1) = alw / ( alw - 1.0D+00 )
            else
              an1(l1) = huge ( an1(l1) )
            end if
            an1(l2) = alt / al2
            an2(l2) = alt / ( alt + 1.0D+00 )
            ic1(i) = l2
            ic2(i) = l1

          end if

        end if

      end if
!
!  If no re-allocation took place in the last M steps, return.
!
      if ( icoun == m ) then
        return
      end if

    end do

  end do

end

