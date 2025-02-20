subroutine clustr ( x, d, dev, b, f, e, i, j, n, nz, k )

!*****************************************************************************80
!
!! clustr() uses the K-means algorithm to cluster data.
!
!  Discussion:
!
!    Given a matrix of I observations on J variables, the
!    observations are allocated to N clusters in such a way that the
!    within-cluster sum of squares is minimised.
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
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58:
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Input:
!
!    real ( kind = rk ) X(I,J), the observed data.
!
!    real ( kind = rk ) D(K,J), the initial cluster centers.
!
!    integer I, the number of observations.
!
!    integer J, the number of variables.
!
!    integer N, the number of clusters.
!
!    integer NZ, the minimum number of observations
!    which any cluster is allowed to have.
!
!    integer K, the maximum number of clusters.
!
!  Output:
!
!    real ( kind = rk ) D(K,J), the updated cluster centers.
!
!    real ( kind = rk ) DEV(K), the sums of squared deviations
!    of observations from their cluster centers.
!
!    integer B(I), indicates the cluster to which
!    each observation has been assigned.
!
!    integer E(K), the number of observations assigned
!    to each cluster.
!
!  Workspace:
!
!    real ( kind = rk ) F(I).
!
  implicit none

  integer, parameter :: rk = kind ( 1.0D+00 )

  integer i
  integer k

  integer b(i)
  real ( kind = rk ), parameter :: big = 1.0D+10
  real ( kind = rk ) d(k,j)
  real ( kind = rk ) da
  real ( kind = rk ) db
  real ( kind = rk ) dc
  real ( kind = rk ) de
  real ( kind = rk ) dev(k)
  integer e(k)
  real ( kind = rk ) f(i)
  real ( kind = rk ) fl
  real ( kind = rk ) fm
  real ( kind = rk ) fq
  integer ic
  integer id
  integer ie
  integer ig
  integer ii
  integer ij
  integer ik
  integer il
  integer in
  integer ip
  integer ir
  integer is
  integer it
  integer iu
  integer iw
  integer j
  integer n
  integer nz
  real ( kind = rk ) x(i,j)

  e(1:n) = 0
!
!  For each observation, calculate the distance from each cluster
!  center, and assign to the nearest.
!
  do ic = 1, i

    f(ic) = 0.0D+00
    da = big

    do id = 1, n

      db = 0.0D+00
      do ie = 1, j
        dc = x(ic,ie) - d(id,ie)
        db = db + dc * dc
      end do

      if ( db < da ) then
        da = db
        b(ic) = id
      end if

    end do

    ig = b(ic)
    e(ig) = e(ig) + 1

  end do
!
!  Calculate the mean and sum of squares for each cluster.
!
  dev(1:n) = 0.0D+00
  d(1:n,1:j) = 0.0D+00

  do ic = 1, i
    ig = b(ic)
    d(ig,1:j) = d(ig,1:j) + x(ic,1:j)
  end do

  do ij = 1, j
    do ii = 1, n
      d(ii,ij) = d(ii,ij) / real ( e(ii), kind = rk )
    end do
  end do

  do ij = 1, j
    do ik = 1, i
      il = b(ik)
      da = x(ik,ij) - d(il,ij)
      db = da * da
      f(ik) = f(ik) + db
      dev(il) = dev(il) + db
    end do
  end do

  do ik = 1, i
    il = b(ik)
    fl = e(il)
    if ( 1 < e(il) ) then
      f(ik) = f(ik) * fl / ( fl - 1.0D+00 )
    end if
  end do
!
!  Examine each observation in turn to see if it should be
!  reassigned to a different cluster.
!
  do

    iw = 0

    do ik = 1, i

      il = b(ik)
      ir = il
!
!  If the number of cluster points is less than or equal to the
!  specified minimum, NZ, then bypass this iteration.
!
      if ( nz < e(il) ) then

        fl = e(il)
        dc = f(ik)

        do in = 1, n

          if ( in /= il ) then

            fm = e(in)
            fm = fm / ( fm + 1.0D+00 )

            de = 0.0D+00
            do ip = 1, j
              da = x(ik,ip) - d(in,ip)
              de = de + da * da * fm
            end do

            if ( de < dc ) then
              dc = de
              ir = in
            end if

          end if

        end do
!
!  Reassignment is made here if necessary.
!
        if ( ir /= il ) then

          fq = e(ir)
          dev(il) = dev(il) - f(ik)
          dev(ir) = dev(ir) + dc
          e(ir) = e(ir) + 1
          e(il) = e(il) - 1

          do is = 1, j
            d(il,is) = ( d(il,is) * fl - x(ik,is) ) / ( fl - 1.0D+00 )
            d(ir,is) = ( d(ir,is) * fq + x(ik,is) ) / ( fq + 1.0D+00 )
          end do

          b(ik) = ir

          do it = 1, i

            ij = b(it)

            if ( ij == il .or. ij == ir ) then
              f(it) = 0.0D+00
              do iu = 1, j
                da = x(it,iu) - d(ij,iu)
                f(it) = f(it) + da * da
              end do
              fl = e(ij)
              f(it) = f(it) * fl / ( fl - 1.0D+00 )
            end if

          end do

          iw = iw + 1

        end if

      end if

    end do
!
!  If any reassignments were made on this pass, then do another pass.
!
    if ( iw == 0 ) then
      exit
    end if

  end do

  return
end

