! Written by Colin Bunner (bunne043@umn.edu) some time in 2018
SUBROUTINE adf(NWAT,BOXLEN,rmin,rmax,nrbin,h1_x,h1_y,h1_z,h2_x,h2_y,h2_z, &
  &o_x,o_y,o_z,m_x,m_y,m_z,hist_g,hist_d,hist_delt)

! Pair correlation function projected onto:
! - dipolar interaction function
! D = 3(mu1^ dot R^)(mu2^ dot R^) - mu1^ dot mu2^,
! - dipolar orientation function
! delta = mu1^ * mu2^
! - identity function 
! I = 1

!!! Variables !!!

! IN
!----------------------------------------------------------
! NWAT       - number of H2O molecules in frame
! nrbin      - number of radial bins for histograms
! BOXLEN     - length of simulation box (cubic implied)
! rmin, rmax - min. and max. distances for histograms
! o_i,m_i,   - x,y,z coordinates of oxygen beads, m sites,
! h1_i,h2_i    and hydrogen beads of TIP4P/2005 water

! INOUT
!----------------------------------------------------------
! hist_g    - center-of-mass RDF
! hist_d    - ensemble average of D
! hist_delt - ensemble average of delta

! DUMMY
!----------------------------------------------------------
! cm_dist       - COM distance between observed water molecules
! dd            - dot product of dipole moment vectors
! d1r,          - dot product of first/second dipole moment with COM
! d2r             with COM separation vector
! dip1_i,       - x,y,z components of first/second dipole moment vector 
! dip2_i
! cm1_i,        - x,y,z coordinats of COM of first/second water molecule
! cm2_i
! cmvec_i       - x,y,z coordinates of COM separation vector
! rbin_width    - width of histogram radial bins
! rbin          - holds calculated bin for observation
! hist_d_sf,    - single-frame histograms for computing frame averages
! hist_delt_sf,
! hist_count

  IMPLICIT NONE
  REAL*8, PARAMETER :: PI=3.14159265359
  REAL*8, PARAMETER :: RTD=180./3.14159265359
  INTEGER :: i, j
  INTEGER, INTENT(IN) :: NWAT, nrbin
  REAL*8, INTENT(IN) :: BOXLEN, rmin, rmax
  REAL*8, DIMENSION(NWAT), INTENT(IN) :: o_x,o_y,o_z, &
   &m_x,m_y,m_z, h1_x,h1_y,h1_z,h2_x,h2_y,h2_z
  REAL*8, DIMENSION(nrbin), INTENT(INOUT) :: hist_g, hist_d, hist_delt
!f2py intent(in,out) :: hist_g,hist_d,hist_delt


  REAL*8 :: cm_dist, dd, d1r, d2r
  REAL*8 :: dip1_x, dip1_y, dip1_z, dip2_x, dip2_y, &
   &dip2_z, cmvec_x, cmvec_y, cmvec_z, cm1_x, cm1_y, cm1_z, &
   &cm2_x, cm2_y, cm2_z
  REAL*8 :: rbin_width
  INTEGER :: rbin
  REAL*8, DIMENSION(nrbin) :: hist_d_sf, hist_delt_sf, hist_count

! Set all values in histogram to zero
  DO i = 1, nrbin
    hist_d_sf(i) = 0.d0
    hist_delt_sf(i) = 0.d0
    hist_count(i) = 0
  END DO

! Compute bin widths 
  rbin_width = (rmax-rmin)/nrbin

! Loop over unique pairs of molecules 
  DO i = 1, NWAT-1

!   Dipole moment points along OM vector
    dip1_x = (m_x(i)-o_x(i))/0.1546
    dip1_y = (m_y(i)-o_y(i))/0.1546
    dip1_z = (m_z(i)-o_z(i))/0.1546

!   COM
    cm1_x = (h1_x(i)*1.0079+h2_x(i)*1.0079+o_x(i)*15.999)/18.0148
    cm1_y = (h1_y(i)*1.0079+h2_y(i)*1.0079+o_y(i)*15.999)/18.0148
    cm1_z = (h1_z(i)*1.0079+h2_z(i)*1.0079+o_z(i)*15.999)/18.0148

    DO j = i+1 , NWAT

!     Dipole moment points along OM vector
!     0.1546 is the length of the OM vector
      dip2_x = (m_x(j)-o_x(j))/0.1546
      dip2_y = (m_y(j)-o_y(j))/0.1546
      dip2_z = (m_z(j)-o_z(j))/0.1546


!     COM of observed water 
      cm2_x = (h1_x(j)*1.0079+h2_x(j)*1.0079+o_x(j)*15.999)/18.0148
      cm2_y = (h1_y(j)*1.0079+h2_y(j)*1.0079+o_y(j)*15.999)/18.0148
      cm2_z = (h1_z(j)*1.0079+h2_z(j)*1.0079+o_z(j)*15.999)/18.0148
             
!     COM separation
      cmvec_x = cm2_x-cm1_x
      cmvec_y = cm2_y-cm1_y
      cmvec_z = cm2_z-cm1_z

!     Periodic boundary conditions
      cmvec_x = cmvec_x - BOXLEN*nint(cmvec_x/boxlen)
      cmvec_y = cmvec_y - BOXLEN*nint(cmvec_y/boxlen)
      cmvec_z = cmvec_z - BOXLEN*nint(cmvec_z/boxlen)

      cm_dist = sqrt(cmvec_x**2+cmvec_y**2+cmvec_z**2)

!     Make COM vector into unit vector
      cmvec_x = cmvec_x/cm_dist
      cmvec_y = cmvec_y/cm_dist
      cmvec_z = cmvec_z/cm_dist

      dd= dip1_x*dip2_x+dip1_y*dip2_y+dip1_z*dip2_z
      d1r= dip1_x*cmvec_x+dip1_y*cmvec_y+dip1_z*cmvec_z
      d2r= dip2_x*cmvec_x+dip2_y*cmvec_y+dip2_z*cmvec_z

            
!     Add observations to histogram if they fit within specified
!     range
      IF (cm_dist .lt. rmax .and. cm_dist .gt. rmin) THEN

        rbin = int((cm_dist-rmin)/rbin_width) + 1

        hist_count(rbin) = hist_count(rbin) + 1
        hist_d_sf(rbin) = hist_d_sf(rbin) + 3*d1r*d2r - dd
        hist_delt_sf(rbin) = hist_delt_sf(rbin) + dd

      END IF

    END DO
  END DO

! Now that we are done adding observations, we compute the average expansion
! coefficient
  DO i = 1, nrbin
    IF (hist_count(i) .gt. 0) THEN
      ! Divide by hist_count(i) because we want average
      hist_d(i) = hist_d(i) + (hist_d_sf(i)/hist_count(i))
      hist_delt(i) = hist_delt(i) + (hist_delt_sf(i)/hist_count(i))
      ! Factor of 2 because we calculated sum over unique pairs (i,j>1),
      ! so each observation for RDF must be doubled by symmetry
      hist_g(i)=hist_g(i)+(hist_count(i)*2*((BOXLEN*BOXLEN*BOXLEN)/(NWAT*NWAT)))
    END IF
  END DO

END SUBROUTINE
