! Written by Colin Bunner (bunne043@umn.edu) some time in 2018
SUBROUTINE ordf(NWAT,BOXLEN,rmin,rfine,rmedium,rmax,nrbin_fine,nrbin_medium,nrbin_coarse&
  ,o_x,o_y,o_z,g_wat_fine,g_wat_medium,g_wat_coarse,nr_wat_fine&
  ,nr_wat_medium,nr_wat_coarse)

! Computes O-O RDF, using 3 variable-width histograms to improve statistics
! in H2O-poor phase

!!! Variables !!!

! PARAMETER
!--------------------------------------------------------------------------
! dp - set double precision

! IN
!--------------------------------------------------------------------------
! NWAT     - number of H2O molecules in frame
! BOXLEN   - length of simulation box (Angstrom). Cubic implied.
! nrbin_x  - x=fine,medium,coarse
!            number of radial bins for fine,medium,coarse portion of 
!            RDF histogram
! rmin     - min radial distances allowed for histogram
! rfine    - cutoff for fine portion of histogram
! rmedium  - cutoff for medium portion of histogram
! rmax     - maximum cutoff for histogram
! o_i      - x,y,z coordinates of O beads

! INOUT
!--------------------------------------------------------------------------
! g_wat_x, - x=fine,medium,coarse
! nr_wat_x   RDF and number integrals for histogram x          

! DUMMY
!--------------------------------------------------------------------------
! cm_dist  - O distance
! cmvec_i  - x,y,z coordinates of O-O separation vector
! cm1_i,   - x,y,z coordinates of O-O of 1st/2nd molecule
! cm2_i
! x_width  - radial bin width for x histogram
! rbin     - holds computed radial bin number
! count_x  - holds number of observations of a water molecule at a
!            distance r-r+bin_width
! vol      - volume of simulation box
! nuwat    - running count of water molecules encountered
  

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=selected_real_kind(15,307)

  INTEGER, INTENT(IN) :: NWAT, nrbin_fine, nrbin_medium, nrbin_coarse
  REAL(dp), INTENT(IN) :: BOXLEN, rmin, rfine, rmedium, rmax
  REAL(dp), DIMENSION(NWAT), INTENT(IN) :: o_x,o_y,o_z
  REAL(dp), DIMENSION(nrbin_fine), INTENT(INOUT) :: g_wat_fine, nr_wat_fine
!f2py intent(in,out) :: g_wat_fine, nr_wat_fine
  REAL(dp), DIMENSION(nrbin_medium), INTENT(INOUT) :: g_wat_medium, nr_wat_medium
!f2py intent(in,out) :: g_wat_medium, nr_wat_medium
  REAL(dp), DIMENSION(nrbin_coarse), INTENT(INOUT) :: g_wat_coarse, nr_wat_coarse
!f2py intent(in,out) :: g_wat_coarse, nr_wat_coarse


  REAL(dp) :: cm_dist, vol, nuwat
  REAL(dp) :: cmvec_x, cmvec_y, cmvec_z, cm1_x, cm1_y, cm1_z, &
   &cm2_x, cm2_y, cm2_z
  REAL(dp) :: fine_width, medium_width, coarse_width
  INTEGER :: rbin
  INTEGER, DIMENSION(nrbin_fine) :: count_fine
  INTEGER, DIMENSION(nrbin_medium) :: count_medium
  INTEGER, DIMENSION(nrbin_coarse) :: count_coarse
  INTEGER :: i, j

  ! single-frame g(r) counters
  count_fine = 0
  count_medium = 0
  count_coarse = 0


! Compute bin widths 
  fine_width = (rfine-rmin)/nrbin_fine
  medium_width = (rmedium-rfine)/nrbin_medium
  coarse_width = (rmax-rmedium)/nrbin_coarse

! Loop over unique pairs of waters
  DO i = 1, NWAT-1

!   O bead (yes, I was too lazy to rename the cm part)
    cm1_x = o_x(i)
    cm1_y = o_y(i)
    cm1_z = o_z(i)

    DO j = i+1 , NWAT

!     O of observed water 
      cm2_x = o_x(j)
      cm2_y = o_y(j)
      cm2_z = o_z(j)
             
!     O separation
      cmvec_x = cm2_x-cm1_x
      cmvec_y = cm2_y-cm1_y
      cmvec_z = cm2_z-cm1_z

!     Periodic boundary conditions
      cmvec_x = cmvec_x - BOXLEN*nint(cmvec_x/boxlen)
      cmvec_y = cmvec_y - BOXLEN*nint(cmvec_y/boxlen)
      cmvec_z = cmvec_z - BOXLEN*nint(cmvec_z/boxlen)

      cm_dist = sqrt(cmvec_x**2+cmvec_y**2+cmvec_z**2)

!     Add observations to histogram if they fit within specified
!     range
      IF (cm_dist .lt. rfine .and. cm_dist .gt. rmin) THEN
        rbin = int((cm_dist-rmin)/fine_width) + 1
        count_fine(rbin) = count_fine(rbin) + 2
      ELSE IF (cm_dist .lt. rmedium .and. cm_dist .gt. rfine) THEN
        rbin = int((cm_dist-rfine)/medium_width) + 1
        count_medium(rbin) = count_medium(rbin) + 2
      ELSE IF (cm_dist .lt. rmax .and. cm_dist .gt. rmedium) THEN
        rbin = int((cm_dist-rmedium)/coarse_width) + 1
        count_coarse(rbin) = count_coarse(rbin) + 2
      END IF

    END DO
  END DO

  nuwat = 0.d0
  vol = BOXLEN*BOXLEN*BOXLEN
! Now that we are done adding observations, perform normalization by the parts
! that fluctuate (volume of bin will be accounted for after)
  DO i = 1, nrbin_fine
    nuwat = nuwat + count_fine(i)
    g_wat_fine(i)=g_wat_fine(i)+(count_fine(i)*(vol/(NWAT*NWAT)))
    nr_wat_fine(i)=nr_wat_fine(i)+(nuwat/NWAT)
  END DO
  DO i = 1, nrbin_medium
    nuwat = nuwat + count_medium(i)
    g_wat_medium(i)=g_wat_medium(i)+(count_medium(i)*(vol/(NWAT*NWAT)))
    nr_wat_medium(i)=nr_wat_medium(i)+(nuwat/NWAT)
  END DO
  DO i = 1, nrbin_coarse
    nuwat = nuwat + count_coarse(i)
    g_wat_coarse(i)=g_wat_coarse(i)+(count_coarse(i)*(vol/(NWAT*NWAT)))
    nr_wat_coarse(i)=nr_wat_coarse(i)+(nuwat/NWAT)
  END DO

END SUBROUTINE
