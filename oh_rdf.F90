! Written by Colin Bunner (bunne043@umn.edu) some time in 2018
SUBROUTINE ohrdf(NWAT,BOXLEN,rmin,rfine,rmedium,rmax,nrbin_fine&
  ,nrbin_medium,nrbin_coarse,o_x,o_y,o_z,h1_x,h1_y,h1_z,h2_x&
  ,h2_y,h2_z,g_fine,g_medium,g_coarse,nr_fine,nr_medium,nr_coarse)

! Computes O-H RDF, using 3 variable-width histograms to improve statistics
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
! h1_i,    - x,y,z coordinates of H beads
! h2_i

! INOUT
!--------------------------------------------------------------------------
! g_x, - x=fine,medium,coarse
! nr_x   RDF and number integrals for histogram x          

! DUMMY
!--------------------------------------------------------------------------
! oh1_dist, - distance of ref O bead from 1st/2nd H bead on observed water
! oh2_dist
! oh1vec_i, - x,y,z coordinates of O-H separation vectors
! oh2vec_i
! o_i       - store O coords of ref molecule to save time
! x_width   - radial bin width for x histogram
! rbin      - holds computed radial bin number
! count_x   - holds number of observations of a water molecule at a
!             distance r-r+bin_width
! vol       - volume of simulation box
! nu        - running count of water molecules for number integral

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=selected_real_kind(15,307)

  INTEGER, INTENT(IN) :: NWAT, nrbin_fine, nrbin_medium, nrbin_coarse
  REAL(dp), INTENT(IN) :: BOXLEN, rmin, rfine, rmedium, rmax
  REAL(dp), DIMENSION(NWAT), INTENT(IN) :: o_x,o_y,o_z&
   ,h1_x,h1_y,h1_z,h2_x,h2_y,h2_z
  REAL(dp), DIMENSION(nrbin_fine), INTENT(INOUT) :: g_fine, nr_fine
!f2py intent(in,out) :: g_fine, nr_fine
  REAL(dp), DIMENSION(nrbin_medium), INTENT(INOUT) :: g_medium, nr_medium
!f2py intent(in,out) :: g_medium, nr_medium
  REAL(dp), DIMENSION(nrbin_coarse), INTENT(INOUT) :: g_coarse, nr_coarse
!f2py intent(in,out) :: g_coarse, nr_coarse


  REAL(dp) :: oh1_dist, oh2_dist, vol, nu
  REAL(dp) :: ohvec1_x, ohvec1_y, ohvec1_z&
   ,ohvec2_x, ohvec2_y, ohvec2_z, ox, oy, oz
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

!   O bead
    ox = o_x(i)
    oy = o_y(i)
    oz = o_z(i)

    DO j = i+1 , NWAT

!     Separation from H1
      ohvec1_x = ox-h1_x(j)
      ohvec1_y = oy-h1_y(j)
      ohvec1_z = oz-h1_z(j)

!     O-H separation 2
      ohvec2_x = ox-h2_x(j)
      ohvec2_y = oy-h2_y(j)
      ohvec2_z = oz-h2_z(j)

!     Periodic boundary conditions
      ohvec1_x = ohvec1_x - BOXLEN*nint(ohvec1_x/BOXLEN)
      ohvec1_y = ohvec1_y - BOXLEN*nint(ohvec1_y/BOXLEN)
      ohvec1_z = ohvec1_z - BOXLEN*nint(ohvec1_z/BOXLEN)

      ohvec2_x = ohvec2_x - BOXLEN*nint(ohvec2_x/BOXLEN)
      ohvec2_y = ohvec2_y - BOXLEN*nint(ohvec2_y/BOXLEN)
      ohvec2_z = ohvec2_z - BOXLEN*nint(ohvec2_z/BOXLEN)

      oh1_dist = sqrt(ohvec1_x**2+ohvec1_y**2+ohvec1_z**2)
      oh2_dist = sqrt(ohvec2_x**2+ohvec2_y**2+ohvec2_z**2)

!     Add observations to histogram if they fit within specified
!     range for H1
      IF (oh1_dist .lt. rfine .and. oh1_dist .gt. rmin) THEN
        rbin = int((oh1_dist-rmin)/fine_width) + 1
        count_fine(rbin) = count_fine(rbin) + 2
      ELSE IF (oh1_dist .lt. rmedium .and. oh1_dist .gt. rfine) THEN
        rbin = int((oh1_dist-rfine)/medium_width) + 1
        count_medium(rbin) = count_medium(rbin) + 2
      ELSE IF (oh1_dist .lt. rmax .and. oh1_dist .gt. rmedium) THEN
        rbin = int((oh1_dist-rmedium)/coarse_width) + 1
        count_coarse(rbin) = count_coarse(rbin) + 2
      END IF

!     Add observations to histogram if they fit within specified
!     range for H2
      IF (oh2_dist .lt. rfine .and. oh2_dist .gt. rmin) THEN
        rbin = int((oh2_dist-rmin)/fine_width) + 1
        count_fine(rbin) = count_fine(rbin) + 2
      ELSE IF (oh2_dist .lt. rmedium .and. oh2_dist .gt. rfine) THEN
        rbin = int((oh2_dist-rfine)/medium_width) + 1
        count_medium(rbin) = count_medium(rbin) + 2
      ELSE IF (oh2_dist .lt. rmax .and. oh2_dist .gt. rmedium) THEN
        rbin = int((oh2_dist-rmedium)/coarse_width) + 1
        count_coarse(rbin) = count_coarse(rbin) + 2
      END IF

    END DO
  END DO

  nu = 0.d0
  vol = BOXLEN*BOXLEN*BOXLEN
! Now that we are done adding observations, perform normalization by the parts
! that fluctuate (volume of bin will be accounted for after)
  DO i = 1, nrbin_fine
    nu = nu + count_fine(i)
    g_fine(i)=g_fine(i)+(count_fine(i)*(vol/(2*NWAT*NWAT)))
    nr_fine(i)=nr_fine(i)+(nu/NWAT)
  END DO
  DO i = 1, nrbin_medium
    nu = nu + count_medium(i)
    g_medium(i)=g_medium(i)+(count_medium(i)*(vol/(2*NWAT*NWAT)))
    nr_medium(i)=nr_medium(i)+(nu/NWAT)
  END DO
  DO i = 1, nrbin_coarse
    nu = nu + count_coarse(i)
    g_coarse(i)=g_coarse(i)+(count_coarse(i)*(vol/(2*NWAT*NWAT)))
    nr_coarse(i)=nr_coarse(i)+(nu/NWAT)
  END DO

END SUBROUTINE
