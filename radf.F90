! Written by Colin Bunner (bunne043@umn.edu) some time in 2018
SUBROUTINE RADF(NWAT,BOXLEN,rmin,rmax,nrbin,amin,amax,nabin,h1_x,&
 &h1_y,h1_z,h2_x,h2_y,h2_z,o_x,o_y,o_z,hist_ga)

!!! Variables !!!

! PARAMETER
!------------------------------------------------------------------------
! RTD - conversion from radians to degrees (rad*RTD ==> degrees)

! IN
!------------------------------------------------------------------------
! NWAT          - number of H2O molecules
! nrbin,        - number of radial/angular bins for 2-d RADF histogram
! nabin
! BOXLEN        - length of simulation box (Angstrom). Cubic implied.
! rmin,rmax,    - min/max radial/angular distances allowed for histogram.
! amin,amax
! o_i,h1_i,h2_i - x,y,z coordinates of o/h1/h2 beads of water model

! INOUT
!------------------------------------------------------------------------
! hist_ga -  2-d RADF histogram

! DUMMY
!------------------------------------------------------------------------
! oo_dist       - O-O distance
! a1,a2,a3,a4   - 1st/2nd/3rd/4th O-O-H angle
! oref_i        - speed ups script by storing x,y,z coords of reference
!                 water molecule
! osep_i        - x,y,z coords of O-O separation vector
! oh1vec1_i,    - two O-H vectors for reference water molecule
! oh2vec1_i,
! oh1vec2_i,    - two O-H vectors for observed water molecule
! oh2vec2_i,
! rbin_width,   - histogram radial/angular bin widths
! abin_width
! rbin, abin    - store computed radial/angular bin number
! hist_ga_sf    - single-frame RADF
! Dens_Norm_Fac - Normalization factor for RADF based on number density,
!                 box volume


  IMPLICIT NONE
  REAL*8, PARAMETER :: RTD=180./3.14159265359

  INTEGER, INTENT(IN) :: NWAT, nrbin, nabin
  REAL*8, INTENT(IN) :: BOXLEN, rmin, rmax, amin, amax
  REAL*8, DIMENSION(NWAT), INTENT(IN) :: o_x,o_y,o_z,h1_x,h1_y, &
   &h1_z,h2_x,h2_y,h2_z
  REAL*8, DIMENSION(nrbin,nabin), INTENT(INOUT) :: hist_ga
!f2py intent(in,out) :: hist_ga

  REAL*8 :: oo_dist, a1, a2, a3, a4, oref_x, oref_y,oref_z, osep_x,&
   &osep_y, osep_z
  REAL*8 :: oh1vec1_x, oh1vec1_y, oh1vec1_z, oh2vec1_x, oh2vec1_y, &
   &oh2vec1_z, oh1vec1_mag, oh2vec1_mag
  REAL*8 :: oh1vec2_x, oh1vec2_y, oh1vec2_z, oh2vec2_x, oh2vec2_y, &
   &oh2vec2_z, oh1vec2_mag, oh2vec2_mag
  REAL*8 :: rbin_width, abin_width
  INTEGER :: rbin, abin
  REAL*8, DIMENSION(nrbin,nabin) :: hist_ga_sf
  DOUBLE PRECISION :: Dens_Norm_Fac
  INTEGER :: i, j

  ! Set all values in histogram to zero
  hist_ga_sf(i,j) = 0.d0

  ! Compute bin widths 
  rbin_width = (rmax-rmin)/nrbin
  abin_width = (amax-amin)/nabin

  ! Loop over unique pairs of molecules 
  DO i = 1, NWAT-1

    oref_x = o_x(i)
    oref_y = o_y(i)
    oref_z = o_z(i)

  ! Compute OH vectors for reference water, then normalize
    oh1vec1_x = h1_x(i)-oref_x
    oh1vec1_y = h1_y(i)-oref_y
    oh1vec1_z = h1_z(i)-oref_z
    oh2vec1_x = h2_x(i)-oref_x
    oh2vec1_y = h2_y(i)-oref_y
    oh2vec1_z = h2_z(i)-oref_z

    oh1vec1_mag = sqrt(oh1vec1_x**2+oh1vec1_y**2+oh1vec1_z**2)
    oh2vec1_mag = sqrt(oh2vec1_x**2+oh2vec1_y**2+oh2vec1_z**2)
 
    oh1vec1_x = oh1vec1_x/oh1vec1_mag
    oh1vec1_y = oh1vec1_y/oh1vec1_mag
    oh1vec1_z = oh1vec1_z/oh1vec1_mag
    oh2vec1_x = oh2vec1_x/oh2vec1_mag
    oh2vec1_y = oh2vec1_y/oh2vec1_mag
    oh2vec1_z = oh2vec1_z/oh2vec1_mag

    DO j = i+1 , NWAT

!     OO Vector between reference and observed waters
      osep_x = o_x(j)-oref_x
      osep_y = o_y(j)-oref_y
      osep_z = o_z(j)-oref_z
          
!     Periodic boundary conditions
      osep_x = osep_x - BOXLEN*nint(osep_x/boxlen)
      osep_y = osep_y - BOXLEN*nint(osep_y/boxlen)
      osep_z = osep_z - BOXLEN*nint(osep_z/boxlen)

      oo_dist = sqrt(osep_x**2+osep_y**2+osep_z**2)

!     Compute OH vectors for observed water, then normalize
      oh1vec2_x = h1_x(j)-o_x(j)
      oh1vec2_y = h1_y(j)-o_y(j)
      oh1vec2_z = h1_z(j)-o_z(j)
      oh2vec2_x = h2_x(j)-o_x(j)
      oh2vec2_y = h2_y(j)-o_y(j)
      oh2vec2_z = h2_z(j)-o_z(j)

      oh1vec2_mag = sqrt(oh1vec2_x**2+oh1vec2_y**2+oh1vec2_z**2)
      oh2vec2_mag = sqrt(oh2vec2_x**2+oh2vec2_y**2+oh2vec2_z**2)

      oh1vec2_x = oh1vec2_x/oh1vec2_mag
      oh1vec2_y = oh1vec2_y/oh1vec2_mag
      oh1vec2_z = oh1vec2_z/oh1vec2_mag
      oh2vec2_x = oh2vec2_x/oh2vec2_mag
      oh2vec2_y = oh2vec2_y/oh2vec2_mag
      oh2vec2_z = oh2vec2_z/oh2vec2_mag

!     Angle reference OH vectors make with observed OO vector
      a1 = acos(oh1vec1_x*(osep_x/oo_dist)+oh1vec1_y*&
       &(osep_y/oo_dist) + oh1vec1_z*(osep_z/oo_dist))
      a2 = acos(oh2vec1_x*(osep_x/oo_dist)+oh2vec1_y*&
       &(osep_y/oo_dist) + oh2vec1_z*(osep_z/oo_dist))

!     Angle observed OH vectors make with reference OO vector
!     Negative because we are really lumping a full iteration
!     over (i,j!=i) into (i,j>i)
      a3 = acos(oh1vec2_x*(-1*osep_x/oo_dist)+oh1vec2_y*&
       &(-1*osep_y/oo_dist) + oh1vec2_z*(-1*osep_z/oo_dist))
      a4 = acos(oh2vec2_x*(-1*osep_x/oo_dist)+oh2vec2_y*&
       &(-1*osep_y/oo_dist) + oh2vec2_z*(-1*osep_z/oo_dist))

      a1 = a1*RTD
      a2 = a2*RTD
      a3 = a3*RTD
      a4 = a4*RTD
              
!   Add observations to histogram if they fit within specified range
      IF (oo_dist .lt. rmax .and. oo_dist .gt. rmin) THEN
        rbin = int((oo_dist-rmin)/rbin_width) + 1
        IF (a1 .lt. amax .and. a1 .gt. amin) THEN
          abin = int((a1-amin)/abin_width) + 1
          hist_ga_sf(rbin,abin) = hist_ga_sf(rbin,abin) + 1
        END IF 
        IF (a2 .lt. amax .and. a2 .gt. amin) THEN
          abin = int((a2-amin)/abin_width) + 1
          hist_ga_sf(rbin,abin) = hist_ga_sf(rbin,abin) + 1
        END IF 
        IF (a3 .lt. amax .and. a3 .gt. amin) THEN
          abin = int((a3-amin)/abin_width) + 1
          hist_ga_sf(rbin,abin) = hist_ga_sf(rbin,abin) + 1
        END IF 
        IF (a4 .lt. amax .and. a4 .gt. amin) THEN
          abin = int((a4-amin)/abin_width) + 1
          hist_ga_sf(rbin,abin) = hist_ga_sf(rbin,abin) + 1
        END IF 
      END IF
    END DO
  END DO

! Now that we are done adding observations, we must normalize w.r.t
! the number density (done here because it varies for every frame)
!
! The normalization by the volume of the bin was done in an external
! python routine that called this subroutine to help compute the 
! ensemble average. This normalization can be done once at the end
! because it doesn't change like the number density.
!
! See the paper for the complete definition.
  DO i = 1, nrbin
    DO j = 1, nabin
      Dens_Norm_Fac = (BOXLEN*BOXLEN*BOXLEN)/(2*NWAT*NWAT)
      hist_ga(i,j) = hist_ga(i,j) + (hist_ga_sf(i,j)*Dens_Norm_Fac)
    END DO
  END DO

END SUBROUTINE
