! Written by Colin Bunner (bunne043@umn.edu) some time in 2018
SUBROUTINE lmfe(NWAT,NH2,BOXLEN,rmin,rmax,nrbin,h1_x,h1_y,h1_z,h2_x,h2_y,h2_z,&
  &o_x,o_y,o_z,hh_x,hh_y,hh_z,lmfe_h2)

!!! Variables !!!

! PARAMETER
!--------------------------------------------------------------------------
! dp - define double precision

! IN
!--------------------------------------------------------------------------
! NWAT          - number of H2O molecules in frame
! NH2           - number of H2 molecules in frame
! BOXLEN        - length of simulation box (Angstrom). Cubic implied.
! nrbin         - number of radial bins for LMFE histogram
! rmin,rmax,    - min/max radial distances allowed for histogram.
! o_i,h1_i,h2_i - x,y,z coordinates of o/h1/h2 beads of water model
! hh_i          - x,y,z coordinates of 1SM H2 beads

! INOUT
!--------------------------------------------------------------------------
! lmfe_h2       - local mole fraction enhancement histogram for H2

! DUMMY
!--------------------------------------------------------------------------
! cm_dist - COM distance
! cmvec_i - x,y,z coordinates of COM separation vector
! cm1_i,  - x,y,z coordinates of COM of 1st/2nd molecule
! cm2_i
! rbin_width - radial bin width for LMFE histogram
! rbin       - holds computed radial bin number
! count_h2,  - holds number of observations of an h2/wat molecule at a
! count_wat    given COM separation distance from an h2 molecule
       


  IMPLICIT NONE
  INTEGER, PARAMETER :: dp=selected_real_kind(15,307)

  INTEGER, INTENT(IN) :: NWAT, NH2, nrbin
  REAL(dp), INTENT(IN) :: BOXLEN, rmin, rmax
  REAL(dp), DIMENSION(NWAT), INTENT(IN) :: o_x,o_y,o_z, &
   &h1_x,h1_y,h1_z,h2_x,h2_y,h2_z
  REAL(dp), DIMENSION(NH2), INTENT(IN) :: hh_x, hh_y, hh_z
  REAL(dp), DIMENSION(nrbin), INTENT(INOUT) :: lmfe_h2
!f2py intent(in,out) :: lmfe_h2


  REAL(dp) :: cm_dist, nuh2, nutot
  REAL(dp) :: cmvec_x, cmvec_y, cmvec_z, cm1_x, cm1_y, cm1_z, &
   &cm2_x, cm2_y, cm2_z
  REAL(dp) :: rbin_width
  INTEGER :: rbin
  REAL(dp), DIMENSION(nrbin) :: count_h2, count_wat
  INTEGER :: i, j

  count_h2 = 0.d0
  count_wat = 0.d0

! Compute bin widths 
  rbin_width = (rmax-rmin)/nrbin

! Loop over all H2 molecules
  DO i = 1, NH2

!   COM (since 1SM, this is just the bead coordinate)
    cm1_x = hh_x(i)
    cm1_y = hh_y(i)
    cm1_z = hh_z(i)

    ! First check against water molecules
    DO j = 1 , NWAT
      ! COM of observed water 
      cm2_x = (h1_x(j)*1.0079+h2_x(j)*1.0079+o_x(j)*15.999)/18.0148
      cm2_y = (h1_y(j)*1.0079+h2_y(j)*1.0079+o_y(j)*15.999)/18.0148
      cm2_z = (h1_z(j)*1.0079+h2_z(j)*1.0079+o_z(j)*15.999)/18.0148
             
      ! COM separation
      cmvec_x = cm2_x-cm1_x
      cmvec_y = cm2_y-cm1_y
      cmvec_z = cm2_z-cm1_z

      ! Periodic boundary conditions
      cmvec_x = cmvec_x - BOXLEN*nint(cmvec_x/BOXLEN)
      cmvec_y = cmvec_y - BOXLEN*nint(cmvec_y/BOXLEN)
      cmvec_z = cmvec_z - BOXLEN*nint(cmvec_z/BOXLEN)

      cm_dist = sqrt(cmvec_x**2+cmvec_y**2+cmvec_z**2)

      IF (cm_dist .lt. rmax .and. cm_dist .gt. rmin) THEN
        rbin = int((cm_dist-rmin)/rbin_width) + 1
        count_wat(rbin) = count_wat(rbin) + 1
      END IF
    END DO

    ! And then other H2 molecules
    DO j = 1 , NH2
      if(j.eq.i) cycle ! Skip self-observations

      cm2_x = hh_x(j)
      cm2_y = hh_y(j)
      cm2_z = hh_z(j)
             
      ! COM separation
      cmvec_x = cm2_x-cm1_x
      cmvec_y = cm2_y-cm1_y
      cmvec_z = cm2_z-cm1_z

      ! Periodic boundary conditions
      cmvec_x = cmvec_x - BOXLEN*nint(cmvec_x/BOXLEN)
      cmvec_y = cmvec_y - BOXLEN*nint(cmvec_y/BOXLEN)
      cmvec_z = cmvec_z - BOXLEN*nint(cmvec_z/BOXLEN)

      cm_dist = sqrt(cmvec_x**2+cmvec_y**2+cmvec_z**2)

      IF (cm_dist .lt. rmax .and. cm_dist .gt. rmin) THEN
        rbin = int((cm_dist-rmin)/rbin_width) + 1
        count_h2(rbin) = count_h2(rbin) + 1
      END IF
    END DO
  END DO

  nuh2 = 0.d0
  nutot = 0.d0

  DO i = 1, nrbin
    nuh2 = nuh2+count_h2(i)
    nutot = nutot+count_h2(i)+count_wat(i)
    ! Avoid NaNs
    ! Also, a previous bug arose from an integer/real confusion, so
    ! cast NH2/(NWAT+NH2) because this global mole fraction is a real
    ! quantity.
    if(nutot.gt.1) then
      lmfe_h2(i)=lmfe_h2(i)+((nuh2/nutot)/(real(NH2)/real(NWAT+NH2)))
    end if
  END DO
END SUBROUTINE
