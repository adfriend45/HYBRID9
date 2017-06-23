!======================================================================!
PROGRAM H9
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global simulation of land surface water and carbon fluxes.
!----------------------------------------------------------------------!
! As far as necessary, uses approach of CESM, based on that code
! directly. Follows as in 'Hydrology2Mod.F90'.
! Calling sequence of CESM:
!  Hydrology2:              surface hydrology driver
!    -> SnowWater:          change of snow mass and snow water onto soil
!    -> SurfaceRunoff:      surface runoff
!    -> Infiltration:       infiltration into surface soil layer
!    -> SoilWater:          soil water movement between layers
!          -> Tridiagonal   tridiagonal matrix solution
!    -> Drainage:           subsurface runoff
!    -> SnowCompaction:     compaction of snow layers
!    -> CombineSnowLayers:  combine snow layers thinner than minimum
!    -> DivideSnowLayers:   subdivide snow layers thicker than maximum
!----------------------------------------------------------------------!
! Uses PGF forcings, 1901-2012.
!----------------------------------------------------------------------!
! Note yet implemented: energy balance, snow, carbon (vegetation and
! and soil), optimisation of chunks.
!----------------------------------------------------------------------!
! Andrew D. Friend
! 19th June, 2017
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE MPI     ! Enable access to the Message Passing Interface library of
            ! parallel routines.
USE NETCDF  ! Enable access to the library of netCDF routines.
USE CONTROl ! Control variables.
USE SHARED  ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

INTEGER :: DOY
REAL, DIMENSION (6) :: LCLIM_array

!----------------------------------------------------------------------!
! Annual sums for diagnostics.
!----------------------------------------------------------------------!
REAL :: npp_sum         ! Sum of NPP over 1-yr            (g[DM]/m^2/yr)
REAL :: plant_mass_sum  ! Sum of plant mass over 1-yr            (g[DM])
REAL :: tas_sum         ! Sum of daily tas over 1-yr                 (K)
REAL :: huss_sum        ! Sum of daily huss over 1-yr            (kg/kg)
REAL :: ps_sum          ! Sum of daily ps over 1-yr                 (Pa)
REAL :: pr_sum          ! Sum of daily pr over 1-yr           (kg/m^2/s)
REAL :: rhs_sum         ! Sum of daily rhs over 1-yr              (%age)
REAL :: h2osoi_sum_total ! Sum of daily total soil water over 1-yr  (mm)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! End of declarations.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL INIT ! Initialise simulation.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Simulation using PGF climate forcing?
!----------------------------------------------------------------------!
IF (PGF) THEN
!----------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Outer main loop: over decades (iDEC 1 = 1901-1910).
  !--------------------------------------------------------------------!
  DO iDEC = iDEC_start, iDEC_end
  !--------------------------------------------------------------------!

    !------------------------------------------------------------------!
    CALL READ_PGF ! Read decade of PGF climate.
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    ! syr is jyear at beginning of decade iDEC.
    !------------------------------------------------------------------!
    syr = (iDEC - 1) * 10 + 1901
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    ! eyr is jyear at end of decade iDEC.
    !------------------------------------------------------------------!
    IF (iDEC < 12) THEN
      eyr = syr + 9
    ELSE
      eyr = syr + 1
    END IF
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    ! Inner main loops: over timepoints in decade at each grid point
    ! with soil in both HWSD texture and theta_s fields.
    !------------------------------------------------------------------!
    DO y = 1, lat_c
      DO x = 1, lon_c
        IF ((soil_tex (x,y) > 0) .AND. (soil_tex (x,y) /= 13) .AND. &
          SUM (theta_s (:,x,y)) > trunc) THEN
          !------------------------------------------------------------!
          nlayers = nsoil_layers_max
          !------------------------------------------------------------!
          ! Loop over years in decade.
          !------------------------------------------------------------!
          DO jyear = syr, eyr
            !----------------------------------------------------------!
            ! Intialise annual diagnostics.
            !----------------------------------------------------------!
          npp_sum        = zero
          plant_mass_sum = zero
          rnf_sum  = zero
          evap_sum = zero
          tas_sum  = zero
          huss_sum = zero
          ps_sum   = zero
          pr_sum   = zero
          rhs_sum  = zero
          h2osoi_sum_total = zero
          theta_sum (:) = zero
          !------------------------------------------------------------!
          ! Innermost main loop: over iTIME points in year.
          !------------------------------------------------------------!
          DO iTIME = time_BOY (jyear-1859), time_BOY (jyear+1-1859) - 1
          !------------------------------------------------------------!

            !----------------------------------------------------------!
            ! Local index within decade for time point.
            !----------------------------------------------------------!
            iT = iTIME-time_BOY(syr-1859) + 1
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Precipitation flux (mm s-1).
            !----------------------------------------------------------!
            prec = 1.0E3 * pr (x,y,iT) / rhow
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Loop over timesteps within day (n).
            !----------------------------------------------------------!
            DO NS = 1, NISURF

              !--------------------------------------------------------!
              ! Integrate hydrology over timestep.
              !--------------------------------------------------------!
              CALL HYDROLOGY
              !--------------------------------------------------------!

            !----------------------------------------------------------!
            END DO ! DO NS = 1, NISURF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Integrate biomass over day.
            !----------------------------------------------------------!
            CALL GROW
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            !IF (INTERACTIVE) THEN
              !WRITE (99,*) iT,theta(1),theta(2),theta(8),prec*dt*48.0
            !END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated over days in year.
            !----------------------------------------------------------!
            tas_sum  = tas_sum  + tas  (x,y,iT) ! K
            huss_sum = huss_sum + huss (x,y,iT) ! kg/kg
            ps_sum   = ps_sum   + ps   (x,y,iT)
            pr_sum   = pr_sum   + pr   (x,y,iT)
            rhs_sum  = rhs_sum  + rhs  (x,y,iT)
            DO K = 1, nplants (x,y)
              plant_mass_sum = plant_mass_sum + plant_mass (K,x,y)
            END DO
            npp_sum = npp_sum + npp
            !----------------------------------------------------------!
            ! Possibly move following to within NS loop.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              theta_sum (I) = theta_sum (I) + theta (I)
              h2osoi_sum_total = h2osoi_sum_total + h2osoi_liq (I,x,y)
              !rnff
            END DO
            !----------------------------------------------------------!

          !------------------------------------------------------------!
          END DO ! over iTIME points in year.
          !------------------------------------------------------------!

          !------------------------------------------------------------!
          ! No. timepoints in year (days).
          !------------------------------------------------------------!
          nt = (time_BOY (jyear + 1 - 1859) - 1) - &
               (time_BOY (jyear - 1859)) + 1
          !------------------------------------------------------------!
          ! Put means for year in global arrays at this year point
          ! relative to beginning of simulation.
          !------------------------------------------------------------!
          iY = jyear-((iDEC_start-1)*10+1901)+1
          !------------------------------------------------------------!
          axy_npp (x,y,iY) = npp_sum                      ! g[DM]/m^2/yr
          !------------------------------------------------------------!
          axy_plant_mass (x,y,iY) = plant_mass_sum / FLOAT (nt)  ! g[DM]
          !------------------------------------------------------------!
          axy_rnf  (x,y,iY) = rnf_sum  / FLOAT (nt * NISURF)      ! mm/s
          axy_evap (x,y,iY) = evap_sum / FLOAT (nt * NISURF)      ! mm/s
          !------------------------------------------------------------!
          axy_tas  (x,y,iY) = tas_sum  / FLOAT (nt)                  ! K
          axy_huss (x,y,iY) = huss_sum / FLOAT (nt)              ! kg/kg
          axy_ps   (x,y,iY) = ps_sum   / FLOAT (nt)                 ! Pa
          axy_pr   (x,y,iY) = pr_sum   / FLOAT (nt)           ! kg/m^2/s
          axy_rhs  (x,y,iY) = rhs_sum  / FLOAT (nt) !             (%age)
          !------------------------------------------------------------!
          DO I = 1, nlayers
            axy_theta (I,x,y,iY) = theta_sum (I) / FLOAT (nt) !(m^3/m^3)
          END DO
          ! 
          axy_theta_total (x,y,iY) = h2osoi_sum_total / FLOAT (nt)
          !------------------------------------------------------------!
        END DO ! jyear = syr, eyr (calendar years in decade)
      END IF ! Soiled grid-box?
    END DO ! x-loop (lons).
  END DO ! y-loop (lats).
  !--------------------------------------------------------------------!

    !------------------------------------------------------------------!
    IF (INTERACTIVE) THEN
      write (*,*) 'theta'
      write (*,*) theta
      write (*,*) 'smp'
      write (*,*) smp
      write (*,*) 'zq'
      write (*,*) zq
      write (*,*) 'nplants'
      write (*,*) nplants
      write (*,*) 'plant_mass'
      write (*,*) plant_mass
      write (*,*) 'NPP'
      write (*,*) axy_npp
      !write (*,*) 'lambda'
      !write (*,*) lambda
      !write (*,*) 'bsw'
      !write (*,*) 1.0/lambda
    END IF
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    ! Free up resources.
    !------------------------------------------------------------------!
    DEALLOCATE (tas)
    DEALLOCATE (huss)
    DEALLOCATE (ps)
    DEALLOCATE (pr)
    DEALLOCATE (rhs)
    !------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  END DO ! iDEC = iDEC_start, iDEC_end
  !--------------------------------------------------------------------!

ELSE ! PGF == .F. so assume use LCLIM.

  !--------------------------------------------------------------------!
  x = 1
  y = 1
  nlayers = nsoil_layers_max
  syr = 2002
  eyr = 2003
  NTIMES = time_BOY (eyr+1-1859) - time_BOY (syr-1859)
  ALLOCATE (pr  (lon_c,lat_c,NTIMES))
  ALLOCATE (tas (lon_c,lat_c,NTIMES))
  ALLOCATE (huss(lon_c,lat_c,NTIMES))
  ALLOCATE (ps  (lon_c,lat_c,NTIMES))
  ALLOCATE (rhs (lon_c,lat_c,NTIMES))
  ! Open file for daily diagnostics.
  OPEN (20,FILE='LCLIM/daily_diag.csv',STATUS='UNKNOWN')
  WRITE (20,*) 'DOY,et,theta_4'
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Open local climate file.
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',LCLIM_filename
  !--------------------------------------------------------------------!
  OPEN (10,FILE=LCLIM_filename,STATUS='OLD')
  !--------------------------------------------------------------------!
  READ (10,*)
  DO jyear = syr, eyr
    DO iTIME = time_BOY (jyear-1859), time_BOY (jyear+1-1859) - 1

      !----------------------------------------------------------------!
      iT = iTIME-time_BOY(syr-1859) + 1
      !----------------------------------------------------------------!

      !----------------------------------------------------------------!
      READ (10,*) DOY,LCLIM_array
      !----------------------------------------------------------------!
      ! Precipitation flux (mm s-1).
      !----------------------------------------------------------------!
      pr   (x,y,iT) = 10.0 * LCLIM_array (2) / sday
      prec = 1.0E3 * pr (x,y,iT) / rhow
      tas  (x,y,iT) = LCLIM_array (3) + tf ! K.
      huss (x,y,iT) = LCLIM_array (5)      ! kg/kg.
      ps   (x,y,iT) = LCLIM_array (6)      ! Pa.
      rhs  (x,y,iT) = LCLIM_array (4)      ! %age.
      !----------------------------------------------------------------!
      ! Daily evaporation (mm/day).
      !----------------------------------------------------------------!
      evap_day = 0.0
      !----------------------------------------------------------------!
      ! Loop over timesteps within day (n).
      !----------------------------------------------------------------!
      DO NS = 1, NISURF
        !--------------------------------------------------------------!
        ! Integrate hydrology over timestep.
        !--------------------------------------------------------------!
        CALL HYDROLOGY
        !--------------------------------------------------------------!
        ! Accumulate daily evaporation                          (mm/day)
        !--------------------------------------------------------------!
        evap_day = evap_day + qflx_evap * dt
        !--------------------------------------------------------------!
      END DO ! DO NS = 1, NISURF
      !----------------------------------------------------------------!
      ! Write daily diagnostics.
      !----------------------------------------------------------------!
      WRITE (20,*) iT,',',evap_day,',',theta(4)
      !----------------------------------------------------------------!
      ! Integrate biomass over day.
      !----------------------------------------------------------------!
      !CALL GROW
      !----------------------------------------------------------------!

    END DO
  END DO
  !--------------------------------------------------------------------!
  CLOSE (10) ! Close file of local climate forcing.
  CLOSE (20) ! Close file for daily diagnostics.
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
END IF ! PGF?
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Loop over decades to write annual diagnostics.
!----------------------------------------------------------------------!
DO iDEC = iDEC_start, iDEC_end
  syr = (iDEC - 1) * 10 + 1901
  IF (iDEC < 12) THEN
    eyr = syr + 9
  ELSE
    eyr = syr + 1
  END IF
  !--------------------------------------------------------------------!
  ! Loop over calendar years in the decade.
  !--------------------------------------------------------------------!
  DO jyear = syr, eyr
  !--------------------------------------------------------------------!

    !------------------------------------------------------------------!
    WRITE (ydate,'(I4)') jyear
    !file_name = '/scratch/adf10/HYBRID9/OUTPUT/axy'//ydate//'.nc'
    file_name = TRIM(PATH_output)//'/axy'//ydate//'.nc'
    WRITE (*,*) 'Writing to ',TRIM(file_name)
    !------------------------------------------------------------------!

    !------------------------------------------------------------------!
    CALL WRITE_NET_CDF_3DR (jyear-((iDEC_start-1)*10+1901)+1)
    !------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  END DO ! Next calendar year.
  !--------------------------------------------------------------------!
END DO ! Next decade.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! DEALLOCATE chunk arrays to free up resources.
!----------------------------------------------------------------------!
DEALLOCATE (plant_mass)
DEALLOCATE (plant_length)
DEALLOCATE (nplants)
DEALLOCATE (data_in_2DI)
DEALLOCATE (theta_s)
DEALLOCATE (hksat)
DEALLOCATE (lambda)
DEALLOCATE (bsw)
DEALLOCATE (psi_s)
DEALLOCATE (theta_m)
DEALLOCATE (lon)
DEALLOCATE (lat)
DEALLOCATE (block_sub)
DEALLOCATE (soil_tex)
DEALLOCATE (topo_slope)
DEALLOCATE (Fmax)
DEALLOCATE (h2osoi_liq)
DEALLOCATE (zwt)
DEALLOCATE (wa)
DEALLOCATE (axy_plant_mass)
DEALLOCATE (axy_rnf)
DEALLOCATE (axy_evap)
DEALLOCATE (axy_tas)
DEALLOCATE (axy_huss)
DEALLOCATE (axy_ps)
DEALLOCATE (axy_pr)
DEALLOCATE (axy_rhs)
DEALLOCATE (axy_theta)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CALL MPI_FINALIZE (ierr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
WRITE (*,*) 'H9 completed successfully on my_id = ',my_id
CALL CPU_TIME (finish(my_id+1))
WRITE (*,*) my_id,'finish-start = ',finish(my_id+1)-start(my_id+1),'s'
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CONTAINS
  SUBROUTINE CHECK (status)
    INTEGER, INTENT (IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM (nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE CHECK
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END PROGRAM H9
!======================================================================!

! Need to check if following really is OK.
! Does making LH a function of TM make a significant difference?
! What happens if make different for TM < 273.15 K (i.e. over ice)?
      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
      USE SHARED, ONLY : MRAT,RVAP,TF
      IMPLICIT NONE
!@var A,B,C   expansion coefficients for QSAT
      REAL, PARAMETER :: A=6.1080*MRAT    !3.797915d0
      REAL, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL, PARAMETER :: C= 1./RVAP        !2.166847d-3
!**** Note that if LH is considered to be a function of temperature, the
!**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
!**** LH = 0.5*(LH(0)+LH(t))
      REAL, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.0,TM)))/PR
      RETURN
      END

