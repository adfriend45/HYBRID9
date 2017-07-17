!======================================================================!
PROGRAM H9
!----------------------------------------------------------------------!
test
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
! Uses PGF forcings, 1901-2012, or local climate if LCLIM option.
!----------------------------------------------------------------------!
! Note yet implemented: energy balance, snow, soil C, N, P,carbon,
! optimisation of chunks, competition.
!----------------------------------------------------------------------!
! Andrew D. Friend
! 17th July, 2017
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

!----------------------------------------------------------------------!
! For spinning up if use local climate.
!----------------------------------------------------------------------!
INTEGER :: iLOOP
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Row input of local daily climate.
!----------------------------------------------------------------------!
REAL, DIMENSION  (6) :: LCLIM_array
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Row input to local sub-daily climate.
!----------------------------------------------------------------------!
REAL, DIMENSION (37) :: LCLIM_array2
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Annual sums for diagnostics.
!----------------------------------------------------------------------!
REAL :: npp_sum         ! Sum of NPP over 1-yr            (g[DM]/m^2/yr)
REAL :: plant_mass_sum  ! Sum of plant mass over 1-yr            (g[DM])
REAL :: tas_sum         ! Sum of daily tas over 1-yr                 (K)
REAL :: rlds_sum        ! Sum of daily rlds over 1-yr            (W/m^2)
REAL :: rsds_sum        ! Sum of daily rsds over 1-yr            (W/m^2)
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
          rlds_sum = zero
          rsds_sum = zero
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
            ! Day of year (d).
            !----------------------------------------------------------!
            DOY = iTIME - time_BOY (jyear-1859) + 1
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Set weather variables.
            !----------------------------------------------------------!
            tak = tas (x,y,iT) ! Surface air temperature (K)
            rh = rhs (x,y,iT)  ! Relative humidity (%age)
            Rnet = 0.92 * rsds (x,y,iT) + rlds (x,y,iT) - &
                   stbo * tas (x,y,iT) ** 4 !***Rneti really!!!! 
                   !assumes 8% albedo across short wavelengths.
            PAR = 0.92 * rsds (x,y,iT) * 2.3
            ppt = pr (x,y,iT) ! Precipitation flux (kg/m^2/s)
            !----------------------------------------------------------!
            ! Precipitation flux (mm s-1).
            !----------------------------------------------------------!
            forc_rain = 1.0E3 * pr (x,y,iT) / rhow
            !----------------------------------------------------------!
            ! Latent heat of vapourisation of water using Eqn. (4) of
            ! Pereira da Silva (2012)                             (J/kg)
            ! Check what CESM does.
            !----------------------------------------------------------!
            lamb = ((2503.0 - 2.386 * (tak - tf))) * 1.0E3
            !----------------------------------------------------------!
            ! Daily evaporation                                 (mm/day)
            !----------------------------------------------------------!
            evap_day = zero
            evap_grnd_day = zero
            !----------------------------------------------------------!
            ! Loop over timesteps within day (n).
            !----------------------------------------------------------!
            DO NS = 1, NISURF

!*****sub-daily climate needed!!!!? try just using an integration of
! light energy absorbed...! timestepping here really for hydrology in
! soil to be stable I think...but transpiration will need it, so do it
! ?

              !--------------------------------------------------------!
              ! Integrate hydrology over timestep.
              !--------------------------------------------------------!
              CALL HYDROLOGY
              !--------------------------------------------------------!        !--------------------------------------------------------------!
              ! Accumulate daily evaporation                    (mm/day)
              !--------------------------------------------------------!
              evap_day = evap_day + &
                         (qflx_evap_grnd + qflx_tran_veg_col) * dt
              evap_grnd_day = evap_grnd_day + qflx_evap_grnd * dt
            !----------------------------------------------------------!
            END DO ! DO NS = 1, NISURF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Integrate biomass over day.
            !----------------------------------------------------------!
            CALL GROW
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            IF (INTERACTIVE) THEN
              !WRITE (99,*) iT,theta(1),theta(2),theta(8),prec*dt*48.0
              WRITE (20,'(2(I5,A1),10(F10.4,A1),F10.4)') &
               jyear,',',DOY,',',&
               evap_day,',',&
               evap_grnd_day,',',theta(1),',',theta(2),',',theta(3),&
               ',',theta(4),',',theta_ma(1),',',LAI,',',LAI_litter,&
               ',',w_i,',',fT
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated over days in year.
            !----------------------------------------------------------!
            tas_sum  = tas_sum  + tas  (x,y,iT) ! K
            rlds_sum = rlds_sum + rlds (x,y,iT) ! W/m^2
            rsds_sum = rsds_sum + rsds (x,y,iT) ! W/m^2
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
          axy_rlds (x,y,iY) = rlds_sum / FLOAT (nt)              ! W/m^2
          axy_rsds (x,y,iY) = rsds_sum / FLOAT (nt)              ! W/m^2
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
    DEALLOCATE (rlds)
    DEALLOCATE (rsds)
    DEALLOCATE (huss)
    DEALLOCATE (ps)
    DEALLOCATE (pr)
    DEALLOCATE (rhs)
    !------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  END DO ! iDEC = iDEC_start, iDEC_end
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  IF (INTERACTIVE) CLOSE (20) ! Daily CSV output.
  !--------------------------------------------------------------------!

ELSE ! PGF == .F. so assume use LCLIM.

  DO iLOOP = 1, NYR_SPIN_UP - (eyr - syr + 1) ! Spin-up

  !--------------------------------------------------------------------!
  ! Open local climate file.
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',TRIM(LCLIM_filename)
  !--------------------------------------------------------------------!
  OPEN (10,FILE=LCLIM_filename,STATUS='OLD')
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  READ (10,*)
  DO jyear = syr, eyr
    IF (jyear == 2002) &
    OPEN (11,FILE='LCLIM/Vaira-grassland-2002-gapfilled-v2015.csv',&
                   STATUS='OLD')
    IF (jyear == 2003) &
    OPEN (11,FILE='LCLIM/Vaira-grassland-2003-gapfilled-v2015.csv',&
                   STATUS='OLD')
    READ (11,*)
    DO iTIME = time_BOY (jyear-1859), time_BOY (jyear+1-1859) - 1

      !----------------------------------------------------------------!
      iT = iTIME-time_BOY(syr-1859) + 1
      !----------------------------------------------------------------!

      !----------------------------------------------------------------!
      READ (10,*) iDOY,LCLIM_array
      !----------------------------------------------------------------!
      evap_obs_day = LCLIM_array (1) / sday ! Read as mm/day.
      !----------------------------------------------------------------!
      ! Precipitation flux (mm s-1).
      !----------------------------------------------------------------!
      pr   (x,y,iT) = LCLIM_array (2) / sday
      forc_rain = 1.0E3 * pr (x,y,iT) / rhow
      tas  (x,y,iT) = LCLIM_array (3) + tf ! K.
      huss (x,y,iT) = LCLIM_array (5)      ! kg/kg.
      ps   (x,y,iT) = LCLIM_array (6)      ! Pa.
      rhs  (x,y,iT) = LCLIM_array (4)      ! %age.
      IF (jyear == 2002) THEN
      IF (iDOY == 1) LAI = 0.88
      IF (iDOY == 59) LAI = 1.17
      IF (iDOY == 79) LAI = 1.87
      IF (iDOY == 94) LAI = 2.23
      IF (iDOY == 108) LAI = 2.55
      IF (iDOY == 122) THEN
        LAI = 1.43
        LAI_litter = LAI_litter + 2.55 - 1.43
      END IF
      IF (iDOY == 136) THEN
        LAI = 0.001
        LAI_litter = LAI_litter + 1.43 - 0.001
      END IF
      IF (iDOY == 357) LAI = 0.61
      END IF
      IF (jyear == 2003) THEN
        IF (iDOY == 29) LAI = 0.96
        IF (iDOY == 52) LAI = 1.58
        IF (iDOY == 76) LAI = 1.82
        IF (iDOY == 95) LAI = 2.63
        IF (iDOY == 106) THEN
          LAI = 2.52
          LAI_litter = LAI_litter + 2.63 - 2.52
        END IF
        IF (iDOY == 120) THEN
          LAI = 1.86
          LAI_litter = LAI_litter + 2.52 - 1.86
        END IF
        IF (iDOY == 141) THEN
          LAI = 0.76
          LAI_litter = LAI_litter + 1.86 - 0.76
        END IF
        IF (iDOY == 158) THEN
          LAI = 0.001
          LAI_litter = LAI_litter + 0.76 - 0.001
        END IF
      END IF
      !----------------------------------------------------------------!
      ! Daily evaporation (mm/day).
      !----------------------------------------------------------------!
      evap_day = zero
      evap_grnd_day = zero
      !----------------------------------------------------------------!
      ! Loop over timesteps within day (n).
      !----------------------------------------------------------------!
      DO NS = 1, NISURF
        !--------------------------------------------------------------!
        READ (11,*) LCLIM_array2
        !--------------------------------------------------------------!
        tak  = LCLIM_array2 (22) + tf
        rh   = LCLIM_array2 (25)
        Rnet = LCLIM_array2 (14)
        PAR  = LCLIM_array2 (16)
        ppt  = LCLIM_array2 (35) / dt ! mm/s
        Gs   = LCLIM_array2 (17)
        !--------------------------------------------------------------!
        ! Use sub-daily forcings.
        !--------------------------------------------------------------!
        forc_rain = ppt
        !--------------------------------------------------------------!
        ! Latent heat of vapourisation of water using Eqn. (4) of
        ! Pereira da Silva (2012)                                 (J/kg)
        ! Check with what CESM does.
        !--------------------------------------------------------------!
        lamb = ((2503.0 - 2.386 * (tak - tf))) * 1.0E3
        !--------------------------------------------------------------!
        ! Obs. ET read as W/m^2, convert to mm/s.
        !--------------------------------------------------------------!
        evap_obs = (LCLIM_array2 (12) / lamb) * 1.0E3 / rhow
        !--------------------------------------------------------------!
        ! Integrate hydrology over timestep.
        !--------------------------------------------------------------!
        CALL HYDROLOGY
        !--------------------------------------------------------------!
        ! Accumulate daily evaporation                          (mm/day)
        !--------------------------------------------------------------!
        evap_day = evap_day + (qflx_evap_grnd + qflx_tran_veg_col) * dt
        evap_grnd_day = evap_grnd_day + qflx_evap_grnd * dt
        !--------------------------------------------------------------!
      END DO ! DO NS = 1, NISURF
      !----------------------------------------------------------------!
      ! Write daily diagnostics.
      !----------------------------------------------------------------!
              WRITE (20,'(2(I5,A1),10(F10.4,A1),F10.4)') &
               jyear,',',DOY,',',&
               evap_day,',',&
               evap_grnd_day,',',theta(1),',',theta(2),',',theta(3),&
               ',',theta(4),',',theta_ma(1),',',LAI,',',LAI_litter,&
               ',',w_i,',',fT
      !----------------------------------------------------------------!
      ! Integrate biomass over day.
      !----------------------------------------------------------------!
!***** also in HYDROLOGY now
!      LAI_litter = 0.99 * LAI_litter ! mostly gone in 1-yr.
      !CALL GROW
      !----------------------------------------------------------------!

    END DO
    CLOSE (11)
  END DO
  CLOSE (10)
  !--------------------------------------------------------------------!
end do
CLOSE (20)
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
DEALLOCATE (axy_rlds)
DEALLOCATE (axy_rsds)
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
IF (.NOT. (PGF)) THEN
  !--------------------------------------------------------------------!
  !CLOSE (10) ! Close file of local climate forcing.
  CLOSE (20) ! Close file for daily diagnostics.
  !--------------------------------------------------------------------!
END IF
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

