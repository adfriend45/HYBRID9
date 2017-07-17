!======================================================================!
SUBROUTINE HYDROLOGY
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Integrate hydrology over timestep.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE CONTROl ! Control variables.
USE SHARED  ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
REAL :: w0,w1 ! For soil water diagnostics.
REAL :: qflx_prec_grnd_rain ! Rain precip. incident on ground     (mm/s)
REAL :: qflx_top_soil  ! Net water input into soil from top       (mm/s)
REAL :: hkdepth     ! Decay factor                                   (m)
REAL :: fff         ! Decay factor                                  (/m)
REAL :: wtfact      ! Max. saturation fraction for a grid cell       (-)
REAL :: fsat        ! Fractional area with water table at surface    (-)
REAL :: fcov        ! Fractional impermable area                     (-)
REAL :: qflx_surf   ! Surface runoff                              (mm/s)
REAL :: frac_h2osfc ! Fraction of ground covered by surface water    (-)
REAL :: Rcan
REAL :: beta_save
REAL :: PMc         ! Closed canopy evaporation                  (W/m^2)
REAL :: PMs         ! Bare substrate evaporation                 (W/m^2)
REAL :: desatdT     ! Rate of change of esat with T           (kg/m^3/K)
REAL :: G           ! Soil heat flux                             (W/m^2)
REAL, PARAMETER :: cp = 1010.0       ! Specific heat of dry air (J/kg/K)
REAL :: VDD         ! Water vapour density deficit              (kg/m^3)
REAL :: VDDs
REAL :: rh_soil
REAL :: gamma       ! Psychrometric parameter                 (kg/m^3/K)
REAL :: rac         ! Bulk boundary layer resistance in canopy     (s/m)
REAL :: Rnets       ! Net radiative flux into substrante         (W/m^2)
REAL :: raa         ! Canopy to air aerodynamic resistance         (s/m)
REAL :: rsc         ! Bulk stomatal resistance of the canopy       (s/m)
REAL :: ras         ! Substrate to canopy aerodynamic resistance   (s/m)
REAL :: rss         ! Surface resistance of the substrate          (s/m)
REAL :: Cc          ! Fraction of PMc to total evaporation           (-)
REAL :: Cs          ! Fraction of PMs to total evaporation           (-)
REAL :: Ra,Rs,Rc    ! Coefficients for Cc and Cs                     (?)
REAL :: LE
REAL :: VDD0
REAL :: LEs
REAL :: LEc
REAL :: LEw
REAL :: decay
REAL :: dchi
REAL :: dsmpdchi
!----------------------------------------------------------------------!
! Soil moisture constraint on evaporative flux                (fraction)
!----------------------------------------------------------------------!
REAL :: beta
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Virtual potential surface temperature                              (K)
!----------------------------------------------------------------------!
REAL :: tsv
!----------------------------------------------------------------------!
! Air density                                                (kg[a]/m^3)
!----------------------------------------------------------------------!
REAL :: rho
!----------------------------------------------------------------------!
! Saturation water vapour mixing ratio                     (kg[w]/kg[a])
!----------------------------------------------------------------------!
REAL :: QSAT
!----------------------------------------------------------------------!
! Water available for evaporation over timestep                   (mm/s)
!----------------------------------------------------------------------!
REAL :: evap_max (nsoil_layers_max)
REAL :: qflx_ev_h2osfc ! Evaporative flux from h2osfc             (mm/s)
REAL :: qflx_in_soil   ! Surface input to soil                    (mm/s)
REAL :: qflx_in_h2osfc ! Surface input to h2osfc                  (mm/s)
REAL :: qinmax    ! Maximum infiltration capacity                 (mm/s)
REAL :: qflx_infl_excess ! Infiltration excess runoff -> h2osfc   (mm/s)
REAL :: qflx_infl ! Infiltration                                  (mm/s)
REAL :: zwtmm    ! Water table depth                                (mm)
REAL :: tempi   ! Temporary variable for calculating vol_eq          (-)
REAL :: temp0   ! Temporary variable for calculating vol_eq          (-)
REAL :: voleq1  ! Temporary variable for calculating vol_eq          (-)
!----------------------------------------------------------------------!
! "s" at interface of layer.
!----------------------------------------------------------------------!
REAL :: s1
!----------------------------------------------------------------------!
! k*s**(2b+2)
!----------------------------------------------------------------------!
REAL :: s2
!----------------------------------------------------------------------!
! Soil wetness                                                       (-)
!----------------------------------------------------------------------!
REAL :: s_node
REAL :: den     ! For trigiagonal solution                          (mm)
REAL :: dzq     ! For trigiagonal solution                          (mm)
REAL :: num     ! For trigiagonal solution                          (mm)
!----------------------------------------------------------------------!
! Soil matric potential in aquifer layer                            (mm)
!----------------------------------------------------------------------!
REAL :: smp1
!----------------------------------------------------------------------!
! d(smp)/d(vol_liq) in aquifer layer                    (mm/(mm^3/mm^3))
!----------------------------------------------------------------------!
REAL :: dsmpdw1
REAL :: BET     ! For trigiagonal solution                           (-)
REAL :: wh_zwt  ! Water head at the water table depth               (mm)
REAL :: ka      ! Hydraulic conductivity of the aquifer           (mm/s)
REAL :: wh      ! smpfz (jwt) - z (jwt)                             (mm)
REAL :: qcharge ! Aquifer recharge rate                           (mm/s)
REAL :: rous    ! Aquifer yield                                      (-)
REAL :: qcharge_tot ! To compute water table depth change in soil (mm/s)
REAL :: s_y      ! Specific yield                                    (-)
REAL :: qcharge_layer ! qcharge in saturated layer                (mm/s)
REAL :: rsub_top_max ! Max. rsub_top                              (mm/s)
REAL :: rsub_top ! Subsurface runoff - topographic control        (mm/s)
REAL :: rsub_top_tot ! To remove water from water table in soil   (mm/s)
REAL :: rsub_top_layer ! To remove water from water table in soil (mm/s)
REAL :: xsi ! Excess soil water above saturation at layer I         (mm)
REAL :: xs1 ! Excess soil water above saturation at layer 1         (mm)
REAL :: qflx_rsub_sat ! Soil saturation excess                    (mm/s)
REAL :: xs  ! Water needed to bring soil moisture to theta_m        (mm)
REAL :: available_h2osoi_liq ! Available soil liq. water in a layer (mm)
!----------------------------------------------------------------------!
! Miscellaneous parameters.
!----------------------------------------------------------------------!
! Minimum soil moisture                                             (mm)
!----------------------------------------------------------------------!
REAL, PARAMETER :: watmin = 0.01
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initial total soil water for diagnostics                          (mm)
!----------------------------------------------------------------------!
w0 = forc_rain * dt + wa (x,y)
DO I = 1, nlayers
  !--------------------------------------------------------------------!
  w0 = w0 + h2osoi_liq (I,x,y)
  !--------------------------------------------------------------------!
  ! Also compute volumetric water                            (mm^3/mm^3)
  !--------------------------------------------------------------------!
  theta    (I) = h2osoi_liq    (I,x,y) / (dz (I) * rhow / 1.0E3)
  theta_ma (I) = h2osoi_liq_ma (I,x,y) / (dz (I) * rhow / 1.0E3)
  !--------------------------------------------------------------------!
END DO
!----------------------------------------------------------------------!

!**********************************************************************!
! Start of CESM 'SLakeHydrology'.
!**********************************************************************!

!----------------------------------------------------------------------!
! Rain precipitation incident on ground                           (mm/s)
!----------------------------------------------------------------------!
qflx_prec_grnd_rain = forc_rain
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Net water input into soil from top                              (mm/s)
!----------------------------------------------------------------------!
qflx_top_soil = qflx_prec_grnd_rain
!----------------------------------------------------------------------!

!**********************************************************************!
! End of CESM 'SLakeHydrology'.
!**********************************************************************!

!**********************************************************************!
! Start of CESM 'SurfaceRunoff'.
!**********************************************************************!

!----------------------------------------------------------------------!
! Decay factor for water table                                       (m)
! In CESM set using line in 'iniTimeConst.F90'.
!----------------------------------------------------------------------!
hkdepth = one / 2.5
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Decay factor for water table                                      (/m)
!----------------------------------------------------------------------!
fff = 1.0 / hkdepth
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Max. saturation fraction for grid cell                             (-)
!----------------------------------------------------------------------!
wtfact = Fmax (x,y)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Fractional area with water table at surface                        (-)
!----------------------------------------------------------------------!
fsat = wtfact * EXP (-0.5 * fff * zwt (x,y))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Fractional impermable area                                         (-)
!----------------------------------------------------------------------!
fcov = fsat
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Surface runoff                                                  (mm/s)
!----------------------------------------------------------------------!
qflx_surf = fcov * qflx_top_soil
!----------------------------------------------------------------------!

!**********************************************************************!
! End of CESM 'SurfaceRunoff'.
!**********************************************************************!

!----------------------------------------------------------------------!
! Fraction of ground covered by surface water
! (full treatment is in H2OSfcMod.F90: do later)                     (-)
! Set to zero for now.
!----------------------------------------------------------------------!
frac_h2osfc = zero
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Virtual potential surface temperature                              (K)
! Correct temperature to use for rho?
! deltx is coeff. of humidity in virtual temperature.
!----------------------------------------------------------------------!
tsv = tak * (one + huss (x,y,iT) * deltx)
!----------------------------------------------------------------------!
! Air density                                                (kg[a] m-3)
!----------------------------------------------------------------------!
rho = ps (x,y,iT) / (rgas * tsv)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Mean rate of change of saturated vapour pressure with temperature
! ('s' of 'Delata' in classic Penman-Monteith equations)      (kg/m^3/K)
! Equation based on  http://www.fao.org/docrep/X0490E/x0490e07.htm
!----------------------------------------------------------------------!
desatdT = (4098.0 * (0.6108 * EXP ((17.27 * (tak - tf)) / &
          (tak - tf + 237.3)))) / &
          ((tak - tf + 237.3) ** 2)
desatDt = desatDt * 18.0 / (gasc * tak)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Saturation vapour density                                     (kg/m^3)
! Eqn from www.fao.org/docrep/X0490E/x0490e07.htm
!----------------------------------------------------------------------!
esat = 0.6108 * EXP (17.27 * (tak - tf) / (tak - tf + 237.3))
esat = esat * 18.0 / (gasc * tak)
!----------------------------------------------------------------------!
! Water vapour density deficit                                  (kg/m^3)
!----------------------------------------------------------------------!
VDD = esat * (one - rh / 100.0)
!----------------------------------------------------------------------!
! Psychrometric parameter                                     (kg/m^3/K)
!----------------------------------------------------------------------!
gamma =  (cp * ps (x,y,iT)/ (lamb * 0.622)) * (18.0E-3 / (gasc * tak))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Moisture control on stomatal conductance of canopy                 (-)
!----------------------------------------------------------------------!
beta_save = zero
DO I = 1, nlayers
  beta = one - (smp (I) - zc (I)) / (-150000.0)
  beta = MIN (one, beta)
  beta = MAX (zero, beta)
  beta_save = beta_save + rootr_col (I,x,y) * beta
END DO
beta = beta_save
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Bulk stomatal resistance of the canopy.
! LAI dependence from Shuttleworth Eqn. (20)                       (s/m)
!----------------------------------------------------------------------!
IF ((LAI (x,y) > zero) .AND. (beta > zero) .AND. (PAR > zero)) &
  THEN
  rsc = (1.0 / (PAR / (PAR + 300.0))) * 400.0 / &
        (2.0 * LAI (x,y) * beta * &
        (2.8 ** (-80.0 * MAX (zero, VDD) / rho)))
ELSE
  rsc = 1.0E6
END IF
!----------------------------------------------------------------------!
! Minimum resistance based on Baldocchi et al. (2004). Is about
! 113 s/m for full leaf area at 0.6 mol/m^s/s.
!----------------------------------------------------------------------!
rsc = MAX (rsc,1.0/((LAI (x,y)/2.7)*0.9/(rho*1.0E3/18.0)))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Bulk boundary layer resistance of the vegetative elements
! in the canopy. Eqn. (20) of SW85                                 (s/m)
!----------------------------------------------------------------------!
IF (LAI (x,y) > zero) THEN
  rac = 25.0 / (2.0 * LAI (x,y))
ELSE
  rac = 1.0E6
END IF
!----------------------------------------------------------------------!
! Aerodynamic resistances between canopy source height and
! reference level, and between the substrate and canopy source height;
! Eqns. (30;31) SW85                                               (s/m)
!----------------------------------------------------------------------!
IF (LAI (x,y) <= 4.0) THEN
  raa = 0.25 * LAI (x,y) *  42.0 + 0.25 * (4.0 - LAI (x,y)) * 34.0
  ras = 0.25 * LAI (x,y) * 128.0 + 0.25 * (4.0 - LAI (x,y)) * 49.0
ELSE
  raa =  42.0
  ras = 128.0
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Surface resistance of the substrate                              (s/m)
! Uses Eqn. (20) of van de Griend & Owe (1994).
!----------------------------------------------------------------------!
IF (theta (1) <= 0.15) THEN
  rss = (10.0 + 1000.0 * LAI_litter (x,y)) * &
        EXP (0.3563 * 100.0 * (0.15 - theta (1)))
ELSE
  rss = (10.0 + 1000.0 * LAI_litter (x,y) * &
        (1.0 - theta(1)/theta_s(1,x,y)))
END IF
!----------------------------------------------------------------------!
! Net radiation into the substrate (Eqn. 21 SW85)                (W/m^2)
!----------------------------------------------------------------------!
Rnets = Rnet * EXP (-0.7 * LAI (x,y))
!----------------------------------------------------------------------!
! Soil heat flux                                                 (W/m^2)
!----------------------------------------------------------------------!
G = 0.2 * Rnets
!----------------------------------------------------------------------!
! Separate Penman-Monteith evaporation between canopy and
! substrate (Eqn. 12 SW85)                                       (W/m^2)
!----------------------------------------------------------------------!
PMc = (desatdT * (Rnet - G) + (rho * cp * VDD - &
      desatdT * rac * (Rnets - G)) / (raa + rac)) / &
      (desatdT + gamma * (one + rsc / (raa + rac)))
!----------------------------------------------------------------------!
! Penman-Monteith bare substrate evaporation (Eqn. 13 SW85)      (W/m^2)
!----------------------------------------------------------------------!
PMs = (desatdT * (Rnet - G) + (rho * cp * VDD - &
      desatdT * ras * (Rnet - Rnets)) / (raa + ras)) / &
      (desatdT + gamma * (one + rss / (raa + ras)))
!----------------------------------------------------------------------!
! Coefficients for Cc and Cr                                       (s/m)
!----------------------------------------------------------------------!
Ra = (desatdT + gamma) * raa
Rs = (desatdT + gamma) * ras + gamma * rss
Rc = (desatdT + gamma) * rac + gamma * rsc
!----------------------------------------------------------------------!
! Fraction of PMc to total evaporation                               (-)
!----------------------------------------------------------------------!
Cc = one / (one + Rc * Ra / (Rs * (Rc + Ra)))
!----------------------------------------------------------------------!
! Fraction of PMs to total evaporation                               (-)
!----------------------------------------------------------------------!
Cs = one / (one + Rs * Ra / (Rc * (Rs + Ra)))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Total latent heat flux                                         (W/m^2)
!----------------------------------------------------------------------!
LE = Cc * PMc + Cs * PMs
!----------------------------------------------------------------------!
! Implied vapour density deficit at the canopy source height    (kg/m^3)
!----------------------------------------------------------------------!
VDD0 = VDD + (desatdT * (Rnet - G) - (desatdT + gamma) * LE) * &
       raa / (rho * cp)
!----------------------------------------------------------------------!
! Implied latent heat fluxes from canopy and substrate           (W/m^2)
!----------------------------------------------------------------------!
LEc = (desatdT * (Rnet - Rnets) + rho * cp * VDD0 / rac) / &
      (desatdT + gamma * (1.0 + rsc / rac))
LEs = (desatdT * (Rnets - G) + rho * cp * VDD0 / ras) / &
      (desatdT + gamma * (1.0 + rss / ras))
!----------------------------------------------------------------------!
! Convert from W/m^2 to mm/s.
!----------------------------------------------------------------------!
qflx_tran_veg_col = LEc * 1.0E3 / (rhow * lamb)
qflx_evap_grnd    = LEs * 1.0E3 / (rhow * lamb)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Limit substrate evaporation to available water in surface layer
! minus canopy transpiration from surface layer  (mm/s)
!----------------------------------------------------------------------!
evap_max (:)= zero
evap_max (1) = dz (1) * (theta (1) - watmin) / dt - &
               qflx_tran_veg_col * rootr_col (1,x,y)
evap_max (1) = MAX (zero, evap_max (1))
qflx_evap_grnd = MIN (evap_max (1), qflx_evap_grnd)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Limit transpiration from each layer to available water in layer (mm/s)
! Fix up...............****************
!----------------------------------------------------------------------!
DO I = 1, nlayers
  evap_max (I) = dz (I) * (theta (I) - watmin) / dt
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Evaporation rate from open water (mm/s). Not used yet.
!----------------------------------------------------------------------!
LEw = (desatdT * (Rnet - Rnets) + rho * cp * VDD0 / rac) &
      / (desatdT + gamma * (one + zero / rac))
qflx_ev_h2osfc = LEw * 1.0E3 / (rhow * lamb)
qflx_ev_h2osfc = zero !****
!----------------------------------------------------------------------!

!**********************************************************************!
! Start of CESM 'Infiltration'.
!**********************************************************************!

!----------------------------------------------------------------------!
DO I = 1, nlayers
  !--------------------------------------------------------------------!
  ! Porosity of soil.
  !--------------------------------------------------------------------!
  eff_porosity (I) = MAX (0.01, theta_s(I,x,y))
  !--------------------------------------------------------------------!
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Local evaporation                                               (mm/s)
!----------------------------------------------------------------------!
qflx_evap = qflx_evap_grnd
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Partition surface inputs between soil and h2osfc                (mm/s)
!----------------------------------------------------------------------!
qflx_in_soil = (one - frac_h2osfc) * (qflx_top_soil - qflx_surf)
qflx_in_h2osfc = frac_h2osfc * (qflx_top_soil - qflx_surf)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Remove evaporation                                              (mm/s)
!----------------------------------------------------------------------!
qflx_in_soil = qflx_in_soil - (one - frac_h2osfc) * qflx_evap
qflx_in_h2osfc = qflx_in_h2osfc - frac_h2osfc * qflx_ev_h2osfc
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Maximum infiltration capacity                                   (mm/s)
!----------------------------------------------------------------------!
qinmax = (one - fsat) * MINVAL (hksat (1:3,x,y))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Infiltration excess runoff -> h2osfc                            (mm/s)
!----------------------------------------------------------------------!
qflx_infl_excess = MAX (zero, qflx_in_soil - &
                   (one - frac_h2osfc) * qinmax)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Soil infiltration                                               (mm/s)
!----------------------------------------------------------------------!
qflx_infl = qflx_in_soil - qflx_infl_excess
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Shift infiltration excess from h2osfc input to surface runoff   (mm/s)
!----------------------------------------------------------------------!
qflx_surf = qflx_surf + qflx_infl_excess
qflx_infl_excess = zero
!----------------------------------------------------------------------!

!**********************************************************************!
! End of CESM 'Infiltration'.
!**********************************************************************!

!**********************************************************************!
! Start of CESM 'SoilWater'.
!**********************************************************************!

!----------------------------------------------------------------------!
! Water table depth                                                 (mm)
!----------------------------------------------------------------------!
zwtmm = 1000.0 * zwt (x,y)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Compute jwt index, the index of the first unsaturated
! layer, i.e. the layer right above the water table (-).
!----------------------------------------------------------------------!
jwt = nlayers
!----------------------------------------------------------------------!
! Allow jwt to equal zero when zwt is in top layer.
!----------------------------------------------------------------------!
DO I = 1, nlayers
 IF (zwt (x,y) <= (zi (I) / 1000.0)) THEN
    jwt = I - 1
    EXIT
  END IF
END DO
IF (INTERACTIVE) write (98,*) 'water table',iTIME,zwt(x,y),jwt
!----------------------------------------------------------------------!

            !----------------------------------------------------------!
            ! Calculate the equilibrium water content based on the water
            ! table depth.
            ! assumes sucsat == -psi_s
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              IF (zwtmm <= zi (I-1)) THEN
                !------------------------------------------------------!
                ! If layer below water table depth, all saturated.
                !------------------------------------------------------!
                vol_eq (I) = theta_s (I,x,y)
                !------------------------------------------------------!
              ELSE IF ((zwtmm < zi (I)) &
                      .AND. (zwtmm > zi (I-1))) THEN
                !------------------------------------------------------!
                ! If water table top in layer, use weighted average from
                ! the saturated part (depth > wtd) and the equilibrium
                ! solution for the rest of the layer.
                !------------------------------------------------------!
                tempi = one
                temp0 = ((((-psi_s (I,x,y)) + zwtmm - zi (I-1)) &
                        / (-psi_s (I,x,y)))) ** (one - one / bsw(I,x,y))
                voleq1 = psi_s (I,x,y) * theta_s (I,x,y) / &
                         (one - one / bsw (I,x,y)) / (zwtmm &
                         -zi (I-1)) * (tempi - temp0)
                vol_eq (I) = (voleq1 * (zwtmm - zi (I-1)) &
                             + theta_s (I,x,y) * (zi (I) - zwtmm)) &
                             / (zi (I) - zi (I-1))
                vol_eq (I) = MIN (theta_s (I,x,y), vol_eq (I))
                vol_eq (I) = MAX (vol_eq (I), zero)
                !------------------------------------------------------!
              ELSE
                !------------------------------------------------------!
                ! Water table below this layer.
                !------------------------------------------------------!
                tempi = (((-psi_s (I,x,y) + zwtmm - zi (I)) &
                        / (-psi_s (I,x,y)))) ** &
                       (1.0 - 1.0 / bsw (I,x,y))
                temp0 = (((-psi_s (I,x,y) + zwtmm - zi (I-1)) &
                        / (-psi_s (I,x,y)))) ** &
                       (1.0 - 1.0 / bsw (I,x,y))
                vol_eq (I) = psi_s (I,x,y) * theta_s (I,x,y) &
                             / (1.0 - 1.0 / bsw (I,x,y)) &
                             / (zi (I) - zi (I-1)) * (tempi - temp0)
                vol_eq (I) = MAX (vol_eq (I), 0.0)
                vol_eq (I) = MIN (theta_s (I,x,y), vol_eq (I))
                !------------------------------------------------------!
              END IF
              !--------------------------------------------------------!
              zq (I) = psi_s (I,x,y) * &
                       (MAX (vol_eq (I) / theta_s (I,x,y), 0.01)) &
                       ** (-bsw (I,x,y))
              zq (I) = MAX (smpmin, zq (I))
              !--------------------------------------------------------!
            END DO ! I = 1, nlayers
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! If water table is below soil column calculate zq for the
            ! 9th layer.
            !----------------------------------------------------------!
            I = nlayers
            !----------------------------------------------------------!
            IF (jwt == nlayers) THEN
              !--------------------------------------------------------!
              tempi = 1.0
              temp0 = (((-psi_s (I,x,y) + zwtmm - zi (I)) / &
                      (-psi_s (I,x,y)))) ** (1.0 - 1.0 / bsw (I,x,y))
              vol_eq (I+1) = psi_s (I,x,y) * theta_s (I,x,y) / &
                             (1.0 - 1.0 / bsw (I,x,y)) &
                             / (zwtmm - zi (I)) * (tempi - temp0)
              vol_eq (I+1) = MAX (vol_eq (I+1), 0.0)
              vol_eq (I+1) = MIN (theta_s (I,x,y), vol_eq (I+1))
              zq (I+1) = psi_s (I,x,y) * (MAX (vol_eq (I+1) &
                         / theta_s (I,x,y), 0.01)) ** (-bsw (I,x,y))
              zq (I+1) = MAX (smpmin, zq (I+1))
              !--------------------------------------------------------!
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Hydraulic conductivity and soil matric potential and their
            ! derivatives.
            ! I refers to bottom of layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              ! Compute saturated hydraulic conductivity at the
              ! interface as mean of the two layers (mm/s).
              !--------------------------------------------------------!
              ! "s" at interface of layer.
              !--------------------------------------------------------!
              s1 =  0.5 * (theta   (I) + &
                           theta   (MIN (nlayers,I+1))) / &
                   (0.5 * (theta_s (I,x,y) + &
                           theta_s (MIN (nlayers,I+1),x,y)))
              s1 = MIN (one, s1)
              s2 = hksat (I,x,y) * s1 ** (2.0 * bsw (I,x,y) + 2.0)
              !--------------------------------------------------------!
              ! Hydraulic conductivity at layer intervace (mm/s).
              ! Need to add ice impedance.
              !--------------------------------------------------------!
              hk (I) = s1 * s2
              !--------------------------------------------------------!
              ! d(hk)/d(vol_liq) (mm s-1) / (mm^3 mm^-3)
              !--------------------------------------------------------!
              dhkdw (I) = (2.0 * bsw (I,x,y) + 3.0) * s2 * &
                          (one / (theta_s (I,x,y) + &
                          theta_s (MIN (nlayers,I+1),x,y)))
              !--------------------------------------------------------!
              ! Compute matric potential and derivative based on liquid
              ! water only.
              !--------------------------------------------------------!
              s_node = MAX (theta (I) / theta_s (I,x,y), 0.01)
              s_node = MIN (one, s_node)
              !--------------------------------------------------------!
              ! Matric potentials and derivatives at each node depth 
              ! (mm).
              ! This is Eqn. 7.94 of Oleson et al. (2013).
              !--------------------------------------------------------!
              smp (I) = psi_s (I,x,y) * s_node ** (-bsw (I,x,y))
              smp (I) = MAX (smpmin, smp (I))
              !--------------------------------------------------------!
              dsmpdw (I) = (-bsw (I,x,y)) * smp (I) / &
                           (s_node * theta_s (I,x,y))
              !--------------------------------------------------------!
            END DO ! I = 1, nlayers
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Aquifer (9th) layer ('zmm' == zc).
            !----------------------------------------------------------!
            zc (nlayers+1) = 0.5 * (zwtmm + zc (nlayers))
            IF (jwt < nlayers) THEN
              dz (nlayers+1) = dz (nlayers)
            ELSE
              dz (nlayers+1) = zwtmm - zc (nlayers)
            ENDIF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Use tridiagonal solution to compute new soil water levels
            ! in each layer given fluxes.
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Set up vectors for top layer based on section 7.4.2.2 O13.
            !----------------------------------------------------------!
            I = 1
            qin (I) = qflx_infl
            den = (zc  (I+1) - zc  (I))
            dzq = (zq  (I+1) - zq  (I))
            num = (smp (I+1) - smp (I)) - dzq
            qout (I) = -hk (I) * num / den     ! -7.116
            dqodw1 (I) = -(-hk (I) * dsmpdw (I)   + num * dhkdw (I)) / &
                         den                   ! -7.119
            dqodw2 (I) = -( hk (I) * dsmpdw (I+1) + num * dhkdw (I)) / &
                         den                   ! -7.120
            rmx (I) = qin (I) - qout (I) - qflx_tran_veg_col * &
                      rootr_col (I,x,y)
            amx (I) = zero                     !  7.136
            bmx (I) = dz (I) / dt + dqodw1 (I) ! -7.137
            cmx (I) = dqodw2 (I)               ! -7.138
            !----------------------------------------------------------!
            ! Nodes I = 2, to I = nlayers-1.
            !----------------------------------------------------------!
            DO I = 2, nlayers - 1
              !--------------------------------------------------------!
              den = zc  (I) - zc  (I-1)
              dzq = zq  (I) - zq  (I-1)
              num = smp (I) - smp (I-1) - dzq
              qin (I) = -hk (I-1) * num / den
              dqidw0 (I) = -(-hk (I-1) * dsmpdw (I-1) + &
                           num * dhkdw (I-1)) / den           ! -7.117
              dqidw1 (I) = -( hk (I-1) * dsmpdw (I)   + &
                           num * dhkdw (I-1)) / den           ! 
              den = zc (I+1) - zc (I)
              dzq = zq (I+1) - zq (I)
              num = (smp (I+1) - smp (I)) - dzq
              qout (I) = -hk (I) * num / den
              dqodw1 (I) = -(-hk (I) * dsmpdw (I)   + &
                           num * dhkdw (I)) / den             ! 
              dqodw2 (I) = -( hk (I) * dsmpdw (I+1) + &
                           num * dhkdw (I)) / den
              rmx (I) = qin (I) - qout (I) - qflx_tran_veg_col * &
                        rootr_col (I,x,y)
              amx (I) = -dqidw0 (I)                           ! -7.140
              bmx (I) = dz (I) / dt - dqidw1 (I) + dqodw1 (I) ! -7.141
              cmx (I) = dqodw2 (I)
              !--------------------------------------------------------!
            END DO ! I = 2, nlayers - 1
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Now compute matrix elements for fluxes at base of soil
            ! column.
            !----------------------------------------------------------!
            I = nlayers
            !----------------------------------------------------------!
            IF (I > jwt) THEN ! Water table is in soil column.
              !--------------------------------------------------------!
              den = zc  (I) - zc  (I-1)
              dzq = zq  (I) - zq  (I-1)
              num   = smp (I) - smp (I-1) - dzq
              qin (I) = -hk (I-1) * num / den
              dqidw0 (I) = -(-hk (I-1) * dsmpdw (I-1) + &
                           num * dhkdw (I-1)) / den
              dqidw1 (I) = -( hk (I-1) * dsmpdw (I) + &
                           num * dhkdw (I-1)) / den
              qout   (I) = zero
              dqodw1 (I) = zero
              rmx (I) =  qin (I) - qout (I) - qflx_tran_veg_col * &
                         rootr_col (I,x,y)
              amx (I) = -dqidw0 (I)
              bmx (I) =  dz (I) / dt - dqidw1 (I) + dqodw1 (I)
              cmx (I) =  zero
              !--------------------------------------------------------!
              ! Next set up aquifer layer, hydrologically inactive.
              !--------------------------------------------------------!
              rmx (I+1) = zero
              amx (I+1) = zero
              bmx (I+1) = dz (I+1) / dt
              cmx (I+1) = zero
              !--------------------------------------------------------!
            ELSE ! Water table is below soil column.
              !--------------------------------------------------------!
              ! Compute aquifer soil moisture as average of layer 9
              ! and saturation.
              !--------------------------------------------------------!
              s_node = MAX (0.5 * (one + theta (I) / &
                       theta_s (I,x,y)), 0.01)
              s_node = MIN (one, s_node)
              !--------------------------------------------------------!
              ! Compute smp for aquifer layer (mm).
              !--------------------------------------------------------!
              smp1 = psi_s (I,x,y) * s_node ** (-bsw (I,x,y))
              smp1 = MAX (smpmin, smp1)
              !--------------------------------------------------------!
              ! Compute dsmpdw for aquifer layer.
              !--------------------------------------------------------!
              dsmpdw1 = -bsw (I,x,y) * smp1 / (s_node * theta_s (I,x,y))
              !--------------------------------------------------------!
              ! First set up bottom layer of soil column.
              !--------------------------------------------------------!
              den = zc  (I) - zc  (I-1)
              dzq = zq  (I) - zq  (I-1)
              num = smp (I) - smp (I-1) - dzq
              qin (I) = -hk (I-1) * num / den
              dqidw0 (I) = -(-hk(I-1) * dsmpdw (I-1) + &
                           num * dhkdw (I-1)) / den
              dqidw1 (I) = -( hk(I-1) * dsmpdw (I)   + &
                           num * dhkdw (I-1)) / den
              den = zc (I+1) - zc (I)
              dzq = zq (I+1) - zq (I)
              num = smp1 - smp (I) - dzq
              qout (I) = -hk (I) * num / den
              dqodw1 (I) = -(-hk (I) * dsmpdw (I) + num * dhkdw (I)) &
                           / den
              dqodw2 (I) = -( hk (I) * dsmpdw1 + num * dhkdw (I)) &
                           / den
              !--------------------------------------------------------!

              !--------------------------------------------------------!
              rmx (I) = qin (I) - qout (I) - qflx_tran_veg_col * &
                        rootr_col (I,x,y)
              amx (I) = -dqidw0 (I)
              bmx (I) = dz (I) / dt - dqidw1 (I) + dqodw1 (I)
              cmx (I) = dqodw2 (I)
              !--------------------------------------------------------!

              !--------------------------------------------------------!
              ! Next set up aquifer layer; den/num unchanged, qin=qout.
              !--------------------------------------------------------!
              qin (I+1) = qout (I)
              dqidw0 (I+1) = -(-hk (I) * dsmpdw (I) + num * dhkdw (I)) &
                             / den
              dqidw1 (I+1) = -( hk (I) * dsmpdw1    + num * dhkdw (I)) &
                             / den
              qout   (I+1) = zero ! Zero-flow bottom boundary condition.
              dqodw1 (I+1) = zero ! Zero-flow bottom boundary condition.
              rmx (I+1) = qin (I+1) - qout (I+1)
              amx (I+1) = -dqidw0 (I+1)
              bmx (I+1) = dz (I+1) / dt - dqidw1 (I+1) + dqodw1 (I+1)
              cmx (I+1) = zero
              !--------------------------------------------------------!

            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Solve for dwat2 using tridiagonal system of equations.
            ! Approach taken from Press et al. (1989), 2.6.
            !----------------------------------------------------------!
            IF (bmx (1) == 0.0) THEN
              WRITE (*,*) 'Problem with tridiagonal 1.'
              WRITE (*,*) 'DiTIME = ',iTIME-time_BOY (jyear-1859)+1
              WRITE (*,*) 'my_id x y ',my_id,x,y
              WRITE (*,*) 'lon, lat, iTIME ',lon(x),lat(y),iTIME
              STOP
            END IF
            BET = bmx (1)
            dwat2 (1) = rmx (1) / BET
            DO I = 2, nlayers + 1
              GAM (I) = cmx (I-1) / BET
              BET = bmx (I) - amx (I) * GAM (I)
              IF (BET == 0.0) THEN
                WRITE (*,*) I,bmx(I),amx(I),GAM(I),amx(I)*GAM(I)
                WRITE (*,*) 'Problem with tridiagonal 2.'
                WRITE (*,*) 'DiTIME = ',iTIME-time_BOY (jyear-1859)+1
                WRITE (*,*) 'my_id x y ',my_id,x,y
                WRITE (*,*) 'lon, lat, iTIME ',lon(x),lat(y),iTIME
                STOP 'Problem with tridiagonal 2.'
              END IF
              dwat2 (I) = (rmx (I) - amx (I) * dwat2 (I-1)) / BET
            END DO
            !----------------------------------------------------------!
            DO I = nlayers, 1, -1
              dwat2 (I) = dwat2 (I) - GAM (I+1) * dwat2 (I+1)
            END DO
            !----------------------------------------------------------!
            ! Set dwat.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              dwat (I) = dwat2 (I)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Renew the mass of liquid water also compute qcharge from
            ! dwat in aquifer layer and update in drainage for case
            ! jwt < nlayers.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + &
                                   dwat2 (I) * dz (I)
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Calculate qcharge for case jwt < nlayers.
            !----------------------------------------------------------!
            IF (jwt < nlayers) THEN
              !--------------------------------------------------------!
              ! Water head at the water table depth                 (mm)
              ! zero since wh_zwt = psi_s - zq_zwt, where zq_zwt = psi_s
              !--------------------------------------------------------!
              wh_zwt = zero
              !--------------------------------------------------------!
              ! Recharge rate qcharge to groundwater (positive to
              ! aquifer).
              !--------------------------------------------------------!
              s_node = MAX (theta (jwt+1) &
                       / theta_s (jwt+1,x,y), 0.01)
              s1 = MIN (one, s_node)
              !--------------------------------------------------------!
              ! Unsaturated hk.
              !--------------------------------------------------------!
              ka = hksat (jwt+1,x,y) * s1 ** &
                  (2.0 * bsw (jwt+1,x,y) + 3.0)
              !--------------------------------------------------------!
              ! Recharge rate qcharge to groundwater (positive to
              ! aquifer).
              !--------------------------------------------------------!
              smp1 = MAX (smpmin, smp (MAX (1, jwt)))
              wh   = smp1 - zq (MAX (1, jwt))
              !--------------------------------------------------------!
              ! scs: original formulation.
              !--------------------------------------------------------!
              IF (jwt == 0) THEN ! Water table at surface.
                qcharge = -ka * (wh_zwt - wh) / (zwtmm + one)
              ELSE
                ! scs: 1/2, assuming flux is at zwt interface,
                ! saturation deeper than zwt.
                qcharge = -ka * (wh_zwt - wh)/ &
                          ((zwtmm - zc (jwt)) * 2.0)
              END IF
              !--------------------------------------------------------!
              ! Limit qcharge (for the first several timesteps).
              !--------------------------------------------------------!
              qcharge = MAX (-10.0 / dt, qcharge)
              qcharge = MIN ( 10.0 / dt, qcharge)
              !--------------------------------------------------------!
            ELSE
              !--------------------------------------------------------!
              ! If water table is below soil column, compute qcharge
              ! from dwat2 (9)                                    (mm/s)
              !--------------------------------------------------------!
              qcharge = dwat2 (nlayers+1) * dz (nlayers+1) / dt
              !--------------------------------------------------------!
            END IF
            !----------------------------------------------------------!

!**********************************************************************!
! End of CESM 'SoilWater'.
!**********************************************************************!

!**********************************************************************!
! Start ofCESM 'Drainage'
!**********************************************************************!

            !----------------------------------------------------------!
            ! Drainage routines.
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! The layer index of the first unsaturated layer, i.e. the
            ! layer right above the water table.
            !----------------------------------------------------------!
            jwt = nlayers
            ! Allow jwt to equal zero when zwt is in top layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              IF (zwt (x,y) <= (zi (I) / 1000.0)) THEN
                jwt = I - 1
                EXIT
              END IF
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Use analytical expression for aquifer specific yield   (-)
            !----------------------------------------------------------!
            rous = theta_s (nlayers,x,y) &
                   * (one - (one + zwtmm / (-psi_s (nlayers,x,y))) &
                   ** (-one / bsw (nlayers,x,y)))
            rous = MAX (rous, 0.02)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water table is below the soil column.
            !----------------------------------------------------------!
            IF (jwt == nlayers) THEN
              !--------------------------------------------------------!
              ! Update water in the unconfined aquifer (mm).
              !--------------------------------------------------------!
              wa (x,y) = wa (x,y) + qcharge * dt
              zwt (x,y) = zwt (x,y) - (qcharge * dt) / 1000.0 / rous
              !--------------------------------------------------------!
            ELSE
              !--------------------------------------------------------!
              ! Water table within soil layers 1-8.
              ! Try to raise water table to account for qcharge.
              !--------------------------------------------------------!
              qcharge_tot = qcharge * dt
              !--------------------------------------------------------!
              IF (qcharge_tot > zero) THEN ! Rising water table.
                DO I = jwt+1, 1, -1
                  ! Use analytical expression for specific yield.
                  s_y = theta_s (I,x,y) * (one - (one + zwtmm &
                        / (-psi_s (I,x,y))) ** (-one / bsw (I,x,y)))
                  s_y = MAX (s_y, 0.02)
                  qcharge_layer = MIN (qcharge_tot, &
                                       s_y * (zwtmm - zi (I-1)))
                  qcharge_layer = MAX (qcharge_layer, zero)
                  IF (s_y > zero) zwt (x,y) = &
                    zwt (x,y) - qcharge_layer / s_y / 1000.0
                  qcharge_tot = qcharge_tot - qcharge_layer
                  IF (qcharge_tot <= zero) EXIT
                END DO
              !--------------------------------------------------------!
              ELSE ! Deepening water table (negative qcharge).
              !--------------------------------------------------------!
                DO I = jwt+1, nlayers
                  ! Use analytical expression for specific yield.
                  s_y = theta_s (I,x,y) * (one - (one + zwtmm &
                        / (-psi_s (I,x,y))) ** (-one / bsw (I,x,y)))
                  s_y = MAX (s_y, 0.02)
                  qcharge_layer = MAX (qcharge_tot, &
                                  -s_y * (zi (I) - zwtmm))
                  qcharge_layer = MIN (qcharge_layer, zero)
                  qcharge_tot = qcharge_tot - qcharge_layer
                  IF (qcharge_tot >= zero) THEN
                    zwt (x,y) = zwt (x,y) - qcharge_layer / s_y / 1000.0
                    EXIT
                  ELSE
                    zwt (x,y) = zi (I) / 1000.0
                  END IF
                END DO
                IF (qcharge_tot > zero) zwt (x,y) = &
                  zwt (x,y) - qcharge_tot / 1000.0 / rous
              END IF
              !--------------------------------------------------------!
              ! Recompute jwt for following calculations.
              ! Allow jwt to equal zero when zwt is in top layer.
              !--------------------------------------------------------!
              jwt = nlayers
              !--------------------------------------------------------!
              DO I = 1, nlayers
                IF (zwt (x,y) <= (zi (I) / 1000.0)) THEN
                  jwt = I - 1
                  EXIT
                END IF
              END DO
              !--------------------------------------------------------!
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Recompute water table height (mm).
            !----------------------------------------------------------!
            zwtmm = 1000.0 * zwt (x,y)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Baseflow.
            ! 'Original formulation' does not use perched water table
            ! drainage, so go with that for now, only computing
            ! topographic runoff.
            !----------------------------------------------------------!
            rsub_top_max = 5.5E-3
            !----------------------------------------------------------!
            ! Subsurface runoff - topographic control (mm/s).
            !----------------------------------------------------------!
            rsub_top = rsub_top_max * EXP (-fff * zwt (x,y))
            !----------------------------------------------------------!
            ! Use analytical expression for aquifer specific yield (-).
            !----------------------------------------------------------!
            rous = theta_s (nlayers,x,y) &
                   * (one - (one + zwtmm / (-psi_s (nlayers,x,y))) &
                   ** (-one / bsw (nlayers,x,y)))
            rous = MAX (rous, 0.02)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Zero-out diagnostic of topographic drainage from each
            ! layer.
            !----------------------------------------------------------!
            rnff (:) = 0.0
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water table is below the soil column.
            !----------------------------------------------------------!
            IF (jwt == nlayers) THEN
              !--------------------------------------------------------!
              wa (x,y) = wa (x,y) - rsub_top * dt
              zwt (x,y) = zwt (x,y) + (rsub_top * dt) / 1000.0 / rous
              h2osoi_liq (nlayers,x,y) = h2osoi_liq (nlayers,x,y) + &
                                         MAX (0.0, (wa (x,y) - 5000.0))
              wa (x,y) = MIN (wa (x,y), 5000.0)
              !--------------------------------------------------------!
              ! Topographic drainage from aquifer.
              !--------------------------------------------------------!
              rnff (nlayers+1) = rsub_top
              !--------------------------------------------------------!
            ELSE ! Water table within soil layers 1-8.
              !--------------------------------------------------------!
              ! Now remove water via rsub_top.
              !--------------------------------------------------------!
              rsub_top_tot = -rsub_top * dt
              !--------------------------------------------------------!
              ! Should never be positive, but include for completeness.
              !--------------------------------------------------------!
              IF (rsub_top_tot > zero) THEN ! Rising water table.
                WRITE (*,*) 'rsub_top_tot is positive in drainage'
                WRITE (*,*) 'HYBRID9 is stopping'
                STOP
              !--------------------------------------------------------!
              ELSE ! Deepening water table.
              !--------------------------------------------------------!
              DO I = jwt + 1, nlayers
                ! Use analytical expression for specific yield.
                s_y = theta_s (I,x,y) &
                      * (one - (one + zwtmm / (-psi_s (I,x,y))) &
                      ** (-one / bsw (I,x,y)))
                s_y = MAX (s_y, 0.02)
                
                rsub_top_layer = MAX (rsub_top_tot, &
                                 -(s_y * (zi (I) - zwtmm)))
                rsub_top_layer = MIN (rsub_top_layer, zero)
                ! N.B. rsub_top_layer is a negative value.
                h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + rsub_top_layer

                rnff (I) = -rsub_top_layer ! Layer drainage.

                rsub_top_tot = rsub_top_tot - rsub_top_layer

                IF (rsub_top_tot >= zero) THEN
                  zwt (x,y) = zwt (x,y) - rsub_top_layer / s_y / 1000.0
                  EXIT
                ELSE
                  zwt (x,y) = zi (I) / 1000.0
                ENDIF
              END DO ! I = jwt + 1, nlayers
              ! Remove residual rsub_top from aquifer.
              zwt (x,y) = zwt (x,y) - rsub_top_tot / 1000.0 / rous
              wa  (x,y) = wa  (x,y) + rsub_top_tot
              rnff (nlayers+1) = rnff (nlayers+1) - rsub_top_tot
              END IF ! (rsub_top_tot > zero)
              !--------------------------------------------------------!

              !--------------------------------------------------------!
              ! Recompute jwt. Allow it to equal zero when zwt is in top
              ! layer.
              !--------------------------------------------------------!
              jwt = nlayers
              DO I = 1, nlayers
                IF (zwt (x,y) <= (zi (I) / 1000.0)) then
                    jwt = I - 1
                    EXIT
                END IF
              END DO
              !--------------------------------------------------------!
            END IF ! End of jwt IF construct.
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            zwt (x,y) = MAX ( 0.0, zwt (x,y))
            zwt (x,y) = MIN (80.0, zwt (x,y))
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Excessive water above saturation added to the above
            ! unsaturated layer like a bucket if column fully saturated,
            ! excess water goes to runoff.
            !----------------------------------------------------------!
            DO I = nlayers, 2, -1
              xsi = MAX (h2osoi_liq (I,x,y) - eff_porosity (I) * &
                    dz (I), zero)
              h2osoi_liq (I,x,y) = MIN (eff_porosity (I) * &
                                   dz (I), h2osoi_liq (I,x,y))
              h2osoi_liq(I-1,x,y) = h2osoi_liq(I-1,x,y) + xsi
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Remove water from top layer if above saturation, adding to
            ! drainage. Not sure what 'pondmx' is, so left out.
            !----------------------------------------------------------!
            xs1 = MAX (MAX (h2osoi_liq (1,x,y), zero) - &
                  MAX (zero, (theta_s (1,x,y) * dz (1))), zero)
            h2osoi_liq (1,x,y) = MIN (MAX (zero, &
                                 theta_s (1,x,y) * dz (1)), &
                                 h2osoi_liq (1,x,y))
            !----------------------------------------------------------!
            ! Send water to drainage.
            !----------------------------------------------------------!
            qflx_rsub_sat = xs1 / dt
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Limit h2osoi_liq to be greater than or equal to watmin.
            ! Get water needed to bring h2osoi_liq equal watmin from
            ! lower layer. If insufficient water in soil layers, get
            ! water from aquifer layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers - 1
              IF (h2osoi_liq (I,x,y) < watmin) THEN
                xs = watmin - h2osoi_liq (I,x,y)
                ! Deepen water table if water is passed from below zwt
                ! layer.
                IF (I == jwt) THEN
                  zwt (x,y) = zwt (x,y) + xs / eff_porosity (I) / 1000.0
                END IF
              ELSE
                xs = zero
              END IF
              h2osoi_liq (I  ,x,y) = h2osoi_liq (I  ,x,y) + xs
              h2osoi_liq (I+1,x,y) = h2osoi_liq (I+1,x,y) - xs
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Get water for bottom layer from layers above if possible.
            !----------------------------------------------------------!
            I = nlayers
            IF (h2osoi_liq (I,x,y) < watmin) THEN
              xs = watmin - h2osoi_liq (I,x,y)
              searchforwater: DO J = nlayers - 1, 1, -1
                available_h2osoi_liq = MAX (h2osoi_liq (J,x,y) -&
                                       watmin - xs, zero)
                IF (available_h2osoi_liq >= xs) THEN
                  h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + xs
                  h2osoi_liq (J,x,y) = h2osoi_liq (J,x,y) - xs
                  xs = zero
                  EXIT searchforwater
                ELSE
                  h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + &
                                       available_h2osoi_liq
                  h2osoi_liq (J,x,y) = h2osoi_liq (J,x,y) - &
                                       available_h2osoi_liq
                  xs = xs - available_h2osoi_liq
                END IF
              END DO searchforwater
            ELSE
               xs = zero
            END IF
            !----------------------------------------------------------!
            ! Needed in case there is no water to be found.
            !----------------------------------------------------------!
            h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + xs
            !----------------------------------------------------------!
            ! Instead of removing water from aquifer where it eventually
            ! shows up as excess drainage to the ocean, take it back out
            ! of drainage.
            !----------------------------------------------------------!
            rsub_top = rsub_top - xs / dt !****
            !----------------------------------------------------------!

!**********************************************************************!
! End ofCESM 'Drainage'
!**********************************************************************!

            !----------------------------------------------------------!
            ! Sum all losses from soil column (mm).
            !----------------------------------------------------------!
            w1 = ((1.0 - frac_h2osfc) * (qflx_surf + qflx_evap_grnd &
                 + qflx_tran_veg_col) + &
                 rsub_top + qflx_rsub_sat) * dt + wa (x,y)
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              ! New total water for diagnostics (mm).
              !--------------------------------------------------------!
              w1 = w1 + h2osoi_liq (I,x,y)
              !--------------------------------------------------------!
              ! Diagnose volumetric water (mm^3/mm^3).
              !--------------------------------------------------------!
              theta (I) = MAX (h2osoi_liq (I,x,y), 1.0E-6) &
                          / (dz (I) * rhow / 1000.0)
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            !IF (INTERACTIVE) THEN
            !----------------------------------------------------------!
            ! Check for conservation of water.
            !----------------------------------------------------------!
            IF (ABS (w1 - w0) > 0.1) THEN
              WRITE (*,*)
              WRITE (*,*) 'Problem in HYDROLOGY'
              WRITE (*,*) 'Water imbalance > 0.1 mm ',w1-w0
              WRITE (*,*) 'DiTIME = ',iTIME-time_BOY (jyear-1859)+1
              WRITE (*,*) 'my_id x y ',my_id,x,y
              WRITE (*,*) 'lon, lat, iTIME ',lon(x),lat(y),iTIME
              WRITE (*,*) 'theta sum = ',SUM (theta(1:nlayers))
              WRITE (*,*) 'fsat      = ',fsat
              WRITE (*,*) '         I      rnff'
              DO I = 1, nlayers + 1
                WRITE (*,*) I,dt*rnff(I)
              END DO
              WRITE (*,*) 'SUM(rnff)  = ',dt*SUM (rnff(1:nlayers+1))
              WRITE (*,*) 'qflx_surf qflx_evap_grnd = ',&
                          dt*qflx_surf,dt*qflx_evap_grnd
              WRITE (*,*) 'qflx_tran_veg_col= ',dt*qflx_tran_veg_col
              WRITE (*,*) 'qout (nlayers) = ',qout(nlayers)*dt
              WRITE (*,*) 'zwt = ',zwt(x,y)
              WRITE (*,*) 'qcharge = ',qcharge*dt
              WRITE (*,*) 'pr = ',dt*pr(x,y,iT)
              WRITE (*,*) 'forc_rain = ',dt*forc_rain
              WRITE (*,*) 'beta = ',beta
              WRITE (*,*) 'wa   =',wa(x,y)
              WRITE (*,*) 'theta theta_m h2osoi_liq h2osoi_liq_ma'
              DO I = 1, nlayers
                WRITE (*,*) I,theta(I),theta_m(I,x,y), &
                            h2osoi_liq (I,x,y),h2osoi_liq_ma (I,x,y)
              END DO
              STOP
            END IF ! (ABS (w1 - w0) > 0.1)
            !----------------------------------------------------------!
            !END IF ! INTERACTIVE
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated within the day.
            !----------------------------------------------------------!
            rnf_sum  = rnf_sum  + qflx_surf * dt !********************
            rnf_sum = rnf_sum + rsub_top * dt !**********************
            !----------------------------------------------------------!

!----------------------------------------------------------------------!
END SUBROUTINE HYDROLOGY
!======================================================================!
