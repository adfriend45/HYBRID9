!======================================================================!
SUBROUTINE GROW
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Increment biomass over day.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!USE NETCDF  ! Enable access to the library of netCDF routines.
USE CONTROl ! Control variables.
USE SHARED  ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
REAL :: grow_plant_mass ! Plant mass growth                  (g[DM]/day)
REAL :: grow_foliage_mass
REAL :: loss_plant_mass ! Plant mass losses                  (g[DM]/day)
REAL :: loss_foliage_mass
REAL :: dplant_mass  ! Plant mass change                     (g[DM]/day)
REAL :: dplant_foliage_mass ! Plant foliage mass change      (g[DM]/day)
REAL :: C_labile_in  ! C input into plant labile pool         (g[C]/day)
REAL :: N_labile_in  ! N input into plant labile pool         (g[N]/day)
REAL :: P_labile_in  ! P input into plant labile pool         (g[P]/day)
REAL :: C_labile_out ! C output from plant labile pool        (g[C]/day)
REAL :: N_labile_out ! N output from plant labile pool        (g[N]/day)
REAL :: P_labile_out ! P output from plant labile pool        (g[P]/day)
REAL :: dC_labile    ! Plant labile C change                  (g[C]/day)
REAL :: dN_labile    ! Plant labile N change                  (g[N]/day)
REAL :: dP_labile    ! Plant labile P change                  (g[P]/day)
REAL :: dwater       ! Change in soil water for plant growth    (mm/day)
REAL :: dLAI         ! Change in LAI                       (m^2/m^2/day)
REAL :: decay        ! For rooting                                   (-)
REAL :: w_i_save
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Grow some plants!
!----------------------------------------------------------------------!
! N.B. Main radial file model is
! /store/MODELS/HXd/SOURCE/RINGS.f90
!----------------------------------------------------------------------!
! .
!----------------------------------------------------------------------!
! Soil moisture control from CESM Eqn. 8.27.
! Currently is a mean value across all vegetation in plot.
! -1.5 MPa is equivalant to about -150,000 mm
!*** inc litter decomp rate
!*** grav on smp
!----------------------------------------------------------------------!
w_i_save = zero
DO I = 1, nlayers
  w_i = (-150000.0 - smp (I)) / (-150000.0 - (-50000.0))
  w_i = MAX (zero, w_i)
  w_i = MIN (one , w_i)
  w_i_save = w_i_save + rootr_col (I,x,y) * w_i
END DO
w_i = w_i_save
!----------------------------------------------------------------------!
! Temperature effect on growth from Hayat et al. (2017), Eqn. 19.
!----------------------------------------------------------------------!
IF ((tas (x,y,iT) - tf) > 18.0) THEN
  fT = one - (ABS (tas (x,y,iT) - tf - 18.0) / 21.0) ** 2
ELSE
  fT = one - (ABS (tas (x,y,iT) - tf - 18.0) / 25.0) ** 2
  fT = MAX (zero, fT)
  fT = MIN (one , fT)
END IF
!----------------------------------------------------------------------!
! Initialise plot relative root distribution over soil layers        (-)
!----------------------------------------------------------------------!
rootr_col (:,x,y) = zero
!----------------------------------------------------------------------!
npp = zero ! For diagnostics.
!----------------------------------------------------------------------!
! Loop over plants in plot.
!----------------------------------------------------------------------!
DO K = 1, nplants (x,y)
  !--------------------------------------------------------------------!
  ! Generalised plant type index                                     (-)
  !--------------------------------------------------------------------!
  iGPT = 1
  !--------------------------------------------------------------------!
  ! Growth currently effectively on ground area basis        (g[DM]/day)
  !--------------------------------------------------------------------!
  grow_plant_mass = (1000.0 / 365.0) * w_i * fT
  grow_foliage_mass = grow_plant_mass / 3.3
  !--------------------------------------------------------------------!
  ! Growth results in a demand for carbon, which is met by
  ! photosynthetic machinery capturing light and using this
  ! energy to reduce CO2 to carbohydrate.
  ! Add a labile carbon pool as a reserve for growth and as
  ! signal to leaves to grow and photosynthesise.
  !--------------------------------------------------------------------!
  ! So, growth creates a demand for C, N, and P from the
  ! labile pools, depending on their concentrations in the
  ! wood. Take these to be fixed for now from:
  !
  !--------------------------------------------------------------------!
  C_labile_in = 0.0
  N_labile_in = 0.0
  P_labile_in = 0.0
  !--------------------------------------------------------------------!
  C_labile_out = 0.47 * grow_plant_mass
  N_labile_out = 0.02 * C_labile_out
  P_labile_out = 0.02 * N_labile_out
  !--------------------------------------------------------------------!
  dC_labile = C_labile_in - C_labile_out
  dN_labile = N_labile_in - N_labile_out
  dP_labile = P_labile_in - P_labile_out
  !--------------------------------------------------------------------!

  ! opening stomata at top of stem section, which results in
  ! demand for water, which is met by pulling water out of
  ! the soil. Less soil water reduces growth, reducing
  ! demand, and hence closing stomata. Need store to
  ! provide communication from growth demand to stomata.
  ! Stomata may also close to reduce chance of cavitation.
  ! Hence they respond to VPD and soil water directly.
  ! N and P also taken up as products of demand
  ! for them in growth and to fix C.
  ! Also have sink for reproduction as a priority.
  ! Objective of growth is to reproduce.
  !--------------------------------------------------------------------!
  ! Simplest: assume WUE of 3 g[DM]/mm, and only root in
  ! top layer.
  !--------------------------------------------------------------------!
  !dwater = MIN (h2osoi_liq (1,x,y), dplant_mass / 3.0)
  !dplant_mass = 3.0 * dwater
  loss_plant_mass = (0.1 / 365.0) * plant_mass (K,x,y)
  !loss_foliage_mass = (8.0 / 365.0) * plant_foliage_mass (K,x,y)
  loss_foliage_mass = (1.0 / 365.0) * plant_foliage_mass (K,x,y) / &
                      MIN(one,MAX (0.01,w_i))
  if (w_i < 0.6) loss_foliage_mass = 0.1 * plant_foliage_mass (K,x,y)
  !IF (w_i < 0.5) THEN
  !  loss_foliage_mass = 0.1 * plant_foliage_mass (K,x,y) / &
  !                      MIN(one,MAX (0.01,w_i*2.0))
  !ELSE
  !  loss_foliage_mass = 0.01 * plant_foliage_mass (K,x,y)
  !END IF
  dplant_mass = grow_plant_mass - loss_plant_mass
  dplant_foliage_mass = grow_foliage_mass - loss_foliage_mass
  !--------------------------------------------------------------------!
  plant_mass (K,x,y) = plant_mass (K,x,y) + dplant_mass
  plant_foliage_mass (K,x,y) = plant_foliage_mass (K,x,y) + &
                               dplant_foliage_mass
  !--------------------------------------------------------------------!
  ! Diagnose plant lengthi assuming all is a cylinder of mass with
  ! length, and with some portion below-ground                      (mm)
  !--------------------------------------------------------------------!
  plant_length (K,x,y) = (400.0 * plant_mass (K,x,y) / &
                         3.142E-3) ** (one / 3.0)
  !--------------------------------------------------------------------!
  ! Sum foliage area over ind's for total LAI of canopy        (m^2/m^2)
  ! For now assume all per unit m^2.
  !--------------------------------------------------------------------!
  dLAI = dplant_foliage_mass * sla (iGPT)
  LAI (x,y) = LAI (x,y) + dLAI
  LAI (x,y) = MAX (0.001, LAI (x,y))
  !--------------------------------------------------------------------!
  ! Update LAI of litter layer (m^2/m^2)
  !--------------------------------------------------------------------!
  LAI_litter (x,y) = LAI_litter (x,y) + MAX (zero, dLAI)
  !--------------------------------------------------------------------!
  ! Rooting depth                                                   (mm)
  !--------------------------------------------------------------------!
  rdepth (K,x,y) = 0.3 * plant_length (K,x,y)
  !--------------------------------------------------------------------!
  ! Rooting distribution (based on Baldocchi et al., 2004).
  ! Assumes 90% of roots are within rdepth.
  !--------------------------------------------------------------------!
  decay = EXP (LOG (0.1) / (rdepth (K,x,y) / 10.0))
  !--------------------------------------------------------------------!
  DO I = 1, nlayers
    rootr_col (I,x,y) = rootr_col (I,x,y) + &
                        (1.0 - decay ** (zi (I) / 10.0)) - &
                        (1.0 - decay ** (zi (I-1) / 10.0))
  END DO
  !--------------------------------------------------------------------!
  ! Sum plant growth to give annual NPP                       (g[DM]/yr)
  !--------------------------------------------------------------------!
  npp = npp + dplant_mass
  !--------------------------------------------------------------------!
END DO ! K = 1, nplants (x,y)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Normalise root fractions over individuals.
!----------------------------------------------------------------------!
IF (nplants (x,y) > 1) rootr_col (:,x,y) = &
                       rootr_col (:,x,y) / FLOAT (nplants (x,y))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Decay of litter layer (m^2/m^2)
!----------------------------------------------------------------------!
LAI_litter (x,y) = LAI_litter (x,y) - 0.02 * LAI_litter (x,y)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END SUBROUTINE GROW
!======================================================================!
