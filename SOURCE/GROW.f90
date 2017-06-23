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

REAL :: w_i ! Soil moisture contraint on plant growth                (-)
REAL :: fT  ! Temperature constraint on plant growth                 (-)
REAL :: dplant_mass         ! Plant mass change              (g[DM]/day)
REAL :: C_labile_in    ! C input into plant labile pool       (g[C]/day)
REAL :: N_labile_in    ! N input into plant labile pool       (g[N]/day)
REAL :: P_labile_in    ! P input into plant labile pool       (g[P]/day)
REAL :: C_labile_out   ! C output from plant labile pool      (g[C]/day)
REAL :: N_labile_out   ! N output from plant labile pool      (g[N]/day)
REAL :: P_labile_out   ! P output from plant labile pool      (g[P]/day)
REAL :: dC_labile           ! Plant labile C change           (g[C]/day)
REAL :: dN_labile           ! Plant labile N change           (g[N]/day)
REAL :: dP_labile           ! Plant labile P change           (g[P]/day)
REAL :: dwater ! Change in soil water for plant growth          (mm/day)

            !----------------------------------------------------------!
            ! Grow some plants!
            !----------------------------------------------------------!
            ! N.B. Main radial file model is
            ! /store/MODELS/HXd/SOURCE/RINGS.f90
            ! -1.5 MPa is equivalant to about -150,000 mm
            !----------------------------------------------------------!
            ! CESM Eqn. 8.27.
            !----------------------------------------------------------!
            w_i = (-150000.0 - smp (1)) / (-150000.0 - (-50000.0))
            w_i = MAX (zero, w_i)
            w_i = MIN (one , w_i)
            !----------------------------------------------------------!
            ! Temperature effect on growth from Hayat et al. (2017),
            ! Eqn. 19.
            !----------------------------------------------------------!
            IF ((tas (x,y,iT) - tf) > 18.0) THEN
              fT = 1.0 - (ABS (tas (x,y,iT) - tf - 18.0) / 21.0) ** 2
            ELSE
              fT = 1.0 - (ABS (tas (x,y,iT) - tf - 18.0) / 25.0) ** 2
              fT = MAX (zero, fT)
              fT = MIN (one , fT)
            END IF
            !----------------------------------------------------------!
            npp = 0.0
            DO K = 1, nplants (x,y)
              !--------------------------------------------------------!
              ! Effectively on ground area basis (g[DM]/day)
              !--------------------------------------------------------!
              dplant_mass = (3000.0 / 365.0) * (w_i ** 0.5) * &
                            (fT ** 0.5)
              !--------------------------------------------------------!
              ! Growth results in a demand for carbon, which is met by
              ! photosynthetic machinery capturing light and using this
              ! energy to reduce CO2 to carbohydrate.
              ! Add a labile carbon pool as a reserve for growth and as
              ! signal to leaves to grow and photosynthesise.
              !--------------------------------------------------------!
              ! So, growth creates a demand for C, N, and P from the
              ! labile pools, depending on their concentrations in the
              ! wood. Take these to be fixed for now from:
              !
              !--------------------------------------------------------!
              C_labile_in = 0.0
              N_labile_in = 0.0
              P_labile_in = 0.0
              !--------------------------------------------------------!
              C_labile_out = 0.47 * dplant_mass
              N_labile_out = 0.02 * C_labile_out
              P_labile_out = 0.02 * N_labile_out
              !--------------------------------------------------------!
              dC_labile = C_labile_in - C_labile_out
              dN_labile = N_labile_in - N_labile_out
              dP_labile = P_labile_in - P_labile_out
              !--------------------------------------------------------!

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
              !--------------------------------------------------------!
              ! Simplest: assume WUE of 3 g[DM]/mm, and only root in
              ! top layer.
              !--------------------------------------------------------!
              dwater = MIN (h2osoi_liq (1,x,y), dplant_mass / 3.0)
              dplant_mass = 3.0 * dwater
              !--------------------------------------------------------!
              plant_mass (K,x,y) = plant_mass (K,x,y) + dplant_mass
              h2osoi_liq (1,x,y) = h2osoi_liq (1,x,y) - dwater
              !--------------------------------------------------------!
              ! Assume all is a cylinder of mass with length, and
              ! with some portion below-ground                      (mm)
              !--------------------------------------------------------!
              plant_length (K,x,y) = (400.0 * plant_mass (K,x,y) / &
                                     3.142E-3) ** (one / 3.0)
              !--------------------------------------------------------!
              ! Sum plant growth to give annual NPP           (g[DM]/yr)
              !--------------------------------------------------------!
              npp = npp + dplant_mass
              !--------------------------------------------------------!
            END DO ! K = 1, nplants (x,y)
            !----------------------------------------------------------!

!----------------------------------------------------------------------!
END SUBROUTINE GROW
!======================================================================!
