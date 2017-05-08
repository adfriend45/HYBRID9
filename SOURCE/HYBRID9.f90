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
! 5th May, 2017
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE MPI    ! Enable access to the Message Passing Interface library of
           ! parallel routines.
USE NETCDF ! Enable access to the library of netCDF routines.
USE SHARED ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For timing the simulation using 'CPU_TIME' on processor 0.
! Results go to screen/slurm-*.out in seconds (see later).
!----------------------------------------------------------------------!
REAL :: start(1025),finish(1025)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CHARACTER (LEN = 200) :: file_in_ts ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_ks ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_lm ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_ps ! Soil input filename.
CHARACTER (LEN = 4) :: ydate ! Character value of jyear for file names.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: ierr              ! MPI error code                        (-)
INTEGER :: ncid              ! netCDF ID                             (-)
INTEGER :: num_procs         ! No. available processors              (-)
INTEGER :: file_free         ! Flag for reading/writing files      (1/0)
INTEGER :: status (MPI_STATUS_SIZE) ! source, tag, and error codes.
INTEGER :: nx_block,ny_block ! x and y dimensions of chunks over Earth.
INTEGER :: i_block,i_block_s ! Chunk indices                         (-)
INTEGER :: x,y,x1,y1         ! Spatial array indices                 (-)
INTEGER :: iTIME             ! Timepoint              (day 1/1/1860 = 1)
INTEGER :: iDEC              ! Decade index              (1 = 1901-1910)
INTEGER :: iDEC_start        ! Index of first decade in simulation   (-)
INTEGER :: iDEC_end          ! Index of last decade in simulation    (-)
INTEGER :: jyear             ! Calendar year                        (CE)
INTEGER :: syr               ! First year in decade                 (CS)
INTEGER :: eyr               ! Last year in decade                  (CE)
INTEGER :: I,J,K             ! Generic loop indices                  (-)
INTEGER :: NISURF            ! No. timepoints in day                 (-)
INTEGER :: NS                ! Timepoints in day index               (-)
INTEGER :: nt                ! No. timesteps in year              (days)
INTEGER :: NYR               ! No. years in simulation           (years)
INTEGER :: nlayers           ! Local no. soil layers                 (-)
!----------------------------------------------------------------------!
! Value of time at beginning of each calendar year from 1860 to 2300.
!----------------------------------------------------------------------!
INTEGER :: time_BOY (2300-1860+1)        !            (day 1/1/1860 = 1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
LOGICAL :: INTERACTIVE       ! Interactive simulation?         (logical)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Longitude of beginning of interactive simulation        (degrees east)
!----------------------------------------------------------------------!
REAL :: lon_w
!----------------------------------------------------------------------!
! Longitude index of beginning of interactive simulation             (-)
!----------------------------------------------------------------------!
INTEGER :: lon_s_w
!----------------------------------------------------------------------!
! Latitude of beginning of interactive simulation        (degrees north)
!----------------------------------------------------------------------!
REAL :: lat_w
!----------------------------------------------------------------------!
! Longitude index of beginning of interactive simulation             (-)
!----------------------------------------------------------------------!
INTEGER :: lat_s_w
!----------------------------------------------------------------------!
! Coefficients to calculate lon_s_w and lat_s_w                      (-)
!----------------------------------------------------------------------!
REAL :: a_lon,b_lon,a_lat,b_lat
!----------------------------------------------------------------------!
! No. longitudes and latitudes in interactive simulation             (-)
!----------------------------------------------------------------------!
REAL :: lon_c_w,lat_c_w
!----------------------------------------------------------------------!
! Local value of iTIME                         (1 = day 1 of iDEC_start)
!----------------------------------------------------------------------!
INTEGER :: iT
!----------------------------------------------------------------------!
! Local value of jyear                        (1 = year 1 of iDEC_start)
!----------------------------------------------------------------------!
INTEGER :: iY
!----------------------------------------------------------------------!
! Global array of processor IDs for diagnostics                      (-)
!----------------------------------------------------------------------!
INTEGER :: chunk (NX,NY)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Climate variables.
!----------------------------------------------------------------------!
! huss:   Specific humidity at time of maximum temperature     (kgw/kga)
! pr:     Precipitation flux                                  (kg/m^2/s)
! ps:     Surface air pressure                                      (Pa)
! rhs:    Relative humidity                                          (-)
! rlds:   Surface downwelling longwave radiation flux in air     (W/m^2)
! rsds:   Surface downwelling shortwave radiation flux in air    (W/m^2)
! tasmax: Maximum surface air temperature                            (K)
! tasmin: Minimum surface air temperature                            (K)
! tas:    Surface air temperature                                    (K)
! wind:   Wind speed at 10 metres                                  (m/s)
!----------------------------------------------------------------------!
! Surface Air Temperature                                            (K)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: tas (:,:,:)
!----------------------------------------------------------------------!
! Specific humidity at time of maximum temperature               (kg/kg)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: huss (:,:,:)
!----------------------------------------------------------------------!
! Surface air pressure                                              (Pa)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: ps (:,:,:)
!----------------------------------------------------------------------!
! Precipitation flux                                          (kg/m^2/s)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: pr (:,:,:)
!----------------------------------------------------------------------!
! Relative humidity                                                  (-)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: rhs (:,:,:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Annual sums for diagnostics.
!----------------------------------------------------------------------!
REAL :: plant_mass_sum  ! Sum of plant mass over 1-yr            (g[DM])
REAL :: rnf_sum         ! Sum of runoff fluxes over 1-yr          (mm/s)
REAL :: evap_sum        ! Sum of evaporation fluxes over 1-yr     (mm/s)
!----------------------------------------------------------------------!
REAL :: tas_sum         ! Sum of daily tas over 1-yr                 (K)
REAL :: huss_sum        ! Sum of daily huss over 1-yr            (kg/kg)
REAL :: ps_sum          ! Sum of daily ps over 1-yr                 (Pa)
REAL :: pr_sum          ! Sum of daily pr over 1-yr           (kg/m^2/s)
REAL :: rhs_sum         ! Sum of daily rhs over 1-yr              (%age)
!----------------------------------------------------------------------!
REAL :: h2osoi_sum_total ! Sum of daily total soil water over 1-yr  (mm)
!----------------------------------------------------------------------!
! Sum of daily volumetric water in each soil layer over 1-yr (mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: theta_sum (:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Miscellaneous variables.
!----------------------------------------------------------------------!
REAL :: w_i ! Soil moisture contraint on plant growth                (-)
REAL :: qflx_prec_grnd_rain ! Rain precip. incident on ground     (mm/s)
REAL :: qflx_top_soil  ! Net water input into soil from top       (mm/s)
REAL :: qflx_in_soil   ! Surface input to soil                    (mm/s)
REAL :: qflx_in_h2osfc ! Surface input to h2osfc                  (mm/s)
REAL :: qflx_surf   ! Surface runoff                              (mm/s)
REAL :: hkdepth     ! Decay factor                                   (m)
REAL :: fff         ! Decay factor                                  (/m)
REAL :: wtfact      ! Max. saturation fraction for a grid cell       (-)
REAL :: fsat        ! Fractional area with water table at surface    (-)
REAL :: fcov        ! Fractional impermable area                     (-)
REAL :: frac_h2osfc ! Fraction of ground covered by surface water    (-)
REAL :: qflx_infl ! Infiltration                                  (mm/s)
REAL :: qinmax    ! Maximum infiltration capacity                 (mm/s)
REAL :: qflx_infl_excess ! Infiltration excess runoff -> h2osfc   (mm/s)
REAL :: qflx_evap      ! Local evaporation                        (mm/s)
REAL :: qflx_ev_h2osfc ! Evaporative flux from h2osfc             (mm/s)
REAL :: qflx_evap_grnd ! Ground surface evaporation rate          (mm/s)
REAL :: dt      ! Integration timestep                               (s)
REAL :: num     ! For trigiagonal solution                          (mm)
REAL :: den     ! For trigiagonal solution                          (mm)
REAL :: dzq     ! For trigiagonal solution                          (mm)
REAL :: BET     ! For trigiagonal solution                           (-)
REAL :: tempi   ! Temporary variable for calculating vol_eq          (-)
REAL :: temp0   ! Temporary variable for calculating vol_eq          (-)
REAL :: voleq1  ! Temporary variable for calculating vol_eq          (-)
REAL :: rous    ! Aquifer yield                                      (-)
REAL :: qcharge ! Aquifer recharge rate                           (mm/s)
REAL :: wh_zwt  ! Water head at the water table depth               (mm)
REAL :: ka      ! Hydraulic conductivity of the aquifer           (mm/s)
REAL :: wh      ! smpfz (jwt) - z (jwt)                             (mm)
REAL :: qcharge_tot ! To compute water table depth change in soil (mm/s)
REAL :: qcharge_layer ! qcharge in saturated layer                (mm/s)
REAL :: s_y      ! Specific yield                                    (-)
REAL :: zwtmm    ! Water table depth                                (mm)
REAL :: rsub_top ! Subsurface runoff - topographic control        (mm/s)
REAL :: rsub_top_max ! Max. rsub_top                              (mm/s)
REAL :: rsub_top_tot ! To remove water from water table in soil   (mm/s)
REAL :: rsub_top_layer ! To remove water from water table in soil (mm/s)
REAL :: xsi ! Excess soil water above saturation at layer I         (mm)
REAL :: xs1 ! Excess soil water above saturation at layer 1         (mm)
REAL :: xs  ! Water needed to bring soil moisture to theta_m        (mm)
REAL :: qflx_rsub_sat ! Soil saturation excess                    (mm/s)
REAL :: available_h2osoi_liq ! Available soil liq. water in a layer (mm)
!----------------------------------------------------------------------!
! Water available for evaporation over timestep                   (mm/s)
!----------------------------------------------------------------------!
REAL :: evap_max
!----------------------------------------------------------------------!
! Soil moisture constraint on evaporative flux                (fraction)
!----------------------------------------------------------------------!
REAL :: beta
!----------------------------------------------------------------------!
! Air density                                                (kg[a]/m^3)
!----------------------------------------------------------------------!
REAL :: rho
!----------------------------------------------------------------------!
! Ratio of air to water density                            (kg[a]/kg[w])
!----------------------------------------------------------------------!
REAL :: rho3
!----------------------------------------------------------------------!
! Virtual potential surface temperature                              (K)
!----------------------------------------------------------------------!
REAL :: tsv
!----------------------------------------------------------------------!
! Specific humidity of ground                              (kg[w]/kg[a])
!----------------------------------------------------------------------!
REAL :: qb
!----------------------------------------------------------------------!
! Saturation water vapour mixing ratio                     (kg[w]/kg[a])
!----------------------------------------------------------------------!
REAL :: QSAT
!----------------------------------------------------------------------!
! Conductance of the atmosphere                                   (mm/s)
!----------------------------------------------------------------------!
REAL :: cna
!----------------------------------------------------------------------!
! Precipitation flux                                              (mm/s)
!----------------------------------------------------------------------!
REAL :: prec
!----------------------------------------------------------------------!
! For soil water diagnostics.
!----------------------------------------------------------------------!
REAL :: w0,w1
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
!----------------------------------------------------------------------!
! Soil matric potential in aquifer layer                            (mm)
!----------------------------------------------------------------------!
REAL :: smp1
!----------------------------------------------------------------------!
! d(smp)/d(vol_liq) in aquifer layer                    (mm/(mm^3/mm^3))
!----------------------------------------------------------------------!
REAL :: dsmpdw1
!----------------------------------------------------------------------!
! Underground drainage from each soil layer and aquifer           (mm/s)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: rnff
!----------------------------------------------------------------------!
! Matric potentials of soil layers                                  (mm)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: smp
!----------------------------------------------------------------------!
! Equilibrium soil matric potentials of soil layers                 (mm)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: zq
!----------------------------------------------------------------------!
! Soil water inflow at top of layer                               (mm/s)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: qin
!----------------------------------------------------------------------!
! Soil water outflow at bottom of layer                           (mm/s)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: qout
!----------------------------------------------------------------------!
! Soil layer interface depths with reference to surface             (mm)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: zi
!----------------------------------------------------------------------!
! Soil layer thicknesses                                            (mm)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dz
!----------------------------------------------------------------------!
! Hydraulic conductivities at soil layer interfaces               (mm/s)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: hk
!----------------------------------------------------------------------!
! d(smp)/d(vol_liq)                                     (mm/(mm^3/mm^3))
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dsmpdw
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I))                              (mm/s1)/(mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqodw1
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I+1))                             (mm/s)/(mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqodw2
!----------------------------------------------------------------------!
! d(hk)/d(vol_liq)                                    (mm/s)/(mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dhkdw
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I-1))                            ((mm/s)/(mm^3/mm^3))
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqidw0
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I))                              ((mm/s)/(mm^3/mm^3))
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqidw1
!----------------------------------------------------------------------!
! "a" left of diagonal term of tridiagonal matrix.
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: amx
!----------------------------------------------------------------------!
! "b" diagonal column for tridiagonal matrix.
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: bmx
!----------------------------------------------------------------------!
! "c" right of diagonal of tridiagonal matrix.
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: cmx
!----------------------------------------------------------------------!
! "r" forcing term of tridiagonal matrix.
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: rmx
!----------------------------------------------------------------------!
! Change in soil water                                         (m^3/m^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dwat,dwat2
!----------------------------------------------------------------------!
! For tridiagonal solution                                           (-)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: GAM
!----------------------------------------------------------------------!
! Volumetric soil water content                              (mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: theta
!----------------------------------------------------------------------!
! Equilibrium volumetric soil water content                  (mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: vol_eq
!----------------------------------------------------------------------!
! Porosity of soil                                           (mm^3/mm^3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: eff_porosity
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Miscellaneous parameters.
!----------------------------------------------------------------------!
! Truncation limit for numerical tests.
!----------------------------------------------------------------------!
REAL, PARAMETER :: trunc = 1.0E-8
!----------------------------------------------------------------------!
! Minimum soil moisture                                             (mm)
!----------------------------------------------------------------------!
REAL, PARAMETER :: watmin = 0.01
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Index of soil layer right above water table                        (-)
!----------------------------------------------------------------------!
INTEGER :: jwt
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! End of declarations.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise MPI.
!----------------------------------------------------------------------!
CALL MPI_INIT (ierr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Get this processor's ID (my_id).
!----------------------------------------------------------------------!
CALL MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Find out how many processors are available (num_procs).
!----------------------------------------------------------------------!
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Some diagnostics.
!----------------------------------------------------------------------!
WRITE (*,*) 'num_procs, my_id ',num_procs,my_id
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Do some diagnostics for establishing run time on processors.
!----------------------------------------------------------------------!
CALL CPU_TIME (start (my_id+1))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set longitudes and latitudes in centres of each grid-box      (degree)
!----------------------------------------------------------------------!
DO x = 1, NX
  lon_all (x) = -179.75 + (x - 1) * 0.5
END DO
DO y = 1, NY
  lat_all (y) = 89.75 - (y - 1) * 0.5
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Volumetric water in each layer, for diagnostics            (mm^3/mm^3)
!----------------------------------------------------------------------!
ALLOCATE (theta_sum (nsoil_layers_max))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate soil arrays over layers.
!----------------------------------------------------------------------!
ALLOCATE (qin    (nsoil_layers_max+1)) ! Inflow at top            (mm/s)
ALLOCATE (qout   (nsoil_layers_max+1)) ! Outflow at bottom        (mm/s)
ALLOCATE (dsmpdw (nsoil_layers_max+1)) ! d(smp)/d(w)    (mm/(mm^3/mm^3))
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I))                              (mm/s1)/(mm^3/mm^3)
!----------------------------------------------------------------------!
ALLOCATE (dqodw1 (nsoil_layers_max+1))
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I+1)) (mm s-1) /(mm^3 mm^-3)
!----------------------------------------------------------------------!
ALLOCATE (dqodw2 (nsoil_layers_max+1))
!----------------------------------------------------------------------!
! d(hk)/d(vol_liq) (mm s-1) / (mm^3 mm^-3)
!----------------------------------------------------------------------!
ALLOCATE (dhkdw  (nsoil_layers_max))
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I-1)) ((mm s-1) / (mm^3 mm^-3)).
!----------------------------------------------------------------------!
ALLOCATE (dqidw0 (nsoil_layers_max+1))
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I)) ((mm s-1) / (mm^3 mm^-3)).
!----------------------------------------------------------------------!
ALLOCATE (dqidw1 (nsoil_layers_max+1))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Terms for tridiagonal matrix.
!----------------------------------------------------------------------!
ALLOCATE (amx    (nsoil_layers_max+1))
ALLOCATE (bmx    (nsoil_layers_max+1))
ALLOCATE (cmx    (nsoil_layers_max+1))
ALLOCATE (rmx    (nsoil_layers_max+1))
!----------------------------------------------------------------------!
ALLOCATE (dwat2  (nsoil_layers_max+1))
ALLOCATE (dwat   (nsoil_layers_max))
ALLOCATE (GAM    (nsoil_layers_max+1))
!----------------------------------------------------------------------!
ALLOCATE (zi     (0:Nlevgrnd))         ! Layer interface depths     (mm)
ALLOCATE (dz     (1:Nlevgrnd))         ! Layer thicknesses          (mm)
ALLOCATE (zc_o   (1:nsoil_layers_max)) ! Soil layer nodes for diag  (mm)
ALLOCATE (zc     (1:Nlevgrnd))         ! Ground layer nodes         (mm)
ALLOCATE (smp    (nsoil_layers_max))   ! Matric potentials          (mm)
ALLOCATE (zq     (nsoil_layers_max+1)) ! Eqm. water potentials      (mm)
!----------------------------------------------------------------------!
! Volumetric soil water content                              (mm^3/mm^3)
!----------------------------------------------------------------------!
ALLOCATE (theta (nsoil_layers_max))
!----------------------------------------------------------------------!
! Equilibrium volumetric soil water (mm^3/mm^3).
!----------------------------------------------------------------------!
ALLOCATE (vol_eq (nsoil_layers_max+1))
!----------------------------------------------------------------------!
! Porosity of soil (mm^3/mm^3).
!----------------------------------------------------------------------!
ALLOCATE (eff_porosity (nsoil_layers_max))
!----------------------------------------------------------------------!
! Hydraulic conductivities                                     (mm s-1).
!----------------------------------------------------------------------!
ALLOCATE (hk  (nsoil_layers_max))
!----------------------------------------------------------------------!
! Underground drainage from each soil layer and aquifer           (mm/s)
!----------------------------------------------------------------------!
ALLOCATE (rnff (nsoil_layers_max+1))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read parameters from driver.txt to control simulation.
! Method for parallel reading based on
! http://rsdavis.mycpanel.princeton.edu/wp/?p=91 method 3.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Master gets permission to read parameter file first, others wait in
! line.
!----------------------------------------------------------------------!
file_free = 0
If (my_id == 0) THEN
  file_free = 1
ELSE
  CALL MPI_RECV (file_free, 1, MPI_INT, my_id-1, 1, MPI_COMM_WORLD, &
                 status, ierr)
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! This process reads the file.
!----------------------------------------------------------------------!
IF (file_free == 1) THEN
  !--------------------------------------------------------------------!
  OPEN (10,FILE='driver.txt',STATUS='OLD')
  READ (10,*) NISURF
  READ (10,*) iDEC_start
  READ (10,*) iDEC_end
  READ (10,*) INTERACTIVE
  READ (10,*) lon_w
  READ (10,*) lat_w
  READ (10,*) lon_c_w
  READ (10,*) lat_c_w
  !--------------------------------------------------------------------!
  ! Read ground layer boundaries, negative as go down from surface,
  ! which is boundary 0 at 0 mm (mm).
  !--------------------------------------------------------------------!
  DO I = 0, Nlevgrnd
    READ (10,*) zi (I)
  END DO
  !--------------------------------------------------------------------!
  CLOSE (10)
!----------------------------------------------------------------------!
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set integration timestep                                           (s)
!----------------------------------------------------------------------!
dt = 86400.0 / FLOAT (NISURF)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For running for particular locations interactively.
!----------------------------------------------------------------------!
IF (INTERACTIVE) THEN
  !--------------------------------------------------------------------!
  ! Indices for a focus point, if used.
  !--------------------------------------------------------------------!
  b_lon = (720.0 - 1.0) / (179.75 - (-179.75))
  a_lon = 1.0 - b_lon * (-179.75)
  b_lat = (360.0 - 1.0) / (-89.75 - 89.75)
  a_lat = 360.0 - b_lat * (-89.75)
  WRITE (*,*) 'Details of focus point'
  write (*,*) a_lon + b_lon * lon_w
  write (*,*) NINT (a_lon + b_lon * lon_w)
  lon_s_w = NINT (a_lon + b_lon * lon_w)
  write (*,*) a_lat + b_lat * lat_w
  write (*,*) NINT (a_lat + b_lat * lat_w)
  lat_s_w = NINT (a_lat + b_lat * lat_w)
  write (*,*) lon_s_w,lat_s_w
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Give read file permission to next process.
!----------------------------------------------------------------------!
IF (my_id /= num_procs-1) CALL MPI_Send (file_free, 1, MPI_INT, &
                                         my_id+1, 1, MPI_COMM_WORLD, &
                                         ierr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Compute ground layer thicknesses and centres (mm).
! Nodes are positive values below surface.
! Layer indices start with 1 at surface.
!----------------------------------------------------------------------!
DO I = 1, Nlevgrnd
  dz (I) = zi (I) - zi (I-1)
END DO
DO I = 1, Nlevgrnd
  zc (I) = zi (I) - dz (I) / 2.0
END DO
!----------------------------------------------------------------------!
! Set soil layer node depths for diagnostics (mm).
!----------------------------------------------------------------------!
DO I = 1, nsoil_layers_max
  zc_o (I) = zc (I)
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Divide Earth's surface into chunks of same size based on number of
! available processors. No. processes must be factoriable into an
! n x n grid, where n is the sqrt of the number of processors.
!----------------------------------------------------------------------!
nx_block = NINT (SQRT (FLOAT (num_procs)))
ny_block = NINT (SQRT (FLOAT (num_procs)))
lon_c = NX / nx_block
lat_c = NY / ny_block
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For running for particular locations interactively.
!----------------------------------------------------------------------!
IF (INTERACTIVE) THEN
  lon_c = lon_c_w
  lat_c = lat_c_w
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set no. years in simulation (yr).
!----------------------------------------------------------------------!
IF (iDEC_end < 12) THEN
  NYR = (iDEC_end - iDEC_start + 1) * 10
ELSE
  NYR = (iDEC_end - iDEC_start + 1 - 1) * 10 + 2
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Allocate array chunks.
!----------------------------------------------------------------------!
! Plant mass       
!----------------------------------------------------------------------!
ALLOCATE (plant_mass (nplants_max,lon_c,lat_c)) ! Plant mass     (g[DM])
ALLOCATE (nplants (lon_c,lat_c)) ! Number plants per grid-box        (-)
!----------------------------------------------------------------------!
ALLOCATE (data_in_2DI (lon_c,lat_c))  ! Generic 2D input integer array.
ALLOCATE (soil_tex    (lon_c,lat_c))  ! HWSD soil textures           (-)
ALLOCATE (topo_slope  (lon_c,lat_c))  ! Gridcell slope             (deg)
ALLOCATE (Fmax        (lon_c,lat_c))  ! Max. sat. fraction    (fraction)
ALLOCATE (lon         (lon_c))        ! Longitudes (degrees east)
ALLOCATE (lat         (lat_c))        ! Latitudes (degrees north)
ALLOCATE (block_sub   (lon_c,lat_c))  ! Processor IDs.
ALLOCATE (axy_plant_mass (lon_c,lat_c,NYR)) ! Mean ann. pl. mass (g[DM])
ALLOCATE (axy_tas     (lon_c,lat_c,NYR)) ! Mean annual tas           (K)
ALLOCATE (axy_huss    (lon_c,lat_c,NYR)) ! Mean annual huss      (kg/kg)
ALLOCATE (axy_ps      (lon_c,lat_c,NYR)) ! Mean annual ps           (Pa)
ALLOCATE (axy_pr      (lon_c,lat_c,NYR)) ! Mean annual pr     (kg/m^2/s)
ALLOCATE (axy_rhs     (lon_c,lat_c,NYR)) ! Mean annual rhs        (%age)
!----------------------------------------------------------------------!
! Mean annual volumetric water content by soil layer            (m3/m^3)
!----------------------------------------------------------------------!
ALLOCATE (axy_theta   (nsoil_layers_max,lon_c,lat_c,NYR))
!----------------------------------------------------------------------!
ALLOCATE (axy_theta_total (lon_c,lat_c,NYR))
!----------------------------------------------------------------------!
ALLOCATE (axy_rnf     (lon_c,lat_c,NYR)) ! Mean annual run-off    (mm/s)
ALLOCATE (axy_evap    (lon_c,lat_c,NYR)) ! Mean annual evap.      (mm/s)
!----------------------------------------------------------------------!
! Liquid soil water mass                                        (kg/m^2)
!----------------------------------------------------------------------!
ALLOCATE (h2osoi_liq (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Water table depth                                                  (m)
!----------------------------------------------------------------------!
ALLOCATE (zwt (lon_c,lat_c))
!----------------------------------------------------------------------!
! Water in the unconfined aquifer (mm).
!----------------------------------------------------------------------!
ALLOCATE (wa (lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunks of soil properties of one soil layer read in at
! 30 arc-seconds.
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1_in(lon_c*60,lat_c*60)) ! 0.001xcm^3cm^3.
ALLOCATE (k_s_l1_in    (lon_c*60,lat_c*60)) ! cm/day.
ALLOCATE (lambda_l1_in (lon_c*60,lat_c*60)) ! 0.001xunitless.
ALLOCATE (psi_s_l1_in  (lon_c*60,lat_c*60)) ! cm.
!----------------------------------------------------------------------!
! Chunks of soil properties of one soil layer gridded to half-degree.
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1   (lon_c,lat_c)) ! 0.001xcm3 cm-3.
ALLOCATE (k_s_l1       (lon_c,lat_c)) ! cm day-1.
ALLOCATE (lambda_l1    (lon_c,lat_c)) ! 0.001xunitless.
ALLOCATE (psi_s_l1     (lon_c,lat_c))                             ! (cm)
!----------------------------------------------------------------------!
! Chunk of saturated volumetric soil water                   (cm^3/cm^3)
!----------------------------------------------------------------------!
ALLOCATE (theta_s (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of saturated soil hydraulic conductivity                  (mm/s)
!----------------------------------------------------------------------!
ALLOCATE (hksat (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of soil pore size distribution index                         (-)
!----------------------------------------------------------------------!
ALLOCATE (lambda (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of Clapp and Hornberger "b"                                  (-)
!----------------------------------------------------------------------!
ALLOCATE (bsw (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of soil saturated capillary/water potential                 (mm)
!----------------------------------------------------------------------!
ALLOCATE (psi_s  (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of minimum volumetric soil matric potential           (cm3/cm^3)
!----------------------------------------------------------------------!
ALLOCATE (theta_m (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise global diagnostic arrays with fill (i.e. NaN and zero)
! values.
!----------------------------------------------------------------------!
axy_plant_mass (:,:,:) = zero / zero
axy_rnf   (:,:,:) = zero / zero
axy_evap  (:,:,:) = zero / zero
axy_tas   (:,:,:) = zero / zero
axy_huss  (:,:,:) = zero / zero
axy_ps    (:,:,:) = zero / zero
axy_pr    (:,:,:) = zero / zero
axy_rhs   (:,:,:) = zero / zero
axy_theta (:,:,:,:) = zero / zero
axy_theta_total (:,:,:) = zero
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IF (my_id == 0) THEN
  WRITE (*,*) 'num_procs nx_block ny_block',num_procs,nx_block,ny_block
  WRITE (*,*) 'num_procs lon_c lat_c',num_procs,lon_c,lat_c
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Assign processor IDs across the entire 'chunk' array.
!----------------------------------------------------------------------!
i_block_s = 0
DO y = 1, NY
  ! First chunk from IDL.
  i_block = i_block_s
  DO x = 1, NX
    ! Assign processor ID to chunk.
    chunk (x,y) = i_block
    ! If reached new chunk along lon increment processor ID.
    IF (MOD (x,lon_c) == 0) THEN
      i_block = i_block + 1
    END IF
  END DO
  ! If reached new chunk along lat increment first chunk from IDL by
  ! no. of chunks in lon dimension.
  IF (MOD (y,lat_c) == 0) THEN
    i_block_s = i_block_s + ny_block
  END IF
END DO
!----------------------------------------------------------------------!
! Find lon_s and lat_s for this processor by working backwards through
! the global array of processor IDs.
!----------------------------------------------------------------------!
DO y = NY, 1, -1
  DO x = NX, 1, -1
    IF (my_id == chunk (x,y)) THEN
      lon_s = x
      lat_s = y
    END IF
  END DO
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For running for particular locations interactively.
!----------------------------------------------------------------------!
IF (INTERACTIVE) THEN
  lon_s = lon_s_w
  lat_s = lat_s_w
  WRITE (*,*) 'my_id lon_s lat_s',my_id,lon_s,lat_s
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read soil textures for this block.
! Currently only used to determine soil grid points.
!----------------------------------------------------------------------!
file_name = '/scratch/adf10/ISIMIP2/INPUT/SOILS/hwsd.final.hlf.nc4'
varid = 3
!----------------------------------------------------------------------!
WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
CALL READ_NET_CDF_2DI
!----------------------------------------------------------------------!
soil_tex (:,:) = data_in_2DI (:,:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read in soil properties for this block.
! Saturated water contents          , theta_s, are 0.001*cm^3/cm^3.
! Saturated hydraulic conductivities, k_s > hksat, are    cm/day.
! Pore size distribution indices    , lambda , are 0.001*unitless.
! Saturated capillary potentials    , psi_s  , are cm.
! Downloaded to /home/adf10/DATA/SOILS
! from http://globalchange.bnu.edu.cn/research/soil4d.jsp#download
! on 22/7/16.
!----------------------------------------------------------------------!
DO I = 1, nsoil_layers_max ! Loop over soil layers
  !--------------------------------------------------------------------!
  SELECT CASE (I)
    CASE (1)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l1.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l1.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l1.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l1.nc'
    CASE (2)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l2.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l2.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l2.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l2.nc'
    CASE (3)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l3.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l3.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l3.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l3.nc'
    CASE (4)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l4.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l4.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l4.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l4.nc'
    CASE (5)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l5.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l5.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l5.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l5.nc'
    CASE (6)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l6.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l6.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l6.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l6.nc'
    CASE (7)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l7.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l7.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l7.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l7.nc'
    CASE (8)
      file_in_ts = '/scratch/adf10/DATA/SOILS/theta_s/theta_s_l8.nc'
      file_in_ks = '/scratch/adf10/DATA/SOILS/k_s/k_s_l8.nc'
      file_in_lm = '/scratch/adf10/DATA/SOILS/lambda/lambda_l8.nc'
      file_in_ps = '/scratch/adf10/DATA/SOILS/psi/psi_s_l8.nc'
  END SELECT
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',TRIM(file_in_ts)
  CALL CHECK(NF90_OPEN (TRIM (file_in_ts), NF90_NOWRITE, ncid))
  WRITE (*,*) 'File ',TRIM(file_in_ts),' opened'
  CALL CHECK(NF90_GET_VAR (ncid, 4, theta_s_l1_in, &
                       start = (/(lon_s-1)*60+1, (lat_s-1)*60+1 /), &
                       count = (/lon_c*60, lat_c*60 /)))
  WRITE (*,*) 'Saturated water contents read for layer',I
  CALL CHECK(NF90_CLOSE (ncid))
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',TRIM(file_in_ks)
  CALL CHECK(NF90_OPEN (TRIM (file_in_ks), NF90_NOWRITE, ncid))
  WRITE (*,*) 'File ',TRIM(file_in_ks),' opened'
  CALL CHECK(NF90_GET_VAR (ncid, 4, k_s_l1_in, &
                       start = (/(lon_s-1)*60+1, (lat_s-1)*60+1 /), &
                       count = (/lon_c*60, lat_c*60 /)))
  WRITE (*,*) 'Saturated hydraulic conductivity read for layer',I
  CALL CHECK(NF90_CLOSE (ncid))
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',TRIM(file_in_lm)
  CALL CHECK(NF90_OPEN (TRIM (file_in_lm), NF90_NOWRITE, ncid))
  WRITE (*,*) 'File ',TRIM(file_in_lm),' opened'
  CALL CHECK(NF90_GET_VAR (ncid, 4, lambda_l1_in, &
                       start = (/(lon_s-1)*60+1, (lat_s-1)*60+1 /), &
                       count = (/lon_c*60, lat_c*60 /)))
  WRITE (*,*) 'Soil pore size distribution read for layer',I
  CALL CHECK(NF90_CLOSE (ncid))
  !--------------------------------------------------------------------!
  WRITE (*,*) 'Opening ',TRIM(file_in_ps)
  CALL CHECK(NF90_OPEN (TRIM (file_in_ps), NF90_NOWRITE, ncid))
  WRITE (*,*) 'File ',TRIM(file_in_ps),' opened'
  CALL CHECK(NF90_GET_VAR (ncid, 4, psi_s_l1_in, &
                       start = (/(lon_s-1)*60+1, (lat_s-1)*60+1 /), &
                       count = (/lon_c*60, lat_c*60 /)))
  WRITE (*,*) 'Soil saturated water potential read for layer',I
  CALL CHECK(NF90_CLOSE (ncid))
  !--------------------------------------------------------------------!
  ! Data is at 30 arc-seconds, so grid to half-degree.
  !--------------------------------------------------------------------!
  theta_s_l1 (:,:) = zero
  k_s_l1     (:,:) = zero
  lambda_l1  (:,:) = zero
  psi_s_l1   (:,:) = zero
  DO y = 1, lat_c
    DO x = 1, lon_c
      j = 0
      DO x1 = (x-1)*60+1, (x-1)*60+60
        do y1 = (y-1)*60+1,(y-1)*60+60
          IF (theta_s_l1_in (x1,y1) >= zero) THEN
            theta_s_l1 (x,y) =  theta_s_l1 (x,y) + theta_s_l1_in (x1,y1)
            k_s_l1     (x,y) =  k_s_l1     (x,y) + k_s_l1_in     (x1,y1)
            lambda_l1  (x,y) =  lambda_l1  (x,y) + lambda_l1_in  (x1,y1)
            psi_s_l1   (x,y) =  psi_s_l1   (x,y) + psi_s_l1_in   (x1,y1)
            j = j + 1
          END IF
        END DO
      END DO
      IF (j > 0) THEN
        theta_s_l1 (x,y) = theta_s_l1 (x,y) / FLOAT (j)
        k_s_l1     (x,y) = k_s_l1     (x,y) / FLOAT (j)
        lambda_l1  (x,y) = lambda_l1  (x,y) / FLOAT (j)
        psi_s_l1   (x,y) = psi_s_l1   (x,y) / FLOAT (j)
      END IF
    END DO
  END DO
  !--------------------------------------------------------------------!
  ! Place saturated soil properties in block arrays for each layer.
  ! Volumetric water contents read in as cm^3/cm^3 x 0.001,
  ! converted to mm.
  ! Saturated hydraulic conductivities read in as cm/day,
  ! converted to mm/s.
  ! Pore size distribution index read in as 0.001 x unitless,
  ! converted to unitless.
  ! Saturated water potential read in as cm, converted to mm.
  !--------------------------------------------------------------------!
  DO y = 1, lat_c
    DO x = 1, lon_c
      theta_s (I,x,y) = theta_s_l1(x,y) / 1.0E3
      hksat   (I,x,y) = 10.0 * k_s_l1 (x,y) / 86400.0
      lambda  (I,x,y) = lambda_l1 (x,y) / 1.0E3
      psi_s   (I,x,y) = 10.0 * psi_s_l1 (x,y)
      !----------------------------------------------------------------!
      ! To stop potential dbz.
      !----------------------------------------------------------------!
      lambda  (I,x,y) = MAX (lambda (I,x,y), trunc)
      !----------------------------------------------------------------!
      ! Clapp and Hornberger "b"                                     (-)
      !----------------------------------------------------------------!
      bsw (I,x,y) = 1.0 / lambda (I,x,y)
      !----------------------------------------------------------------!
    END DO
  END DO
  !--------------------------------------------------------------------!
END DO ! I = 1, nsoil_layers_max ! Loop over soil layers.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Free up resources.
!----------------------------------------------------------------------!
DEALLOCATE (theta_s_l1_in)
DEALLOCATE (k_s_l1_in)
DEALLOCATE (lambda_l1_in)
DEALLOCATE (psi_s_l1_in)
DEALLOCATE (theta_s_l1)
DEALLOCATE (k_s_l1)
DEALLOCATE (lambda_l1)
DEALLOCATE (psi_s_l1)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read maximum soil saturated fractions for this block (fraction).
!----------------------------------------------------------------------!
file_name = 'Fmax_half.nc'
varid = 1
!----------------------------------------------------------------------!
WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
!----------------------------------------------------------------------!
! Read Fmax values into data_in_2DI.
!----------------------------------------------------------------------!
CALL READ_NET_CDF_2DI
!----------------------------------------------------------------------!
Fmax (:,:) = zero / zero
!----------------------------------------------------------------------!
DO y = 1, lat_c
  DO x = 1, lon_c
    IF ((soil_tex (x,y) > 0) .AND. (soil_tex (x,y) /= 13) .AND. &
      SUM (theta_s (:,x,y)) > trunc) THEN
      !----------------------------------------------------------------!
      ! Set Fmax for soiled points.
      !----------------------------------------------------------------!
      ! Set any missing Fmax to global arithmetic half-degree mean.
      !----------------------------------------------------------------!
      IF (data_in_2DI (x,y) == -9999) THEN
        data_in_2DI (x,y) = 3809
      END IF
      !----------------------------------------------------------------!
      Fmax (x,y) = FLOAT (data_in_2DI (x,y)) / 10000.0
      !----------------------------------------------------------------!
    END IF
  END DO
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write out saturated water potentials for diagnostics (mm).
!----------------------------------------------------------------------!
file_name = 'psi_s.nc'
WRITE (*,*) 'Writing to ',TRIM(file_name)
CALL WRITE_NET_CDF_3DR_soils
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Minimum soil volumetric water contents                     (mm^3/mm^3)
! Following assumes theta_m is hygroscopic water at -31 bar.
! https://en.wikipedia.org/wiki/Water_content defines the remaining
! water at high tension as 'Residual water content'. Assume this is the
! same as theta_m. Assume is hygroscopic water at -31 bar as in
! http://www.napavalley.edu/people/ajohnson/Documents/VWT%20132%
! 20Spring%202015/13%20-%20Basic%20Properties%20of%20Soil%20Water%20
! (171%20-%20203).pdf
! They use:
! http://lawr.ucdavis.edu/classes/ssc107/SSC107Syllabus/chapter2-00.pdf
! to compute -31 bar in mm tension, and then invert equation for psi to
! give theta_m at -31 bar.
! Not sure if this is the amount that would stop flow - too high?
! Also set water table depth and aquifer water as in CESM.
!----------------------------------------------------------------------!
theta_m    (:,:,:) = zero
h2osoi_liq (:,:,:) = zero
plant_mass (:,:,:) = zero
zwt (:,:) = zero
wa  (:,:) = zero
nlayers = nsoil_layers_max
DO y = 1, lat_c
  DO x = 1, lon_c
    IF ((soil_tex (x,y) > 0) .AND. (soil_tex (x,y) /= 13) .AND. &
      SUM (theta_s (:,x,y)) > trunc) THEN
      DO I = 1, nsoil_layers_max
        !--------------------------------------------------------------!
        ! In CESM following is watmin = 0.01 mm?
        !--------------------------------------------------------------!
        theta_m (I,x,y) = theta_s (I,x,y) * ((-3.1E9 / (1000.0 * 9.8)) &
                          / psi_s (I,x,y)) ** (-lambda (I,x,y))
        !--------------------------------------------------------------!
        ! Initial liquid water in each soil layer (kg/m^2).
        !--------------------------------------------------------------!
        h2osoi_liq (I,x,y) = 0.5 * theta_s (I,x,y) * dz (I) * rhow &
                             / 1000.0
        !--------------------------------------------------------------!
      END DO
      !----------------------------------------------------------------!
      ! Initial water table depth from p. 28 of O13                  (m)
      !----------------------------------------------------------------!
      zwt (x,y) = (zi (nlayers) + 5000.0) / 1000.0
      !----------------------------------------------------------------!
      ! Initial water stored in the unconfined aquifer and unsaturated
      ! soil from p. 28 of O13 (mm).
      !----------------------------------------------------------------!
      wa (x,y) = 4000.0
      !----------------------------------------------------------------!
      ! Initial number of plants in grid-box                         (-)
      !----------------------------------------------------------------!
      nplants (x,y) = 1
      !----------------------------------------------------------------!
      ! Initial plant masses                                     (g[DM])
      !----------------------------------------------------------------!
      DO K = 1, nplants (x,y)
        plant_mass (K,x,y) = 1.0
      END DO
      !----------------------------------------------------------------!
    END IF
  END DO
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Put processor ID into block_sub for diagnostics.
!----------------------------------------------------------------------!
block_sub (:,:) = my_id
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write chunk numbers for diagnostics.
!----------------------------------------------------------------------!
file_name = 'chunk.nc'
var_name = 'chunk'
!----------------------------------------------------------------------!
CALL WRITE_NET_CDF_2DI (block_sub)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write HWSD soil textures for diagnostics.
!----------------------------------------------------------------------!
file_name = 'soil_tex.nc'
var_name = 'soil_tex'
!----------------------------------------------------------------------!
CALL WRITE_NET_CDF_2DI (soil_tex)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! time is days since 1860-01-01, so first time is 0. Is actually 1.
! Climate input 1901 starts Time = 14975.
!----------------------------------------------------------------------!
! time at beginning of each calendar year starting in 1860.
!----------------------------------------------------------------------!
time_BOY (1) = 1 ! By testing.
DO jyear = 1861, 2300
  IF (MOD (jyear-1,4) .NE. 0) THEN
    time_BOY (jyear-1859) = time_BOY (jyear-1859-1) + 365
  ELSE
    IF (MOD (jyear-1, 100) .NE. 0) THEN
      time_BOY (jyear-1859) = time_BOY (jyear-1859-1) + 366
    ELSE
      IF (MOD (jyear-1, 400) .NE. 0) THEN
        time_BOY (jyear-1859) = time_BOY (jyear-1859-1) + 365
      ELSE
        time_BOY (jyear-1859) = time_BOY (jyear-1859-1) + 366
      END IF
    END IF
  END IF
END DO
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Get the latitudes of the grid-boxes in this block.
!----------------------------------------------------------------------!
file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
            &'tas_pgfv2.1_1901_1910.nc4'
CALL READ_NET_CDF_1DR (lat,2,lat_s,lat_c) ! Get latitudes.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Get the longitudes of the grid-boxes in this block.
!----------------------------------------------------------------------!
file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
            &'tas_pgfv2.1_1901_1910.nc4'
CALL READ_NET_CDF_1DR (lon,1,lon_s,lon_c) ! Get latitudes.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IF (INTERACTIVE) THEN
  WRITE (*,*) 'longitudes'
  WRITE (*,*) lon
  WRITE (*,*) 'latitudes'
  WRITE (*,*) lat
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Outer main loop: over decades (iDEC 1 = 1901-1910).
!----------------------------------------------------------------------!
DO iDEC = iDEC_start, iDEC_end
!----------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Get the number of timesteps in the tas file for this decade.
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'tas_pgfv2.1_'//decade(iDEC)//'.nc4'
  CALL READ_NET_CDF_0D (NTIMES) ! Get NTIMES.
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Surface Air Temperature (K) for this decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (tas(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (tas)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Specific Humidity at time of Maximum Temperature (kg kg-1)
  ! for this decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'huss_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (huss(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (huss)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Surface pressure (Pa) for this decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'ps_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (ps(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (ps)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Precipitation Flux (kg m-2 s-1) for this decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'pr_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (pr(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (pr)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Relative Humidity (no units) for this decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'rhs_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (rhs(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (rhs)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Loop over days.
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! syr is jyear at beginning of decade iDEC.
  !--------------------------------------------------------------------!
  syr = (iDEC - 1) * 10 + 1901
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! eyr is jyear at end of decade iDEC.
  !--------------------------------------------------------------------!
  IF (iDEC < 12) THEN
    eyr = syr + 9
  ELSE
    eyr = syr + 1
  END IF
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Inner main loops: over timepoints in decade at each grid point
  ! with soil in both HWSD texture and theta_s fields.
  !--------------------------------------------------------------------!
  DO y = 1, lat_c
    DO x = 1, lon_c
      IF ((soil_tex (x,y) > 0) .AND. (soil_tex (x,y) /= 13) .AND. &
        SUM (theta_s (:,x,y)) > trunc) THEN
        !--------------------------------------------------------------!
        nlayers = nsoil_layers_max
        !--------------------------------------------------------------!
        ! Loop over years in decade.
        !--------------------------------------------------------------!
        DO jyear = syr, eyr
          !------------------------------------------------------------!
          ! Intialise annual diagnostics.
          !------------------------------------------------------------!
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
            iT = iTIME-time_BOY(syr-1859)+1
            !----------------------------------------------------------!

            !CALL HYDROLOGY

            !----------------------------------------------------------!
            ! Precipitation flux (mm s-1).
            !----------------------------------------------------------!
            prec = 1.0E3 * pr (x,y,iT) / rhow
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Loop over timesteps within day (n).
            !----------------------------------------------------------!
            DO NS = 1, NISURF

            !----------------------------------------------------------!
            ! Initial total soil water for INTERACTIVE diagnostics (mm).
            ! Also compute volumetric water (mm^3/mm^3).
            !----------------------------------------------------------!
            w0 = prec * dt + wa (x,y)
            DO I = 1, nlayers
              !--------------------------------------------------------!
              w0 = w0 + h2osoi_liq (I,x,y)
              !--------------------------------------------------------!
              theta (I) = h2osoi_liq (I,x,y) / (dz (I) * rhow / 1000.0)
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

!**********************************************************************!
! Start of CESM 'SLakeHydrology'.
!**********************************************************************!

            !----------------------------------------------------------!
            ! Precipitation onto ground, assuming constant during
            ! day. Need to test if pulses changes results         (mm/s)
            !----------------------------------------------------------!
            qflx_prec_grnd_rain = prec
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Net water input into soil from top                  (mm/s)
            !----------------------------------------------------------!
            qflx_top_soil = qflx_prec_grnd_rain
            !----------------------------------------------------------!

!**********************************************************************!
! Start of CESM 'SurfaceRunoff'.
!**********************************************************************!

            !----------------------------------------------------------!
            ! Decay factor                                           (m)
            ! In CESM set using line in 'iniTimeConst.F90'.
            !----------------------------------------------------------!
            hkdepth = 1.0 / 2.5
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Decay factor                                          (/m)
            !----------------------------------------------------------!
            fff = 1.0 / hkdepth
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Max. saturation fraction for grid cell                 (-)
            !----------------------------------------------------------!
            wtfact = Fmax (x,y)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Fractional area with water table at surface            (-)
            !----------------------------------------------------------!
            fsat = wtfact * EXP (-0.5 * fff * zwt (x,y))
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Fractional impermable area                             (-)
            !----------------------------------------------------------!
            fcov = fsat
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Surface runoff                                      (mm/s)
            !----------------------------------------------------------!
            qflx_surf = fcov * qflx_top_soil
            !----------------------------------------------------------!

!**********************************************************************!
! End of CESM 'SurfaceRunoff'.
!**********************************************************************!

            !----------------------------------------------------------!
            ! Fraction of ground covered by surface water
            ! (full treatment is in H2OSfcMod.F90: do later) (-).
            ! Set to zero for now.
            !----------------------------------------------------------!
            !frac_h2osfc = fsat
            frac_h2osfc = zero !***************!
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Soil moisture constraint on evaporative flux (0-1).
            !----------------------------------------------------------!
            IF (theta_s (1,x,y) > trunc) THEN
              beta = MIN (one, theta (1) / theta_s (1,x,y))
            ELSE
              beta = 0.0
            END IF
            !----------------------------------------------------------!
            ! Virtual potential surface temperature                 (K).
            ! Correct temperature to use for rho?
            ! deltx is coeff. of humidity in virtual temperature.
            !----------------------------------------------------------!
            tsv = tas (x,y,iT) * (one + huss (x,y,iT) * deltx)
            !----------------------------------------------------------!
            ! Air density                                   (kg[a] m-3).
            !----------------------------------------------------------!
            rho = ps (x,y,iT) / (rgas * tsv)
            !----------------------------------------------------------!
            ! Ratio of air to water density             (kg[a] kg[w]-1).
            !----------------------------------------------------------!
            rho3 = rho / rhow
            !----------------------------------------------------------!
            ! Specific humidity of ground               (kg[w] kg[a]-1).
            !----------------------------------------------------------!
            qb = QSAT (tas (x,y,iT), lhe, ps (x,y,iT) / 100.0)
            !----------------------------------------------------------!
            ! Conductance of the atmosphere (mm s-1). Hack for now.
            !----------------------------------------------------------!
            cna = 0.02E3
            !----------------------------------------------------------!
            ! Water available for evaporation over timestep (mm s-1).
            !----------------------------------------------------------!
            evap_max = dz (1) * (theta (1) - theta_m (1,x,y)) / dt
            !----------------------------------------------------------!
            ! Ground surface evaporation rate (mm s-1).
            !----------------------------------------------------------!
            huss (x,y,iT) = qb * rhs (x,y,iT) / 100.0
            qflx_evap_grnd = beta * rho3 * cna * &
                          (qb - qb * rhs (x,y,iT) / 100.0)
            qflx_evap_grnd = MAX (zero, qflx_evap_grnd)
            qflx_evap_grnd = MIN (evap_max,qflx_evap_grnd)
            !----------------------------------------------------------!
            ! Evaporative flux from h2osfc (mm/s).
            !----------------------------------------------------------!
            qflx_ev_h2osfc = rho3 * cna * &
                          (qb - qb * rhs (x,y,iT) / 100.0)
            qflx_ev_h2osfc = MAX (zero, qflx_ev_h2osfc)
            qflx_ev_h2osfc = MIN (evap_max,qflx_ev_h2osfc)
            !----------------------------------------------------------!

!**********************************************************************!
! Start of CESM 'Infiltration'.
!**********************************************************************!

            !----------------------------------------------------------!
            DO I = 1, nlayers
              ! Porosity of soil.
              eff_porosity (I) = MAX (0.01, theta_s(I,x,y))
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Local evaporation                                   (mm/s)
            !----------------------------------------------------------!
            qflx_evap = qflx_evap_grnd
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Partition surface inputs between soil and h2osfc    (mm/s)
            !----------------------------------------------------------!
            qflx_in_soil = (one - frac_h2osfc) * &
                           (qflx_top_soil - qflx_surf)
            qflx_in_h2osfc = frac_h2osfc * (qflx_top_soil - qflx_surf)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Remove evaporation                                  (mm/s)
            !----------------------------------------------------------!
            qflx_in_soil = qflx_in_soil - &
                           (one - frac_h2osfc) * qflx_evap
            qflx_in_h2osfc = qflx_in_h2osfc - &
                             frac_h2osfc * qflx_ev_h2osfc
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Maximum infiltration capacity                       (mm/s)
            !----------------------------------------------------------!
            qinmax = (one - fsat) * MINVAL (hksat (1:3,x,y))
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Infiltration excess runoff -> h2osfc                (mm/s)
            !----------------------------------------------------------!
            qflx_infl_excess = MAX (0.0, qflx_in_soil - &
                               (1.0 - frac_h2osfc) * qinmax)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Soil infiltration                                   (mm/s)
            !----------------------------------------------------------!
            qflx_infl = qflx_in_soil - qflx_infl_excess
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Shift infiltration excess from h2osfc input to surface
            ! runoff                                              (mm/s)
            !----------------------------------------------------------!
            qflx_surf = qflx_surf + qflx_infl_excess
            qflx_infl_excess = zero
            !----------------------------------------------------------!

!**********************************************************************!
! End of CESM 'Infiltration'.
!**********************************************************************!

!**********************************************************************!
! Start of CESM 'SoilWater'.
!**********************************************************************!

            !----------------------------------------------------------!
            ! Water table depth                                     (mm)
            !----------------------------------------------------------!
            zwtmm = 1000.0 * zwt (x,y)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Compute jwt index, the index of the first unsaturated
            ! layer, i.e. the layer right above the water table (-).
            !----------------------------------------------------------!
            jwt = nlayers
            !----------------------------------------------------------!
            ! Allow jwt to equal zero when zwt is in top layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              IF (zwt (x,y) <= (zi (I) / 1000.0)) THEN
                jwt = I - 1
                EXIT
              END IF
            END DO
            IF (INTERACTIVE) &
              write (98,*) 'water table',iTIME,zwt(x,y),jwt
            !----------------------------------------------------------!

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
            rmx (I) = qin (I) - qout (I)
            amx (I) = 0.0                      !  7.136
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
              rmx (I) = qin (I) - qout (I)
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
              rmx (I) =  qin (I) - qout (I)
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
              den   = zc  (I) - zc  (I-1)
              dzq = zq  (I) - zq  (I-1)
              num   = smp (I) - smp (I-1) - dzq
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
              rmx (I) = qin (I) - qout (I)
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
              IF (qcharge_tot > 0.0) THEN ! Rising water table.
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
            rsub_top = rsub_top - xs / dt
            !----------------------------------------------------------!

!**********************************************************************!
! End ofCESM 'Drainage'
!**********************************************************************!

            !----------------------------------------------------------!
            ! Sum all losses from soil column (mm).
            !----------------------------------------------------------!
            w1 = ((1.0 - frac_h2osfc) * (qflx_surf + qflx_evap_grnd) + &
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
              WRITE (*,*) 'qout (nlayers) = ',qout(nlayers)*dt
              WRITE (*,*) 'zwt = ',zwt(x,y)
              WRITE (*,*) 'qcharge = ',qcharge*dt
              WRITE (*,*) 'pr = ',dt*pr(x,y,iT)
              WRITE (*,*) 'beta = ',beta
              WRITE (*,*) 'wa   =',wa(x,y)
              WRITE (*,*) 'theta theta_m h2osoi_liq'
              DO I = 1, nlayers
                WRITE (*,*) I,theta(I),theta_m(I,x,y),h2osoi_liq (I,x,y)
              END DO
              STOP
            END IF ! (ABS (w1 - w0) > 0.1)
            !----------------------------------------------------------!
            !END IF ! INTERACTIVE
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated within the day.
            !----------------------------------------------------------!
            rnf_sum  = rnf_sum  + qflx_surf
            evap_sum = evap_sum + qflx_evap_grnd
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            END DO ! DO NS = 1, NISURF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Grow some plants!
            ! main file model is /store/MODELS/HXd/SOURCE/RINGS.f90
            ! -1.5 MPa is equivalant to about -150,000 mm
            !----------------------------------------------------------!
            DO K = 1, nplants (x,y)
              w_i = (-150000.0 - smp (1)) / (-150000.0 - (-50000.0))
              w_i = MAX (0.0, w_i)
              w_i = MIN (1.0, w_i)
              plant_mass (K,x,y) = plant_mass (K,x,y) +  &
                                   (1500.0 / 365.0) * w_i
            END DO
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
            !----------------------------------------------------------!
            ! Possibly move following to within NS loop.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              theta_sum (I) = theta_sum (I) + theta (I)
              h2osoi_sum_total = h2osoi_sum_total + h2osoi_liq (I,x,y)
              !rnff
            END DO
            !----------------------------------------------------------!

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
        END DO
      END IF ! Soiled grid-box?
    END DO
  END DO
  !--------------------------------------------------------------------!

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
  !write (*,*) 'lambda'
  !write (*,*) lambda
  !write (*,*) 'bsw'
  !write (*,*) 1.0/lambda
END IF

  !--------------------------------------------------------------------!
  ! Free up resources.
  !--------------------------------------------------------------------!
  DEALLOCATE (tas)
  DEALLOCATE (huss)
  DEALLOCATE (ps)
  DEALLOCATE (pr)
  DEALLOCATE (rhs)
  !--------------------------------------------------------------------!

!----------------------------------------------------------------------!
END DO ! iDEC = iDEC_start, iDEC_end
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
    file_name = 'axy'//ydate//'.nc'
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
      USE SHARED
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

