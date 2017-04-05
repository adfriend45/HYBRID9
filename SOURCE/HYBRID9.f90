!======================================================================!
PROGRAM H9
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global simulation of land surface water and carbon fluxes.
! Uses PGF forcings, 1901-2012.
! Andrew D. Friend
! 5th April, 2017
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
! Results go to 'fort.81' in seconds (see later).
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
INTEGER :: ierr              ! MPI error code                       (n).
INTEGER :: ncid              ! netCDF ID                            (n).
INTEGER :: num_procs         ! No. available processors             (n).
INTEGER :: file_free         ! Flag for reading/writing files     (1/0).
INTEGER :: status (MPI_STATUS_SIZE) ! source, tag, and error codes.
INTEGER :: nx_block,ny_block ! x and y dimensions of chunks over Earth.
INTEGER :: i_block,i_block_s ! Chunk indices                        (n).
INTEGER :: x,y,x1,y1         ! Spatial array indices                (n).
INTEGER :: iTIME             ! Timepoint             (day 1/1/1860 = 1).
INTEGER :: iDEC              ! Decade index             (1 = 1901-1910).
INTEGER :: iDEC_start        ! Index of first decade in simulation  (n).
INTEGER :: iDEC_end          ! Index of last decade in simulation   (n).
INTEGER :: jyear             ! Calendar year                       (CE).
INTEGER :: syr               ! First year in decade                (CS).
INTEGER :: eyr               ! Last year in decade                 (CE).
INTEGER :: nlayers           ! Local no. soil layers                (n).
INTEGER :: I,J               ! Generic loop indices                 (n).
INTEGER :: NISURF            ! No. timepoints in day                (n).
INTEGER :: NS                ! Timepoints in day index              (n).
INTEGER :: nt                ! No. timesteps in year             (days).
INTEGER :: NYR               ! No. years in simulation          (years).
!----------------------------------------------------------------------!
! Value of time at beginning of each calendar year from 1860 to 2300.
!----------------------------------------------------------------------!
INTEGER :: time_BOY (2300-1860+1)        !           (day 1/1/1860 = 1).
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
LOGICAL :: INTERACTIVE       ! Interactive simulation? (logical).
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Longitude of beginning of interactive simulation (degrees east).
!----------------------------------------------------------------------!
REAL :: lon_w
!----------------------------------------------------------------------!
! Longitude index of beginning of interactive simulation (n).
!----------------------------------------------------------------------!
INTEGER :: lon_s_w
!----------------------------------------------------------------------!
! Latitude of beginning of interactive simulation (degrees north).
!----------------------------------------------------------------------!
REAL :: lat_w
!----------------------------------------------------------------------!
! Longitude index of beginning of interactive simulation (n).
!----------------------------------------------------------------------!
INTEGER :: lat_s_w
!----------------------------------------------------------------------!
! Coefficients to calculate lon_s_w and lat_s_w.
!----------------------------------------------------------------------!
REAL :: a_lon,b_lon,a_lat,b_lat
!----------------------------------------------------------------------!
! No. longitudes and latitudes in interactive simulation (n).
!----------------------------------------------------------------------!
REAL :: lon_c_w,lat_c_w
!----------------------------------------------------------------------!
! Local value of iTIME (1 = day 1 of iDEC_start).
!----------------------------------------------------------------------!
INTEGER :: iT
!----------------------------------------------------------------------!
! Local value of jyear (1 = year 1 of iDEC_start).               
!----------------------------------------------------------------------!
INTEGER :: iY
!----------------------------------------------------------------------!
! Global array of processor IDs for diagnostics.
!----------------------------------------------------------------------!
INTEGER :: chunk (NX,NY)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Climate variables.
!----------------------------------------------------------------------!
! huss:   Specific humidity at time of maximum temperature  (kgw kga-1).
! pr:     Precipitation flux                               (kg m-2 s-1).
! ps:     Surface air pressure                                     (Pa).
! rhs:    Relative humidity                                  (no units).
! rlds:   Surface downwelling longwave radiation flux in air    (W m-2).
! rsds:   Surface downwelling shortwave radiation flux in air   (W m-2).
! tasmax: Maximum surface air temperature                           (K).
! tasmin: Minimum surface air temperature                           (K).
! tas:    Surface air temperature                                   (K).
! wind:   Wind speed at 10 meter                                (m s-1).
!----------------------------------------------------------------------!
! Surface Air Temperature (K).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: tas (:,:,:)
!----------------------------------------------------------------------!
! Specific humidity at time of Maximum Temperature (kg kg-1).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: huss (:,:,:)
!----------------------------------------------------------------------!
! Surface air pressure (Pa).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: ps (:,:,:)
!----------------------------------------------------------------------!
! Precipitation flux (kg m-2 s-1).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: pr (:,:,:)
!----------------------------------------------------------------------!
! Relative humidity (no units).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: rhs (:,:,:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Annual sums for diagnostics.
!----------------------------------------------------------------------!
REAL :: rnf_sum         ! Sum of runoff over 1 yr              (mm s-1).
REAL :: evap_sum        ! Sum of evaporation over 1 yr         (mm s-1).
!----------------------------------------------------------------------!
REAL :: tas_sum         ! Sum of daily tas over 1 yr                (K).
REAL :: huss_sum        ! Sum of daily huss over 1 yr         (kg kg-1).
REAL :: ps_sum          ! Sum of daily ps over 1 yr                (Pa).
REAL :: pr_sum          ! Sum of daily pr over 1 yr        (kg m-2 s-1).
REAL :: rhs_sum         ! Sum of daily rhs over 1 yr             (%age).
!----------------------------------------------------------------------!
REAL :: theta_sum_total ! Sum of daily total soil water over 1 yr  (mm).
!----------------------------------------------------------------------!
! Water in each soil layer                                   (mm3 mm-3).
!----------------------------------------------------------------------!
REAL, ALLOCATABLE :: theta_sum (:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Miscellaneous variables.
!----------------------------------------------------------------------!
REAL :: thetan ! New theta for testing over/undersaturation   (mm3 mm-3)
REAL :: dtheta ! Time derivative of theta                 (mm3 mm-3 s-1)
REAL :: q_infl ! Water infiltration into surface soil layer    (mm s^-1)
REAL :: qflx_evap_soi ! Ground surface evaporation rate           (mm/s)
REAL :: dt     ! Integration timestep                                (s)
REAL :: dflux  ! Negative runoff removal                        (mm s-1)
REAL :: drnf   ! Correction to runoff                           (mm s-1)
REAL :: slope  ! Slope                                           (ratio)
REAL :: num    ! 
REAL :: den    ! 
REAL :: dpsie  ! 
REAL :: zwt    ! Water table depth, positive downwards              (mm)
REAL :: BET    ! 
!----------------------------------------------------------------------!
! Water available for evaporation over timestep                (mm s-1).
!----------------------------------------------------------------------!
REAL :: evap_max
!----------------------------------------------------------------------!
! Soil moisture constraint on evaporative flux                    (0-1).
!----------------------------------------------------------------------!
REAL :: beta
!----------------------------------------------------------------------!
! Air density                                               (kg[a] m-3).
!----------------------------------------------------------------------!
REAL :: rho
!----------------------------------------------------------------------!
! Ratio of air to water density                         (kg[a] kg[w]-1).
!----------------------------------------------------------------------!
REAL :: rho3
!----------------------------------------------------------------------!
! Virtual potential surface temperature                             (K).
!----------------------------------------------------------------------!
REAL :: tsv
!----------------------------------------------------------------------!
! Specific humidity of ground                           (kg[w] kg[a]-1).
!----------------------------------------------------------------------!
REAL :: qb
!----------------------------------------------------------------------!
! Saturation water vapour mixing ratio                  (kg[w] kg[a]-1).
!----------------------------------------------------------------------!
REAL :: QSAT
!----------------------------------------------------------------------!
! Conductance of the atmosphere                                (mm s-1).
!----------------------------------------------------------------------!
REAL :: cna
!----------------------------------------------------------------------!
! Precipitation flux                                           (mm s-1).
!----------------------------------------------------------------------!
REAL :: prec
!----------------------------------------------------------------------!
! Surface runoff                                               (mm s-1).
!----------------------------------------------------------------------!
REAL :: rnf
!----------------------------------------------------------------------!
! For soil water diagnostics.
!----------------------------------------------------------------------!
REAL :: w0,w1
!----------------------------------------------------------------------!
! Saturated hydraulic conductivity at soil layer interface      (mm s-1).
!----------------------------------------------------------------------!
REAL :: k_s_int
!----------------------------------------------------------------------!
! "s" at interface of layer.
!----------------------------------------------------------------------!
REAL :: s1
!----------------------------------------------------------------------!
! k*s**(2b+2)
!----------------------------------------------------------------------!
REAL :: s2
!----------------------------------------------------------------------!
! Soil size distribution at layer interface                   (unitless)
!----------------------------------------------------------------------!
REAL :: lambda_int
!----------------------------------------------------------------------!
! Soil wetness (-).
!----------------------------------------------------------------------!
REAL :: s_node
!----------------------------------------------------------------------!
! Soil matric potential in aquifer layer (mm).
!----------------------------------------------------------------------!
REAL :: smp1
!----------------------------------------------------------------------!
! d(smp)/d(vol_liq) in aquifer layer (mm/(mm^3 mm^-3))
!----------------------------------------------------------------------!
REAL :: dsmpdw1
!----------------------------------------------------------------------!
! Underground runoff from each soil layer                      (mm s-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: rnff
!----------------------------------------------------------------------!
! Matric potentials of soil layers                                 (mm).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: smp
!----------------------------------------------------------------------!
! Equilibrium soil matric potentials of soil layers                (mm).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: zq
!----------------------------------------------------------------------!
! Soil water fluxes upwards                                    (mm s-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: f
!----------------------------------------------------------------------!
! Soil water inflow at top of layer                            (mm s^-1)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: qin
!----------------------------------------------------------------------!
! Soil water outflow at bottom of layer                        (mm s^-1)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: qout
!----------------------------------------------------------------------!
! Soil layer boundaries with reference to surface                  (mm).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: zb
!----------------------------------------------------------------------!
! Soil layer thicknesses                                           (mm).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dz
!----------------------------------------------------------------------!
! Hydraulic conductivities at soil layer interfaces            (mm s-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: xk
!----------------------------------------------------------------------!
! d(smp)/d(vol_liq) (mm/(mm^3 mm^-3))
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dsmpdw
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I)) (mm s-1) /(mm^3 mm^-3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqodw1
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I+1)) (mm s-1) /(mm^3 mm^-3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqodw2
!----------------------------------------------------------------------!
! d(hk)/d(vol_liq) (mm s-1) / (mm^3 mm^-3)
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dhkdw
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I-1)) ((mm s-1) / (mm^3 mm^-3)).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dqidw0
!----------------------------------------------------------------------!
! d(qin)/d(vol_liq(I)) ((mm s-1) / (mm^3 mm^-3)).
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
! Change in soil water (m^3/m^3).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: dwat,dwat2
!----------------------------------------------------------------------!
! For tridiagonal solution (-).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: GAM
!----------------------------------------------------------------------!
! Volumetric soil water (mm^3/mm^3).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: theta
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Miscellaneous parameters.
!----------------------------------------------------------------------!
! Truncation limit for numerical tests.
!----------------------------------------------------------------------!
REAL, PARAMETER :: trunc = 1.0E-8
!----------------------------------------------------------------------!
! Conductivity for underground runoff, xkud from GHY           (mm s-1).
!----------------------------------------------------------------------!
REAL, PARAMETER :: xkud = 2.78E-5 * 1.0E3
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Index of soil layer right above water table (-).
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

!***************************************************!
! Attempt to get lats and lons into axy files.
! lat (1) is 89.75 in the PGF files.
DO y = 1, NY
  lat_all (y) = 89.75 - (y - 1) * 0.5
END DO
DO x = 1, NX
  lon_all (x) = -179.75 + (x - 1) * 0.5
END DO
!***************************************************!

!----------------------------------------------------------------------!
! Allocate soil arrays over layers.
!----------------------------------------------------------------------!
ALLOCATE (f      (nsoil_layers_max+1)) ! Upwards water fluxes  (mm s-1).
ALLOCATE (qin    (nsoil_layers_max+1)) ! Inflow at top         (mm s^-1)
ALLOCATE (qout   (nsoil_layers_max+1)) ! Outflow at top        (mm s^-1)
ALLOCATE (dsmpdw (nsoil_layers_max+1)) ! d(smp)/d(w)   (mm/(mm^3 mm^-3))
!----------------------------------------------------------------------!
! d(qout)/d(vol_liq(I)) (mm s-1) / (mm^3 mm^-3)
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
ALLOCATE (zb     (0:Nlevgrnd))         ! Layer boundaries           (mm)
ALLOCATE (dz     (1:Nlevgrnd))         ! Layer thicknesses          (mm)
ALLOCATE (zsoi   (1:Nlevgrnd))         ! Soil layer nodes           (mm)
ALLOCATE (zc_o   (1:nsoil_layers_max)) ! Soil layer nodes           (mm)
ALLOCATE (zc     (1:Nlevgrnd))         ! Ground layer nodes         (mm)
ALLOCATE (smp    (nsoil_layers_max))   ! Matric potentials          (mm)
ALLOCATE (zq     (nsoil_layers_max+1)) ! Eqm. water potentials      (mm)
!----------------------------------------------------------------------!
! Volumetric soil water (mm^3 mm^-3).
!----------------------------------------------------------------------!
ALLOCATE (theta (nsoil_layers_max))
!----------------------------------------------------------------------!
ALLOCATE (theta_sum(nsoil_layers_max)) ! For diagnostics   (mm^3 mm^-3).
!----------------------------------------------------------------------!
! Hydraulic conductivities                                     (mm s-1).
!----------------------------------------------------------------------!
ALLOCATE (xk  (nsoil_layers_max+1))
!----------------------------------------------------------------------!
! Underground runoff from each soil layer                      (mm s-1).
!----------------------------------------------------------------------!
ALLOCATE (rnff (nsoil_layers_max))
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
  ! Read soil layer boundaries, negative as go down from surface,
  ! which is boundary 1 at 0 mm (mm).
  !--------------------------------------------------------------------!
  DO I = 0, Nlevgrnd
    READ (10,*) zb (I)
  END DO
  !--------------------------------------------------------------------!
  CLOSE (10)
!----------------------------------------------------------------------!
END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set integration timestep                                          (s).
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
! Compute soil layer thicknesses and centres (mm).
! Nodes are positive  values below surface.
! Layer indices start with 1 at surface.
!----------------------------------------------------------------------!
DO I = 1, Nlevgrnd
  dz (I) = zb (I) - zb (I-1)
END DO
DO I = 1, Nlevgrnd
  zc (I) = zb (I) - dz (I) / 2.0
END DO
!----------------------------------------------------------------------!
! Set soil layer node depths for diagnostics (mm).
!----------------------------------------------------------------------!
DO I = 1, nsoil_layers_max
  zc_o (I) = zc (I)
END DO
!----------------------------------------------------------------------!
DO I = 1, Nlevgrnd
  zsoi (I) = 0.025 * (EXP (0.5 * (I - 0.5)) - 1.0) ! Node depths (mm).
END DO
DO I = 1, Nlevgrnd
  write (*,*) I,dz(I),zc(I),1.0E3*zsoi(I)
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
ALLOCATE (data_in_2DI (lon_c,lat_c))  ! Generic 2D input integer array.
ALLOCATE (soil_tex    (lon_c,lat_c))  ! Soil textures.
ALLOCATE (lon         (lon_c))        ! Longitudes.
ALLOCATE (lat         (lat_c))        ! Latitudes.
ALLOCATE (block_sub   (lon_c,lat_c))  ! Processor IDs.
ALLOCATE (axy_tas     (lon_c,lat_c,NYR)) ! Accumulated tas.
ALLOCATE (axy_huss    (lon_c,lat_c,NYR)) ! Accumulated huss.
ALLOCATE (axy_ps      (lon_c,lat_c,NYR)) ! Accumulated ps.
ALLOCATE (axy_pr      (lon_c,lat_c,NYR)) ! Accumulated pr.
ALLOCATE (axy_rhs     (lon_c,lat_c,NYR)) ! Accumulated rhs.
!----------------------------------------------------------------------!
! Accumulated theta.
!----------------------------------------------------------------------!
ALLOCATE (axy_theta   (nsoil_layers_max,lon_c,lat_c,NYR))
ALLOCATE (axy_theta_total (lon_c,lat_c,NYR))
!----------------------------------------------------------------------!
ALLOCATE (axy_rnf     (lon_c,lat_c,NYR)) ! Accumulated runoff.
ALLOCATE (axy_evap    (lon_c,lat_c,NYR)) ! Accumulated evap.
!----------------------------------------------------------------------!
! Liquid soil water mass (kg/m^2).
!----------------------------------------------------------------------!
ALLOCATE (h2osoi_liq (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunks of soil properties of one soil layer read in at
! 30 arc-seconds.
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1_in(lon_c*60,lat_c*60)) ! 0.001xcm3 cm-3.
ALLOCATE (k_s_l1_in    (lon_c*60,lat_c*60)) ! cm day-1.
ALLOCATE (lambda_l1_in (lon_c*60,lat_c*60)) ! 0.001xunitless.
ALLOCATE (psi_s_l1_in  (lon_c*60,lat_c*60)) ! cm.
!----------------------------------------------------------------------!
! Chunks of soil properties of one soil layer gridded to half-degree.
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1   (lon_c,lat_c)) ! 0.001xcm3 cm-3.
ALLOCATE (k_s_l1       (lon_c,lat_c)) ! cm day-1.
ALLOCATE (lambda_l1    (lon_c,lat_c)) ! 0.001xunitless.
ALLOCATE (psi_s_l1     (lon_c,lat_c)) ! cm.
!----------------------------------------------------------------------!
! Chunk of saturated volumetric soil water (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_s (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of saturated soil hydraulic conductivity (mm s-1).
!----------------------------------------------------------------------!
ALLOCATE (k_s (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of soil pore size distribution index (unitless).
!----------------------------------------------------------------------!
ALLOCATE (lambda (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of Clapp and Hornberger "b".
!----------------------------------------------------------------------!
ALLOCATE (bsw (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of soil saturated capillary/water potential (mm).
!----------------------------------------------------------------------!
ALLOCATE (psi_s  (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of minimum volumetric soil matric potential (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_m (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise global diagnostic arrays with fill (i.e. NaN and zero)
! values.
!----------------------------------------------------------------------!
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
CALL READ_NET_CDF_2D
!----------------------------------------------------------------------!
soil_tex (:,:) = data_in_2DI (:,:)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read in soil properties for this block.
! Saturated water contents          , theta_s, are 0.001*cm^3/cm^3.
! Saturated hydraulic conductivities, k_s    , are    cm/day.
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
  ! Saturated water potential read in as cm,
  ! converted to mm.
  !--------------------------------------------------------------------!
  DO y = 1, lat_c
    DO x = 1, lon_c
      theta_s (I,x,y) = theta_s_l1(x,y) / 1.0E3
      k_s     (I,x,y) = 10.0 * k_s_l1 (x,y) / 86400.0
      lambda  (I,x,y) = lambda_l1 (x,y) / 1.0E3
      psi_s   (I,x,y) = 10.0 * psi_s_l1 (x,y)
      !----------------------------------------------------------------!
      ! To stop potential dbz.
      !----------------------------------------------------------------!
      lambda  (I,x,y) = MAX (lambda (I,x,y),trunc)
      !----------------------------------------------------------------!
      ! Clapp and Hornberger "b".
      !----------------------------------------------------------------!
      bsw (I,x,y) = 1.0 / lambda (I,x,y)
      !----------------------------------------------------------------!
    END DO
  END DO
  !--------------------------------------------------------------------!
END DO ! I = 1, nsoil_layers_max ! Loop over soil layers.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write out saturated water potentials for diagnostics (mm).
!----------------------------------------------------------------------!
file_name = 'psi_s.nc'
WRITE (*,*) 'Writing to ',TRIM(file_name)
CALL WRITE_NET_CDF_3DR_soils
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Minimum soil volumetric water contents (mm^3 mm^-3).
! Following assumes theta_m is hygroscopic water at -31 bar,
! Not sure if this is the amount that would stop flow - too high?
!----------------------------------------------------------------------!
DO y = 1, lat_c
  DO x = 1, lon_c
    DO I = 1, 8
      !theta_m (I,x,y) = theta_s (I,x,y) / 100.0
      ! 
      theta_m (I,x,y) = theta_s (I,x,y) * ((-3.1E9 / (1000.0 * 9.8)) / &
                        psi_s (I,x,y)) ** (-lambda (I,x,y))
      !----------------------------------------------------------------!
      ! Initial liquid water in each soil layer (kg/m^2).
      !----------------------------------------------------------------!
      h2osoi_liq (I,x,y) = 0.5 * theta_s (I,x,y) * dz (I) * rhow &
                           / 1000.0
      !----------------------------------------------------------------!
    END DO
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
write (*,*) 'mid = ',my_id,lat_c,lon_c,lat_s,lon_s
!----------------------------------------------------------------------!
! Write soil textures for diagnostics.
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
  ! with soil in both texture and theta_s fields.
  !--------------------------------------------------------------------!
  DO y = 1, lat_c
    DO x = 1, lon_c
      IF ((soil_tex (x,y) > 0) .AND. (soil_tex (x,y) /= 13) .AND. &
        SUM (theta_s (:,x,y)) > trunc) THEN
        nlayers = nsoil_layers_max
        !--------------------------------------------------------------!
        ! Loop over years in decade.
        !--------------------------------------------------------------!
        DO jyear = syr, eyr
          !------------------------------------------------------------!
          ! Intialise annual diagnostics.
          !------------------------------------------------------------!
          rnf_sum  = zero
          evap_sum = zero
          tas_sum  = zero
          huss_sum = zero
          ps_sum   = zero
          pr_sum   = zero
          rhs_sum  = zero
          theta_sum_total = zero
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
            ! Initial total soil water for diagnostics (mm).
            ! Also compute volumetric water (mm^3/mm^3).
            !----------------------------------------------------------!
            w0 = prec * dt
            DO I = 1, nlayers
              !--------------------------------------------------------!
              w0 = w0 + h2osoi_liq (I,x,y)
              !--------------------------------------------------------!
              theta (I) = MAX (h2osoi_liq (I,x,y), 1.0E-6) &
                          / (dz (I) * rhow / 1000.0)
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Holder for now.
            ! Water table depth (mm).
            !----------------------------------------------------------!
            zwt = 5000.0
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
              IF (zwt <= zc (I)) THEN
                jwt = I - 1
                EXIT
              END IF
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Equilibrium matric potentials of each layer (mm).
            ! Assumed relationshop from Eqn. 7.86 of
            ! Oleson et al. (2013), with fixed water table depth of 4 m.
            ! Requires psi_s at WTD, assume is same as bottom layer for
            ! now. The WTD depth cancels out for the flux equation, but
            ! is really needed here for the value of psi_s, but just
            ! assume is always in bottom layer for now.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              zq (I) = psi_s (nlayers,x,y) + zwt - zc (I)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! If water table is below soil column calculate zq for 
            ! next layer below.
            !----------------------------------------------------------!
            I = nlayers
            !----------------------------------------------------------!
            IF (jwt == nlayers) THEN
              zq (I+1) = psi_s (nlayers,x,y) + zwt - zc (I)
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              ! Water saturation fraction of layers (0-1).
              !--------------------------------------------------------!
              IF (theta_s (I,x,y) > trunc) THEN
                s_node = theta (I) / theta_s (I,x,y)
              ELSE
                s_node = 0.01
              END IF
              s_node = MIN (one, s_node)
              !--------------------------------------------------------!
              ! Matric potentials and derivatives at each node depth 
              ! (mm).
              ! This is Eqn. 7.94 of Oleson et al. (2013).
              !--------------------------------------------------------!
              smp (I) = psi_s (I,x,y) * s_node ** (-bsw (I,x,y))
              smp (I) = MAX (-1.0E8, smp (I))
              !--------------------------------------------------------!
              dsmpdw (I) = (-bsw (I,x,y)) * smp (I) / &
                           (s_node * theta_s (I,x,y))
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Hydraulic conductivities at interfaces between layers
            ! based on Eqn. 7.89 of Oleson et al. (2013) (mm s-1).
            ! Refers to bottom of layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              ! Compute saturated hydraulic conductivity at the
              ! interface as mean of the two layers (mm/s).
              !--------------------------------------------------------!
              ! "s" at interface of layer.
              !--------------------------------------------------------!
              s1 =  0.5 * (theta   (I) + &
                           theta   (MIN(nlayers-1,I+1))) / &
                   (0.5 * (theta_s (I,x,y) + &
                           theta_s (MIN(nlayers-1,I+1),x,y)))
              s1 = MIN (1.0, s1)
              s2 = k_s (I,x,y) * s1 ** (2.0 * bsw (I,x,y) + 2.0)
              !--------------------------------------------------------!
              ! Hydraulic conductivity at layer intervace (mm/s).
              ! Need to add ice impedance.
              !--------------------------------------------------------!
              xk (I) = s1 * s2
              !--------------------------------------------------------!
              ! d(hk)/d(vol_liq) (mm s-1) / (mm^3 mm^-3)
              !--------------------------------------------------------!
              dhkdw (I) = (2.0 * bsw (I,x,y) + 3.0) * s2 * &
                          (1.0 / (theta_s (I,x,y) + &
                                  theta_s (MIN (nlayers-1,I+1),x,y)))
              !--------------------------------------------------------!
              k_s_int = 0.5 * (k_s (I,x,y) + k_s (I+1,x,y))
              !--------------------------------------------------------!
              ! Compute soil size distribution index at the
              ! interface as mean of the two layers (unitless).
              !--------------------------------------------------------!
              lambda_int = 0.5 * (lambda (I,x,y) + lambda (I+1,x,y))
              !--------------------------------------------------------!
            END DO

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
            ! Evaporative flux (mm s-1).
            !----------------------------------------------------------!
            huss (x,y,iT) = qb * rhs (x,y,iT) / 100.0
            !evap = beta * rho3 * cna * (qb - huss (x,y,iT))
            qflx_evap_soi = beta * rho3 * cna * &
                            (qb - qb * rhs (x,y,iT) / 100.0)
            qflx_evap_soi = MAX (zero, qflx_evap_soi)
            qflx_evap_soi = MIN (evap_max,qflx_evap_soi)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Surface runoff (mm s-1).
            !----------------------------------------------------------!
            rnf = (theta (1) * dz (1) + prec * dt - &
                  theta_s (1,x,y) * dz (1)) / dt
            rnf = MAX (zero, rnf)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water infiltration into surface soil layer (mm s^-1).
            !----------------------------------------------------------!
            q_infl = prec - rnf - qflx_evap_soi
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Underground runoff from each layer               (mm s-1).
            ! Hack for now.
            ! Drainage based on Eqn. 7.168 of Oleson et al. (2013).
            ! (mm s-1).
            !----------------------------------------------------------!
            slope = 0.035 ! Radians of slope (e.g. 2 degrees = 0.035).
            !fdrai = 2.5E-3 ! Decay factor (mm-1).
            ! Need to diagnose WTD.
            !WTD = 4000.0 ! Water table depth (mm).
            ! Maximum drainage when the water table depth is at the
            ! surface (mm s-1).
            !qdrai_max = 10.0 * SIN (slope)
            ! Drainage (mm s-1).
            !qdrai = qdrai_max * EXP (-fdrai * WTD)
            rnff (:) = 0.0
            !rnff (nlayers) = drai
            !DO I = 1, nlayers
            !  rnff (I) = xk (I) * slope * dz (I) * &
            !             theta (I) / theta_s (I,x,y)
            !END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water fluxes at each layer interface (mm s-1).
            ! Uses Darcy's Law as in Eqn. 7.116 of Oleson et al. (2013).
            ! Also, to make it positive upwards, denominator is reversed
            ! as z is positive downwards in CLM but negative downwards
            ! in HYBRID9.
            ! f (I) is the interface at the bottom of layer I, so
            ! f (1) is at the bottom of the first layer.
            !----------------------------------------------------------!
            !DO I = 1, nlayers - 1
              !--------------------------------------------------------!
            !  f (I) = -xk (I) * ((smp (I) - smp (I+1)) + &
            !          (zq (I+1) - zq (I))) / &
            !          (zc (I) - zc (I+1))
              !--------------------------------------------------------!
            !END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Bottom boundary flux (mm s-1).
            !----------------------------------------------------------!
            !f (nlayers) = zero
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Use tridiagonal solution to compute new soil water levels
            ! in each layer given fluxes.
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Set up vectors for top layer based on section 7.4.2.2 O13.
            !----------------------------------------------------------!
            I = 1
            qin (I) = q_infl
            den   = (zc  (I+1) - zc  (I))
            dpsie = (zq  (I+1) - zq  (I))
            num   = (smp (I+1) - smp (I)) - dpsie
            qout (I) = -xk (I) * num / den
            dqodw1 (I) = -(-xk (I) * dsmpdw (I)   + num * dhkdw (I)) / &
                         den
            dqodw2 (I) = -(-xk (I) * dsmpdw (I+1) + num * dhkdw (I)) / &
                         den
            rmx (I) = qin (I) - qout (I)
            amx (I) = 0.0
            bmx (I) = dz (I) / dt + dqodw1 (I)
            cmx (I) = dqodw2 (I)
            !----------------------------------------------------------!
            ! Nodes I = 2, to I = nlayers-1.
            !----------------------------------------------------------!
            DO I = 2, nlayers - 1
              !--------------------------------------------------------!
              den   = zc  (I) - zc  (I-1)
              dpsie = zq  (I) - zq  (I-1)
              num   = smp (I) - smp (I-1) - dpsie
              qin (I) = -xk (I-1) * num / den
              dqidw0 (I) = -(-xk (I-1) * dsmpdw (I-1) + &
                           num * dhkdw (I-1)) / den
              dqidw1 (I) = -( xk (I-1) * dsmpdw (I)   + &
                           num * dhkdw (I-1)) / den
              den    = zc (I+1) - zc (I)
              dpsie  = zq (I+1) - zq (I)
              num   = (smp (I+1) - smp (I)) - dpsie
              qout (I) = -xk (I) * num / den
              dqodw1 (I) = -(-xk (I) * dsmpdw (I)   + &
                           num * dhkdw (I)) / den
              dqodw2 (I) = -(-xk (I) * dsmpdw (I+1) + &
                           num * dhkdw (I)) / den
              rmx (I) = qin (I) - qout (I)
              amx (I) = -dqidw0 (I)
              bmx (I) = dz (I) / dt - dqidw1 (I) + dqodw1 (I)
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
              den   = zc  (I) - zc  (I-1)
              dpsie = zq  (I) - zq  (I-1)
              num   = smp (I) - smp (I-1) - dpsie
              qin (I) = -xk (I-1) * num / den
              dqidw0 (I) = -(-xk (I-1) * dsmpdw (I-1) + &
                           num * dhkdw (I-1)) / den
              dqidw1 (I) = -(-xk (I-1) * dsmpdw (I) + &
                           num * dhkdw (I-1)) / den
              qout   (I) = 0.0
              dqodw1 (I) = 0.0
              rmx (I) =  qin (I) - qout (I)
              amx (I) = -dqidw0 (I)
              bmx (I) =  dz (I) / dt - dqidw1 (I) + dqodw1 (I)
              cmx (I) =  0.0
              !--------------------------------------------------------!
              ! Next set up aquifer layer, hydrologically inactive.
              !--------------------------------------------------------!
              rmx (I+1) = 0.0
              amx (I+1) = 0.0
              bmx (I+1) = dz (I+1) /dt
              cmx (I+1) = 0.0
              !--------------------------------------------------------!
            ELSE ! Water table is below soil column.
              !--------------------------------------------------------!
              ! Compute aquifer soil moisture as average of layer 10
              ! and saturation.
              !--------------------------------------------------------!
              s_node = MAX (0.5 * (1.0 + theta (I) / &
                       theta_s (I,x,y)), 0.01)
              s_node = MIN (1.0, s_node)
              !--------------------------------------------------------!
              ! Compute smp for aquifer layer (mm).
              !--------------------------------------------------------!
              smp1 = psi_s (I,x,y) * s_node ** (-bsw (I,x,y))
              smp1 = MAX (-1.0E8, smp1)
              !--------------------------------------------------------!
              ! Compute dsmpdw for aquifer layer.
              !--------------------------------------------------------!
              dsmpdw1 = -bsw (I,x,y) * smp1 / (s_node * theta_s (I,x,y))
              !--------------------------------------------------------!
              ! First set up bottom layer of soil column.
              !--------------------------------------------------------!
              den   = zc  (I) - zc  (I-1)
              dpsie = zq  (I) - zq  (I-1)
              num   = smp (I) - smp (I-1) - dpsie
              qin (I) = -xk (I) * num / den
              dqidw0 = -(-xk(I-1) * dsmpdw (I-1) + num * dhkdw (I-1)) &
                       / den
              dqidw1 = -(-xk(I-1) * dsmpdw (I)   + num * dhkdw (I-1)) &
                       / den
              den   = zc (I+1) - zc (I)
              dpsie = zq (I+1) - zq (I)
              num = smp1 - smp (I) - dpsie
              qout (I) = -xk (I) * num / den
              dqodw1 (I) = -(-xk (I) * dsmpdw (I) + num * dhkdw (I)) &
                           / den
              dqodw2 (I) = -(-xk (I) * dsmpdw1 + num * dhkdw (I)) &
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
              dqidw0 (I+1) = -(-xk (I) * dsmpdw (I) + num * dhkdw (I)) &
                             / den
              dqidw1 (I+1) = -(-xk (I) * dsmpdw1    + num * dhkdw (I)) &
                             / den
              qout   (I+1) = 0.0 ! Zero-flow bottom boundary condition.
              dqodw1 (I+1) = 0.0 ! Zero-flow bottom boundary condition.
              rmx (I+1) = qin (I+1) - qout (I+1)
              amx (I+1) = -dqidw0 (I+1)
              bmx (I+1) = dz (I+1) / dt - dqidw1 (I+1) + dqodw1 (I+1)
              cmx (I+1) = 0.0
              !--------------------------------------------------------!

            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Solve for dwat2 using tridiagonal system of equations.
            !----------------------------------------------------------!
            IF (bmx (1) == 0.0) STOP 'Problem with tridiagonal 1.'
            BET = bmx (1)
            dwat2 (1) = rmx (1) / BET
            DO I = 2, nlayers + 1
              GAM (I) = cmx (I-1) / BET
              BET = bmx (I) - amx (I) * GAM (I)
              IF (BET == 0.0) THEN
                WRITE (*,*) I,bmx(I),amx(I),GAM(I)
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
            ! Renew the mass of liquid water also compute qcharge from
            ! dwat in aquifer layer and update in drainage for case
            ! iwt < nlayers.
            DO I = 1, nlayers
              h2osoi_liq (I,x,y) = h2osoi_liq (I,x,y) + &
                                   dwat2 (I) * dz (I)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            w1 = (rnf + qflx_evap_soi + SUM (rnff (:))) * dt
            !----------------------------------------------------------!
            DO I = 1, nlayers
              !--------------------------------------------------------!
              ! New liquid soil water mass (kg/m^2).
              !--------------------------------------------------------!
              !h2osoi_liq (I,x,y) = theta (I) * dz (I) * rhow / 1000.0
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
            ! Check for conservation of water.
            !----------------------------------------------------------!
            !IF (ABS (w1 - w0) > 0.001) THEN
            !  WRITE (*,*)
            !  WRITE (*,*) 'Water imbalance > 0.001 mm ',w1-w0
            !  WRITE (*,*) 'my_id x y ',my_id,x,y
            !  WRITE (*,*) 'lon, lat, iTIME ',lon(x),lat(y),iTIME
            !  WRITE (*,*) 'theta sum = ',SUM (theta(1:nlayers))
            !  WRITE (*,*) 'rnff  = ',dt*SUM (rnff(1:nlayers))
            !  WRITE (*,*) 'rnf qflx_evap_soi = ',dt*rnf,dt*qflx_evap_soi
            !  WRITE (*,*) 'pr = ',dt*pr(x,y,iTIME)
            !  WRITE (*,*) 'beta = ',beta
            !  DO I = 1, nlayers
            !    WRITE (*,*) I,theta(I),theta_m(I,x,y)
            !  END DO
              !STOP
            !END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated within the day.
            !----------------------------------------------------------!
            rnf_sum  = rnf_sum  + rnf
            evap_sum = evap_sum + qflx_evap_soi
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            END DO ! DO NS = 1, NISURF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            IF (INTERACTIVE) THEN
              WRITE (99,*) iT,theta(1),theta(2),theta(8),prec*dt*48.0
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated over days in year.
            !----------------------------------------------------------!
            tas_sum  = tas_sum  + tas  (x,y,iT)
            huss_sum = huss_sum + huss (x,y,iT)
            ps_sum   = ps_sum   + ps   (x,y,iT)
            pr_sum   = pr_sum   + pr   (x,y,iT)
            rhs_sum  = rhs_sum  + rhs  (x,y,iT)
            !----------------------------------------------------------!
            ! Possibly move following to within NS loop.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              theta_sum (I) = theta_sum (I) + theta (I)
              theta_sum_total = theta_sum_total + theta (I) * dz (I)
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
          axy_rnf  (x,y,iY) = rnf_sum  / FLOAT (nt * NISURF) ! mm s-1
          axy_evap (x,y,iY) = evap_sum / FLOAT (nt * NISURF) ! mm s-1
          !------------------------------------------------------------!
          axy_tas  (x,y,iY) = tas_sum  / FLOAT (nt) ! 
          axy_huss (x,y,iY) = huss_sum / FLOAT (nt)
          axy_ps   (x,y,iY) = ps_sum   / FLOAT (nt)
          axy_pr   (x,y,iY) = pr_sum   / FLOAT (nt)
          axy_rhs  (x,y,iY) = rhs_sum  / FLOAT (nt)
          !------------------------------------------------------------!
          DO I = 1, nlayers
            axy_theta (I,x,y,iY) = theta_sum (I) / FLOAT (nt)
          END DO
          axy_theta_total (x,y,iY) = theta_sum_total / FLOAT (nt)
          !------------------------------------------------------------!
        END DO
      END IF ! Soiled grid-box?
    END DO
  END DO
  !--------------------------------------------------------------------!

IF (INTERACTIVE) THEN
  write (*,*) theta
  write (*,*) smp
  write (*,*) zq
  !write (*,*) lambda
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
DEALLOCATE (soil_tex)
DEALLOCATE (block_sub)
DEALLOCATE (data_in_2DI)
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

