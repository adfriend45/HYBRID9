!======================================================================!
PROGRAM H9
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global simulation of land surface water and carbon fluxes.
! Uses PGF forcings, 1901-2012.
! Andrew D. Friend
! 28th November, 2016
! Just some edits.
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
CHARACTER (LEN = 200) :: file_in ! Generic filename.
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
REAL, ALLOCATABLE ::pr (:,:,:)
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
REAL :: wz_sum          ! Sum of daily wz over 1 yr                (mm).
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
REAL :: thetan ! New theta for testing over/undersaturation  (mm3 mm-3).
REAL :: dtheta ! Derivative of theta                     (mm3 mm-3 s-1).
REAL :: evap   ! Evaporative flux                              (mm s-1).
REAL :: dt     ! Integration timestep                               (s).
REAL :: dflux  ! Negative runoff removal                       (mm s-1).
REAL :: drnf   ! Correction to runoff                          (mm s-1).
REAL :: slope  ! Slope                                          (ratio).
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
! Water table depth	                                           (mm).
!----------------------------------------------------------------------!
REAL :: wz  
!----------------------------------------------------------------------!
! soil matric potential for calculating water table                (mm).
!----------------------------------------------------------------------!
REAL :: hmat 
!----------------------------------------------------------------------!



!----------------------------------------------------------------------!
! Underground runoff from each soil layer                      (mm s-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: rnff
!----------------------------------------------------------------------!
! Water saturation fraction of layers                             (0-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: wv
!----------------------------------------------------------------------!
! Water potentials of soil layers                                  (mm).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: h
!----------------------------------------------------------------------!
! Soil water fluxes upwards                                    (mm s-1).
!----------------------------------------------------------------------!
REAL, DIMENSION (:), ALLOCATABLE :: f
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

!----------------------------------------------------------------------!
! Miscellaneous parameters.
!----------------------------------------------------------------------!
! Truncation limit for numerical tests.
!----------------------------------------------------------------------!
REAL, PARAMETER :: trunc = 1.0E-8
!----------------------------------------------------------------------!
! Saturated hydraulic conductivity from Zeng and Decker        (mm s-1).
!----------------------------------------------------------------------!
REAL, PARAMETER :: xks = 0.0038
!----------------------------------------------------------------------!
! Conductivity for underground runoff, xkud from GHY           (mm s-1).
!----------------------------------------------------------------------!
REAL, PARAMETER :: xkud = 2.78E-5 * 1.0E3
!----------------------------------------------------------------------!
! Textural parameter taken from Zeng and Decker (2009).
!----------------------------------------------------------------------!
REAL, PARAMETER :: B = 9.3
!----------------------------------------------------------------------!
! Saturated soil water potential from Zeng and Decker (2009)       (mm).
!----------------------------------------------------------------------!
REAL, PARAMETER :: swp_s = -227.0
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
ALLOCATE (wv  (nsoil_layers_max))      ! Saturation fraction      (0-1).
ALLOCATE (f   (nsoil_layers_max+1))    ! Upwards water fluxes  (mm s-1).
ALLOCATE (zb  (nsoil_layers_max+1))    ! Layer boundaries          (mm).
ALLOCATE (dz  (nsoil_layers_max))      ! Layer thicknesses         (mm).
ALLOCATE (zc  (nsoil_layers_max))      ! Layer centres             (mm).
ALLOCATE (h   (nsoil_layers_max))      ! Water potentials          (mm).
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
  DO I = 1, nsoil_layers_max + 1
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
! Centres are negative  values below surface.
! Layers start with 1 at surface.
!----------------------------------------------------------------------!
DO I = 1, nsoil_layers_max
  dz (I) = zb (I) - zb (I+1)
  zc (I) = zb (I) - dz (I) / 2.0
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
ALLOCATE (axy_wz      (lon_c,lat_c,NYR))! Mean yearly water table depth.
!----------------------------------------------------------------------!
! Accumulated theta.
!----------------------------------------------------------------------!
ALLOCATE (axy_theta   (nsoil_layers_max,lon_c,lat_c,NYR))
ALLOCATE (axy_theta_total (lon_c,lat_c,NYR))
!----------------------------------------------------------------------!
ALLOCATE (axy_rnf     (lon_c,lat_c,NYR)) ! Accumulated runoff.
ALLOCATE (axy_evap    (lon_c,lat_c,NYR)) ! Accumulated evap.
!----------------------------------------------------------------------!
! Soil water (mm^3 mm^-3).
!----------------------------------------------------------------------!
ALLOCATE (theta       (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
!  Chunk of theta_s of one soil layer read in at
! 30 arc-seconds (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1_in(lon_c*60,lat_c*60))
!----------------------------------------------------------------------!
! Chunk of theta_s of one soil layer gridded to half-degree (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_s_l1   (lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of saturated volumetric soil water (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_s (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!
! Chunk of volumetric soil matric potential (cm3 cm-3).
!----------------------------------------------------------------------!
ALLOCATE (theta_m (nsoil_layers_max,lon_c,lat_c))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise global diagnostic arrays with fill (i.e. NaN) values.
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
axy_wz (:,:,:) = zero / zero
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
! Read in soil saturated water contents for this block.
! Units are cm^3/cm^3. Downloaded to
! 
! from http://globalchange.bnu.edu.cn/research/soil4d.jsp#download
! on 22/7/16.
!----------------------------------------------------------------------!
DO I = 1, nsoil_layers_max ! Loop over soil layers
  IF (I == 1) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l1.nc'
  IF (I == 2) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l2.nc'
  IF (I == 3) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l3.nc'
  IF (I == 4) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l4.nc'
  IF (I == 5) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l5.nc'
  IF (I == 6) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l6.nc'
  IF (I == 7) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l7.nc'
  IF (I == 8) file_in = '/scratch/adf10/DATA/SOILS/theta_s/&
                          theta_s_l8.nc'
  WRITE (*,*) 'Opening ',TRIM(file_in)
  CALL CHECK(NF90_OPEN (TRIM (file_in), NF90_NOWRITE, ncid))
  WRITE (*,*) 'File ',TRIM(file_in),' opened'
  CALL CHECK(NF90_GET_VAR (ncid, 4, theta_s_l1_in, &
                       start = (/(lon_s-1)*60+1, (lat_s-1)*60+1 /), &
                       count = (/lon_c*60, lat_c*60 /)))
  WRITE (*,*) 'Saturated water contents read for layer',I
  CALL CHECK(NF90_CLOSE (ncid))
  !--------------------------------------------------------------------!
  ! Data is at 30 arc-seconds, so grid to half-degree.
  !--------------------------------------------------------------------!
  theta_s_l1 (:,:) = zero
  DO y = 1, lat_c
    DO x = 1, lon_c
      j = 0
      DO x1 = (x-1)*60+1, (x-1)*60+60
        do y1 = (y-1)*60+1,(y-1)*60+60
          IF (theta_s_l1_in (x1,y1) >= zero) THEN
            theta_s_l1 (x,y) =  theta_s_l1 (x,y) + theta_s_l1_in (x1,y1)
            j = j + 1
          END IF
        END DO
      END DO
      IF (j > 0) THEN
        theta_s_l1 (x,y) = theta_s_l1 (x,y) / FLOAT (j)
      END IF
    END DO
  END DO
  !--------------------------------------------------------------------!
  ! Place saturated soil volumetric water contents in block array
  ! for each layer.
  ! Read in as cm^3/cm^3 x 0.001, converted to mm.
  !--------------------------------------------------------------------!
  DO y = 1, lat_c
    DO x = 1, lon_c
      theta_s (I,x,y) = theta_s_l1(x,y) / 1.0E3
    END DO
  END DO
  !--------------------------------------------------------------------!
END DO ! I = 1, nsoil_layers_max ! Loop over soil layers.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Minimum soil volumetric water contents (mm^3 mm^-3).
! Following is hack for now.
!----------------------------------------------------------------------!
DO y = 1, lat_c
  DO x = 1, lon_c
    DO I = 1, 8
      theta_m (I,x,y) = theta_s (I,x,y) / 100.0
      ! Initial water in each soil layer (mm^3 mm^-3).
      !theta (I,x,y) = theta_m (I,x,y)
      theta (I,x,y) = theta_s (I,x,y)
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
!stop
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
          wz_sum   = zero
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
            !----------------------------------------------------------!
            w0 = prec * dt
            DO I = 1, nlayers
              w0 = w0 + theta (I,x,y) * dz (I)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water saturation fraction of layers (0-1).
            !----------------------------------------------------------!
            DO I = 1, nlayers
              IF (theta_s (I,x,y) > trunc) THEN
                wv (I) = theta (I,x,y) / theta_s (I,x,y)
              ELSE
                wv (I) = zero
              END IF
              wv (I) = MIN (one, wv (I))
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Hydraulic conductivities at interfaces between layers
            ! based on Eqn. 7.89 of Oleson et al. (2013) (mm s-1).
            ! Refers to bottom of layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers - 1
              xk (I) = xks * (&
                (0.5 * (theta   (I,x,y) + theta   (I+1,x,y))) / &
                (0.5 * (theta_s (I,x,y) + theta_s (I+1,x,y))) &
                             ) &
                       ** (2.0 * B + 3.0)
              !--------------------------------------------------------!
              ! Impose limits just to be sure.
              !--------------------------------------------------------!
              xk (I) = MIN (xk (I), xks)
              xk (I) = MAX (zero, xk (I))
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water potentials in each layer (mm).
            ! This is Eqn. 7.94 of Oleson et al. (2013) for the matric
            ! potential component and the gravitational potential
            ! added as in Eqn. 7.81. zc is the elevation of the centre
            ! of the soil layer with reference to the surface (mm).
            !----------------------------------------------------------!
            DO I = 1, nlayers
              IF (wv (I) < 0.01) THEN
                h (I) = swp_s * 0.01 ** (-B) + zc (I)
              ELSE
                h (I) = swp_s * wv (I) ** (-B) + zc (I)
              END IF
              h (I) = MAX (-1.0E8, h(I))
            END DO
            !----------------------------------------------------------!
 
            !----------------------------------------------------------!
            ! Soil moisture constraint on evaporative flux (0-1).
            !----------------------------------------------------------!
            beta = MIN (one, wv (1))
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
            evap_max = dz (1) * (theta (1,x,y) - theta_m (1,x,y)) / dt
            !----------------------------------------------------------!
            ! Evaporative flux (mm s-1).
            !----------------------------------------------------------!
            huss (x,y,iT) = qb * rhs (x,y,iT) / 100.0
            !evap = beta * rho3 * cna * (qb - huss (x,y,iT))
            evap = beta * rho3 * cna * (qb - qb * rhs (x,y,iT) / 100.0)
            evap = MAX (zero, evap)
            evap = MIN (evap_max,evap)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Surface runoff (mm s-1).
            !----------------------------------------------------------!
            rnf = (theta (1,x,y) * dz (1) + prec * dt - &
                  theta_s (1,x,y) * dz (1)) / dt
            rnf = MAX (zero, rnf)
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Underground runoff from each layer               (mm s-1).
            ! Hack for now.
            !----------------------------------------------------------!
            slope = 0.05
            DO I = 1, nlayers
              rnff (I) = xk (I) * slope * dz (I) * &
                         theta (I,x,y) / theta_s (I,x,y)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Water fluxes at each layer interface (mm s-1).
            ! Uses Darcy's Law as in Eqn. 7.116 of Oleson et al. (2013),
            ! except not (yet) including the Zeng and Decker (2009)
            ! modification.
            ! Also, to make it positive upwards, denominator is reversed
            ! as z is positive downwards in CLM.
            ! f (I) is the interface at the bottom of layer I, so
            ! f (1) is at the bottom of the first layer.
            !----------------------------------------------------------!
            DO I = 1, nlayers - 1
              f (I) = -xk (I) * (h (I) - h (I+1)) / (zc (I) - zc (I+1))
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Bottom boundary flux (mm s-1).
            !----------------------------------------------------------!
            f (nlayers) = zero
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Prevent over/undersaturation of layers 1 - nlayers.
            ! Uses method of GISS modelE (N.B. f refers to bottom of
            ! layer here rather than top as in modelE).
            !----------------------------------------------------------!
            DO I = nlayers, 2, -1
              !--------------------------------------------------------!
              dtheta = (f (I) - f (I-1) - rnff (I)) / dz (I)
              thetan = theta (I,x,y) + dt * dtheta
              !--------------------------------------------------------!
              ! Compensate oversaturation by increasing flux up.
              !--------------------------------------------------------!
              IF (thetan - theta_s (I,x,y) > trunc) THEN
                f (I-1) = f (I-1) + dz (I) * &
                          (thetan - theta_s (I,x,y) + trunc) / dt
              END IF
              !--------------------------------------------------------!
              ! Compensate undersaturation by decreasing runoff.
              !--------------------------------------------------------!
              IF (thetan - theta_m (I,x,y) < trunc) THEN
                rnff (I) = rnff (I) + &
                  dz (I) * (thetan - theta_m (I,x,y) - trunc) / dt
              END IF
              !--------------------------------------------------------!
              ! See if have to compensate with f.
              !--------------------------------------------------------!
              IF (rnff (I) < zero) THEN
                f (I-1) = f (I-1) + rnff (I)
                rnff (I) = zero
              END IF
              !--------------------------------------------------------!
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Prevent over/undersaturation of first layer.
            ! f (1) is the flux upwards into layer 1 (mm/s).
            ! Uses method of GISS modelE.
            !----------------------------------------------------------!
            dtheta = (f(1) + prec - rnf - rnff (1) - evap) / dz (1)
            thetan = theta (1,x,y) + dt * dtheta
            !----------------------------------------------------------!
            ! Compensate oversaturation by increasing surface runoff.
            !----------------------------------------------------------!
            IF (thetan - theta_s (1,x,y) > trunc) rnf = &
              rnf + dz (1) * (thetan - theta_s (1,x,y) + trunc) / dt
            !----------------------------------------------------------!
            ! Compensate undersaturation by decreasing surface runoff.
            !----------------------------------------------------------!
            IF (thetan - theta_m (1,x,y) < trunc) rnf = &
              rnf + dz (1) * (thetan - theta_m (1,x,y) - trunc) / dt
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Now trying to remove any negative runoff.
            !----------------------------------------------------------!
            I = 1
            !----------------------------------------------------------!
            DO WHILE ((rnf < zero) .AND. (I <= nlayers))
              IF (I > 1) THEN
                !------------------------------------------------------!
                ! This is how much water we can take from layer I
                ! (mm s-1).
                !------------------------------------------------------!
                dflux = f (I) + &
                  dz (I) * (theta (I,x,y) - theta_m (I,x,y)) / dt &
                  - f (I-1) - rnff (I)
                f (I-1) = f (I-1) - rnf
                rnf = rnf + MIN (-rnf, dflux)
                !------------------------------------------------------!
              END IF
              IF (rnff (I) < zero) THEN
                !------------------------------------------------------!
                WRITE (*,*) 'Negative underground runoff',my_id,x,y
                STOP
                !------------------------------------------------------!
              END IF
              drnf = MIN (-rnf, rnff (I))
              rnf = rnf + drnf
              rnff (I) = rnff (I) - drnf
              I = I + 1
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Apply water fluxes (positive upwards) to layers to give
            ! new thetas, volumetric water (mm^3 mm^-3).
            !----------------------------------------------------------!
            dtheta = (f (1) + prec - rnf - rnff (1) - evap) / dz (1)
            theta (1,x,y) = theta (1,x,y) + dt * dtheta
            !----------------------------------------------------------!
            DO I = 2, nlayers - 1
              dtheta = (-f (I-1) + f (I) - rnff (I)) / dz (I)
              theta (I,x,y) = theta (I,x,y) + dt * dtheta
            END DO
            !----------------------------------------------------------!
            dtheta = (-f (nlayers-1) + rnff (nlayers)) / dz (nlayers)
            theta (nlayers,x,y) = theta (nlayers,x,y) + dt * dtheta
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Check no theta is negative.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              IF (theta (I,x,y) < zero) THEN
                WRITE (*,*)
                WRITE (*,*) my_id ,x,y,'theta (I,x,y) < zero'
                WRITE (*,*) 'lon, lat ',lon(x),lat(y)
                WRITE (*,*) i,theta(I,x,y)
                WRITE (*,*) 'f ',f(I-1)*dt/dz(I),f(I)*dt/dz(I)
                WRITE (*,*) 'pr rnf ',pr(x,y,iTIME)*dt/dz(I),rnf
                WRITE (*,*) 'theta_m ',theta_m(I,x,y)
                WRITE (*,*) 'evap beta ',evap*dt/dz(I),beta
                WRITE (*,*) 'rnff drnf dflux ',rnff(I),drnf,dflux
                STOP
              END IF
            END DO
            !----------------------------------------------------------!
          
            !----------------------------------------------------------!
            ! Calculation of water table depth wz (mm) based on 
            ! Abramopoulos et al 1988 and as implemented in modelE. 
            ! Testing for first non-saturated layer from bottom up
            !----------------------------------------------------------!
            DO I = nlayers, 1, -1
               IF((wv(I) < 1) .AND. (wv(I-1)==1)) THEN
                hmat = h(I) - zc(I)
                 WRITE (*,*) 'wv(I) ', wv(I)
                 WRITE (*,*) 'wv(I-1) ',wv(I-1)
                 WRITE (*,*) 'lon, lat ',lon(x),lat(y)
                 WRITE (*,*) 'hmat ,NS ', hmat, NS
                 WRITE (*,*) 'layer ', I
                wz = zb(I) - sqrt( (-2.E0 * hmat * dz(I)) / (f(I-1) / xk(I)) )
                 
               END IF
            
            END DO
            !----------------------------------------------------------!


            !----------------------------------------------------------!
            ! New total water for diagnostics (mm).
            !----------------------------------------------------------!
            w1 = (rnf + evap + SUM (rnff (:))) * dt
            DO I = 1, nlayers
              w1 = w1 + theta (I,x,y) * dz (I)
            END DO
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Check for conservation of water.
            !----------------------------------------------------------!
            IF (ABS (w1 - w0) > 0.001) THEN
              WRITE (*,*)
              WRITE (*,*) 'Water imbalance > 0.001 mm ',w1-w0
              WRITE (*,*) 'my_id x y ',my_id,x,y
              WRITE (*,*) 'lon, lat ',lon(x),lat(y)
              WRITE (*,*) 'theta = ',SUM (theta(1:nlayers,x,y))
              WRITE (*,*) 'rnff  = ',dt*SUM (rnff(1:nlayers))
              WRITE (*,*) 'rnf evap = ',dt*rnf,dt*evap
              WRITE (*,*) 'pr = ',dt*pr(x,y,iTIME)
              WRITE (*,*) 'beta = ',beta
              DO I = 1, nlayers
                WRITE (*,*) I,theta(I,x,y),theta_m(I,x,y)
              END DO
              STOP
            END IF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            ! Diagnostics accumulated within the day.
            !----------------------------------------------------------!
            rnf_sum  = rnf_sum  + rnf
            evap_sum = evap_sum + evap
            !wz_sum = wz_sum + wz
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            END DO ! DO NS = 1, NISURF
            !----------------------------------------------------------!

            !----------------------------------------------------------!
            IF (INTERACTIVE) THEN
              WRITE (*,*) iT,theta (1,x,y),prec*dt*48.0,prec,dt
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
            wz_sum   = wz_sum   + wz   
            !----------------------------------------------------------!
            ! Possibly move following to within NS loop.
            !----------------------------------------------------------!
            DO I = 1, nlayers
              theta_sum (I) = theta_sum (I) + theta (I,x,y)
              theta_sum_total = theta_sum_total + theta (I,x,y) * dz (I)
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
          !axy_wz  (x,y,iY) = wz_sum   / FLOAT (nt * NISURF) ! mm
          !------------------------------------------------------------!
          axy_tas  (x,y,iY) = tas_sum  / FLOAT (nt)          ! K 
          axy_huss (x,y,iY) = huss_sum / FLOAT (nt)          ! kg kg-1
          axy_ps   (x,y,iY) = ps_sum   / FLOAT (nt)          ! Pa
          axy_pr   (x,y,iY) = pr_sum   / FLOAT (nt)          ! kg m-2 s-1
          axy_rhs  (x,y,iY) = rhs_sum  / FLOAT (nt)          ! no units
          axy_wz   (x,y,iY) = wz_sum   / FLOAT (nt)          ! mm
          !------------------------------------------------------------!
          DO I = 1, nlayers
            axy_theta (I,x,y,iY) = theta_sum (I) / FLOAT (nt)
          END DO
          axy_theta_total (x,y,iY) = theta_sum_total / FLOAT (nt)
          !------------------------------------------------------------!
        END DO ! Loop over years in decade
      END IF ! Soiled grid-box?
    END DO ! Loop over lon_c 
  END DO ! Loop over lat_c 
  !--------------------------------------------------------------------!

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

