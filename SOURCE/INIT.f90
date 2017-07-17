!======================================================================!
SUBROUTINE INIT
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Initialise simulation.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE NETCDF  ! Enable access to the library of netCDF routines.
USE CONTROl ! Control variables.
USE SHARED  ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
REAL :: decay ! For rooting.
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
! Specific leaf area for each GPT                           (m^2/kg[DM])
!----------------------------------------------------------------------!
ALLOCATE (sla (nGPTs))
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
! Volumetric soil water content in macropores                (mm^3/mm^3)
!----------------------------------------------------------------------!
ALLOCATE (theta_ma (nsoil_layers_max))
!----------------------------------------------------------------------!
ALLOCATE (S (nsoil_layers_max))
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
! Set GPT parameters.
!----------------------------------------------------------------------!
! Following for beech from Bouriaud et al., 2003, CJRS 29, 371-380.
!----------------------------------------------------------------------!
sla (1) = 23.0E-3 ! Specific leaf area (m^2/g[DM])
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
  READ (10,*) PATH_output
  READ (10,*) NISURF
  READ (10,*) PGF
  READ (10,*) iDEC_start
  READ (10,*) iDEC_end
  READ (10,*) INTERACTIVE
  READ (10,*) LCLIM
  READ (10,*) LCLIM_filename
  READ (10,*) LSOIL_filename
  READ (10,*) syr
  READ (10,*) eyr
  READ (10,*) NYR_SPIN_UP
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
! Set sub-daily integration timestep                                 (s)
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
  WRITE (*,*) 'Details of focus point',lon_w,lat_w
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
ALLOCATE (plant_mass   (nplants_max,lon_c,lat_c)) ! Plant mass   (g[DM])
! Plant foliage mass (g[DM]).
ALLOCATE (plant_foliage_mass (nplants_max,lon_c,lat_c))
ALLOCATE (plant_length (nplants_max,lon_c,lat_c)) ! Plant length    (mm)
ALLOCATE (rdepth       (nplants_max,lon_c,lat_c)) ! Rooting depth   (mm)
ALLOCATE (C_labile     (nplants_max,lon_c,lat_c)) ! Labile C      (g[C])
ALLOCATE (N_labile     (nplants_max,lon_c,lat_c)) ! Labile N      (g[N])
ALLOCATE (P_labile     (nplants_max,lon_c,lat_c)) ! Labile P      (g[P])
!----------------------------------------------------------------------!
ALLOCATE (nplants    (lon_c,lat_c)) ! Number plants per grid-box     (-)
ALLOCATE (LAI        (lon_c,lat_c)) ! Leaf area index          (m^2/m^2)
ALLOCATE (LAI_litter (lon_c,lat_c)) ! LAI of litter            (m^2/m^2)
!----------------------------------------------------------------------!
! Effective root fraction in each soil layer                         (-)
!----------------------------------------------------------------------!
ALLOCATE (rootr_col (1:Nlevgrnd,lon_c,lat_c))
!----------------------------------------------------------------------!
ALLOCATE (data_in_2DI (lon_c,lat_c))  ! Generic 2D input integer array.
ALLOCATE (soil_tex    (lon_c,lat_c))  ! HWSD soil textures           (-)
ALLOCATE (topo_slope  (lon_c,lat_c))  ! Gridcell slope             (deg)
ALLOCATE (Fmax        (lon_c,lat_c))  ! Max. sat. fraction    (fraction)
ALLOCATE (lon         (lon_c))        ! Longitudes (degrees east)
ALLOCATE (lat         (lat_c))        ! Latitudes (degrees north)
ALLOCATE (block_sub   (lon_c,lat_c))  ! Processor IDs.
ALLOCATE (axy_npp     (lon_c,lat_c,NYR)) ! Total ann. NPP (g[DM]/m^2/yr)
ALLOCATE (axy_plant_mass (lon_c,lat_c,NYR)) ! Mean ann. pl. mass (g[DM])
ALLOCATE (axy_tas     (lon_c,lat_c,NYR)) ! Mean annual tas           (K)
ALLOCATE (axy_rlds    (lon_c,lat_c,NYR)) ! Mean annual rlds      (W/m^2)
ALLOCATE (axy_rsds    (lon_c,lat_c,NYR)) ! Mean annual rsds      (W/m^2)
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
ALLOCATE (h2osoi_liq    (nsoil_layers_max,lon_c,lat_c)) ! Micropores.
ALLOCATE (h2osoi_liq_ma (nsoil_layers_max,lon_c,lat_c)) ! Macropores.
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
ALLOCATE (theta_s    (nsoil_layers_max,lon_c,lat_c)) ! Micropores
ALLOCATE (theta_ma_s (nsoil_layers_max,lon_c,lat_c)) ! Macropores.
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
! Initialise global diagnostic arrays with fill values (i.e. NaN and
! zero).
!----------------------------------------------------------------------!
axy_npp        (:,:,:) = zero / zero
axy_plant_mass (:,:,:) = zero / zero
axy_rnf   (:,:,:) = zero / zero
axy_evap  (:,:,:) = zero / zero
axy_tas   (:,:,:) = zero / zero
axy_rlds  (:,:,:) = zero / zero
axy_rsds  (:,:,:) = zero / zero
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
      ! Set macropore thetas from BG81.
      !----------------------------------------------------------------!
      theta_ma_s (I,x,y) = 0.1 ! Macroporosity.
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
file_name = '/scratch/adf10/DATA/SOILS/Fmax_half.nc'
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
theta_m       (:,:,:) = zero
h2osoi_liq    (:,:,:) = zero
h2osoi_liq_ma (:,:,:) = zero
plant_mass    (:,:,:) = zero
C_labile      (:,:,:) = zero
N_labile      (:,:,:) = zero
P_labile      (:,:,:) = zero
zwt           (:,:) = zero
wa            (:,:) = zero
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
        h2osoi_liq (I,x,y) = 0.4 * theta_s (I,x,y) * dz (I) * rhow &
                             / 1000.0
        h2osoi_liq_ma (I,x,y) = 0.4 * theta_ma_s (I,x,y) * dz (I) * &
                                rhow / 1000.0
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
      ! Initial LAI of litter layer                            (m^2/m^2)
      !----------------------------------------------------------------!
      LAI_litter (x,y) = 0.001
      !----------------------------------------------------------------!
      ! Initial number of plants in grid-box                         (-)
      !----------------------------------------------------------------!
      nplants (x,y) = 1
      !----------------------------------------------------------------!
      ! Initialise LAI of canopy                               (m^2/m^2)
      !----------------------------------------------------------------!
      LAI (x,y) = zero
      !----------------------------------------------------------------!
      ! Initialise plot relative root distribution over soil layers  (-)
      !----------------------------------------------------------------!
      rootr_col (:,x,y) = zero
      !----------------------------------------------------------------!
      DO K = 1, nplants (x,y)
        !--------------------------------------------------------------!
        ! Generalised plant type index                               (-)
        !--------------------------------------------------------------!
        iGPT = 1
        !--------------------------------------------------------------!
        ! Initial plant masses                                   (g[DM])
        !--------------------------------------------------------------!
        plant_mass         (K,x,y) = 1.0
        plant_foliage_mass (K,x,y) = 0.0435 ! LAI = 0.001.
        !--------------------------------------------------------------!
        ! Assume all is a cylinder of mass with length, and
        ! with some portion below-ground                            (mm)
        !--------------------------------------------------------------!
        plant_length (K,x,y) = (400.0 * plant_mass (K,x,y) / &
                               3.142E-3) ** (one / 3.0)
        !--------------------------------------------------------------!
        ! Sum foliage area over ind's for total LAI of canopy  (m^2/m^2)
        !--------------------------------------------------------------!
        LAI (x,y) = LAI (x,y) + &
                    plant_foliage_mass (K,x,y) * sla (iGPT) / plot_area
        !--------------------------------------------------------------!
        ! Rooting depth                                             (mm)
        !--------------------------------------------------------------!
        rdepth (K,x,y) = 0.3 * plant_length (K,x,y)
        !--------------------------------------------------------------!
        ! Rooting distribution (based on Baldocchi et al., 2004).
        ! Assumes 90% of roots are within rdepth.
        !--------------------------------------------------------------!
        decay = EXP (LOG (0.1) / (rdepth (K,x,y) / 10.0))
        !--------------------------------------------------------------!
        DO I = 1, nlayers
          rootr_col (I,x,y) = rootr_col (I,x,y) + &
                      (1.0 - decay ** (zi (I) / 10.0)) - &
                      (1.0 - decay ** (zi (I-1) / 10.0))
        END DO
        !--------------------------------------------------------------!
        ! Initial plant labile contents                              (g)
        ! Values as for leaf in Figure 4 of Bell et al., 2013
        ! (10.1111/nph.12531), eye-balled means.
        !--------------------------------------------------------------!
        C_labile (K,x,y) = plant_mass (K,x,y) * 0.5 * 0.1
        N_labile (K,x,y) = C_labile (K,x,y) * 0.035 ! N:C = 3.5%
        P_labile (K,x,y) = N_labile (K,x,y) * 0.025 ! P:N = 2.5%
        !--------------------------------------------------------------!
      END DO
      !----------------------------------------------------------------!
    END IF ! Soiled with some water capacity.
  END DO ! x = 1, lon_c
END DO ! y = 1, lat_c
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
!----------------------------------------------------------------------!
! Open file for daily diagnostics.
!----------------------------------------------------------------------!
OPEN (20,FILE='LCLIM/daily_diag.csv',STATUS='UNKNOWN')
WRITE (20,'(A200)') 'jyear,  DOY,  et,  et_g,  theta_1,  theta_2,  &
              &theta_3,  theta_4,&
              &  theta_ma_1,  lai,  lai_litter, w_i, fT'
!----------------------------------------------------------------------!
END IF

!----------------------------------------------------------------------!
IF (.NOT. (PGF)) THEN
  !--------------------------------------------------------------------!
  x = 1
  y = 1
  nlayers = nsoil_layers_max
  NTIMES = time_BOY (eyr+1-1859) - time_BOY (syr-1859)
  ALLOCATE (pr  (lon_c,lat_c,NTIMES))
  ALLOCATE (tas (lon_c,lat_c,NTIMES))
  ALLOCATE (rlds(lon_c,lat_c,NTIMES))
  ALLOCATE (rsds(lon_c,lat_c,NTIMES))
  ALLOCATE (huss(lon_c,lat_c,NTIMES))
  ALLOCATE (ps  (lon_c,lat_c,NTIMES))
  ALLOCATE (rhs (lon_c,lat_c,NTIMES))
  OPEN (10,FILE=TRIM(LSOIL_filename),STATUS='OLD')
  READ (10,*)
  DO I = 1, nlayers
    READ (10,*) J,hksat(I,x,y),lambda(I,x,y),psi_s(I,x,y),theta_s(I,x,y)
    bsw (I,x,y) = 1.0E3 / lambda (I,x,y)
!************************
! Vaira (very rocky silt loam)
! work up best esimates
!bsw (I,x,y) = 3.3
!bsw (I,x,y) = 6.0
!bsw (I,x,y) = 2.013! Baldocchi et al. (2004).
!lambda (I,x,y) = 1.0E3 / bsw (I,x,y)
!psi_s (I,x,y) = -0.00786 * 10000.0
!write (*,*) hksat(I,x,y)
!hksat (I,x,y) = 2.59*10.0E3/3600.0
!psi_s (I,x,y) = psi_s (I,x,y)
!theta_s (I,x,y) = theta_s (I,x,y) / 1.4
!hksat (I,x,y) = hksat (I,x,y) * 0.1
!theta_s (I,x,y) = 300.0
!************************
  END DO
!theta_s (1,x,y) = 290.0
!theta_s (2,x,y) = 300.0
!theta_s (3,x,y) = 300.0
!theta_s (4,x,y) = 350.0
  CLOSE (10)
  hksat   (:,:,:) = hksat   (:,:,:) * 10.0 / 86400.0
  lambda  (:,:,:) = lambda  (:,:,:) / 1.0E3
  psi_s   (:,:,:) = psi_s   (:,:,:) * 10.0
  theta_s (:,:,:) = theta_s (:,:,:) / 1.0E3
  DO I = 1, nlayers
    theta_m (I,x,y) = theta_s (I,x,y) * ((-3.1E9 / (1000.0 * 9.8)) &
                      / psi_s (I,x,y)) ** (-lambda (I,x,y))
!***************
!write (*,*) psi_s(I,x,y)
!psi_s (I,x,y) = -1780.0
!if (zi (i) > 1000.0) hksat (I,x,y) = hksat (I,x,y) / 100000.0
!if (zi (i) > 1000.0) hksat (I,x,y) = 0.0
!hksat (6,x,y) = hksat (6,x,y) / 1000.0 ! rock layer
!***************
  END DO
  LAI = 0.001 !*******************
  LAI_litter = 0.0 !**********************
  !--------------------------------------------------------------------!

END IF
!----------------------------------------------------------------------!

!======================================================================!
CONTAINS
  SUBROUTINE CHECK (status)
    INTEGER, INTENT (IN) :: status

    IF (status /= nf90_noerr) THEN
      PRINT *, TRIM (nf90_strerror(status))
      STOP "Stopped"
    END IF
  END SUBROUTINE CHECK
!======================================================================!

!----------------------------------------------------------------------!
END SUBROUTINE INIT
!======================================================================!
