!======================================================================!
SUBROUTINE WRITE_NET_CDF_3DR (iyr)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write a 3D real netCDF file.
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
INTEGER, INTENT (IN) :: iyr
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: ncid                    ! netCDF ID.
!INTEGER :: x_dimid,y_dimid,z_dimid ! Dimension IDs.
INTEGER :: lon_dimid,lat_dimid,z_dimid ! Dimension IDs.
INTEGER :: lon_varid
INTEGER :: lat_varid
INTEGER :: z_varid
INTEGER :: rnf_varid
INTEGER :: evap_varid
INTEGER :: tas_varid
INTEGER :: huss_varid
INTEGER :: ps_varid
INTEGER :: pr_varid
INTEGER :: rhs_varid
INTEGER :: theta_varid
INTEGER :: thetaz_varid
INTEGER :: dimids_two   (2) ! Array of two dimids.
INTEGER :: dimids_three (3) ! Array of three dimids.
INTEGER :: wz_varid
REAL :: fillval
!----------------------------------------------------------------------!
CHARACTER (LEN = *), PARAMETER :: LON_NAME    = "longitude"
CHARACTER (LEN = *), PARAMETER :: LAT_NAME    = "latitude"
CHARACTER (LEN = *), PARAMETER :: DEPTH_NAME  = "layer_centre_depth"
CHARACTER (LEN = *), PARAMETER :: RNF_NAME    = "runoff"
CHARACTER (LEN = *), PARAMETER :: EVAP_NAME   = "evaporation"
CHARACTER (LEN = *), PARAMETER :: TAS_NAME    = "temperature"
CHARACTER (LEN = *), PARAMETER :: HUSS_NAME   = "specific_humidity"
CHARACTER (LEN = *), PARAMETER :: PS_NAME     = "air_pressure"
CHARACTER (LEN = *), PARAMETER :: PR_NAME     = "precipitation"
CHARACTER (LEN = *), PARAMETER :: RHS_NAME    = "relative_humidity"
CHARACTER (LEN = *), PARAMETER :: THETA_NAME  = "soil_water"
CHARACTER (LEN = *), PARAMETER :: THETAS_NAME = "soil_water_layers"
CHARACTER (LEN = *), PARAMETER :: WZ_NAME     = "water_table_depth"
!----------------------------------------------------------------------!
CHARACTER (LEN = *), PARAMETER :: LON_UNITS        = "degrees_east"
CHARACTER (LEN = *), PARAMETER :: LAT_UNITS        = "degrees_north"
CHARACTER (LEN = *), PARAMETER :: RNF_UNITS        = "mm/s"
CHARACTER (LEN = *), PARAMETER :: EVAP_UNITS       = "mm/s"
CHARACTER (LEN = *), PARAMETER :: TAS_UNITS        = "K"
CHARACTER (LEN = *), PARAMETER :: HUSS_UNITS       = "kg[water]/kg[air]"
CHARACTER (LEN = *), PARAMETER :: PS_UNITS         = "Pa"
CHARACTER (LEN = *), PARAMETER :: PR_UNITS         = "kg/m^2/s"
CHARACTER (LEN = *), PARAMETER :: RHS_UNITS        = "percent"
CHARACTER (LEN = *), PARAMETER :: THETA_TOTAL_UNITS = "mm"
CHARACTER (LEN = *), PARAMETER :: THETA_UNITS       = "mm^3/mm^3"
CHARACTER (LEN = *), PARAMETER :: DEPTH_UNITS       = "mm"
!----------------------------------------------------------------------!
CHARACTER (LEN = *), PARAMETER :: UNITS       = "units"
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Create a netCDF file. The 'CHECK' wrapping checks the
! return code to make sure any return not equal to nf90noerr (0) will
! print a netCDF error message and exit.
! nf90_create creates a new netCDF dataset, returing a netCDF ID.
! The new dataset is opened for write access and placed in define mode.
! The NF90_NETCDF4 flag causes a HDF5/netCDF-4 file to be created.
! The NF90_MPIIO flag selects MPI/IO (rather than MPI/POSIX).
! The comm and info parameters cause parallel I/O to be enabled.
!----------------------------------------------------------------------!
CALL CHECK (nf90_create(TRIM(file_name),IOR(NF90_NETCDF4, NF90_MPIIO), &
  ncid, comm = MPI_COMM_WORLD, info = MPI_INFO_NULL))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the dimensions. netCDF hands back an ID for each.
! Metadata operations must take place on all processors.
!----------------------------------------------------------------------!
CALL CHECK (nf90_def_dim(ncid, LAT_NAME, NY, lat_dimid))
CALL CHECK (nf90_def_dim(ncid, LON_NAME, NX, lon_dimid))
CALL CHECK (nf90_def_dim(ncid, DEPTH_NAME, nsoil_layers_max, z_dimid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the coordinate variables.
!----------------------------------------------------------------------!
CALL CHECK (NF90_DEF_VAR(ncid, LAT_NAME, NF90_FLOAT, &
            lat_dimid, lat_varid))
CALL CHECK (NF90_DEF_VAR(ncid, LON_NAME, NF90_FLOAT, &
            lon_dimid, lon_varid))
CALL CHECK (NF90_DEF_VAR(ncid, DEPTH_NAME, NF90_FLOAT, &
            z_dimid, z_varid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Assign units attributes to coordinate var data.
!----------------------------------------------------------------------!
CALL CHECK (NF90_PUT_ATT(ncid, lat_varid, UNITS, LAT_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, lon_varid, UNITS, LON_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, z_varid  , UNITS, DEPTH_UNITS))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! The dimids array is used to pass the IDs of the dimensions of the
! variables. Note that in Fortran arrays are stored in column-major
! format (i.e. leftmost indices vary fastest).
!----------------------------------------------------------------------!
dimids_two   = (/lon_dimid,lat_dimid/)
dimids_three = (/z_dimid,lon_dimid,lat_dimid/)
!dimids_two (1) = lon_dimid
!dimids_two (2) = lat_dimid
!dimids_three (1) = z_dimid
!dimids_three (2) = lon_dimid
!dimids_three (3) = lat_dimid
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the data variables.
!----------------------------------------------------------------------!
CALL CHECK (nf90_def_var(ncid, RNF_NAME  , NF90_float, &
  dimids_two, rnf_varid))
CALL CHECK (nf90_def_var(ncid, EVAP_NAME , NF90_float, &
  dimids_two, evap_varid))
CALL CHECK (nf90_def_var(ncid, TAS_NAME  , NF90_float, &
  dimids_two, tas_varid))
CALL CHECK (nf90_def_var(ncid, HUSS_NAME , NF90_float, &
  dimids_two, huss_varid))
CALL CHECK (nf90_def_var(ncid, PS_NAME   , NF90_float, &
  dimids_two, ps_varid))
CALL CHECK (nf90_def_var(ncid, PR_NAME   , NF90_float, &
  dimids_two, pr_varid))
CALL CHECK (nf90_def_var(ncid, RHS_NAME  , NF90_float, &
  dimids_two, rhs_varid))
CALL CHECK (nf90_def_var(ncid, THETA_NAME, NF90_float, &
  dimids_two, theta_varid))
CALL CHECK (nf90_def_var(ncid, THETAS_NAME, NF90_float, &
  dimids_three, thetas_varid ))
CALL CHECK (nf90_def_var(ncid, WZ_NAME, NF90_float, &
 dimids_two, wz_varid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Assign units attributed to variables.
!----------------------------------------------------------------------!
CALL CHECK (NF90_PUT_ATT(ncid, rnf_varid, UNITS, RNF_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, evap_varid,UNITS, EVAP_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, tas_varid, UNITS, TAS_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, huss_varid,UNITS, HUSS_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, ps_varid,  UNITS, PS_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, pr_varid,  UNITS, PR_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, rhs_varid, UNITS, RHS_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, theta_varid, &
                         UNITS,THETA_TOTAL_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, thetaz_varid,UNITS, THETA_UNITS))
CALL CHECK (NF90_PUT_ATT(ncid, wz_varid, UNITS,DEPTH_UNITS))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the fill values to be used where there are no data.
!----------------------------------------------------------------------!
fillval = zero / zero
CALL CHECK (nf90_PUT_ATT(ncid, rnf_varid   , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, evap_varid  , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, tas_varid   , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, huss_varid  , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, ps_varid    , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, pr_varid    , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, rhs_varid   , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, theta_varid , "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, thetaz_varid, "_FillValue", fillval))
CALL CHECK (nf90_PUT_ATT(ncid, wz_varid    , "_FillValue", fillval))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! End define mode. This tells netCDF we are done defining metadata.
! This operation is collective and all processors will write their
! metadata to disk.
!----------------------------------------------------------------------!
CALL CHECK (nf90_enddef(ncid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Set up the position and dimensions of the chunk to write.
!----------------------------------------------------------------------!
start_two = (/lon_s,lat_s/)
count_two = (/lon_c,lat_c/)
start_three = (/1,lon_s,lat_s/)
count_three = (/nsoil_layers_max,lon_c,lat_c/)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write the coordinate variable data.
!----------------------------------------------------------------------!
!IF (my_id == 0) THEN
  CALL CHECK (NF90_PUT_VAR(ncid, lat_varid, lat_all))
  CALL CHECK (NF90_PUT_VAR(ncid, lon_varid, lon_all))
  ! Following does not work. Not sure why.
  !CALL CHECK (NF90_PUT_VAR(ncid, z_varid, zc))
  !CALL CHECK (NF90_PUT_VAR(ncid, lat_varid, lat, &
  !  start = (/lat_s/), count = (/lat_c/)))
  !CALL CHECK (NF90_PUT_VAR(ncid, lon_varid, lon, &
  !  start = (/lon_s/), count = (/lon_c/)))
  CALL CHECK (NF90_PUT_VAR(ncid, z_varid, zc))
!END IF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write the data to the netCDF file. Each processor writes one chunk.
!----------------------------------------------------------------------!
CALL CHECK(nf90_put_var(ncid, rnf_varid  , axy_rnf  (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, evap_varid , axy_evap (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, tas_varid  , axy_tas  (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, huss_varid , axy_huss (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, ps_varid   , axy_ps   (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, pr_varid   , axy_pr   (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, rhs_varid  , axy_rhs  (:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, theta_varid, axy_theta_total(:,:,iyr), &
  start = start_two, count = count_two))
CALL CHECK(nf90_put_var(ncid, thetaz_varid, axy_theta (:,:,:,iyr), &
  start = start_three, count = count_three))
CALL CHECK(nf90_put_var(ncid, wz_varid   , axy_wz   (:,:,iyr), &
  start = start_two, count = count_two)) 
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Close the netCDF file.
!----------------------------------------------------------------------!
CALL CHECK (nf90_close(ncid))
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
END SUBROUTINE WRITE_NET_CDF_3DR
!======================================================================!

