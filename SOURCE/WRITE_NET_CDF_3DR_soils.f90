!======================================================================!
SUBROUTINE WRITE_NET_CDF_3DR_soils
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write a 3D real netCDF file.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE MPI    ! Enable access to the Message Passing Interface library of
           ! parallel routines.
USE NETCDF ! Enable access to the library of netCDF routines.
! Control parameters and variables.
USE CONTROL, ONLY : file_name,NY,NX,lon_s,lat_s,start_three, &
                    count_three,lon_c,lat_c,lat_all,lon_all
USE SHARED ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: ncid                    ! netCDF ID.
INTEGER :: lon_dimid,lat_dimid,z_dimid ! Dimension IDs.
INTEGER :: lon_varid
INTEGER :: lat_varid
INTEGER :: z_varid
INTEGER :: psi_s_varid
INTEGER :: dimids_three (3) ! Array of three dimids.
REAL :: fillval
!----------------------------------------------------------------------!
CHARACTER (LEN = *), PARAMETER :: LON_NAME    = "longitude"
CHARACTER (LEN = *), PARAMETER :: LAT_NAME    = "latitude"
CHARACTER (LEN = *), PARAMETER :: DEPTH_NAME  = "layer_centre_depth"
CHARACTER (LEN = *), PARAMETER :: PSI_S_NAME  = "saturated_psi"
!----------------------------------------------------------------------!
CHARACTER (LEN = *), PARAMETER :: LON_UNITS        = "degrees_east"
CHARACTER (LEN = *), PARAMETER :: LAT_UNITS        = "degrees_north"
CHARACTER (LEN = *), PARAMETER :: PSI_S_UNITS      = "mm"
CHARACTER (LEN = *), PARAMETER :: DEPTH_UNITS      = "mm"
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
dimids_three = (/z_dimid,lon_dimid,lat_dimid/)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the data variables.
!----------------------------------------------------------------------!
CALL CHECK (nf90_def_var(ncid, PSI_S_NAME, NF90_float, &
  dimids_three, psi_s_varid))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Assign units attributed to variables.
!----------------------------------------------------------------------!
CALL CHECK (NF90_PUT_ATT(ncid, psi_s_varid, UNITS, PSI_S_UNITS))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Define the fill values to be used where there are no data.
!----------------------------------------------------------------------!
fillval = zero / zero
CALL CHECK (nf90_PUT_ATT(ncid, psi_s_varid   , "_FillValue", fillval))
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
start_three = (/1,lon_s,lat_s/)
count_three = (/nsoil_layers_max,lon_c,lat_c/)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write the coordinate variable data.
!----------------------------------------------------------------------!
CALL CHECK (NF90_PUT_VAR(ncid, lat_varid, lat_all))
CALL CHECK (NF90_PUT_VAR(ncid, lon_varid, lon_all))
CALL CHECK (NF90_PUT_VAR(ncid, z_varid, zc_o))
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Write the data to the netCDF file. Each processor writes one chunk.
!----------------------------------------------------------------------!
CALL CHECK(nf90_put_var(ncid, psi_s_varid, psi_s (:,:,:), &
  start = start_three, count = count_three))
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
END SUBROUTINE WRITE_NET_CDF_3DR_soils
!======================================================================!

