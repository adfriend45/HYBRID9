!======================================================================!
SUBROUTINE READ_NET_CDF_2DI
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read a 2D netCDF integer variable file.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE MPI    ! Enable access to the Message Passing Interface library of
           ! parallel routines.
USE NETCDF ! Enable access to the library of netCDF routines.
! Control parameters and variables.
USE CONTROL, ONLY : file_name,varid,lon_s,lat_s,lon_c,lat_c
USE SHARED ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: ncid ! netCDF ID.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Open existing netCDF dataset for access.
! NF90_NOWRITE specifies open with read-only access.
!----------------------------------------------------------------------!
CALL CHECK (NF90_OPEN (TRIM(file_name), NF90_NOWRITE, ncid))
!----------------------------------------------------------------------!
CALL CHECK (NF90_GET_VAR (ncid, varid, data_in_2DI, &
                     start = (/lon_s, lat_s /), &
                     count = (/lon_c, lat_c /)))
CALL CHECK (NF90_CLOSE (ncid))
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
END SUBROUTINE READ_NET_CDF_2DI
!======================================================================!

