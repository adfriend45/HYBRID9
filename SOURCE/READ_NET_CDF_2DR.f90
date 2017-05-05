!======================================================================!
SUBROUTINE READ_NET_CDF_2DR (data_in,var_id)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read from a 2D netCDF file.
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
INTEGER, INTENT ( IN) :: var_id
REAL, INTENT (OUT) :: data_in (lon_c,lat_c)

!----------------------------------------------------------------------!
INTEGER :: ncid ! netCDF ID.
INTEGER :: nDims
INTEGER :: nVars
INTEGER :: nGlobalAtts
INTEGER :: unlimdimid
INTEGER :: RecordDimId
INTEGER :: len
INTEGER :: Records
INTEGER :: xtype
INTEGER :: attnum
REAL :: r_value
CHARACTER (LEN = NF90_MAX_NAME) :: RecordDimName,c_value,name
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Open existing netCDF dataset for access.
! NF90_NOWRITE specifies open with read-only access.
!----------------------------------------------------------------------!
CALL CHECK (NF90_OPEN (TRIM(file_name), NF90_NOWRITE, ncid))
!----------------------------------------------------------------------!
CALL CHECK (NF90_INQUIRE (ncid,nDims,nVars,nGlobalAtts,unlimdimid))
!write (*,*) ncid,nDims,nVars,nGlobalAtts,unlimdimid
!----------------------------------------------------------------------!
DO RecordDimID = 1, nDims
  CALL CHECK (NF90_INQUIRE_DIMENSION (ncid,RecordDimID, &
              name = RecordDimName, len = Records))
  IF (RecordDimName == 'time') NTIMES = Records
END DO
!----------------------------------------------------------------------!
! Have a look at the attributes of the time variable.
! varid = 1 is longitude (degrees_east)
! varid = 2 is latitude (degrees_north)
! varid = 3 is time (days since 1860-01-01 00:00:00)
! varid = 4 is main variable
!----------------------------------------------------------------------!
!varid = 3
!----------------------------------------------------------------------!
!DO attnum = 1, nGlobalAtts-5
!  CALL CHECK (NF90_INQ_ATTNAME (ncid,varid,attnum,name))
!  WRITE (*,*) TRIM(name)
!  CALL CHECK (NF90_INQUIRE_ATTRIBUTE (ncid,varid,name,xtype,len))
!  WRITE (*,*) xtype,len,attnum
!  IF (xtype == 2) THEN
!    CALL CHECK (NF90_GET_ATT (ncid,varid,name,c_value))
!    WRITE (*,*) TRIM(name),' ',TRIM(c_value)
!  ENDIF
!  IF (xtype == 5) THEN
!    CALL CHECK (NF90_GET_ATT (ncid,varid,name,r_value))
!    WRITE (*,*) TRIM(name),' ',r_value
!  ENDIF
!END DO
!----------------------------------------------------------------------!
! Have a look at the attributes of the main variable.
!varid = 4
!----------------------------------------------------------------------!
!DO attnum = 1, nGlobalAtts-5
!  CALL CHECK (NF90_INQ_ATTNAME (ncid,varid,attnum,name))
!  WRITE (*,*) TRIM(name)
!  CALL CHECK (NF90_INQUIRE_ATTRIBUTE (ncid,varid,name,xtype,len))
!  WRITE (*,*) xtype,len,attnum
!  IF (xtype == 2) THEN
!    CALL CHECK (NF90_GET_ATT (ncid,varid,name,c_value))
!    WRITE (*,*) TRIM(name),' ',TRIM(c_value)
!  ENDIF
!  IF (xtype == 5) THEN
!    CALL CHECK (NF90_GET_ATT (ncid,varid,name,r_value))
!    WRITE (*,*) TRIM(name),' ',r_value
!  ENDIF
!END DO
!----------------------------------------------------------------------!
CALL CHECK (NF90_GET_VAR (ncid, var_id, data_in, &
                     start = (/lon_s, lat_s/), &
                     count = (/lon_c, lat_c/)))
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
END SUBROUTINE READ_NET_CDF_2DR
!======================================================================!

