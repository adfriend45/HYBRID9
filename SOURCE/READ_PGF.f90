!======================================================================!
SUBROUTINE READ_PGF
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Read PGF climate forcings for decade.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
!USE NETCDF  ! Enable access to the library of netCDF routines.
USE CONTROl ! Control variables.
USE SHARED  ! Shared variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
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
  ! Read Surface Downwelling Longwave Radiation Flux (W/m^2) for this
  ! decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'rlds_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (rlds(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (rlds)
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  ! Read Surface Downwelling Shortwave Radiation Flux (W/m^2) for this
  ! decade and block.
  !--------------------------------------------------------------------!
  varid = 4
  !--------------------------------------------------------------------!
  file_name = '/scratch/adf10/ISIMIP2/Input_Hist_obs/PGFv2.1/'//&
              &'rsds_pgfv2.1_'//decade(iDEC)//'.nc4'
  WRITE (*,*) my_id,'Reading from ',TRIM(file_name)
  ALLOCATE (rsds(lon_c,lat_c,NTIMES))
  CALL READ_NET_CDF_3DR (rsds)
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

!----------------------------------------------------------------------!
END SUBROUTINE READ_PGF
!======================================================================!
