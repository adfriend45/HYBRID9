!======================================================================!
MODULE CONTROL
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Declare control parameters and variables.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
USE MPI ! Enable access to the Message Passing Interface library of
        ! parallel routines.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
IMPLICIT NONE
!----------------------------------------------------------------------!
SAVE
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER, PARAMETER :: sday = 86400 ! s/day.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Global dimensions (no. half-degree grid-boxes).
!----------------------------------------------------------------------!
INTEGER, PARAMETER :: NX = 720
INTEGER, PARAMETER :: NY = 360
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Longitudes and latitudes at centre of grid-box)               (degree)
!----------------------------------------------------------------------!
REAL :: lon_all (NX)
REAL :: lat_all (NY)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! No. timesteps in input climate files                               (n)
!----------------------------------------------------------------------!
INTEGER :: NTIMES
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: start_two (2)    ! Array position for start of netCDF write.
INTEGER :: count_two (2)    ! Array count for netCDF write.
INTEGER :: start_three (3)  ! Array position for start of netCDF write.
INTEGER :: count_three (3)  ! Array count for netCDF write.
INTEGER :: lon_s,lat_s      ! Array positions for start of chunk.
INTEGER :: lon_c,lat_c      ! Array counts for chunk.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
INTEGER :: varid ! Variable ID.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CHARACTER (LEN = 200) :: file_name ! Generic filename.
CHARACTER (LEN = 200) :: var_name  ! Generic variable name.
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! For timing the simulation using 'CPU_TIME' on processor 0.
! Results go to screen/slurm-*.out in seconds (see later).
!----------------------------------------------------------------------!
REAL :: start(1025),finish(1025)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
CHARACTER (LEN = 200) :: file_in_ts  ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_ks  ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_lm  ! Soil input filename.
CHARACTER (LEN = 200) :: file_in_ps  ! Soil input filename.
CHARACTER (LEN = 200) :: PATH_output ! Path for axy output files.
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
INTEGER :: syr               ! First year in decade/simulation      (CE)
INTEGER :: eyr               ! Last year in decade/simulation       (CE)
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
LOGICAL :: PGF               ! Use PGF climate forcing?        (logical)
LOGICAL :: INTERACTIVE       ! Interactive simulation?         (logical)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
LOGICAL :: LCLIM ! Climate from local climate file? (logical)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Name of local climate file, if used.
!----------------------------------------------------------------------!
CHARACTER (LEN = 200) :: LCLIM_filename
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Name of local soil properties file, if used.
!----------------------------------------------------------------------!
CHARACTER (LEN = 200) :: LSOIL_filename
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Number of years to spin-up if using local climate                 (yr)
!----------------------------------------------------------------------!
INTEGER :: NYR_SPIN_UP
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
! Names of decades.
!----------------------------------------------------------------------!
CHARACTER (LEN =  9) :: decade   (12) &
&      = (/ '1901_1910','1911_1920','1921_1930','1931_1940', &
&           '1941_1950','1951_1960','1961_1970','1971_1980', &
&           '1981_1990','1991_2000','2001_2010','2011_2012'/)
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
! Index of soil layer right above water table                        (-)
!----------------------------------------------------------------------!
INTEGER :: jwt
!----------------------------------------------------------------------!

!----------------------------------------------------------------------!
END MODULE CONTROL
!======================================================================!
