F90 = mpif90
FFLAGS = -O3
LDFLAGS = -lnetcdff -lnetcdf

SRC = \
	HYBRID9.f90 \
	CONTROL.f90 \
	SHARED.f90 \
	INIT.f90 \
        READ_PGF.f90 \
        HYDROLOGY.f90 \
        GROW.f90 \
	READ_NET_CDF_0D.f90 \
	READ_NET_CDF_1DI.f90 \
	READ_NET_CDF_1DR.f90 \
	READ_NET_CDF_2DI.f90 \
	READ_NET_CDF_3DR.f90 \
	WRITE_NET_CDF_2DI.f90 \
	WRITE_NET_CDF_3DR.f90 \
	WRITE_NET_CDF_3DR_soils.f90

OBJ = $(SRC:.f90=.o)

HYBRID9.exe : $(OBJ)
	$(F90) $(FFLAGS) -o HYBRID9.exe $(OBJ) $(LDFLAGS)

# Main routine.
HYBRID9.o : CONTROL.o SHARED.o INIT.o READ_PGF.o HYDROLOGY.o GROW.o READ_NET_CDF_0D.o READ_NET_CDF_1DI.o READ_NET_CDF_1DR.o READ_NET_CDF_2DI.o READ_NET_CDF_3DR.o WRITE_NET_CDF_2DI.o WRITE_NET_CDF_3DR.o WRITE_NET_CDF_3DR_soils.o HYBRID9.f90
	$(F90) $(FFLAGS) -c HYBRID9.f90

# Subroutines.
INIT.o : CONTROL.o SHARED.o INIT.f90
	$(F90) $(FFLAGS) -c INIT.f90

READ_PGF.o : CONTROL.o SHARED.o READ_PGF.f90
	$(F90) $(FFLAGS) -c READ_PGF.f90

HYDROLOGY.o : HYDROLOGY.f90
	$(F90) $(FFLAGS) -c HYDROLOGY.f90

GROW.o : GROW.f90
	$(F90) $(FFLAGS) -c GROW.f90

READ_NET_CDF_0D.o : SHARED.o READ_NET_CDF_0D.f90
	$(F90) $(FFLAGS) -c READ_NET_CDF_0D.f90

READ_NET_CDF_1DI.o : SHARED.o READ_NET_CDF_1DI.f90
	$(F90) $(FFLAGS) -c READ_NET_CDF_1DI.f90

READ_NET_CDF_1DR.o : SHARED.o READ_NET_CDF_1DR.f90
	$(F90) $(FFLAGS) -c READ_NET_CDF_1DR.f90

READ_NET_CDF_2DI.o : SHARED.o READ_NET_CDF_2DI.f90
	$(F90) $(FFLAGS) -c READ_NET_CDF_2DI.f90

READ_NET_CDF_3DR.o : SHARED.o READ_NET_CDF_3DR.f90
	$(F90) $(FFLAGS) -c READ_NET_CDF_3DR.f90

WRITE_NET_CDF_2DI.o : SHARED.o WRITE_NET_CDF_2DI.f90
	$(F90) $(FFLAGS) -c WRITE_NET_CDF_2DI.f90

WRITE_NET_CDF_3DR.o : SHARED.o WRITE_NET_CDF_3DR.f90
	$(F90) $(FFLAGS) -c WRITE_NET_CDF_3DR.f90

WRITE_NET_CDF_3DR_soils.o : SHARED.o WRITE_NET_CDF_3DR_soils.f90
	$(F90) $(FFLAGS) -c WRITE_NET_CDF_3DR_soils.f90

# Modules
CONTROL.o : CONTROL.f90
	$(F90) $(FFLAGS) -c CONTROL.f90

SHARED.o : SHARED.f90
	$(F90) $(FFLAGS) -c SHARED.f90

