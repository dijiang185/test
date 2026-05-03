FC = mpif90
opt = -O3 -xHost -ip -assume byterecl -fp-model precise -ftz -no-fma -fc=ifort
lib = -L ~/netcdf/lib -l netcdf -l netcdff
inc = -I ~/netcdf/include
src = $(wildcard *.f90)
obj = $(patsubst %.f90, %.o, $(src))

run: $(obj)
	$(FC) $(obj) $(lib) $(opt) -o run

%.o: 
	$(FC) $(opt) $(inc) -c $*.f90

clean:
	rm -rf *.o *.mod run
base.o: jd.o
mod_io.o: base.o
mpi_routine.o:base.o jd.o mod_io.o
solver.o: base.o jd.o mod_io.o mpi_routine.o
baro_vort.o: base.o mod_io.o solver.o mpi_routine.o


