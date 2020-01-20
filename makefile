fcflags = -c -O3 -cpp -ffree-form #compile
flflags = -O3 #link

f90=gfortran

mod = xyz
build = /home/moogmt/PHYSIX_Utils/GPlib/Fortran
src = /home/moogmt/PHYSIX_Utils/

objects = $(mod:%=$(build)/%.o)

build_module:
	$(f90) $(fcflags) -o $(build)/xyz.o $(build)/xyz.f90

test_module: build_module
	$(build_module)
	$(f90) $(flflags) -o $(src)/test.o $(src)/test_module.f90 $(build)/xyz.o

clean:
	rm -f  $(objects) $(mod:%=%.mod)

cclean:
	rm -f  $(objects) $(exec) $(mod:%=%.mod)
