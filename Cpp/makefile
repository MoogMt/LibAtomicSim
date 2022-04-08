# C++ Compilator
CPP= g++

# C++ compilation options
CFLAGS= -std=gnu++11 -O3 -W -Wall -ansi -pedantic #-I.

DEPS= linear_algebra.h complex_numbers.h general_io.h read_xyz.h utils.h frame.h traj.h 
OBJ=  linear_algebra.o complex_numbers.o general_io.o read_xyz.o utils.o frame.o traj.o

# Makes everything
#-------------------------------------
%.o: %.c $(DEPS)
	$(CPP) -o $@ -c $< $(CFLAGS)
lib: $(OBJ)
	$(CPP) -o $@ -c $< $(CFLAGS)
#-------------------------------------

# Make test scripts
#-------------------------------------
# General testing (matrix, file i/o)
OBJtest=linear_algebra.o complex_numbers.o general_io.o utils.o traj.o frame.o test.o 
test: $(OBJtest)
		$(CPP) -o $@ $^ $(CFLAGS)
# File I/O testing
OBJfile= linear_algebra.o complex_numbers.o general_io.o utils.o traj.o frame.o test_file.o
test_file: $(OBJfile)
		$(CPP) -o $@ $^ $(CFLAGS)
#-------------------------------------

# Clean
#-----------------------------------
clean:
	rm -rf *.o
#-----------------------------------




