# C++ Compilator
CPP=g++

# C++ compilation options
CFLAGS=-std=gnu++11 -O3 -W -Wall -ansi -pedantic #-I.

DEPS=complex_numbers.h general_io.h read_xyz.h utils.h frame.h traj.h g
OBJ=complex_numbers.o general_io.o read_xyz.o utils.o frame.o traj.o

# Makes everything
#-------------------------------------
%.o: %.c $(DEPS)
	$(CPP) -o $@ -c $< $(CFLAGS)
lib: $(OBJ)
	$(CPP) -o $@ -c $< $(CFLAGS)
#-------------------------------------

# Make test scripts
#-------------------------------------
OBJtest = complex_numbers.o general_io.o test.o utils.o traj.o frame.o
test: $(OBJtest)
	 $(CPP) -o $@ $^ $(CFLAGS)
#-------------------------------------

# Clean
#-----------------------------------
clean:
	rm -rf *.o
#-----------------------------------




