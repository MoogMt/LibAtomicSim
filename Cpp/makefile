# C++ Compilator
CPP=g++

# C++ compilation options
CFLAGS=-std=gnu++11 -O3 -W -Wall -ansi -pedantic #-I.

DEPS=utils.h frame.h traj.h 
OBJ=utils.o frame.o traj.o

# Makes everything
#-------------------------------------
%.o: %.c $(DEPS)
	$(CPP) -o $@ -c $< $(CFLAGS)
lib: $(OBJ)
	$(CPP) -o $@ -c $< $(CFLAGS)
#-------------------------------------

# Make test scripts
#-------------------------------------
OBJtest = test.o utils.o traj.o frame.o
test: $(OBJtest)
	$(CPP) -o $@ $^ $(CFLAGS)
#-------------------------------------

# Clean
#-----------------------------------
clean:
	rm -rf *.o
#-----------------------------------




