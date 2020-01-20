#------------
# Makefile
#------------

#-------
CPP=g++
#-------

#-------------------------------------
CFLAGS=-std=gnu++11 -O3 -W -Wall -ansi -pedantic -I.
#-------------------------------------

#-------------------------------------
DEPS = utils.h atom.h cell.h cutoff.h contact_matrix.h xyz.h  histogram.h sim.h pdb.h molecules.h lut.h bonds.h
#-------------------------------------

#-------------------------------------
OBJP = pressure.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJNN = NN.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJmols = make_molecules.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJangles = angles.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJhistratio = ratio_counting.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJdiff = diffusion.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJwrap = wrap.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJcomp = compress.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJmolid = mol_id.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJcHist = computeHist.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJhistangles = hist_angles.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJcoord = softCoord.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJcorr = correlation.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJnormalize = normalize.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJdune = dynamics.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
OBJcolvar = readcolvar.o cutoff.o contact_matrix.o cell.o atom.o xyz.o utils.o histogram.o sim.o pdb.o molecules.o lut.o bonds.o
#-------------------------------------

#-------------------------------------
%.o: %.c $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)
#-------------------------------------

#-------------------------------------
pressure: $(OBJP)
	$(CPP) -o $@ $^ $(CFLAGS)
NN: $(OBJNN)
	$(CPP) -o $@ $^ $(CFLAGS)
molecule: $(OBJmols)
	$(CPP) -o $@ $^ $(CFLAGS)	
angles: $(OBJangles)
	$(CPP) -o $@ $^ $(CFLAGS)
ratio: $(OBJhistratio)
	$(CPP) -o $@ $^ $(CFLAGS)
diffusion: $(OBJdiff)
	$(CPP) -o $@ $^ $(CFLAGS)
wrap: $(OBJdyn)
	$(CPP) -o $@ $^ $(CFLAGS)
dune: $(OBJdune)
	$(CPP) -o $@ $^ $(CFLAGS)
compress: $(OBJcomp)
	$(CPP) -o $@ $^ $(CFLAGS)
molid: $(OBJmolid)
	$(CPP) -o $@ $^ $(CFLAGS)
hist_p: $(OBJcHist)
	$(CPP) -o $@ $^ $(CFLAGS)
hist_angles: $(OBJhistangles)
	$(CPP) -o $@ $^ $(CFLAGS)
coord: $(OBJcoord)
	$(CPP) -o $@ $^ $(CFLAGS)
correlation: $(OBJcorr)
	$(CPP) -o $@ $^ $(CFLAGS)
normalize: $(OBJnormalize)
	$(CPP) -o $@ $^ $(CFLAGS)
readcol: $(OBJcolvar)
	$(CPP) -o $@ $^ $(CFLAGS)
#-------------------------------------

#----------------
# Cleaning files
#-----------------------------------
clean:
	rm -rf *.o
#-----------------------------------



