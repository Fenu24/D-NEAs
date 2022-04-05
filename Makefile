include make.flags
#FC=gfortran
#FFLAGS=-J.mod -g -O 
#LINK=-llapack -fopenmp# links towards external libraries


#########################################################
# executable names and corresponding sources

SRC0=gamma_est_mc.f90 used_const.f90 yarko_force.f90
EXE0=gamma_est_mc.x

#########################################################
# liste of object files with /obj prefix

OBJ_NOPREFIX0=$(SRC0:.f90=.o)
OBJ0=$(addprefix .obj/,$(OBJ_NOPREFIX0))

BIN=bin

################################################################
# creation of executable programs

all: $(BIN)/$(EXE0) links_ex

$(BIN)/$(EXE0): $(OBJ0)
	$(FC) $^ -o $@ $(LINK)

links_ex:
	rm -rf *.x test/*.x 2>/dev/null 
	ln -s bin/$(EXE0)

#########################################################
# compilation of sources

#obj/%.o: src/%.f90
#	$(FC) $(FCFLAGS) -c $< -o $@

.obj/gamma_est_mc.o: src/gamma_est_mc.f90 .obj/yarko_force.o .obj/used_const.o
	$(FC) $(FFLAGS) -c $< -o $@

.obj/yarko_force.o: src/yarko_force.f90 .obj/used_const.o
	$(FC) $(FFLAGS) -c $< -o $@

.obj/used_const.o: src/used_const.f90 
	$(FC) $(FFLAGS) -c $< -o $@

#########################################################
# suppression of compilation/output files

clean:
	rm -f *.x .obj/*.o .mod/*.mod bin/* 

#########################################################
