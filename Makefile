include make.flags

#########################################################
# executable names and corresponding sources

SRC0=gamma_est_mc.f90 used_const.f90 yarko_force.f90
EXE0=gamma_est_mc.x

SRC1=gamma_est_mc_threadsafe.f90 used_const.f90 yarko_force.f90
EXE1=gamma_est_mc_threadsafe.x

SRC2=yarko_grid_rho_D.f90 used_const.f90 yarko_force.f90
EXE2=yarko_grid_rho_D.x

SRC3=yarko_compare.f90 used_const.f90 yarko_force.f90
EXE3=yarko_compare.x

#########################################################
# liste of object files with /obj prefix

OBJ_NOPREFIX0=$(SRC0:.f90=.o)
OBJ0=$(addprefix .obj/,$(OBJ_NOPREFIX0))

OBJ_NOPREFIX1=$(SRC1:.f90=.o)
OBJ1=$(addprefix .obj/,$(OBJ_NOPREFIX1))

OBJ_NOPREFIX2=$(SRC2:.f90=.o)
OBJ2=$(addprefix .obj/,$(OBJ_NOPREFIX2))

OBJ_NOPREFIX3=$(SRC3:.f90=.o)
OBJ3=$(addprefix .obj/,$(OBJ_NOPREFIX3))

BIN=bin

################################################################
# creation of executable programs

all: $(BIN)/$(EXE0) $(BIN)/$(EXE1) $(BIN)/$(EXE2) $(BIN)/$(EXE3) links_ex

$(BIN)/$(EXE0): $(OBJ0)
	$(FC) $(FFLAGS) $^ -o $@ $(LINK)

$(BIN)/$(EXE1): $(OBJ1)
	$(FC) $(FFLAGS) $^ -o $@ $(LINK)

$(BIN)/$(EXE2): $(OBJ2)
	$(FC) $(FFLAGS) $^ -o $@ $(LINK)

$(BIN)/$(EXE3): $(OBJ3)
	$(FC) $(FFLAGS) $^ -o $@ $(LINK)

links_ex:
	rm -rf *.x 2>/dev/null 
	ln -s bin/$(EXE0)
	ln -s bin/$(EXE1)
	ln -s bin/$(EXE2)
	ln -s bin/$(EXE3)

#########################################################
# compilation of sources

#obj/%.o: src/%.f90
#	$(FC) $(FCFLAGS) -c $< -o $@

.obj/yarko_compare.o: src/yarko_compare.f90 .obj/yarko_force.o .obj/used_const.o
	$(FC) $(FFLAGS) -c $< -o $@

.obj/yarko_grid_rho_D.o: src/yarko_grid_rho_D.f90 .obj/yarko_force.o .obj/used_const.o
	$(FC) $(FFLAGS) -c $< -o $@

.obj/gamma_est_mc.o: src/gamma_est_mc.f90 .obj/yarko_force.o .obj/used_const.o
	$(FC) $(FFLAGS) -c $< -o $@

.obj/gamma_est_mc_threadsafe.o: src/gamma_est_mc_threadsafe.f90 .obj/yarko_force.o .obj/used_const.o
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
