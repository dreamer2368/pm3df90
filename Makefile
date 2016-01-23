### Compilers & flags
F90=gfortran

FFTWLIBS=~/bin/FFTW/lib/libfftw3.a
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(LAPACKLIB) $(FFTWLIBS)


EXE = exec
F90SRC = constants.f90 random.f90 MatrixSolver.f90 modPlasma.f90 modMesh.f90 modAssign.f90 modRecord.f90 modPM3D.f90 timeStep.f90 init.f90 testmodule.f90 main.f90
F90OBJ = constants.o random.o MatrixSolver.o modPlasma.o modMesh.o modAssign.o modRecord.o modPM3D.o timeStep.o init.o testmodule.o main.o

### Targets
all: $(EXE)
run: $(EXE) 
	./$(EXE)

# Link object files to executables
$(EXE): $(F90OBJ)
	$(F90) -o $(EXE) $(F90OBJ) $(LIBS)

# All .o files depend on the corresponding .f90 file
%.o: %.f90
	$(F90) -c $<

# Dependencies
MatrixSolver.o : constants.o
modPlasma.o : constants.o
modMesh.o : MatrixSolver.o
modAssign.o : modPlasma.o modMesh.o
modRecord.o : modPlasma.o modMesh.o
modPM3D.o : modPlasma.o modMesh.o modAssign.o modRecord.o
timeStep.o : modPM3D.o
init.o: modPM3D.o random.o
testmodule.o : init.o timeStep.o
main.o : testmodule.o
#assignFunctions.o : declaration.o
#main.o : testmodule.o
#convergence.o : init.o timeStep.o
#main.o: convergence.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


