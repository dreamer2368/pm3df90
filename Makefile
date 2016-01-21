### Compilers & flags
F90=gfortran

FFTWLIBS=~/bin/FFTW/lib/libfftw3.a
VECLIBSMACOSX=
LAPACKLIB=-L/opt/local/lib/lapack-3.5.0 -llapack -lblas
BLASLIB=/opt/local/lib/lapack-3.5.0/librefblas.a
PNETCDFLIBS=

LIBS    = $(LAPACKLIB) $(FFTWLIBS)


EXE = exec
F90SRC = main.f90 random.f90 declaration.f90 init.f90 modPlasma.f90 timeStep.f90 assignFunctions.f90 convergence.f90 MatrixSolver.f90 constants.f90 testmodule.f90
F90OBJ = main.o random.o declaration.o init.o modPlasma.o timeStep.o assignFunctions.o convergence.o MatrixSolver.o constants.o testmodule.o

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
declaration.o : modPlasma.o
init.o: declaration.o random.o MatrixSolver.o
assignFunctions.o : declaration.o
timeStep.o : assignFunctions.o MatrixSolver.o
testmodule.o : init.o timeStep.o
main.o : testmodule.o
#convergence.o : init.o timeStep.o
#main.o: convergence.o

clean:
	rm *.o *.mod $(EXE)

.PHONY: all run clean


