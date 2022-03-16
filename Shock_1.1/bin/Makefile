# path to the math routines
# external
# MATH=/pathtollapacklib
# local
MATH=../mtx

# compilo='intel'
compilo='gfortran'
# compilo='pgf90'

# Intel fortran
ifeq ($(compilo), 'intel')
 FC = ifort
 DEBUG= #  -g  -traceback -check bounds#-g # -s
 FFLAGS =  -fpe0 -r8 -O2  #$(DEBUG)
endif
# gfortran
ifeq ($(compilo), 'gfortran')
FC=gfortran
DEBUG= -g -C -fbacktrace -fcheck=bounds
#-fcheck=bounds
# FULL DEBUG MODE
# FFLAGS = -g3 -O0 -fno-second-underscore -Wall -Wextra -Wno-compare-reals -fcheck=all -fimplicit-none -std=f2003 -pedantic -ffpe-trap=zero,invalid,overflow -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function -framework CoreFoundation
FFLAGS = -g3 -O0 -fno-second-underscore -Wall -Wextra -Wno-compare-reals -fcheck=all -fimplicit-none -std=f2003 -pedantic -ffpe-trap=zero,invalid,overflow -fbacktrace -gdwarf-2 -fall-intrinsics -Wno-unused-function
# OPTIMIZED MODE
# FFLAGS= -O2  -fdefault-double-8 -fdefault-real-8 -ffpe-trap=invalid,zero,overflow  -finit-real=nan
endif
# pgf90
ifeq ($(compilo), 'pgf90')
FC=pgf90
DEBUG= -g -C -Miomutex -Minform,inform 
FFLAGS= -r8 -O2
endif

 OBJ =\
       techconfig.o \
       numerical_recipes.o \
       tools.o \
       constants.o \
       var_VODE.o \
       variables_mhd.o \
       read_fe_data.o \
       chemical_species.o \
       chemical_reactions.o \
       refh2.o \
       variables_thermal_balance.o \
       molecular_cooling.o \
       H2.o \
       CO.o \
       SiO.o \
       line_excit.o \
       shield.o \
       dust.o \
       radiation.o \
       evolution.o \
       integrateur_vode.o \
       energetics.o \
       profil_tables.o \
       initialize.o \
       outputs.o \
       outputs_hdf5.o \
       switch_cj.o \
       drive.o \
       mhd_vode.o

# Executable

mhd_vode: $(MATH)/libmtx.a $(OBJ)
	$(FC) $(FFLAGS) $(DEBUG) $(OBJ) -o ../mhd_vode -L$(MATH) -lmtx -L/usr/lib


# Build Local Lapack subroutines.
$(MATH)/libmtx.a:
	cd $(MATH); make

# Standard conversion rule
%.o:../src/%.f90
	$(FC) $(FFLAGS) $(DEBUG) -c $^ -o $@

%.o:../src/%.f
	$(FC) $(FFLAGS) $(DEBUG) -c $^ -o $@

clean:
	rm -f *.o *.mod ../mhd_vode

cleanall:
	make clean ; cd $(MATH) ; make clean
