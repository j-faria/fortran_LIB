# no console output during compilation:
#.SILENT:

# *******************************************************
# ***          Comecar por limpar os sufixos          ***
# *******************************************************
.SUFFIXES:

# *******************************************************
# ***     Especificar os sufixos para .f .o .do       ***
# *******************************************************
#.SUFFIXES: .f90 .o

#
# *******************************************************
# ***                       Macros                    ***
# *******************************************************
#FC = /usr/bin/f77
#FC	= g77 -Wall
# compiler:
FC	= gfortran
FFLAGS = -O2 -ffree-form -fbacktrace -fno-range-check # -Wuninitialized
LINK = -llapack -L/usr/lib
# archiver:
ARQ = ar r

#dbx     = -O5 -r8 -g
#profil  = -p -O5 -r8 
#samedir = .
#FTN     = ftnchek

lib = /home/joao/Programs/fortran/lib

#
# *******************************************************
# *** Regra que por defeito produz os ficheiros .o **
# *******************************************************
%.o: %.f
	$(FC) $(FFLAGS) -c -o $@ $*.f $(LINK)
%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $*.f90

#
# *******************************************************
# ***   Especificar as directorias com as subrotinas  ***
# *******************************************************


OBJS =\
lib_algebra.o\
lib_array.o\
lib_matrix.o \
lib_assert.o\
lib_messages.o\
lib_random.o\
lib_statistics.o\
lib_opt.o \
lib_simplex.o \
lib_pikaia.o \
lib_pikaia12.o \
lib_periodogram.o \
lib_plot.o \
lib_splines.o \
lib_sort.o \
lib_io.o \
lib_utils.o \
lib_mcmc.o

OBJ_test = test_libs.o

# **********************************************************
# ***             Compilar os programas                  *** 
# **********************************************************

libmodules.a: $(OBJS)
	$(ARQ) $@ $^


test: $(OBJ_test)
	gfortran -fbacktrace -o $@ $^ -L$(lib) -I$(lib) $(LINK) -lmodules


clean:
	rm -f $(OBJS) libmodules.a *~
	
cleant:
	rm -f $(OBJ_test)
