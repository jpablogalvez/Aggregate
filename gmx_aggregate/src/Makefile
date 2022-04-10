# Paths definition
IDIR =../include
BDIR =../bin
ODIR =./obj
MDIR =./mod
SDIR =.
#libpath=/home/users/quicuan/jpablo/lib/xdrfile-1.1b/build/lib/libxdrfile.so
#libpath=/home/jpablo/lib/xdrfile-1.1b/build/lib/libxdrfile.so
libpath=/home/jpablo/lib/xdrfile/lib/libxdrfile.so
#libpath=/usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

# Compiler/Linker settings
FC = gfortran -I$(IDIR) -J$(MDIR)  # -Wno-argument-mismatch 
FLFLAGS = -g
FCFLAGS = -g -c -Wall -Wextra -Wconversion -Og -fcheck=bounds -fmax-errors=5
PROGRAM = aggregate
SRCF90  = gmx_aggregate.f90
SRCMOD  = datatypes.f90      \
          utils.f90          \
          input_section.f90  \
          geometry.f90       \
          graphtools.f90     \
          sorting.f90        \
          xdr.F90 

# Find all source files, create a list of corresponding object files
OBJ     = $(patsubst %.f90,$(ODIR)/%.o,$(SRCF90))
OBJ    := $(OBJ) $(patsubst %.F90,$(ODIR)/%.o,$(SRC))
MODOBJ  = $(patsubst %.f90,$(ODIR)/%.o,$(SRCMOD)) 
MODOBJ := $(patsubst %.F90,$(ODIR)/%.o,$(MODOBJ))
MODMOD  = $(patsubst $(ODIR)/%.o,$(MDIR)/%.mod,$(MODOBJ))

default: $(MODOBJ) $(OBJ)
	$(FC) $(MODOBJ) $(OBJ) -o $(BDIR)/$(PROGRAM) $(libpath)

# Compiler steps for all objects
$(ODIR)/%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ $^ 

$(ODIR)/%.o: %.f90
	$(FC) $(FCFLAGS) -o $@ $^ 

# Linker
$(BDIR)/$(PROGRAM) : $(OBJ)
	$(FC) $(FLFLAGS) -o $@ $< $(libpath) 

# Dependencies : Main program depends on all modules
$(ODIR)/$(PRGOBJ) : $(MODOBJ)

debug: 

	@echo "SRCF90  = $(SRCF90)"
	@echo "OBJ     = $(OBJ)"
	@echo "SRCMOD  = $(MOD)"
	@echo "MODOBJ  = $(MODOBJ)"
	@echo "PROGRAM = $(PROGRAM)"
	@echo "PRGOBJ  = $(PRGOBJ)"
	@echo "MODMOD  = $(MODMOD)"
#	DEBUG=1 make clean $(BDIR)/$(PROGRAM)

install:
	cp $(BDIR)/$(PROGRAM) ~/mis_bin

clean:
	rm -f $(ODIR)/*.o $(MDIR)/*mod $(BDIR)/$(PROGRAM) $(BDIR)/debug

