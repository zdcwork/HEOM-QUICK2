.SILENT:
SHELL = /bin/sh
NAME=/path/to/HEOM-QUICK2.x.x.x/bin/HEOM-QUICK2.x

DFT=$(NAME)

F77     = ifort
FFLAGS  = -qopenmp -O2 -module objmod
LIBDIR  = -L/MKLPATH -L/opt/intel/lib/intel64/ -L/usr/lib64/ -L/usr/lib/ -I/MKLINCLUDE
LIBS    = $(LIBDIR)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -qopenmp -lpthread

all:;
	make project

project:;
	make DIR=module           obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=main             obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=td               obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=bicg             obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=tfqmr_double     obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=jacobi           obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=utility          obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=fermion          obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=openmp           obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=corrfunc         obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=qubits           obj     SHELL=/bin/sh  FF=$(F77)
	make DIR=heat_bath        obj     SHELL=/bin/sh  FF=$(F77)
	make $(DFT)


$(DFT): objects/*.o 
	rm -f $(DFT)
	$(F77) $(FOPTS) -module objmod objects/*.o $(LIBS) -o $(DFT)
	echo Compiling --$(DFT)-- done !

obj:;
	ext=".f90";\
	for d in  $$DIR/*.f90; \
	do \
	  filename=`basename $$d $$ext`; \
	  make FILE=$$filename DIR=$$DIR FF=$$FF objects/$$filename.o; \
	done 

objects/$(FILE).o: $(DIR)/$(FILE).f90;
	$(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90;\
	mv $(FILE).o objects;\
	echo $(FF) $(FFLAGS) -c $(DIR)/$(FILE).f90
	echo update $(DIR)/$(FILE)

clean:
	rm -f $(DFT) objects/*.o objmod/*.mod
	echo delete objects/*.o objmod/*.mod $(DFT)

