# This is the Makefile for the ODEExamples of StoichPack.
# Make sure that ../Makefile.inc is set up correctly.

include ../Makefile.inc

.PHONY = clean

# libs
LIBDIR = ../MyBiochemicalSystems
LIB = MyBiochemicalSystems

default: ODEExample

# build libMyBiochemicalSystems.a if not available
$(LIBDIR)/lib$(LIB).a: $(LIBDIR)/*.h $(LIBDIR)/*.cpp $(STOICHPACK)/*.h
	cd $(LIBDIR); $(MAKE)

# extract eigen files if not available. remove this if eigen is already installed.
$(EIGEN):
	cd ..; $(MAKE) eigen

# target
ODEExample: $(LIBDIR)/lib$(LIB).a ODEExample.o ODEConstants.o
	$(CC) $(DEFAULT_CC_FLAGS) -L$(LIBDIR) -o ODEExample ODEExample.o ODEConstants.o -l$(LIB) $(DEFAULT_LINK)

# sub targets in this directory
ODEExample.o: $(EIGEN) ODEExample.cpp ODEConstants.h Makefile ../Makefile.inc ../StoichPackEigen.h $(STOICHPACK)/*.h
	$(CC) $(DEFAULT_CC_FLAGS) -I$(LIBDIR) -c ODEExample.cpp

ODEConstants.o: ODEConstants.h ODEConstants.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c ODEConstants.cpp

# clean up
clean:
	rm -f *.o *~ out_* ODEExample
