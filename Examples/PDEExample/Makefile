# This is the Makefile for the PDEExamples for StoichPack.
# Make sure, that ../Makefile.inc is set up correctly.

include ../Makefile.inc

.PHONY = clean

# libs
LIBDIR1 = ../MyBiochemicalSystems
LIB1 = MyBiochemicalSystems
LIBDIR2 = ../Utility
LIB2 = Utility

default: PDEExample

# make libMyBiochemicalSystems.a if not available
$(LIBDIR1)/lib$(LIB1).a: $(LIBDIR1)/*.h $(LIBDIR1)/*.cpp $(STOICHPACK)/*.h
	cd $(LIBDIR1); $(MAKE)

# make libUtility.a if not available
$(LIBDIR2)/lib$(LIB2).a: $(LIBDIR2)/*.h $(LIBDIR2)/*.cpp
	cd $(LIBDIR2); $(MAKE)

# extract eigen if not available. remove this if eigen is already installed.
$(EIGEN):
	cd ..; $(MAKE) eigen

# target
PDEExample: PDEExample.o PDEConstants.o $(LIBDIR2)/lib$(LIB2).a $(LIBDIR1)/lib$(LIB1).a
	$(CC) $(DEFAULT_CC_FLAGS) -L$(LIBDIR1) -L$(LIBDIR2) -o PDEExample PDEExample.o PDEConstants.o -l$(LIB1) -l$(LIB2) $(DEFAULT_LINK)

# subtargets in this directory
PDEExample.o: $(EIGEN) PDEExample.cpp PDEConstants.h Makefile ../Makefile.inc ../StoichPackEigen.h $(STOICHPACK)/*.h $(LIBDIR1)/*.h $(LIBDIR2)/*.h
	$(CC) $(DEFAULT_CC_FLAGS) -I$(LIBDIR1) -I$(LIBDIR2) -c PDEExample.cpp

PDEConstants.o: PDEConstants.h PDEConstants.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c PDEConstants.cpp

# clean up
clean:
	rm -f *.o *~ out_* PDEExample
