# This is the Makefile for the Utility directory of the StoichPack examples.
# Make sure thate ../Makefile.inc is set up correctly.

include ../Makefile.inc

.PHONY = clean

# target
LIB = libUtility.a

default: $(LIB)

# extract eigen files if not available. remove this if eigen is already installed.
$(EIGEN):
	cd ..; $(MAKE) eigen

# sub target
FDMMesh.o: FDMMesh.cpp FDMMesh.h $(EIGEN)
	$(CC) $(DEFAULT_CC_FLAGS) -c FDMMesh.cpp

# target
$(LIB): FDMMesh.o
	$(ARCH) $(LIB) FDMMesh.o

#clean up
clean:
	rm -f $(LIB) *.o *~
