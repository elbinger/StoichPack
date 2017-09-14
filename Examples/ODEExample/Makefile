include ../Makefile.inc

.PHONY = clean

LIBDIR = ../MyStoichiometries
LIB = MyStoichiometries

default: ODEExample

ODEExample: $(LIBDIR)/lib$(LIB).a ODEExample.o ODEConstants.o
	$(CC) $(DEFAULT_CC_FLAGS) -L$(LIBDIR) -o ODEExample ODEExample.o ODEConstants.o -l$(LIB) $(DEFAULT_LINK)

ODEExample.o: ODEExample.cpp ODEConstants.h Makefile ../Makefile.inc ../StoichPackEigen.h $(STOICHPACK)/*.h
	$(CC) $(DEFAULT_CC_FLAGS) -I$(LIBDIR) -c ODEExample.cpp

ODEConstants.o: ODEConstants.h ODEConstants.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c ODEConstants.cpp

clean:
	rm -f *.o *~ out_* ODEExample
