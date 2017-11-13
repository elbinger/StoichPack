# This is the Makefile for the MyBiochemicalSystems example of StoichPack.
# Make sure thate ../Makefile.inc is set up correctly.

include ../Makefile.inc

.PHONY = clean
STDDEP = ../Makefile.inc ../../StoichPack/*.h Makefile
LIB = MyBiochemicalSystems
OBJECTS = MyReactions.o MyBiochemicalSystems.o

default: lib$(LIB).a

# sub targets
MyReactions.o: $(STDDEP) MyReactions.h MyReactions.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c MyReactions.cpp

MyBiochemicalSystems.o: $(STDDEP) MyReactions.h MyBiochemicalSystems.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c MyBiochemicalSystems.cpp

# target
lib$(LIB).a: $(OBJECTS)
	$(ARCH) lib$(LIB).a $(OBJECTS)

# clean up
clean:
	rm -f *~ *.o *.a
