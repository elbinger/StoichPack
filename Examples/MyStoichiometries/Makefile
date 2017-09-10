include ../Makefile.inc

.PHONY = clean
STDDEP = ../Makefile.inc ../../StoichPack/*.h Makefile
LIB = MyStoichiometries

default: lib$(LIB).a

MyReactions.o: $(STDDEP) MyReactions.h MyReactions.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c MyReactions.cpp

MyStoichiometries.o: $(STDDEP) MyReactions.h MyStoichiometries.h MyStoichiometries.cpp
	$(CC) $(DEFAULT_CC_FLAGS) -c MyStoichiometries.cpp

lib$(LIB).a: MyReactions.o MyStoichiometries.o
	$(ARCH) lib$(LIB).a MyReactions.o MyStoichiometries.o

clean:
	rm -f *~ *.o *.a
