include ../Makefile.inc

.PHONY = clean
LIB = libUtility.a

default: $(LIB)

$(EIGEN):
	cd ..; $(MAKE) eigen

FDMMesh.o: FDMMesh.cpp FDMMesh.h $(EIGEN)
	$(CC) $(DEFAULT_CC_FLAGS) -c FDMMesh.cpp

$(LIB): FDMMesh.o
	$(ARCH) $(LIB) FDMMesh.o

clean:
	rm -f $(LIB) *.o *~
