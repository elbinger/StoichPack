# This is the global Makefile. If the Makefile and the Makefile.inc of each subdirectory is set up correctly, this Makefile will also work.

.PHONY = default clean

# build: cd to relevant directories and build there
default:
	cd StoichPack; $(MAKE)
	cd Examples; $(MAKE)

# clean up: cd to relevant directories and clean up there
clean:
	cd Examples; $(MAKE) clean
	cd StoichPack; $(MAKE) clean
	rm -f *~
