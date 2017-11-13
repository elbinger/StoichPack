# This is the Makefile for the Examples directory of the StoichPack library.
# Make sure that Makefile.inc is set up correctly.

.PHONY = clean distclean

# build: cd to relevant directories and build there.
default: eigen
	cd Utility; $(MAKE)
	cd MyBiochemicalSystems; $(MAKE)
	cd ODEExample; $(MAKE)
	cd PDEExample; $(MAKE)

# extract eigen. remove this if eigen is already installed.
eigen:
	tar -xjf ../eigen-eigen-5a0156e40feb.tar.bz2

# clean up: cd to relevant directories and clean up there.
clean:
	rm -f *~
	cd Utility; $(MAKE) clean
	cd MyBiochemicalSystems; $(MAKE) clean
	cd ODEExample; $(MAKE) clean
	cd PDEExample; $(MAKE) clean
	rm -rf eigen

