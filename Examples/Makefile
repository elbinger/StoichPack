.PHONY = clean distclean

default: eigen
	cd Utility; $(MAKE)
	cd MyBiochemicalSystems; $(MAKE)
	cd ODEExample; $(MAKE)
	cd PDEExample; $(MAKE)

eigen:
	tar -xjf ../eigen-eigen-5a0156e40feb.tar.bz2

clean:
	rm -f *~
	cd Utility; $(MAKE) clean
	cd MyBiochemicalSystems; $(MAKE) clean
	cd ODEExample; $(MAKE) clean
	cd PDEExample; $(MAKE) clean
	rm -rf eigen

