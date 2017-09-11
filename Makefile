.PHONY = default distclean

default:
	cd StoichPack; $(MAKE)
	cd Examples; $(MAKE)

clean:
	cd Examples; $(MAKE) clean
	cd StoichPack; $(MAKE) clean
	rm -f *~
