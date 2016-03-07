
default: epw

all: epw

epw: 
	(cd src ; make )
	(cd bin ; ln -fs ../src/epw.x . )

clean:
	cd src ; rm -f *.o *.mod  *~ 

release:
	cd ../ ; cp -r EPW EPW-release; cd EPW-release ; \
	rm -f src/*.o src/*.mod src/*.F90 src/*~ ; \
	rm bin/*.x ; \
	rm -rf examples/*/epw/out/* examples/*/epw/tmp/* \
	 examples/*/phonons/out/* examples/*/phonons/tmp/* \
	 examples/*/phonons/save/* ; \
	rm -rf .svn */.svn */*/*.svn */*/*/*.svn */*/*/*/*.svn 
	cd .. ; tar cfz EPW/EPW-release.tgz EPW-release ; \
	rm -rf EPW-release ; cd EPW
	
