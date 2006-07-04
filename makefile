# Intel fortran compiler case:
#LIBRARIES = -L/usr/X11R6/lib/ -L/home/hutch/accis/ifc/ -laccisX -lXt -lX11 -lPEPCF90 -lSM -lICE
#G77=ifc
#COMPILE-SWITCHES = -O3 -w90 -cm

# g77 default case:
LIBRARIES = -L/usr/X11R6/lib/ -L./accis/ -laccisX -lXt -lX11 
G77=mpif77
COMPILE-SWITCHES = -Wall -O2  -I. 
#not necessary -fno-backslash
# For debugging.
#  -g  -ffortran-bounds-check
# For profiling
#COMPILE-SWITCHES = -Wall -O2 -pg

# Old, and obsolete reinjection schemes.
#REINJECT=newinject.o
#REINJECT=reinject.o
#REINJECT=orbitinject.o extint.o
# Current: To allow both collisional and collision free:
REINJECT=fvinject.o orbitinject.o extint.o maxreinject.o

# MPI version needs the beowulf libraries. Edit the MPILIBS to point to your
# Local versions. 

MPICOMPILE-SWITCHES = -Wall  -Wno-globals -DMPI -O2 -I.
OBJECTS = initiate.o advancing.o randc.o randf.o diags.o outputs.o chargefield.o  $(REINJECT) damp.o stringsnames.o
MPIOBJECTS=sor2dmpi.o mpibbdy.o shielding.o

sceptic : sceptic.F piccom.f  ./accis/libaccisX.a $(OBJECTS) makefile
	$(G77) $(COMPILE-SWITCHES) -o sceptic sceptic.F $(OBJECTS) $(LIBRARIES)

scepticmpi : sceptic.F piccom.f bbdydecl.f ./accis/libaccisX.a $(OBJECTS) $(MPIOBJECTS) makefile
	$(G77) $(MPICOMPILE-SWITCHES) -o scepticmpi  sceptic.F $(OBJECTS) $(MPIOBJECTS) $(LIBRARIES)

./accis/libaccisX.a : ./accis/*.f
	make -C accis

fvinjecttest : fvinjecttest.F makefile fvinject.o initiate.o advancing.o chargefield.o randf.o fvcom.f
	$(G77)  -o fvinjecttest $(COMPILE-SWITCHES) fvinjecttest.F fvinject.o initiate.o advancing.o chargefield.o randf.o  $(LIBRARIES)

fvinject.o : fvinject.f fvcom.f piccom.f
	$(G77) -c $(COMPILE-SWITCHES) fvinject.f

#pattern rule
%.o : %.f piccom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)

sceptic.tar.gz : ./accis/libaccisX.a sceptic
	make -C accis/ifc clean
	make -C accis mproper
	make -C tools clean
	make clean
	./copyattach.sh
	tar chzf sceptic.tar.gz -C .. SCEPTIC

clean :
	rm -f *.o
	rm -f *.ps
	rm -f *.dat
	rm -f *.orb
	rm -f *.frc
	rm -f *.html
	rm -f sihh snew Orbits.txt
	rm -f *~


ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic.F *.f

