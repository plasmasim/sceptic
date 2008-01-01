
LIBRARIES =  -L./accis/ -L/usr/X11R6/lib/ -laccisX -lXt -lX11 
G77=mpif77

# To silence warnings when compiling with g77 uncomment the next statement:
#NOGLOBALS=-Wno-globals

COMPILE-SWITCHES = -Wall -Wno-unused-variable -Wno-unused-labels -O2  -I.
# For debugging.
#  -g  -ffortran-bounds-check
# For profiling
#COMPILE-SWITCHES = -Wall -O2 -pg

# Old, and obsolete reinjection schemes.
#REINJECT=newinject.o
#REINJECT=reinject.o
#REINJECT=orbitinject.o extint.o
# Current: To allow both collisional and collision free:
REINJECT=fvinject.o orbitinject.o extint.o maxreinject.o ogeninject.o

MPICOMPILE-SWITCHES = $(NOGLOBALS) -DMPI $(COMPILE-SWITCHES)

OBJECTS = initiate.o advancing.o randc.o randf.o diags.o outputs.o outputlive.o chargefield.o  $(REINJECT) damp.o stringsnames.o rhoinfcalc.o

OBJECTSO = initiate.o advancingo.o randc.o randf.o diags.o outputs.o chargefield.o  $(REINJECT) damp.o stringsnames.o rhoinfcalc.o

MPIOBJECTS=sor2dmpi.o mpibbdy.o shielding.o

sceptic : sceptic.F piccom.f  ./accis/libaccisX.a $(OBJECTS) makefile
	$(G77) $(COMPILE-SWITCHES) -o sceptic sceptic.F $(OBJECTS) $(LIBRARIES)

sceptico : sceptic.F piccom.f  ./accis/libaccisX.a $(OBJECTSO) makefile
	$(G77) $(COMPILE-SWITCHES) -o sceptico sceptic.F $(OBJECTSO) $(LIBRARIES)

scepticmpi : sceptic.F piccom.f bbdydecl.f ./accis/libaccisX.a $(OBJECTS) $(MPIOBJECTS) makefile
	$(G77) $(MPICOMPILE-SWITCHES) -o scepticmpi  sceptic.F $(OBJECTS) $(MPIOBJECTS) $(LIBRARIES)

./accis/libaccisX.a : ./accis/*.f
	make -C accis

orbitint : orbitint.f coulflux.o $(OBJECTS) ./accis/libaccisX.a makefile
	$(G77) $(COMPILE-SWITCHES) -o orbitint orbitint.f $(OBJECTS) coulflux.o $(LIBRARIES)

coulflux.o : tools/coulflux.f
	$(G77) -c $(COMPILE-SWITCHES) tools/coulflux.f

fvinjecttest : fvinjecttest.F makefile fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o fvcom.f
	$(G77)  -o fvinjecttest $(COMPILE-SWITCHES) fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o  $(LIBRARIES)

fvinject.o : fvinject.f fvcom.f piccom.f
	$(G77) -c $(COMPILE-SWITCHES) fvinject.f

#pattern rule
%.o : %.f piccom.f fvcom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)

sceptic.tar.gz : ./accis/libaccisX.a sceptic scepticmpi
	make -C accis mproper
	make -C tools clean
	make clean
	./copyattach.sh
	tar chzf sceptic.tar.gz -C .. sceptic
	./copyremove.sh

clean :
	rm -f *.o
	rm -f *.ps
	rm -f *.dat
	rm -f *.orb
	rm -f *.frc
	rm -f *.html
	rm -f sihh snew Orbits.txt
	rm -f *~
	rm -f *.liv

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic.F *.f

