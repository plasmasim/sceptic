#__________________________________________________________________________
#
#     This code is copyright (c)
#              Ian H Hutchinson    hutch@psfc.mit.edu.
#              Leonardo Patacchini patacchi@mit.edu
#
#     It may be used freely with the stipulation that any scientific or
#     scholarly publication concerning work that uses the code must give
#     an acknowledgement referring to the relevant papers
#
#     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
#     1953 (2002), vol 45, p 1477 (2003).
#
#     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
#     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
#
#     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
#     p013505 (2007)
#
#     The code may not be redistributed except in its original package.
#
#     No warranty, explicit or implied, is given. If you choose to build
#     or run the code, you do so at your own risk.
#___________________________________________________________________________

# Universal Makefile for sceptic

#Defaults compiler (mpif77 compiler)
ifeq ("$(G77)","")
	G77=mpif77
endif
#Default Xlib (32 bit)
ifeq ("$(XLIB)","")
	XLIB=/usr/X11R6/lib
endif
#Default Accis lib
ifeq ("$(ACCISLIB)","")
	ACCISLIB=./accis
endif

LIBRARIES =  -L$(XLIB) -L$(ACCISLIB) -laccisX -lXt -lX11 

#Default No Warnings
ifeq ("$(NOWARN)","")
	NOWARN=
endif

COMPILE-SWITCHES =-Wall -Wno-unused-variable  $(NOWARN) -O2  -I.
# For debugging.
#  -g  -ffortran-bounds-check
# For profiling
#COMPILE-SWITCHES = -Wall -O2 -pg

REINJECT=fvinject.o orbitinject.o extint.o maxreinject.o ogeninject.o

MPICOMPILE-SWITCHES = -DMPI $(COMPILE-SWITCHES)

OBJECTS = initiate.o advancing.o randc.o randf.o diags.o outputs.o	\
outputlive.o chargefield.o $(REINJECT) damp.o stringsnames.o		\
rhoinfcalc.o shielding.o

OBJECTSO = initiate.o advancingo.o randc.o randf.o diags.o outputs.o	\
chargefield.o $(REINJECT) damp.o stringsnames.o rhoinfcalc.o		\
shielding.o

MPIOBJECTS=mpibbdy.o sor2dmpi.o shielding_par.o 

# So make does not do multiple tries. The default target is makefile first.
all : makefile sceptic

sceptic : sceptic.F  piccom.f  ./accis/libaccisX.a $(OBJECTS) makefile
	$(G77) $(COMPILE-SWITCHES) -o sceptic sceptic.F  $(OBJECTS) $(LIBRARIES)

# The real Makefile
MAKEFILE=makefile

makefile : Makefile MFSconfigure
	@echo Configuring the Makefile for this platform.
	rm -f ./accis/makefile
	rm -f *.o
	rm -f *./accis/*.o
	make -C accis
	./MFSconfigure
	@echo Now running make again using the new Makefile
	make -f $(MAKEFILE)

sceptico : sceptic.F  piccom.f  ./accis/libaccisX.a $(OBJECTSO) makefile
	$(G77) $(COMPILE-SWITCHES) -o sceptico sceptic.F  $(OBJECTSO) $(LIBRARIES)

scepticmpi : sceptic.F  piccom.f piccomsor.f ./accis/libaccisX.a $(OBJECTS) $(MPIOBJECTS) makefile
	$(G77) $(MPICOMPILE-SWITCHES) -o scepticmpi  sceptic.F   $(OBJECTS) $(MPIOBJECTS) $(LIBRARIES)

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

sceptic.tar.gz : ./accis/libaccisX.a sceptic scepticmpi
	make -C accis mproper
	make -C tools clean
	make makefile
	make clean
	./copyattach.sh
	tar chzf sceptic.tar.gz  --exclude *.tar.gz -C .. sceptic
	./copyremove.sh

clean :
	rm -f makefile
	rm -f *.o
	rm -f *.ps
	rm -f *.orb
	rm -f *.html
	rm -f Orbits.txt
	rm -f *~
	rm -f *.liv
	rm -f sceptic.tar.gz

cleanall :
	make clean
	rm -f *.dat
	rm -f *.frc

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic.F *.f


#pattern rules need to be at end not to override specific rules
%.o : %.f piccom.f fvcom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F makefile
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)
