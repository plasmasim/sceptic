
libraries = -L/usr/X11R6/lib/ -L../accis/ -laccisX -lXt -lX11 
G77=g77
COMPILE-SWITCHES = -Wall -O2
#COMPILE-SWITCHES = -Wall -O2 -pg

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : %.f piccompost.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(libraries)

% : %.F
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(libraries)

all : ../accis/libaccisX.a postprocB postproc parrot fluxVp newsumplot fcol 

../accis/libaccisX.a : ../accis/*.f
	make -C ../accis

postproc : postproc.f piccompost.f ../accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g  -o postproc postproc.f coulflux.f $(libraries)

parrot : parrot.f piccompost.f ../accis/libaccisX.a
	$(G77) $(COMPILE-SWITCHES) -o parrot -fno-backslash parrot.f -laccisX $(libraries)

fluxVp : fluxVp.f ../accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g -o fluxVp fluxVp.f coulflux.f $(libraries)

postprocB : postprocB.f piccompost.f ../accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g  -o postprocB postprocB.f coulflux.f $(libraries)

fcol : fcol.f kennedyallen.o
	$(G77)  -o fcol $(COMPILE-SWITCHES) fcol.f  kennedyallen.o $(libraries)

clean :

	rm -f *.o
	rm -f postproc postprocB parrot fluxVp newsumplot fcol phioutplot rhopicplot pdiag
	rm -f *.dat *.frc
	rm -f *.ps *.txt
