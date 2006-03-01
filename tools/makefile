
libraries = -L/usr/X11R6/lib/ -L/home/hutch/accis/ -laccisX -lXt -lX11 
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

all : postprocB postproc parrot fluxVp newsumplot

postproc : postproc.f piccompost.f /home/hutch/accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g  -o postproc postproc.f coulflux.f $(libraries)

parrot : parrot.f piccompost.f /home/hutch/accis/libaccisX.a
	$(G77) $(COMPILE-SWITCHES) -o parrot -fno-backslash parrot.f -laccisX $(libraries)

fluxVp : fluxVp.f /home/hutch/accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g -o fluxVp fluxVp.f coulflux.f $(libraries)

postprocB : postprocB.f piccompost.f /home/hutch/accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g  -o postprocB postprocB.f coulflux.f $(libraries)


clean :

	rm -f *.o
	rm -f postproc parrot fluxVp newsumplot
	rm -f *.dat *.frc
	rm -f *.ps *.txt
