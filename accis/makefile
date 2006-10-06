# Attempts to get a general make file going.
# Primary source is now in linux.
ifeq ($(G77),gfortran)
	G77no=gfortran -fno-backslash
	#G77=gfortran
	#NOBACKSLASH=
else
	G77=g77 -L.
	G77no=g77 -fno-backslash
	NOBACKSLASH= -fno-backslash
# -g  -ffortran-bounds-check 
endif

object_files= \
autopl.o \
contour.o \
fitrange.o \
axis.o \
axlabels.o \
axregion.o \
charsize.o \
lautopl.o \
pltinit.o \
vecw.o \
winset.o \
getfont.o \
polyline.o \
labeline.o \
minmax.o \
folwcont.o \
scale.o \
vec4014.o \
ticset.o \
hid.o \
vec3.o \
examine.o \
axon.o \
drwstr.o \
sysfont.o 

#vecnp needs no -E switch.

# Directory of source files and headers.
src = /home/hutch/accis/

libraries = -L/usr/X11R6/lib/ -lXt -lX11

#header include files
headers= hidcom.h plotcom.h world3.h

all : libaccis.a accissum.html

accissum.html : accissum.tex
	tth -e2 accissum
	purge.tex

#update the library.
libaccis.a : $(object_files) vecnp.o vecx.o
	echo "Updating libaccis."
	ar -rs libaccis.a $(object_files) vecnp.o
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecx.o


#  Tell we depend on .f files, so keep them.
.PRECIOUS : %.f

# Headers must be updated, but this section seems to override the
# implicit rule unless the header dependency is :: double.
#Objects depend on headers
$(object_files) :: $(headers) makefile

#vecnp needs no -E switch. Jan 97
vecnp.o : vecnp.f makefile
	$(G77no) -c vecnp.f

#vecx is the C linkage to X-window.
vecx.o : vecx.c makefile
	gcc -c vecx.c

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : %.f ;
	$(G77no) -c $*.f

plottest : plottest.f libaccisX.a
	$(G77) -o plottest $(NOBACKSLASH) plottest.f -laccisX $(libraries)

hidtest : hidtest.f libaccisX.a
	$(G77) -o hidtest $(NOBACKSLASH) hidtest.f -laccisX  $(libraries)

hidtest2 : hidtest2.f libaccisX.a
	$(G77) -o hidtest2 $(NOBACKSLASH) hidtest2.f -laccisX  $(libraries)

fontshow : fontshow.f libaccisX.a
	$(G77) -o fontshow $(NOBACKSLASH) fontshow.f -laccisX  $(libraries)

sysfontshow : sysfontshow.f libaccisX.a
	$(G77) -o sysfontshow $(NOBACKSLASH) sysfontshow.f -laccisX  $(libraries)

eshow : eshow.f libaccisX.a
	$(G77) -o eshow $(NOBACKSLASH) eshow.f -laccisX  $(libraries)

contest : contest.f libaccisX.a
	$(G77) -o contest $(NOBACKSLASH) contest.f -laccisX  $(libraries)

box : box.f libaccisX.a
	$(G77) -o box $(NOBACKSLASH) box.f -laccisX  $(libraries)

colortest : colortest.f libaccisX.a
	$(G77) -o colortest $(NOBACKSLASH) colortest.f -laccisX  $(libraries)

libaccisf2c : $(object_files) vecnp.o vecxf2c.o
	echo "Updating libaccisf2c."
	ar -rs libaccis.a $(object_files) vecnp.o
	cp libaccis.a libaccisf2c.a
	ar -d libaccisf2c.a vec4014.o
	ar -q libaccisf2c.a vecxf2c.o

% : %.f
	$(G77) -o $*  $*.f -laccisX  $(libraries)

mproper : clean
	rm -f *.a

clean :
	rm -f *.o
	rm -f *.ida
	rm -f plot*.ps
	rm -f contest
	rm -f plottest
	rm -f hidtest
	rm -f fontshow
	g77 -fno-backslash -c drwstr.f          

% : %.f libaccis.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  -laccisX $(libraries)

