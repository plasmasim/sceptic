# Attempts to get a general make file going.
# Primary source is now in linux.
G77=g77
NOBACKSLASH=-fno-backslash
#G77=gfortran
#NOBACKSLASH=

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
drwstr.o 
#vecnp.o \
#vecnp needs no -E switch.

# Directory of source files and headers.
src = /home/hutch/accis/

libraries = -L/usr/X11R6/lib/ -lXt -lX11

#header include files
headers= hidcom.h plotcom.h world3.h

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
vecnp.o : vecnp.f
	$(G77) -c $(NOBACKSLASH) vecnp.f

#vecx is the C linkage to X-window.
vecx.o : vecx.c
	gcc -c vecx.c

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : %.f ;
	$(G77) -c $(NOBACKSLASH) $*.f

#Problems 17 May 96. 
# vecnp causes f2c segmentation fault. Data statements at fault.
# 18 May 96 vecnp can't be fixed, but dos f2c works on it. Also used the
# ca definition fix so we can use -E switch.
#
# 24 May recompiled f2c from source using the -DCRAY switch on malloc.c
# it now works on these two files.
# We still have a problem with commons in vecnp unless this is compiled
# WITHOUT the -E switch, because its commons seem to be declared external.
# In xterm does not autoswitch to Tek mode. Fixed by getting the PFC
# version of vec4014. However, now we have lost the alternate fonts.
# Old version has the same problem. Perhaps we never had them.
# Ah! need to compile plottest with '-!bs' switch to get the backslash
# sequences. Then greek letters are there.
# Discovered the problem to lie with the fact that ichar(char(x)) is 
# negative if x > 128. Fixed in drwstr.for.
#
# 5 Jun 99 rewrote makefile for $(G77). '-!bs' becomes -fno-backslash.
# Had some problems with parameter char() in vecnp, vec4014, drwstr.
# Intrinsic not allowed in parameter statement. Used ' ' strings instead.
# xargc/v had to become f__xargc/v for $(G77) command line access.

# 10 Jul 99 found things were not fixed properly. There were some additional
# char() in parameters that had to be fixed. Also vecnp ' ' approach was
# not working. In the end had to use printable character ' \' to denote line
# end, and old style in lnswrt. Also there was a deep bug in wstr that it
# was not being properly terminated and so the string length was wrong.

# 26 July 99 got f2c going by recompiling the library from source and removing
# the .so. version which was being used by default. Then created an f2c
# version of the accis library by replacing the vecx with vecxf2c which
# accesses the command line arguments as xargv/c. Other stuff left the same.
# It seems to work ok with geqread.make


plottest : plottest.f libaccisX.a
	$(G77) -o plottest $(NOBACKSLASH) plottest.f -laccisX $(libraries)

hidtest : hidtest.f libaccisX.a
	$(G77) -o hidtest $(NOBACKSLASH) hidtest.f -laccisX  $(libraries)

hidtest2 : hidtest2.f libaccisX.a
	$(G77) -o hidtest2 $(NOBACKSLASH) hidtest2.f -laccisX  $(libraries)

fontshow : fontshow.f libaccisX.a
	$(G77) -o fontshow $(NOBACKSLASH) fontshow.f -laccisX  $(libraries)

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

mproper : clean
	rm -f *.a

clean :
	rm *.o
	rm -f *.ida
	rm -f plot*.ps
	rm -f contest
	rm -f plottest
	rm -f hidtest
	rm -f fontshow
	g77 -fno-backslash -c drwstr.f		

% : %.f
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  -laccisX $(libraries)

