#G77 compiler version.
G77=g77 -fno-backslash

# Intel fortran Compiler:
#G77=ifc
#C_COMPILE_DEFINES=-DIFC
#COMPILE_SWITCHES= -nbs -w90 -cm
#EXTRA_LIBS= -lPEPCF90 -lSM -lICE 

LIBRARIES = -L/usr/X11R6/lib/ -lXt -lX11 $(EXTRA_LIBS)

# Directory of source files and headers.
SRC=

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
vecnp.o \
drwstr.o 

fortran_files=$(addprefix $(SRC), $(addsuffix .f,$(basename $(object_files) )))

#header include files
headers= $(SRC)hidcom.h $(SRC)plotcom.h $(SRC)world3.h

c_files=$(SRC)vecx.c

#update the accis libraries.
libaccis.a : $(object_files) vecx.o
	echo "Updating libaccis."
	ar -rs libaccis.a $(object_files)
	cp libaccis.a libaccisX.a
	ar -d libaccisX.a vec4014.o
	ar -q libaccisX.a vecx.o

# Headers must be updated, but this section seems to override the
# implicit rule unless the header dependency is :: double.
# Objects depend on headers
$(object_files) :: $(headers)


#vecx is the C linkage to X-window.
vecx.o : $(SRC)vecx.c
	gcc $(C_COMPILE_DEFINES) -o vecx.o -c $(SRC)vecx.c

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : $(SRC)%.f ;
	$(G77) $(COMPILE_SWITCHES) -c -o $*.o $(SRC)$*.f

#This is how to compile any fortran code that calls accis, that resides here.
% : $(SRC)%.f libaccisX.a
	$(G77)  -o $* $(COMPILE_SWITCHES) $(SRC)$*.f  ./libaccisX.a $(LIBRARIES)

clean :
	rm *.o
	rm *.a

