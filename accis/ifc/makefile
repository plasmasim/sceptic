# This is the real makefile for this directory.
# Do a make to create the intel fortran compiler version of accis libraries.
G77=ifc
C_COMPILE_DEFINES=-DIFC
COMPILE_SWITCHES= -nbs -w90 -cm
libraries = -L/usr/X11R6/lib/ -lXt -lX11 -lPEPCF90 -lSM -lICE 
# Directory of source files and headers.
#src=accis/
src=../

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

fortran_files=$(addprefix $(src), $(addsuffix .f,$(basename $(object_files) )))

#header include files
headers= $(src)hidcom.h $(src)plotcom.h $(src)world3.h

c_files=$(src)vecx.c

#update the libraries.
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
vecx.o : $(src)vecx.c
	gcc $(C_COMPILE_DEFINES) -o vecx.o -c $(src)vecx.c

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : $(src)%.f ;
	$(G77) $(COMPILE_SWITCHES) -c -o $*.o $(src)$*.f

% : $(src)%.f libaccisX.a
	$(G77)  -o $* $(COMPILE_SWITCHES) $(src)$*.f  ./libaccisX.a $(libraries)

%.f : $(src)%.f ;
	cp $(src)$*.f .

clean :
	rm -f *.o
	rm -f *.a

mproper :
	rm -f *.o
	rm -f *.a

accis.tar.gz : $(fortran_files) $(headers) $(c_files) $(src)Makefile $(src)accissum.tex
	tar czf accis.tar.gz $(fortran_files) $(headers) $(c_files) $(src)Makefile $(src)accissum.tex
	@echo Created Archive
