/* Xwindow driver for accis plotting */
/* Fortran callable (f2c) routines */
/* ********************************************************************** */
/*
Refreshing version.
*/
  #include  <X11/StringDefs.h> 
  #include  <X11/Intrinsic.h> 
  #include  <X11/Core.h> 

/* Globals for these routines*/
Widget accis_wshell;
Widget accis_drawing;
Display *accis_display;
Drawable accis_window;
Pixmap accis_pixmap;

GC accis_gc;
int accis_depth;
Colormap accis_colormap;

struct Screen_Size {  
    Dimension width; /* unsigned short */
    Dimension height;
  };
    static struct Screen_Size s_s;
/* Static maximum number of points in path */
#define accis_path_max 4000
XPoint accis_path[accis_path_max];
int accis_pathlen=0;

#define a_gradPixno 256
unsigned long a_gradPix[a_gradPixno];

#define a_maxPixels 16
unsigned long accis_pixels[a_maxPixels];

char *accis_colornames[a_maxPixels]=
{
  "White",
  "MediumBlue",
  "SeaGreen",
  "MediumTurquoise",
  "Firebrick",
  "Orchid",
  "Brown",
  "LightGrey",
  "SlateGrey",
  "Blue",
  "Green",
  "Cyan",
  "Red",
  "Magenta",
  "Yellow",
  "Black"
};



/* Subroutine */ 
int svga_(scrxpix, scrypix, vmode, ncolor)
short *scrxpix, *scrypix, *vmode, *ncolor;
{  
  static int second=0;
  extern int f__xargc;
  extern char **f__xargv;
  int svga_argc=0;
  char *svga_argv[1]; 
  static int n;
  static Arg wargs[10];
  int theDepth;
  Colormap theColormap;
  extern void config_handler();
  extern void accis_refresh();
  XWindowAttributes attributes_return;
  long zero=0,s4k=65535;

  if(second == 0){
    accis_wshell = XtInitialize("accis","Accis", NULL, 0,
      &f__xargc, f__xargv);  /* This refers to the f2c command line args.
			It may only work with f2c, therefore. */
      /* &svga_argc, svga_argv); alternate */

    accis_drawing = XtCreateManagedWidget("drawing",coreWidgetClass,
		 accis_wshell, NULL, 0);

    *vmode=88;
    *ncolor=15;
    /* Set up a default size of the drawing window widget. 
       This is overruled by Accis*geometry resources setting wshell size.
       */
    n = 0;
    XtSetArg(wargs[n], XtNheight, 480); n++;
    XtSetArg(wargs[n], XtNwidth, 640); n++;
    XtSetValues(accis_drawing, wargs, n);
    XtRealizeWidget(accis_wshell);
    accis_display = XtDisplay(accis_drawing);
    accis_window = XtWindow(accis_drawing);
    accis_gc = XCreateGC(accis_display, accis_window, 0, NULL);
    accis_depth=DefaultDepth(accis_display,0);
    accis_colormap=DefaultColormap(accis_display,0);
    initDefaultColors();
    /* Leave setup for resizing. */
    n = 0;
    XtSetArg(wargs[n], XtNheight, &s_s.height); n++;
    XtSetArg(wargs[n], XtNwidth, &s_s.width); n++;
    /* Pixmap setup */
    XtGetValues(accis_wshell, wargs, n);
    accis_pixmap=XCreatePixmap(accis_display,accis_window,
			       s_s.width,s_s.height,accis_depth);
    XtAddEventHandler(accis_drawing,ExposureMask,FALSE,
		      accis_refresh,NULL);
    XSelectInput(accis_display,accis_window,
		 KeyPressMask | ExposureMask | ButtonPress
		 | FocusChangeMask | EnterWindowMask |
		 ButtonReleaseMask | ButtonMotionMask ); /* Tell events */
    /*Default gradient is a gray-scale.*/
    accisgradinit_(&zero,&zero,&zero,&s4k,&s4k,&s4k);
    second++;
  }else{
    /* Need to have this to make the get correct. */
    ManageEvents();
  } 
  /* This doesn't need a configuration event handler to work correctly. */
  XtGetValues(accis_wshell, wargs, n);
  *scrxpix=s_s.width;
  *scrypix=s_s.height;
  XFlush(accis_display);
  XClearWindow(accis_display,accis_window);
  XSetForeground(accis_display,accis_gc,accis_pixels[0]);  
  /* Clear the window to background color. VMS doesn't do correctly.
     But also this seems to fix the bad match error. */
  XFillRectangle(accis_display,accis_window,accis_gc,0,0,
		 s_s.width,s_s.height);
  /* Clear the pixmap  */
  XFillRectangle(accis_display,accis_pixmap,accis_gc,0,0,
		 s_s.width,s_s.height);
  XSetForeground(accis_display,accis_gc,accis_pixels[15]); 
  XRaiseWindow(accis_display, accis_window);
  XGetWindowAttributes (accis_display, accis_window,&attributes_return);
  /* Then set the focus. We should not get a bad match that way*/
  if(attributes_return.map_state == IsViewable){
    /* XSetInputFocus(accis_display, accis_window, RevertToPointerRoot,
       CurrentTime);
    Sometimes it is easier not to do this. */
  }
  return 0;
}

/* ******************************************************************** */
initDefaultColors()
{
  XColor theRGBColor,theHardColor;
  int status;
  int i;
  for(i=0;i<a_maxPixels;i++){
    status=XLookupColor(accis_display,accis_colormap,accis_colornames[i],
			&theRGBColor,&theHardColor);
    if(status !=0){
      status=XAllocColor(accis_display,accis_colormap,&theHardColor);
      if(status !=0){
	accis_pixels[i]=theHardColor.pixel;
      }else{
	accis_pixels[i]=BlackPixel(accis_display,0);
      }
    }else{
      accis_pixels[i]=BlackPixel(accis_display,0);
    }
  }
}

/* ******************************************************************** */
/* End plotting and return to text editing. */
/* Subroutine */ 
/* #include <curses.h>*/
int txtmode_()
{
  XEvent event; 
  XFlush(accis_display);
  do{
    /*    printf("Executing XtNextEvent"); */
    XtNextEvent(&event);
    /* XNextEvent(accis_display,&event); is equivalent */
    XtDispatchEvent(&event);  
    /*    printf("The event type: %d\n",event); */
  }while(event.type != ButtonPress && event.type != KeyPress );
  /* Here we should give the focus back to parent, but I don't see how.
    XSetInputFocus(accis_display, ??, PointerRoot,
		   CurrentTime); */
}

/* ******************************************************************** */
/* Flush the plot buffer
/* Subroutine */ 
int accisflush_()
{
  XFlush(accis_display);
}

/* ********************************************************************* */
/* Subroutine */ int scolor_(li)
long *li;
{
  /* *ncolor=*li; */
  if((*li < a_maxPixels) && (*li >= 0)){
    XSetForeground(accis_display,accis_gc,accis_pixels[(int) *li]);
    return 1;
  }else{    
    return 0;
  }
} /* scolor_ */

/* ******************************************************************** */
/* Subroutine */ int vec_(px, py, ud)
long *px, *py, *ud;
{ /*  Draw vector on screen, with pen up or down. */
    static int px1=0,py1=0,px2=0,py2=0;
    extern XPoint accis_path[];
    extern int accis_pathlen;

    px1=px2;
    py1=py2;
    px2 = *px;
    py2 = *py;
    if( *ud != 0) {
      XDrawLine(XtDisplay(accis_drawing),XtWindow(accis_drawing), accis_gc,
		  px1,py1,px2,py2);
      XDrawLine(XtDisplay(accis_drawing),accis_pixmap, accis_gc,
		  px1,py1,px2,py2);
      if(accis_pathlen<accis_path_max){      /* Add point to path */
	accis_pathlen++;
      }
    }else{ /* Restart path */
      accis_pathlen=0;
    }
    accis_path[accis_pathlen].x=*px;
    accis_path[accis_pathlen].y=*py;
/*    XFlush(accis_display);
 Flush removed here. Now relies on txtmode to flush display.  */
    return 0;
} /* vec_ */
/* ******************************************************************** */
int vecfill_()
{
    extern XPoint accis_path[];
    extern int accis_pathlen;
    if(accis_pathlen>1){ /* If path is more than 2 points, fill. */
      XFillPolygon(accis_display,accis_window,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
      XFillPolygon(accis_display,accis_pixmap,accis_gc,
		   accis_path,accis_pathlen+1,Nonconvex,CoordModeOrigin);
    }
}
/* ******************************************************************** */
/* Not currently in use. */
void config_handler(w,cs_s,event)
Widget w;
struct Screen_Size cs_s;
XEvent *event;
{
cs_s.width=event->xconfigure.width;
cs_s.height=event->xconfigure.height;
printf("config_handler width= %d, height= %d\n",cs_s.width,cs_s.height);
}
/* ******************************************************************** */
/* ******************************************************************** */
void accis_refresh(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  XCopyArea(XtDisplay(w),accis_pixmap,accis_window,accis_gc,0,0,
	    s_s.width,s_s.height,0,0);
  XFlush(accis_display);
}
/* ******************************************************************** */

ManageEvents()
{
  XEvent event; 
  while(XtPending()){
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
  }
}
/* ******************************************************************** */
/* Testing only */
/*
main()
 {
   long li=1;
   int x=200,y=150,pen=1;
   int a,b,c,d;
   int ch;
   while(ch != 'q'){
     printf("Calling svga\n");
     svga_(&a,&b,&c,&d);
     printf("Returned, x= %d, y= %d, mode= %d, ncolor= %d,\n",
	    a,b,c,d);
     printf("Drawing stuff directly\n");
     draw_graphics(accis_drawing);
     printf("Drawing line using vec_\n");
     for(a=0;a<16;a++){
       li=a;
       scolor_(&li);
       for(b=0;b<10;b++){
	 pen=0;x=200;y=150+10*a+b;
	 vec_(&x,&y,&pen); 
	 pen=1;x=400;y=150+10*a+b;
	 vec_(&x,&y,&pen); 
       }
     }
     ch=getchar();
   }
 }

  draw_graphics(w)
  Widget w; {
       Display *display;
       Drawable window;
       GC gc;
       int store;
       XWindowAttributes attributes;
       XSetWindowAttributes sattributes;
       unsigned long valuemask;

       display = XtDisplay(w);
       window = XtWindow(w);
       store= DoesBackingStore(DefaultScreenOfDisplay(display));
       printf("DoesBackingstore return, %d of 
               NotUseful,WhenMapped,Always %d,%d,%d\n",
	      store,NotUseful,WhenMapped,Always);
       XGetWindowAttributes(display,window,&attributes);
       printf("Got attributes. backing_store,planes,pixel,%d %x %u\n"
	      ,attributes.backing_store,attributes.backing_planes,
	      attributes.backing_pixel);


       gc = XCreateGC(display, window, 0, NULL);
       XSetForeground(display, gc, 50);
       XSetBackground(display, gc, 0);

       XDrawLine(display, window, gc, 10, 10, 400, 400);
       XDrawRectangle(display, window, gc, 75, 110, 150, 100);
       XDrawArc(display, window, gc, 75, 110, 150, 100, 45*64, 120*64);

       XFreeGC(display, gc);
  }


*/

/* ******************************************************************** */
/* Routines for interactive examination of the plot. Accomplished
 by calling eye3d(ival) from the fortran code, prior to pltend.
 If ival is 1 on return then we just moved the view. If zero we hit a key,
 implying we should carry on.
*/

float xeye,yeye,zeye;
float xeye0,yeye0,zeye0;
float accis_x0,accis_y0;

void accis_butdown(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  accis_x0=event->xbutton.x;
  accis_y0=event->xbutton.y;
  butdown_(&xeye0,&yeye0,&zeye0);
  xeye=xeye0; yeye=yeye0; zeye=zeye0;
}

void accis_butup(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  butup_(&xeye,&yeye,&zeye);
}


void accis_moved(w,data,event)
Widget w;
caddr_t data;
XEvent *event;
{
  float xmoved,ymoved;
  /* need screen to normal scaling not done here at present*/ 
  xmoved=event->xbutton.x-accis_x0;
  ymoved=event->xbutton.y-accis_y0;
/*    printf("eye: %f %f %f\n",xeye,yeye,zeye); */
  viewrot_(&xmoved,&ymoved,&xeye0,&yeye0,&zeye0,&xeye,&yeye,&zeye);
  cubeupd_(&xeye,&yeye,&zeye);
}
/* ******************************************************************** */
int eye3d_(value)
     int *value;
{
  extern void accis_butdown();
  extern void accis_butup();
  extern void accis_moved();
  XEvent event; 
  XFlush(accis_display);  
  XtAddEventHandler(accis_drawing,ButtonPressMask,FALSE,
		      accis_butdown,NULL);
  XtAddEventHandler(accis_drawing,ButtonReleaseMask,FALSE,
		      accis_butup,NULL);
  XtAddEventHandler(accis_drawing,ButtonMotionMask,FALSE,
		      accis_moved,NULL);
  do{
/*      printf("Executing XtNextEvent "); */
    XtNextEvent(&event);
    XtDispatchEvent(&event);  
/*      printf("The event type: %d\n",event); */
  }while(event.type != KeyPress && event.type != ButtonRelease  );
  if( event.type != ButtonRelease || 
      ( xeye==xeye0 && yeye==yeye0 && zeye==zeye0) ) *value=0; else *value=1;
}
/* ******************************************************************** */
/* Routines for using 256 color gradients in preference to 16 fixed colors.*/
/************** Setup The Gradient **********************/
int accisgradinit_(r1,g1,b1,r2,g2,b2)
     long *r1,*g1,*b1,*r2,*g2,*b2;
     /* RGB are specified in the range 0 to 65535 */
{
  int i,j,status;
  XColor theRGBcolor;
  for (i=0;i<a_gradPixno;i++){
    j=(i* *r2+(a_gradPixno-1-i)* *r1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    theRGBcolor.red=j;
    j=(i* *g2+(a_gradPixno-1-i)* *g1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    theRGBcolor.green=j;
    j=(i* *b2+(a_gradPixno-1-i)* *b1)/(a_gradPixno-1.) ;
    if(j<0)j=0; else if(j>65535)j=65535;
    theRGBcolor.blue=j;
    /*      printf("theRGBcolor %d,%d,%d\n",theRGBcolor.red,theRGBcolor.green,theRGBcolor.blue); */
    if(XAllocColor(accis_display,accis_colormap,&theRGBcolor)){
      a_gradPix[i]=theRGBcolor.pixel;
    }else{
      a_gradPix[i]=BlackPixel(accis_display,0);
    }
    }
}
/********** Use a gradient color out of 256 *********************************/
/* Subroutine */ int acgradcolor_(li)
long *li;
{
  /* *ncolor=*li; */
  if((*li < a_gradPixno) && (*li >= 0)){
    XSetForeground(accis_display,accis_gc,a_gradPix[(int) *li]);
    return 1;
  }else{    
    return 0;
  }
}
/********** Tell the current rgb color *********************************/
/* Subroutine */ int getrgbcolor_(ipixel,red,green,blue)
long *ipixel,*red,*blue,*green;
{
  XColor theRGBcolor;
  theRGBcolor.pixel=a_gradPix[(int) *ipixel];
  XQueryColor(accis_display,accis_colormap,&theRGBcolor);
  *red=theRGBcolor.red;
  *green=theRGBcolor.green;
  *blue=theRGBcolor.blue;
  return 0;
}

