C compare results with those of Parrot et al.
c Common storage:
c      include 'piccom.f'

      character*200 charin
      character*50 filename,string
      integer nr
      real phipic(1000),rhopic(1000)
      real rpic(1000)
      
      integer nti0
      parameter (nti0=100)
      real rti0(nti0)
      real phiti0(nti0)
      
      logical fullrange,zoom,debyezoom
      integer nsmax
      parameter (nsmax=2000)
      real fluxprobe(nsmax)

      integer nparrot
      parameter(nparrot=17)
      real rparrot(nparrot),phiparrot(nparrot)
c      real rhoparrot(nparrot)
c      data rhoparrot/.611,.782,.856,.895,.919,.936,.957,.969,.976
c     $     ,.982,.985,.988,.990,.991,.993,.993,.995/
      data rparrot/1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2,3.6,4.0,4.4,
     $     4.8,5.2,5.6,6.0,7.0/
      data phiparrot/-.493,-.246,-.156,-.111,-.084,
     $     -.066,-.044,-.032,-.024
     $     ,-.019,-.015,-.012,-.010,-.009,-.008,-.007,-.005/
      integer nlfm
      parameter (nlfm=17)
      real lfmz(nlfm),lfmr(nlfm),lfmphi01(nlfm),lfmn01(nlfm)
      data lfmz/0,0.5,1,1.5,2,4,6,8,10,12,14,16,18,20,22,24,26/
c Laframboise data for Rp/lambda=10 figs 13/14. 
c Probe potential -25 Te. Ti=Te
      data lfmn01/0.362069,0.350575,0.345977,0.343678,0.33908,0
     $     .345977,0.406897,0.511494,0.632184,0.712644,0.771264,0.817241
     $     ,0.848276,0.866667,0.886207,0.902299,0.913793/
      data lfmphi01/24.98848,21.64901,18.655,15.89129,13.58821,6.909258
     $     ,2.982497,1.140028,0.575772,0.356978,0,0,0,0,0,0,0/
      integer nalpert
      parameter (nalpert=16)
      real ralpert(nalpert),phialpert(nalpert)
      data ralpert/1.,1.001,1.002,1.004,1.006,
     $     1.01,1.02,1.04,1.06,
     $     1.1,1.2,1.4,1.6,2.0,2.5,3.0/
      data phialpert/-.690,-.658,-.638,-.610,-.591,
     $     -.569,-.536,-.480,-.430,
     $     -.362,-.261,-.164,-.117,-.068,-.041,-.028/

      fullrange=.true.
      zoom=.true.
      debyezoom=.true.
      filename(1:1)=' '
c Deal with arguments
      do 1 i=1,20
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            if(filename(1:1).eq.' ')goto 51
            goto 3
         elseif(string(1:2) .eq. '-f')then
            fullrange=.false.
         elseif(string(1:2) .eq. '-z')then
            zoom=.false.
         elseif(string(1:2) .eq. '-?')then
            goto 51
         else
            filename=string
         endif
 1    continue
 3    continue
c      call getarg(1,filename)
      if(filename(1:1) .eq. ' ') then
         open(10,file='picout.dat',status='old',err=101)
      else
         open(10,file=filename,status='old',err=101)
      endif
c Two lines for nothing.
      read(10,*)charin
      read(10,*,err=201)dt,vd,Ti,isteps,rhoinf,phiinf,fave,
     $     debyelen,vprobe
 201  continue
      write(*,*)'dt,vd,Ti,isteps,rhoinf,phiinf,fave'
      write(*,*)dt,vd,Ti,isteps,rhoinf,phiinf,fave
      read(10,*,err=202)nr
      do i=1,nr
         read(10,*,err=203)rpic(i),phipic(i)
c eventually we'll read the new variable ,rhopic(i)
         goto 213
 203     write(*,*)"rhopic error"
 213     continue
      enddo
      read(10,*)charin
      read(10,*)nsteps
      if(nsteps.gt.nsmax)nsteps=nsmax
      read(10,*)(fluxprobe(j),j=1,nsteps)

      close(10)
      do i=1,nr
         phipic(i)=phipic(i)-phipic(nr)
      enddo

      ir1=.6*nr
c Obtain the slope of the potential at the outer edge.
      call slopegen(phipic,rpic,nr,ir1,nr,slope,average)
c Now take the variation of perturbation to be \propto 1/r^2
c So slope is -2/r times the perturbation. Or perturbation is -(r/2)*slope.
      rave=(rpic(nr)+rpic(ir1))*0.5
      phiinf=average+rave*0.5*slope
      write(*,*)'average=',average,' slope=',slope,' phiinf=',phiinf   

c Ti=0 quasineutral case:
      do j=1,nti0
         phiti0(j)=-0.5*j/float(nti0)
         rti0(j)=sqrt(exp(-0.5-phiti0(j))/sqrt(-2.*phiti0(j)))
      enddo

c Laframboise data.
      write(*,*)'Laframboise Data'
      do i=1,nlfm
         lfmr(i)=1.+debyelen*lfmz(i)
         lfmphi01(i)=-lfmphi01(i)
         if(lfmphi01(i).eq.0.)then
            nlfmactual=i-1
            goto 61
         endif
         write(*,*)i,lfmr(i),lfmphi01(i)
      enddo
 61   continue

      call pfset(2)
c      call autoplot(rpic,phipic,nr)
c      call winset(.true.)
      do j=1,nr
         phipic(j)=phipic(j)-phiinf
      enddo
c      if(debyelen.gt..01)debyezoom=.false.
      nf=0
      if(fullrange)nf=nf+1
      if(zoom)nf=nf+1
      if(debyezoom)nf=nf+1
      if(nf.gt.1) call multiframe(nf,1,3)
      if(fullrange) then
         call autoplot(rpic,phipic,nr)
         call axlabels('r','phi')
c         call polymark(rpic,phipic,nr,4)
         call winset(.true.)
         call legendline(.6,.1,0," PIC")
         if(.not.debyezoom)then
            call color(2)
            call dashset(2)
c     call polyline(rparrot,phiparrot,nparrot)
            call polymark(rparrot,phiparrot,nparrot,1)
            call legendline(.6,.2,1," Parrot")
            call dashset(0)
            call color(10)
            call vecw(rpic(ir1),average+slope*(rpic(ir1)-rave),0)
            call vecw(rpic(nr),average+slope*(rpic(nr)-rave),1)
            call color(13)
            call polyline(ralpert,phialpert,nalpert)
            call legendline(.6,.3,0," Alpert")
         endif
c     Laframboise data here
         if(debyelen.eq.0.1)then
            call color(5)
            call polyline(lfmr,lfmphi01,nlfmactual)
            call polymark(lfmr,lfmphi01,nlfmactual,6)
            zoom=.false.
         endif
         call color(15)
         call winset(.false.)
      endif
      if(debyezoom)then
         call pltinit(1.,rpic(nr),-4.,0.)
         call axis()
         call axlabels('r','phi')
         call winset(.true.)
         call polyline(rpic,phipic,nr)
         call color(5)
         call polyline(lfmr,lfmphi01,nlfmactual)
         call polymark(lfmr,lfmphi01,nlfmactual,6)
         call color(15)
         call winset(.false.)
      endif
      if(zoom)then
         call pltinit(1.,1.25,-.8,0.1)
         call axis()
         call winset(.true.)
         call polyline(rpic,phipic,nr)
         call legendline(.6,.1,0," PIC")
         call color(2)
         call dashset(2)
         call polyline(rparrot,phiparrot,nparrot)
         call polymark(rparrot,phiparrot,nparrot,1)
         call legendline(.6,.2,-1," Parrot")
         call dashset(0)
         call color(10)
         call vecw(rpic(ir1),average+slope*(rpic(ir1)-rave),0)
         call vecw(rpic(nr),average+slope*(rpic(nr)-rave),1)
         call color(13)
         call polyline(ralpert,phialpert,nalpert)
         call polymark(ralpert,phialpert,nalpert,2)
         call legendline(.6,.3,-2," Alpert")
         call winset(.false.)
      endif
      call pltend()

      call multiframe(0,0,0)
      call yautoplot(fluxprobe,nsteps)
      call axlabels('step','Particles to probe')
      call pltend()
      return
 101  continue
      write(*,*)'Error opening ',filename
      return
 202  write(*,*)"nr error"
      return
 51   write(*,*)'Usage parrot [-f -z] filename'
      end
c***************************************************************************
c Get the average and slope over the rmesh range i1,i2.
      subroutine slopegen(phi,r,nr,i1,i2,slope,average)
      integer nr
      real phi(nr),r(nr)

c Assume r-mesh is linear
      rmom0=0.
      rmom1=0.
      rmid=(r(i2)+r(i1))/2.
      do i=i1,i2
         rmom0=rmom0+phi(i)
         rmom1=rmom1+(r(i)-rmid)*phi(i)
      enddo
      average=rmom0/(i2-i1+1)
      rave=rmom1/(i2-i1+1)
      rlen=r(i2)-r(i1)
      slope=12.*(rmom1)/(rlen*rlen)/(i2-i1+1)
c      write(*,*)rmom0,rmom1,r(i1),r(i2),rmid,rlen
      end
