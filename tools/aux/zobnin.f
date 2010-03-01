c     Floating potential for spherical collisional collection of Zobnin
c     JETP 91, 483 (2000). For neon, Ti/Te=.01, versus \lambda_s/mfp.
c
      function potlzob(r,it)
c The ratio of probe radius to electron Debye length:
      real r
c Curve number
      integer it
c
      integer nkapts,nz2,nz3,nzmax
      parameter (nkapts=14,nz2=11,nz3=17,nzmax=40)
      integer xyarr(2,nkapts),xya2(2,nz2),xya3(2,nz3)
      integer nlen(3)
      real roverl(nzmax,3),rlog(nzmax,3),fltptl(nzmax,3)
      logical linit
c Traced data:
c a/R_D=0.06
      data xyarr/
     $   1967,2955,2497,3285,2977,3585,3677,4015,4257,4345,4757,4625,
     $   5227,4865,5617,5015,5937,5095,6277,5125,6557,5065,6857,4925,
     $   7067,4775,7357,4515/          
c 0.012
      data xya2/
     $	 1987,3135,2907,3695,4007,4375,4767,4835,5407,5185,5917,5425,
     $	 6287,5565,6597,5645,6957,5715,7327,5745,8107,5735/
c 0.24
      data xya3/
     $   1997,2595,2457,2935,2947,3225,3357,3435,3777,3615,4307,3805,
     $	 4757,3925,5087,3975,5327,3995,5617,3955,5837,3875,6137,3715,
     $ 	 6367,3555,6597,3355,6807,3135,7067,2835,7377,2445/
      data nlen/nkapts,nz2,nz3/
      data linit/.false./
      save
c Default curve to use
      if(it.eq.0)it=1
c Initialize if needed.
      if(.not.linit)then
c axis locator points of tracing
         xm1=3630
         xp1=8340
         y0=7350
         y2=2565
         do i=1,nkapts
            rlog(i,1)=(-1.+2.*(xyarr(1,i)-xm1)/(xp1-xm1))
            roverl(i,1)=10**rlog(i,1)
            fltptl(i,1)=2.*(xyarr(2,i)-y0)/(y2-y0)
         enddo
         do i=1,nz2
            rlog(i,2)=(-1.+2.*(xya2(1,i)-xm1)/(xp1-xm1))
            roverl(i,2)=10**rlog(i,2)
            fltptl(i,2)=2.*(xya2(2,i)-y0)/(y2-y0)
         enddo
         do i=1,nz3
            rlog(i,3)=(-1.+2.*(xya3(1,i)-xm1)/(xp1-xm1))
            roverl(i,3)=10**rlog(i,3)
            fltptl(i,3)=2.*(xya3(2,i)-y0)/(y2-y0)
         enddo
         linit=.true.
c Testing only, first time, do a plot.
         call lautoplot(roverl(1,1),fltptl(1,1),nkapts,.true.,.false.)
         do i=2,3
            call polyline(roverl(1,i),fltptl(1,i),nlen(i))
         enddo
      endif
      xlog=log10(r)
      if(xlog.lt.rlog(1,it))then 
         potlzob=0.
      elseif(xlog.gt.rlog(nlen(it),it))then
         potlzob=fltptl(nlen(it),it)
      else
         potlzob=finterp(rlog(1,it),fltptl(1,it),nlen(it),xlog,xx)
         if(xx.eq.0)then
            write(*,*)'****Error: potlzob argument out of range'
            potlzob=-1.e-30
         endif
      endif
      end
c**************************************************************
      program testpotlzob
c Assumes the lautoplot is uncommented.
      nmark=100
      icurve=3
      do k=1,nmark
         xlog=(-2.+2.6*k/nmark)
         x=10**xlog
         y=potlzob(x,icurve)
         write(*,*)x,y
         call polymark(x,y,1,10)
      enddo
      call pltend()
      end
c********************************************************************
c Given two monotonic functions Q(x), R(x) 
c on a 1-D grid x=1..nq, solve Q(x)=y for x and interpolate R(x)
c That is, invert Q to give x=Q^-1(y) and give R(x), or simply
c interpolate the value of R when Q=y.
c If x is returned as 1, then the y is outside Q's range.
      function finterp(Q,R,nq,y,x)
      real Q(nq),R(nq)
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql

      finterp=0.
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         xpart=(y-Ql)/(Qr-Ql)
         x=xpart+iql
         finterp=R(iql)+xpart*(R(iqr)-R(iql))
      else
         x=iql
         write(*,*)'****** Error!: finterp coincident points'
      endif
      end
c**********************************************************************
