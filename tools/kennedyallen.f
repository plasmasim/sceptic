c     Floating potential for spherical ABR traced from Kennedy and Allen
c     J.Plasma Phys 67 243 (2002). This curve is for argon. Given as +ve.
c
      function potlka(r)
c The ratio of probe radius to electron Debye length:
      real r
c
      integer nkapts
      parameter (nkapts=23)
      integer xyarr(2,nkapts)
      real roverl(nkapts),rlog(nkapts),fltptl(nkapts)
      logical linit
c Traced data:
      data xyarr/
     &     227,7635,507,7615,947,7595,1447,7515,1887,7435,2337,7305,
     &     2687,7175,3037,6975,3317,6785,3827,6305,4467,5635,5187,4795,
     &     6307,3435,6657,3075,7057,2715,7427,2435,7817,2215,8097,2075,
     &     8497,1945,9047,1835,9617,1775,10297,1735,11137,1725/
      data linit/.false./
      save

      if(.not.linit)then
c axis locator points of tracing
         xm4=237
         xp4=11187
         y0=7635
         y5=1965
         do i=1,nkapts
            rlog(i)=(-4.+8.*(xyarr(1,i)-xm4)/(xp4-xm4))
            roverl(i)=10**rlog(i)
            fltptl(i)=5.*(xyarr(2,i)-y0)/(y5-y0)
         enddo
         linit=.true.
c Testing only, first time, do a plot.
c         call lautoplot(roverl,fltptl,nkapts,.true.,.false.)
      endif
      xlog=log10(r)
      if(xlog.lt.rlog(1))then 
         potlka=0.
      elseif(xlog.gt.rlog(nkapts))then
         potlka=fltptl(nkapts)
      else
         potlka=finterp(rlog,fltptl,nkapts,xlog,xx)
         if(xx.eq.0)then
            write(*,*)'****Error: potlka argument out of range'
            potlka=-1.e-30
         endif
      endif
      end
c**************************************************************
c      program testpotlka
c Assumes the lautoplot is uncommented.
c      it=100
c      do k=1,it
c         xlog=(-4.5+9.*k/it)
c         x=10**xlog
c         y=potlka(x)
c         write(*,*)x,y
c         call polymark(x,y,1,10)
c      enddo
c      call pltend()
c      end
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
