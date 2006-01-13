c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c This code is copyright(c) (2003-5) Ian H Hutchinson hutch@psfc.mit.edu
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build
c or run the code, you do so at your own risk.
c
c Version 2.6 Aug 2005.
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c Version 2.5; Aug 2004
c***********************************************************************
      subroutine chargediag(dt,istep)
      real dt
      integer istep
c Common data:
      include 'piccom.f'
      real rhoplot(nrsize),rho1theta(nthsize),rhomidtheta(nthsize)
      real rhomidave(nthsize),rho1ave(nthsize)
      real phiave(nrsize)
      real riave
      save

c This allows us to restart with nstepsave .ne. 1 if rhoinf is set.
      if(riave.eq.0)riave=rhoinf

c Calculate rhoplot,diagphi,diagrho,rho1theta,rhomidtheta
      do i=1,nr
 510     format(10f8.1)
         rhoplot(i)=0.
         phiave(i)=0.
         nrp=0
         do j=1,NTHUSED
            nrp=nrp+1
c     rhoplot is unnormalized. All others are normalized.
            rhoplot(i)=rhoplot(i)+rho(i,j)*rhoinf
            phiave(i)=phiave(i)+phi(i,j)
         enddo
         rhoplot(i)=rhoplot(i)/nrp
         phiave(i)=phiave(i)/nrp
c This needs to be fixed for Debye code.
         diagphi(i)=(diagphi(i)*(nstepsave-1)+phiave(i))/nstepsave
c old quasineutral way:
c         log(rhoplot(i)/rhoinf))/nstepsave
         diagrho(i)=(diagrho(i)*(nstepsave-1) + (rhoplot(i)))/nstepsave
      enddo

c Calculate diagchi (outer potential as a function of nth normalized
c to the ion thermal velocity) Necessay for the reinjection (Not here)
      do j=1,NTHUSED
        diagchi(j)=(diagchi(j)*(nstepsave-1)+phi(NRUSED,j)/Ti)/nstepsave
      enddo
      

c*******
c New rhoinf calculation.
      if(nrein .gt. 0) then
c estimate of the rhoinf based on flux from this step.

c Trial of different scheme. Combination equivalent to phihere usage.
         averein=(diagphi(NRFULL)+diagphi(NRUSED))*.5

         if(averein.gt.0.5*Ti)then
c This is necessary to prevent smaxflux errors. smaxflux is not correct
c for repulsive potentials.
c            write(*,*)'Excessive averein',averein,' capped'
            averein=0.5*Ti
         endif

c     if nbc, we use the linearized distribution out of the boundary
c     (valid for vd<<vTi and -averein<<Ti)
         if (.not.nbc) then
            riest=(nrein/dt) /
     $           (sqrt(2*Ti)*
     $           smaxflux(vd/sqrt(2.*Ti),(-averein/Ti))
     $           *r(NRFULL)**2 )
            
         else
            riest=(nrein/dt) /
     $           (sqrt(2*Ti)*
     $           smaxflux2(-averein/Ti)
     $           *r(NRFULL)**2 )
         endif


  
c         write(*,*)'nrein=',nrein,'  psum=',psu,
c     $        '  averein=',averein,' riest=',riest
      else
c If no valid estimate, just keep it the same.
c         write(*,*) 'nrein=0'
         if(.not. rhoinf.gt.1.e-4)then
            write(*,*)'rhoinf error in chargediag:',rhoinf,
     $           ' approximated.'
c estimate rhoinf approximately:
            rhoinf=numprocs*npart/(4.*pi*r(NRFULL)**3/3.)
         endif
         averein=0.
         riest=rhoinf
      endif
c      write(*,*)'riest,riave',riest,riave
      riave=(riave*(nstepsave-1) + riest)/nstepsave 
c*******
 540  format(2f12.4,i6,3f12.4)
      rhoinf=riave
      rhoplotmax=1.
      do j=1,NTHUSED
         rho1theta(j)=rho(1,j)
         rhomidtheta(j)=0.
         icount=0
         do k=nr/2,NRUSED
            icount=icount+1
            rhomidtheta(j)=rhomidtheta(j)+rho(k,j)
c     $           psum(k,j)*volinv(k)*(nth-1.)*np/rhoinf
         enddo
c Average flux, should use ninthstep which is partreduced.
c         finthave(j)=((nstepsave-1)*finthave(j) + ninth(j))/nstepsave
c Gives same result, in single processor as:
         finthave(j)=((nstepsave-1)*finthave(j)
     $        + ninthstep(j,istep-1))/nstepsave
         rhomidtheta(j)=rhomidtheta(j)/float(icount)
         rhomidave(j)=((nstepsave-1)*rhomidave(j)+rhomidtheta(j))
     $        /nstepsave
         rho1ave(j)=((nstepsave-1)*rho1ave(j) + rho1theta(j))/nstepsave
         if(rhomidave(j).gt.rhoplotmax)rhoplotmax=rhomidave(j)
         if(rho1ave(j).gt.rhoplotmax)rhoplotmax=rho1ave(j)
      enddo

      if(lplot)then
c Plot density.
c Sideways theta cross-sections.
         call pltinit(.12*nint(10.*rhoplotmax),0.,th(1),th(nth))
         call axis()
         call vecw(1.,th(1),0)
         call vecw(1.,th(nth),1)
         call axlabels('density','cos!Aq!@')
         call winset(.true.)
         call polyline(rho1theta(1),tcc(1),NTHUSED)
         call polyline(rhomidtheta(1),tcc(1),NTHUSED)
         call dashset(2)
         call polyline(rhomidave(1),tcc(1),NTHUSED)
         call polyline(rho1ave(1),tcc(1),NTHUSED)
         call dashset(0)
         call winset(.false.)
      endif
c      call pltend()
      end
c***********************************************************************
c***********************************************************************
c Plot phi and phi averaged over angles.
      subroutine phidiag()
c Common data:
      include 'piccom.f'
      real phiplot(nrsize)


c      write(*,*)'phi='
      do i=1,nr
c         write(*,510)(phi(i,j),j=1,nth)
 510      format(10f8.3)
          phiplot(i)=0.
          do j=1,NTHUSED
             phiplot(i)=phiplot(i)+phi(i,j)
          enddo
          phiplot(i)=phiplot(i)/(nth)
      enddo
 511  format(10f8.4)
      phimax=log(1.1*float(npart)/((4.*pi/3.)*(r(nr)**3-r(1)**2)))
      phimax=phimax+log(float(numprocs))
      call pltinit(r(1),r(nr),phimax-1.,phimax)
      call axis()
      call axlabels('radius','!Af!@/T!de!d')
      call winset(.true.)
      call polyline(r,phiplot,nr)
      call dashset(2)
      call polyline(r,diagphi,nr)
      call dashset(0)
      call winset(.false.)
c      call pltend()
      end
c***********************************************************************
c***********************************************************************
      subroutine vrdiag()
c Common data:
      include 'piccom.f'
      real vrplot(nrsize),psumtot(nrsize),vrsmooth(nrsize)
      real vtplot(nrsize),vpplot(nrsize),vtsmooth(nrsize),
     $     vpsmooth(nrsize),
     $     trsmooth(nrsize),ttpsmooth(nrsize)
      real tempr(nrsize,nthsize),temprplot(nrsize),
     $     temptp(nrsize,nthsize),temptpplot(nrsize)
      character*30 char1
      save

c      call pltinit(0.,1.,0.,1.)
c      return

      do i=1,NRUSED
          vrplot(i)=0.
          vtplot(i)=0.
          vpplot(i)=0.
          temprplot(i)=0.
          temptpplot(i)=0.
          psumtot(i)=0.
          do j=1,NTHUSED
             vrplot(i)=vrplot(i)+vrsum(i,j)
             vtplot(i)=vtplot(i)+vtsum(i,j)
             vpplot(i)=vpplot(i)+vpsum(i,j)
             psumtot(i)=psumtot(i)+psum(i,j)
             tempr(i,j)=vr2sum(i,j)*psum(i,j)-vrsum(i,j)**2
             if(psum(i,j).gt.0.)tempr(i,j)=tempr(i,j)/psum(i,j)**2
             temprplot(i)=temprplot(i)+tempr(i,j)*psum(i,j)
             temptp(i,j)=0.5*(vtp2sum(i,j)*psum(i,j)
     $            - (vtsum(i,j)**2+ vpsum(i,j)**2))
             if(psum(i,j).gt.0.)temptp(i,j)=temptp(i,j)/psum(i,j)**2
             temptpplot(i)=temptpplot(i)+temptp(i,j)*psum(i,j)
             if(psum(i,j).gt.0.)
     $            diagvr(i,j)=(diagvr(i,j)*(nstepsave-1) +
     $            vrsum(i,j)/(psum(i,j)+.001))/nstepsave
          enddo
          vrplot(i)=vrplot(i)/psumtot(i)
          vtplot(i)=vtplot(i)/psumtot(i)
          vpplot(i)=vpplot(i)/psumtot(i)
          temprplot(i)=temprplot(i)/psumtot(i)
          temptpplot(i)=temptpplot(i)/psumtot(i)
          vrsmooth(i)=(vrsmooth(i)*(nstepsave-1) + vrplot(i))/nstepsave
          vtsmooth(i)=(vtsmooth(i)*(nstepsave-1) + vtplot(i))/nstepsave
          vpsmooth(i)=(vpsmooth(i)*(nstepsave-1) + vpplot(i))/nstepsave
          trsmooth(i)=(trsmooth(i)*(nstepsave-1) + temprplot(i))
     $         /nstepsave
          ttpsmooth(i)=(ttpsmooth(i)*(nstepsave-1) + temptpplot(i))
     $         /nstepsave
      enddo
      call minmax(vrplot,NRUSED,vrmin,vrmax)
      vrmin=int(vrmin-.7)
      vrmin=min(vrmin,-1.1*sqrt(1.+0.5*Ti))
      call pltinit(r(1),r(nr),vrmin,0.1)
      call axis()
      call axlabels('radius','')
      call jdrwstr(wx2nx(r(1))-.04,wy2ny(-0.65),'v!dr!d',2.)
      call winset(.true.)
      call polyline(rcc(1),vrplot,NRUSED)
c      call polyline(rcc(1),vtplot,NRUSED)
      call vecw(1.,0.,0)
      call vecw(r(nr),0.,1)
      call dashset(2)
      call polyline(rcc(1),vrsmooth,NRUSED)
c      call polyline(rcc(1),vtsmooth,NRUSED)
      call dashset(0)
c Temperature diagnostics.
c      write(*,*)trsmooth
      call minmax(temptpplot,NRUSED,tpmin,tpmax)
      call fitrange(0.,tpmax,2,nxfac,xfac,xtic,x1st,xlast)
      tpmax=xlast
      if(Ti.ge.1.)tpmax=nint(max(tpmax,2.*Ti+.07))
      call scalewn(1.,rcc(NRUSED),0.,tpmax,.false.,.false.)
      call axptset(1.,0.)
      call ticrev()
      if(.not.ldist) call winset(.false.)
      call altyaxis(0.,1.)
      call ticrev()
c      call axlabels('',' T!di!d')
      call axptset(0.,0.)
      call winset(.true.)
      call color(1)
      call polyline(rcc(1),temprplot,NRUSED)
      call jdrwstr(wx2nx(r(nr)),wy2ny(0.8*Ti),' T!dr!d',-2.)
      call color(3)
      call polyline(rcc(1),temptpplot,NRUSED) 
      call jdrwstr(wx2nx(r(nr)),wy2ny(1.2*Ti),' T!d!Aq!@!d',-2.)
      call dashset(2)
      call polyline(rcc(1),ttpsmooth,NRUSED)
      call color(1)
      call polyline(rcc(1),trsmooth,NRUSED)
      call color(15)
      call dashset(0)
      call fwrite(tpmin,iwidth,4,char1)
      call jdrwstr(wx2nx(r(nr)),wy2ny(0.3*tpmax),
     $     ' T!d!Aq!@min!d='//char1(1:iwidth),-2.)
      call winset(.false.)
c      call pltend()
      end
c***********************************************************************
c Contouring of the charge density, rho, on distorted mesh.
      subroutine rhodisplay(ir,rhomin,rhomax)
      integer ir
c Common data:
      include 'piccom.f'
c      save
      character cworka(nrsize*nthsize+13)
      integer ncont
      parameter (ncont=9)
      real zclv(ncont)
      real xrho(nrsize+1,0:nthsize),zrho(nrsize+1,0:nthsize)
      logical lfirst
      data lfirst/.true./
      save xrho,zrho,lfirst,first
      
      rmi=rhomin
      rma=rhomax
      do i=1,NRUSED
         do j=1,NTHUSED
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
         if(.not.LCIC)then
            zrho(i,0)=rcc(i)
            xrho(i,0)=0.
            zrho(i,nth)=-rcc(i)
            xrho(i,nth)=0.
         endif
      enddo

      if(rma.eq.0)then
         call minmax(rho(1,1),NRUSED,rmi,rma)
         call minmax(rho(1,NTHUSED),NRUSED,rn2,rx2)
         rmi=min(rmi,rn2)
         rma=max(rma,rx2)
         if(rma.lt.1.2)rma=1.
         if(rma.gt.5.)then
            rma=5.
            rmi=0.
         endif
         if(Ezext.ne.0. .and. rma.gt.2.)then
            rma=2.
            rmi=0.
         endif
         if(rmi.lt.0.7)rmi=0.
      endif
      call fitrange(rmi,rma,ncont-1,ipow,fac10,delta,first,xlast)
c      write(*,*)'rhomin/max=',rmi,rma,first,xlast,delta
      ncl=0
      do j=1,ncont
         zclv(j)=(j-1)*delta + first
         if(zclv(j).gt. xlast*1.0001) goto 440
         ncl=ncl+1
      enddo
 440  continue
      icl=-ncl
      call pltinaspect(-r(nr),r(nr),0.,r(nr))
      if(lfirst)then
         call accisgradinit(-50000,10000,-10000,130000,65000,65000)
c         call accisgradinit(-25000,00000,25000,130000,65000,130000)
         lfirst=.false.
      endif
c     There are signs that the following is overflowing bounds somehow.
      ntype=2+16+32
c      call contourl(rho(1,1),cworka,nr,nr,nth,zclv,icl,xrho,zrho,ntype)
      if(LCIC)then
         call contourl(rho(1,1),cworka,nrsize+1,NRUSED,nth,
     $        zclv,icl,zrho(1,1),xrho(1,1),ntype)
      else
         call contourl(rho(1,0),cworka,nrsize+1,NRUSED,nth+1,
     $        zclv,icl,zrho,xrho,ntype)
      endif
c      call color(15)
      call gradlegend(zclv(1),zclv(ncl),
     $     .1,1.25,.9,1.25,0.02,.true.)
      call axis()
      call axlabels('z     density color, velocity arrows','r sin!Aq!@')
c      call legendline(.1,-.3,258,'density color-contours')
      call color(12)
      if(ir.le.0.or.ir.ge.100) ir=10
      j0=1
      if(Lcic)j0=2
      do j=j0,NTHUSED-1,NTHUSED/5+1
         do i=1,NRUSED,max(nr/ir,1)
            vri=vrsum(i,j)/(psum(i,j)+1.e-5)
            vti=vtsum(i,j)/(psum(i,j)+1.e-5)
            size=0.02*sqrt(vri**2+vti**2)
            angle=atan2(vti,vri)+acos(tcc(j))
            call charsize(size,0.5*size)
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $           char(ichar('_')+128)//char(0),0.)
         enddo
      enddo
      call charangl(0.)
      call charsize(0.,0.)
      call color(15)
c      call pltend()
      end
c********************************************************************
      subroutine vdiagsout(lbig)
      logical lbig
c Plot relaunched distribution. nvdiag calculated in padvnc. 
c Common storage:
      include 'piccom.f'
      call autoinit(vdiag,nvdiag,nvmax)
      if(lbig) call charsize(0.02,0.02)
      call axis()
      call polyline(vdiag,nvdiag,nvmax)
      call legendline(.7,.9,258,'!Bf!@!dr!d(r!dmax!d)')
      call axlabels('velocity','')
      if(lbig) call charsize(0.0,0.0)
      end
c********************************************************************
      subroutine vdiagsin(lbig)
      logical lbig
c Plot distributions at probe edge, calculated in padvnc. 
c Common storage:
      include 'piccom.f'
      character*40 rval,tval
c      call autoinit(vdiag,vrdiagin,nvmax)
c Do this instead to ensure both plots fit.
      call minmax(vdiag,nvmax,xmin,xmax)
      call fitrange(xmin,xmax,6,nxfac,xfac,xdelta,xmin,xmax)
      call minmax(vrdiagin,nvmax,ymin,ymax)
      call minmax(vtdiagin,nvmax,ytmin,ytmax)
      if(ytmax.gt.ymax) ymax=ytmax
      if(ytmin.lt.ymin) ymin=ytmin
      call fitrange(ymin,ymax,6,nxfac,xfac,xdelta,ymin,ymax)
      call pltinit(xmin,xmax,ymin,ymax)

      call polyline(vdiag,vrdiagin,nvmax)
      if(lbig) call charsize(0.02,0.02)
      call axis()
      call axlabels('velocity','')
      call fwrite(rcc(ircell),irwidth,1,rval)
c      call fwrite(acos(tcc(itcell))*180./3.14159,itwidth,1,tval)
      call fwrite(tcc(itcell),itwidth,2,tval)
      call legendline(.8,.85,258,'!Bf!@!dr!d')
      call legendline(.05,.95,258,'Position r='//rval(1:irwidth)//
     $     ', cos!Aq!@='//tval(1:itwidth)//':')
      call color(3)
      call polyline(vdiag,vtdiagin,nvmax)
      call legendline(.8,.75,258,'!Bf!@!d!Aq!@!d')
      call color(15)
      if(lbig) call charsize(0.02,0.02)
      end
c*********************************************************************
      subroutine slices(jstepth,rhomin,rhomax,time)
c
      include 'piccom.f'
      character charin*100,ch2*100
      real rccleft(nrsize),phip1(nrsize)
      logical sfirst
      data sfirst/.true./
      save rccleft,sfirst

      if(sfirst)then
         do j=1,NRUSED
            rccleft(j)=-rcc(j)
         enddo
         sfirst=.false.
      endif

      jstepth=NTHUSED/2-1
      if(rhomax.eq.0.)then
         call minmax2(rho(1,1),NRFULL,NRUSED,NTHUSED,rhomin,rhopm)
         rhopm=max(1.2,.1*nint(10.*rhopm))
      else
         rhopm=rhomax
      endif
      if(rhopm.gt.4.)then
         phi0rho=rhopm*.9
      else
         phi0rho=1.2
      endif
      if(vprobe.ne.0.) then
         phiscale=phi0rho*abs(1./vprobe)
      else
         phiscale=phi0rho*0.25
      endif
      call pltinit(-rcc(NRUSED),rcc(NRUSED),
     $     -phi0rho/phiscale,(rhopm-phi0rho)/phiscale)
      call axptset(1.,1.)
      call jdrwstr(wx2nx(rcc(NRUSED)),wy2ny(-0.05/phiscale),
     $     '!Af!@',-1.5)
      call ticrev()
      call yaxis(0.,0.)
      call ticrev()
      call axptset(0.,0.)
      call scalewn(-rcc(NRUSED),rcc(NRUSED),0.,rhopm,
     $     .false.,.false.)
      call axis()
      call boxtitle('Density and Potential')
      call axlabels('Radial position','n/n!A!d;!d!@')
c     call polyline(rcc,rhoright,NRUSED)
c     call polyline(rccleft,rholeft,NRUSED)
      do j=1,NTHUSED/2,jstepth
         call color(mod(j,15)+1)
         call winset(.true.)
         call polyline(rccleft,rho(1,NTHUSED-j+1),NRUSED)
         call polyline(rcc(1),rho(1,j),NRUSED)
         call dashset(2)
         do k=1,NRUSED
            phip1(k)=phiscale*phi(k,NTHUSED-j+1)+phi0rho
         enddo
         call polyline(rccleft,phip1,NRUSED)
         do k=1,NRUSED
            phip1(k)=phiscale*phi(k,j)+phi0rho
         enddo
         call polyline(rcc(1),phip1,NRUSED)
         call dashset(0)
         call winset(.false.)
         write(charin,'(f4.0)')thang(j)*180./3.1415927
         call legendline(0.7,0.08*(j-1)/jstepth+.1,
     $        0,charin(1:4)//char(0))
      enddo
      call color(15)
      call vecw(-rcc(NRUSED),1.,0)
      call vecw(rcc(NRUSED),1.,1)
      call vecw(-rcc(NRUSED),phi0rho,0)
      call vecw(rcc(NRUSED),phi0rho,1)
      call fwrite(vd,iwdth,2,charin)
      call fwrite(time,iw2,2,ch2)
      call jdrwstr(0.03,0.37,'v!dd!d='//charin(1:iwdth)//' t='
     $     //ch2(1:iw2),1.)

      end
c********************************************************************
      real function smaxflux(uc,chi)
c     Return the total flux to a unit radius sphere from a unit density
c     maxwellian distribution shifted by velocity
      real uc
c     normalized to sqrt(2T/m), in a spherically symmetric potential
c     having a value on the sphere normalized to Ti of minus
      real chi

      real eps,pi
      data eps/1.e-3/pi/3.1415927/


      erf=1.-erfcc(uc)
      sqpi=sqrt(pi)
      if(abs(uc).lt.eps) then
         erfbyu=(2./sqpi)*(1.-uc**2 /3.)
      else
         erfbyu=erf/uc
      endif

      smaxflux=pi*(uc*erf + (0.5+chi)*erfbyu + exp(-uc**2)/sqpi)
      end
c*******************************************************************
      real function smaxflux2(chi)
c     Return the total flux to a unit radius sphere from a unit density
c     slightly shifted maxwellian,
c     normalized to sqrt(2T/m)=vtion, in a low potential non necessary
c     symmetric, with whatever B field, having a value on the sphere
c     normalized to Ti of minus
      real chi
      real pi
      data pi/3.1415927/
      sqpi=sqrt(pi)

      smaxflux2=sqpi*2*(1+chi)

      end





c*******************************************************************

c Overplot orbits on existing plot.
      subroutine plotorbits
      include 'piccom.f'
      call winset(.true.)
      do k=1,norbits
         call color(7)
         call polyline(zorbit(1,k),rorbit(1,k),iorbitlen(k))
         call color(15)
         call charsize(0.01,0.01)
         call accircle(wx2nx(zorbit(iorbitlen(k),k)),
     $        wy2ny(rorbit(iorbitlen(k),k)))
         call charsize(0.0,0.0)
      enddo
      call winset(.false.)
      end
