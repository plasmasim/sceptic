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
c Version of Mar 07 includes icolntypes with bit 3 set to do all
c collisions at the end of the step.
c Version 2.5; Jan 2005: subcycling of padvnc.
c Version 2.5; Jan 2005: fixed reinjection flux option.
c Advance the particles
      subroutine padvnc(dtin,icolntype,step)
      integer step
      real dtin
c Common data:
      include 'piccom.f'
      include 'colncom.f'

      real accel(3)
      real dt
c moved to piccom.f      logical lsubcycle
      real cosomdt,sinomdt
c temp data:
      real temp
      integer idum

c Data for the domain sub
      real rn0

c Place to read/write on the add particule injection
      integer nread
      
      nread=mod(step,addhist)+1
      idum=1
c Reset the number of particles that enter the inner domain
      xpstonum(nread)=0
      rsp=r(rsplit)
      rn0=1.
c Xp is the three x-coordinates followed by the 3 v coordinates.
c Use a leapfrog scheme, so interpret the v-coords as half a step
c behind the x-coords. 
c If lsubcycle, use multiple fractional steps near inner boundary.
      dt=dtin
      rp2=r(1)**2
c Zero the sums now these are assigned here.
      nrein=0
      nreintry=0
      spotrein=0.
      ninner=0
      fluxrein=0.
      ntrapre=0
      zmomprobe=0.
      zmout=0.
      iocthis=0.
      do j=1,nth
         ninth(j)=0
      enddo
      do i=0,NRFULL
         do j=0,NTHFULL
            psum(i,j)=0.
            vrsum(i,j)=0.
            vtsum(i,j)=0.
            vpsum(i,j)=0.
            v2sum(i,j)=0.
            vr2sum(i,j)=0.
            vtp2sum(i,j)=0.
            vzsum(i,j)=0.
         enddo
      enddo
c      if(lfixedn)then
cc  Fixed number of particles. The ninjcomp is already=npartmax.
c         ipmax=npart
c      else
cc Fixed flux. Rely on ninjcomp to terminate the particle slot treatment.
c         ipmax=npartmax
c      endif
c      write(*,*)'Starting cycle',ipmax,ninjcomp


c End of setup. Start Cycling through particles.
      if (dsub) then
         ido=npart+npartadd
      else
         ido=npart
      endif


      do i=1,ido
         if(ipf(i).gt.0) then
c     Is this an active slot?
c     Find the mesh position and the trigonometry.
c     Here we do need half quantities.
            ih=1
            hf=88.
 
            call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $           ,zetap,ih,hf)

c  Now we know where we are in radius rp. 
c  We decide the level of subcycling.
            if(lsubcycle) then
               isubcycle=r(nrfull)/rp
c     if(mod(i,1000).eq.0) write(*,*)isubcycle
               dt=dtin/isubcycle
            else
               isubcycle=1
            endif

            do 81 ic=1,isubcycle

c     Except for the first time, find new position.
               if(ic.ne.1)then 
                  ih=1
                  hf=77.
                  call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)
               
               endif
               call getaccel(i,accel,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)
c getaccel returns the accel based on the charge-field calculation.
c We then add on the acceleration due to the neutral-collisions-implied
c electric field.
c Trying adding on the Eneutral partially at the end.
               accel(3)=accel(3)+Eneutral*.5

c     write(*,501)accel,(xp(j,i),j=1,3)

               if(dsub) then
                  rn2=0.
                  do j=1,3
                     rn2=rn2+xp(j,i)**2
                  enddo
                  rn0=sqrt(rn2)
               endif
               
c     Fixings for the subcycling
               if(dtprec(i).eq.0.)dtprec(i)=dt
               dtnow=0.5*(dt+dtprec(i))
c     Parameters for the Lorentz force
               cosomdt=cos(Bz*dtnow)
               sinomdt=sin(Bz*dtnow)

c First half of velocity advance:    AccelPhi/2+AccelBz+AccelPhi/2
               do j=4,6
                  xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
               enddo
c B-field rotation
               temp=xp(4,i)
               xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
               xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt
c Second half of velocity advance
               do j=4,6
                  xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
               enddo
               dtprec(i)=dt
               rn2=0.
               xdv=0.
               v2=0.
c Position advance (and accumulate coordinate terms).
               do j=1,3
                  xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                  rn2=rn2+xp(j,i)**2
                  xdv=xdv+xp(j,i)*xp(j+3,i)
                  v2=v2+xp(j+3,i)**2
               enddo
               tm=xdv/v2
               rn=sqrt(rn2)
c     Test if we went through the probe and came back out.
               if((0..lt.tm .and. tm.lt.dt .and.
     $              (rn2 - tm**2/v2).lt.rp2))then
                  if(rn.gt.r(1))then
c     write(*,*)'Through probe',tm,(rn2 - tm**2/v2)
                     rn=0.
                  endif
               endif
c Trying adding on the Eneutral partially at the end.
               xp(6,i)=xp(6,i)+Eneutral*dt*0.5
c     Handling boundaries for 'real particles' :
            if(i.le.npart) then
               if(rn.le.r(1)) then
                  ninner=ninner+1
                  nrealin=nrealin-1
c     Solve for sphere crossing step fraction, s.
                  a=0.
                  b=0.
                  c=0.
                  do j=1,3
                     a=a+(dt*xp(j+3,i))**2
                     b=b-2.*xp(j,i)*(dt*xp(j+3,i))
                     c=c+xp(j,i)**2
                  enddo
                  c=c-r(1)**2
                  s=(-b+sqrt(b**2-4.*a*c))/(2.*a)
                  xc=xp(1,i)-s*dt*xp(4,i)
                  yc=xp(2,i)-s*dt*xp(5,i)
                  zc=xp(3,i)-s*dt*xp(6,i)
                  ctc=zc/sqrt(xc**2+yc**2+zc**2)
c     Interpolate onto the theta mesh as in ptomesh               
                  ithc=interpth(ctc,tfc)
c     
                  if(LCIC)then
                     icell=nint(ithc+tfc)
                  else
                     icell=ithc
                  endif
                  ninth(icell)=ninth(icell)+1
                  zmomprobe=zmomprobe+xp(6,i)
               elseif(rn.ge.r(nr))then
                  zmout=zmout-xp(6,i)
c     We did not leave the grid inside.
               elseif(dsub) then
c     Check if the particle entered the inner domain, and if yes
c     store it. We also update nrealin
                  if ((rn0.ge.rsp) .and. (rn.lt.rsp)) then
                     xpstonum(nread)=xpstonum(nread)+1
                     do k=1,6
                        xpstorage(k,xpstonum(nread),nread)=xp(k,i)
                     enddo  
                     nrealin=nrealin+1
                  elseif ((rn.gt.rsp) .and. (rn0.le.rsp)) then
                     nrealin=nrealin-1
                  endif
                  goto 81
               else
                  goto 81
               endif
c We left. 
c If we haven't exhausted complement, restart particle i.
               if(nrein.lt.ninjcomp) then
                  call reinject(i,dtin,icolntype,bcr)
                  ipf(i)=1
                  zmout=zmout+xp(6,i)
                  if(i.le.norbits) then
                     if (.not.(orbinit))
     $                    iorbitlen(i)=0
                  endif
               else
                  ipf(i)=0
c                  if(i.gt.190000) write(*,*)'Leaving empty slot',i
               endif
c     Break from subcycles.
               goto 82
c     Now we care about the add particles. If they leave the domain, we
c     just reinject them, since we don't use them for diagnostics
c     as fluxes
            else
               if((rn.le.r(1)).or.(rn.ge.rsp)) then
                  a=nint(ran0(idum)*(addhist-1))+1
                  b=nint(ran0(idum)*(xpstonum(a)-1))+1
                  do k=1,6
                     xp(k,i)=xpstorage(k,b,a)
                  enddo
               endif
            endif
 81         continue
 82         continue
            if(ldist) then
               rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
c     Diagnostics of f_r(rmax):
               if(rn.gt.r(nr-1))then
                  v=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i)+xp(6,i)*xp(3,i))/rn
                  ivdiag=1+max(0,nint(nvmax*(v/vrange + .499)))
                  if(ivdiag.gt.nvmax) ivdiag=nvmax
                  nvdiag(ivdiag)=nvdiag(ivdiag)+1
               elseif(rn.gt.r(ircell).and.rn.le.r(ircell+1))then
c     Inner distribution Diagnostics: Assumes reinject never gets here.
                  ctc=xp(3,i)/rn
                  ithc=interpth(ctc,thc)
                  if(ithc.eq.itcell)then
                     vz=xp(6,i)
                     vxy=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i))/
     $                    sqrt(xp(1,i)**2+ xp(2,i)**2)
                     vr=vz*ct+vxy*st
                     vt=-vz*st+vxy*ct
c     Radial
                     ivdiag=1+max(0,nint(nvmax*(vr/vrange + .499)))
                     if(ivdiag.gt.nvmax) ivdiag=nvmax
                     vrdiagin(ivdiag)=vrdiagin(ivdiag)+1
c     Angular
                     ivdiag=1+max(0,nint(nvmax*(vt/vrange + .499)))
                     if(ivdiag.gt.nvmax) ivdiag=nvmax
                     vtdiagin(ivdiag)=vtdiagin(ivdiag)+1
c     write(*,502)rn,ithc,vr
                  endif
               endif
            endif
c     Orbit diagnostics
            if(i.le.norbits) then
               iorbitlen(i)=iorbitlen(i)+1
               xorbit(iorbitlen(i),i)=xp(1,i)
               yorbit(iorbitlen(i),i)=xp(2,i)
               rorbit(iorbitlen(i),i)=sqrt(xp(1,i)**2+xp(2,i)**2)
               zorbit(iorbitlen(i),i)=xp(3,i)
               vxorbit(iorbitlen(i),i)=xp(4,i)
               vyorbit(iorbitlen(i),i)=xp(5,i)
               vzorbit(iorbitlen(i),i)=xp(6,i)
c     write(*,503)i,iorbitlen(i),xorbit(iorbitlen(i),i)
c     $           ,yorbit(iorbitlen(i),i),zorbit(iorbitlen(i),i)
c     $           ,rorbit(iorbitlen(i),i)
            endif
            if(ipf(i).gt.0)iocthis=i
         elseif(nrein.lt.ninjcomp)then
c Case for ipf(i) le 0 (empty slot) but still wanting to inject. 
c We should not come here unless .not.lfixedn.
c            write(*,*)'Reinjecting empty slot',i
            call reinject(i,dtin,icolntype,bcr)
            ipf(i)=1
            iocthis=i
         elseif(i.ge.iocprev)then
c Break if inactive slot and we have exhausted the complement of injections.
c And we have reached the maximum occupied slot of previous run.
            goto 401
         endif
      enddo
 401  continue
c We just want the diagnostics with the true particles for now
c      iocthis=min(iocthis,npartmax)
      iocprev=iocthis
c      if(.not.lfixedn)write(*,504)ninjcomp,nrein,i,iocprev
 504  format('  ninjcomp=',i6,'  nrein=',i6,'  i=',i6,
     $     '  iocprev=',i6)
 503  format('Orbit',i3,' length=',i5,' position=',4f7.3)
 502  format('Distrib. rn=',f6.3,' ithc=',i4,' vr=',f6.3)
 501  format('accel=',3f11.4,' xp=',3f11.4)

      end
c***********************************************************************
c Version using precalculated functions. About 30% faster.
      subroutine ptomesh(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp
     $     ,zetap,ih,hf)
c Return the left hand mesh point and the fractional mesh distance of the
c position of particle i, in irl,rf,itl,tf,ipl,pf 
c Return the sines and cosines of theta and phi in st,ct,sp,cp
c If ih.ne.0 on entry, calculate the half-mesh postion, zetap,ih,hf.
      implicit none
      integer i
      integer irl,ithl,ipl
      real rf,thf,pf
      real ct,st,cp,sp,rp
      real zetap,hf
      integer ih
      real hp
c      common /angles/ct,st,cp,sp,rp
c Common data:
      include 'piccom.f'
      real rsp
      real x,y,z
      external interpth
      integer interpth

c Testing
      if(.not. xp(1,i).le.400.)then
         write(*,*)'Ptomesh particle overflow on entry'
         write(*,*)i,(xp(ipl,i),ipl=1,6)
         write(*,*)i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp,zetap,ih,hf
         stop
      endif

C Find the cell and cell fraction we are at.
      x=xp(1,i)
      y=xp(2,i)
      z=xp(3,i)

      rsp=x**2+y**2
c The square roots here cost perhaps 1/3 of this routine. 
      rp=sqrt(rsp+z**2)
c 
      if(.not. rp.le.r(nr))then
         write(*,*)'Ptomesh particle outside on entry'
         write(*,*)'xp:',(xp(ipl,i),ipl=1,6)
         write(*,*)'i,r(nr),rp,zetap,ih,hf'
         write(*,*)i,r(nr),rp,zetap,ih,hf
         write(*,*) 'x:',x,' y:',y,' z:',z
         stop
      endif
      rsp=sqrt(rsp)
      if(rsp .gt. 1.e-9) then
         cp=x/rsp
         sp=y/rsp
      else
         cp=1.
         sp=0.
      endif
      st=rsp/rp
      ct=z/rp
c      pf=atan2(sp,cp)
      pf=1
c Not using pf at present.
c      ipl=pf*dpinv +1
c      pf=pf*dpinv -ipl
c
      if(abs(1+int((ct-th(1))*tfac)).gt.ntpre)then
         write(*,*)'ptomesh overflow. Probably particle NAN'
         write(*,*)'i,irl,rf,ithl,thf',i,irl,rf,ithl,thf
         write(*,*)'ct,th(1),tfac,z,rp',ct,th(1),tfac,z,rp
         write(*,*)'xp',xp(1,i),xp(2,i),xp(3,i),xp(4,i),xp(5,i),xp(6,i)
         write(*,*)'x,y,z',x,y,z
         stop
      endif
      ithl=interpth(ct,thf)

      irl=irpre(1+int((rp-r(1))*rfac))
      rf=(rp-r(irl))/(r(irl+1)-r(irl))
      if(rf.lt.0.) then
         write(*,*)"Negative rf from irpre. i,ih,irl,rf,rp="
     $        ,i,ih,irl,rf,rp
         write(*,*) 'r: ',xp(1,i)**2+xp(2,i)**2+xp(3,i)**2
      endif
c "While not"      
 402  if(rf.le.1.)goto 401
      if(irl.eq.nr)then
         write(*,*)'ptomesh rf gt 1 error:',rf,irl
         stop
      else
         irl=irl+1
         rf=(rp-r(irl))/(r(irl+1)-r(irl))
      endif
      goto 402
 401  continue
c      return
c New section for halfmesh quantities. Adds about 10% to time.
c Now we have identified the whole mesh position. The half mesh is very
c near it, either irl or irl+1.
      if(ih.ne.0)then
         ih=irl+1
         hp=rp-r(1)
         zetap=sqrt(2.*hp)
         hf=zetap-zetahalf(ih)
         if(hf.lt.0.)ih=ih-1
c     This is the halfmesh fraction 'x'
         hf=(zetap-zetahalf(ih))/(zetahalf(ih+1)-zetahalf(ih))
         
         if(hf.gt.1.or.hf.lt.0.or.zetap.lt.0..or.ih.le.0
     $        )then
c     $        .or. ih.eq.NRFULL)then
            write(*,*)'hf error, ih,irl,rf,zetahalf',ih,irl,rf,
     $           zetahalf(ih),zetahalf(ih+1)
            write(*,*)'zetap,zeta(ih),zeta(ih+1),hf',
     $           zetap,zeta(ih),zeta(ih+1),hf
         endif
      endif
      end
c***********************************************************************
c Calculate potential phi from rho.
      subroutine fcalc_lambda(dt,icolntype,colnwt)
      real dt
c Common data:
      include 'piccom.f'
      real relax
      real phislopeconst(nth),phislopefac(nth)
c Chebychev acceleration. Wild guess at the Jacoby convergence radius.
      rjac=1.-4./max(10,NRUSED)**2
      omega=1.
      maxits=2.5*NRUSED
      dconverge=1.e-5
c cic boundary is at i=1, ngp at 0+(1/2) (sort of).
      if(LCIC)then
         imin=1
      else
         imin=0
      endif

c Do SOR iteration with boundary set at probe potential.
c      do j=1,NTHUSED
c         phi(imin,j)=vprobe
c      enddo
      call innerbc(imin,dt)
c Outer boundary
c Rindex is - the factor multiplying phi.
c dp/dr = - rindex p/r + blfac
c so p_N = p_(N-1)*psislopefac + psislopeconst
c Simplistic inverse square BC.
c      rindex=2.
c Simplistic inverse r BC.
c      rindex=1.
c      phislopefac=(redge-delredge*rindex*0.5)/
c     $     (redge+delredge*rindex*0.5)
c      phislopeconst=0.
c
      redge= (rcc(NRFULL)+rcc(NRFULL-1))*0.5
      delredge=rcc(NRFULL)-rcc(NRFULL-1)
c Screening k-number combines electrons and ions.
      if(debyelen.gt. 1.e-10) then
         el2=(1.+1./Ti)/debyelen**2
      else
 2       el2=2.e20
      endif
      el=sqrt(el2)
      afactor=0.02
      alpha=1./(1.+(afactor*redge/debyelen)**2)
      rxl=el*redge
      expE1=(alog(1.+1./rxl) - 0.56/(1.+4.1*rxl+0.9*rxl**2))
      rindex=alpha*(redge*el+1.)+ (1.-alpha)*2.
c At high collisionality reduce the debye gradient term
      if(icolntype.eq.1 .or. icolntype.eq.2)then
         rindex=(rindex-1.)/(1.+(colnwt*redge)**2/Ti)+1.
      endif
c      if(icolntype.eq.2)then
c Remove the deficit term
c         expE1=0.
c Simplistic 1/r trial.
c         rindex=1.
c      endif
      adeficit=0
c Boundary slope factor calculations:
      do j=1,NTHUSED
c Current fractional ion deficit due to collection.
c Coefficient of 1/r^2 in modified shielding equation is
c a = deficitj * r_edge^2 / \lambda_De^2
         deficitj=1-phi(NRUSED,j)/Ti -rho(NRUSED,j)
c         write(*,*)rho(NRUSED,j),phi(NRUSED,j),deficitj
         blfac1=(deficitj/debyelen**2) * redge
         adeficit=adeficit+blfac1
c BC modification is (a/r_edge)[exp(EL*r) E_1(El*r)] given by approx.
         blfac=blfac1*expE1
         blfac=alpha*blfac
         phislopeconst(j)=blfac*redge*delredge/
     $        (redge+delredge*rindex*0.5)
         phislopefac(j)=(redge-delredge*rindex*0.5)/
     $        (redge+delredge*rindex*0.5)
      enddo
c Actual a factor averaged over angles:
      adeficit=adeficit*redge/NTHUSED
      if(adeficit.lt.0.)then
c         write(*,*)'Negative adeficit',adeficit,' set to zero'
         adeficit=0.
      endif
c SOR iteration.
      do k=1,maxits
c Use over-relaxation if debyelen is large, or straight newton otherwise.
         relax=(omega*debyelen**2+1.)/(debyelen**2+1.)
         deltamax=0.
c Alternate iteration directions
         if(mod(k/2,2).eq.0)then
         do j=1,NTHUSED
            do i=imin+1,NRFULL-1
               expphi=exp(phi(i,j))
               dnum= apc(i)*phi(i+1,j)+bpc(i)*phi(i-1,j) + cpc(i,j)
     $              *phi(i,j+1)+dpc(i,j)*phi(i,j-1) -fpc(i,j)*phi(i,j)
     $              + rho(i,j) - expphi
               dden=fpc(i,j) + expphi
               delta=relax*dnum/dden
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               phi(i,j)=phi(i,j)+delta
            enddo
c     Outer boundary.
            if(Ezext.eq.0)then
               delta=phi(NRFULL-1,j)*phislopefac(j)-phislopeconst(j)
     $              -phi(NRFULL,j)
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               phi(NRFULL,j)=phi(NRFULL,j)+relax*delta
            else
               phi(NRFULL,j)=Ezext*tcc(j)*r(NRFULL)
            endif
         enddo
         else
         do j=NTHUSED,1,-1
            do i=NRFULL-1,imin+1,-1
               expphi=exp(phi(i,j))
               dnum= apc(i)*phi(i+1,j)+bpc(i)*phi(i-1,j) + cpc(i,j)
     $              *phi(i,j+1)+dpc(i,j)*phi(i,j-1) -fpc(i,j)*phi(i,j)
     $              + rho(i,j) - expphi
               dden=fpc(i,j) + expphi
               delta=relax*dnum/dden
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               phi(i,j)=phi(i,j)+delta
            enddo
c     Outer boundary.
            if(Ezext.eq.0)then
               delta=phi(NRFULL-1,j)*phislopefac(j)-phislopeconst(j)
     $              -phi(NRFULL,j)
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               phi(NRFULL,j)=phi(NRFULL,j)+relax*delta
            else
               phi(NRFULL,j)=Ezext*tcc(j)*r(NRFULL)
            endif
         enddo
         endif
         if(abs(deltamax).lt.dconverge.and.k.ge.2)goto 11
         if(k.eq.1)then
            omega=1./(1.-0.5*rjac**2)
         else
            omega=1./(1.-0.25*rjac**2*omega)
        endif
      enddo
c      write(*,*)'SOR not converged. deltamax=',deltamax
 11   continue
      write(*,'('':'',i3,$)')k
c      write(*,201)k,deltamax,relax
 201  format(' SOR iteration',I4,' delta:',f10.6,' relax=',f8.4)
c Calculate electric force on probe. Moved to main.
c Inner Boundary values
      do j=1,NTHUSED
         phi(0,j)=2.*phi(imin,j)-phi(imin+1,j)
      enddo
      do i=1,NRUSED
         phi(i,0)=phi(i,imin+1)
         phi(i,NTHUSED+1)=phi(i,NTHUSED-imin)
      enddo
c      write(*,*)'phi(rmax)=',phi(NRFULL,NTHUSED/2)
      end
c*******************************************************************
      subroutine innerbc(imin,dt)
      include 'piccom.f'
      real flogfac
      real fluxofangle(nthsize)
      if(linsulate.or.lfloat) then
         flogfac=0.5*alog(2.*pi/(rmtoz*1837.))
         totflux=0.
         do j=1,nthused
            if(lcic)then
               fluxofangle(j)=finthave(j)*(nthused-1.)/
     $              (4.*pi*rhoinf*dt*r(1)**2)
               if(j.eq.1 .or. j.eq.nthused)
     $              fluxofangle(j)=fluxofangle(j)*2.
            else
               fluxofangle(j)=finthave(j)*(nthused)/
     $              (4.*pi*rhoinf*dt*r(1)**2)
            endif
            totflux=totflux+fluxofangle(j)
            if(linsulate)then
               if(fluxofangle(j).gt.0.)then
                  phi(imin,j)=alog(fluxofangle(j))+flogfac
               else
                  phi(imin,j)=vprobe
               endif
            endif
         enddo
         if(totflux.gt.0.)then
            vprobe=alog(totflux/nthused)+flogfac
         endif
         if(lfloat)then
            do j=1,nthused
               phi(imin,j)=vprobe
            enddo
         endif
c         write(*,*)
c         write(*,*)'fluxofangle=',(fluxofangle(j),j=1,NTHUSED)
c         write(*,*)'phi=',(phi(imin,j),j=1,NTHUSED)
      else
         do j=1,NTHUSED
            phi(imin,j)=vprobe+Ezext*tcc(j)
         enddo     
      endif
      end
c***********************************************************************
c Initialization for Collisions
      subroutine colninit(colnwt,icolntype)
      include 'piccom.f'
      include 'colncom.f'
c      write(*,*)'Initialized collisions',colnwt,icolntype
      if(icolntype.eq.1 .or. icolntype.eq.2
     $     .or. icolntype.eq.5 .or. icolntype.eq.6)then
c Constant nu collisions. The Eneutral must be consistent with vd:
         Eneutral=colnwt*(vd-vneutral)
c Testing
c         Eneutral=0.
         if(myid .eq.0) write(*,*)'colnwt,vd,vneutral=',colnwt,vd
     $        ,vneutral,' Eneutral=',Eneutral
     $        
      elseif(icolntype.eq.0)then
c Need more code here for other types. Not yet implemented.
c Must set Eneutral to zero by default.
         Eneutral=0.
      else
         write(*,*)'Incorrect icolntype',icolntype
         stop
      endif
      end
c***********************************************************************
c Master Collision subroutine.
      subroutine collide(dt,colnwt,icolntype)
      real dt,colnwt
      integer icolntype
c Call the appropriate collision routine.
      if(icolntype.eq.1 .or. icolntype.eq.2
     $     .or. icolntype.eq.5 .or. icolntype.eq.6)then
         call nucollide(dt,colnwt,icolntype)
      elseif(icolntype.gt.8)then
         call mfpcollide(dt,colnwt,icolntype)
      else
         write(*,*)'Incorrect collision type:',icolntype
         stop
      endif
      end
c******************************************************************
      subroutine nucollide(dt,cnu,icolntype)
c Collision type with specified collision frequency.
c Input dt=timestep, colnwt=cnu=collision freq (\propto target density)
c icolntype bit 3 (4) if set says do all collisions at step end.
c That is less accurate for treating multiple collisions but should
c avoid errors at the boundary from reinjection etc.
      include 'piccom.f'
      include 'colncom.f'
      real accel(3)

c Don't attempt if collision freq is zero.
      if(cnu.le.0.) return
      ncollide=0
      idum=1
      tisq=sqrt(Ti)

c     Here we need to invert a poisson distribution to tell if we had
c     a collision. The cumulative poisson distribution is 
c         P(<t)=1-exp(-nu.t)
c     Consequently, the time after the start of the last step,
c     of a collision is t=-(1/nu)ln(1-y) for a random number y, 
c     and if this is less than dt we collided.
c     If nu.dt is small, this is approximately t=y/nu, but if not, we
c     still get a correct answer with the full solution, by repeating
c     until we have used up the whole time-step.
c     When nu.dt is very small we gain efficiency by only treating a
c     subset of the particles. And multiplying the nu by icycle.
      icycle=1./(20.*cnu*dt)
c     We keep icycle.nu.dt < 1/20  to avoid bias by the cycle.
      if(.not.icycle.ge.1) icycle=1
      ichoose=ran0(idum)*icycle+1
      if(ichoose.gt.icycle)ichoose=icycle
      ytestfull=1.-exp(-cnu*icycle*dt)
c      write(*,*)'nucollide: icycle,ichse,dt,ytstf='
c     $     ,icycle,ichoose,dt,ytestfull
c Cycle through all the particles, skipping if necessary.
      do i=ichoose,iocprev,icycle
         if(ipf(i).gt.0)then
c This is an active particle slot.
            dtd=dt
            ytest=ytestfull 
            if(mod(icolntype/4,2).eq.1)then
c Simple: all collisions at end of step.
               y=ran0(idum)
               if(y.lt.ytest)then
c Get new velocity; reflects neutral maxwellian shifted by vneutral.
                  xp(4,i)=tisq*gasdev(idum)
                  xp(5,i)=tisq*gasdev(idum)
                  xp(6,i)=tisq*gasdev(idum)+ vneutral
               endif
            else
c Start of multiple collision loop           
 1             y=ran0(idum)
               if(y.lt.ytest)then
c A collision occured at a time dtc during the last [partial] step
                  ncollide=ncollide+1
c at a time dtc after its start.
                  dtc=-alog(1.-y)/(cnu*icycle)
c Adjust the step duration to that remaining after collision.
                  dtd=dtd-dtc
c               write(*,*)i,ytest,y,ytest,dtc,dtd
c Calculate the probability of a collision in remaining partial step.
                  ytest=1.-exp(-cnu*icycle*dtd)
c Get the current acceleration, needed below.
c We don't bother about the place we calculate this being correct,
c because that would need more elaborate testing for being inside
c the computational region, and is a second order effect.
                  ih=1
                  hf=66.
                  call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $                 ,zetap,ih,hf)
                  call getaccel(i,accel,il,rf,ith,tf,ipl,pf,st,ct,
     $                 sp,cp,rp,zetap,ih,hf)
                  accel(3)=accel(3)+Eneutral
c Save the position at end of step (inside the grid)
                  xpi1=xp(1,i)
                  xpi2=xp(2,i)
                  xpi3=xp(3,i)
c Back track the position to the point of last collision
c [None of this is correct for finite magnetic field. For that we would
c need a rotation of the perpendicular position/velocity.]
                  xp(1,i)=xp(1,i)-dtd*xp(4,i)
                  xp(2,i)=xp(2,i)-dtd*xp(5,i)
                  xp(3,i)=xp(3,i)-dtd*xp(6,i)
c Get new velocity; reflects neutral maxwellian shifted by vneutral.
                  xp(4,i)=tisq*gasdev(idum)
                  xp(5,i)=tisq*gasdev(idum)
                  xp(6,i)=tisq*gasdev(idum)+ vneutral
c Forward track fractional position with new velocity
                  xp(1,i)=xp(1,i)+dtd*xp(4,i)
                  xp(2,i)=xp(2,i)+dtd*xp(5,i)
                  xp(3,i)=xp(3,i)+dtd*xp(6,i)
                  rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
c Make sure we are inside the grid. If not, it was either because the forward
c tracking was too much or because this was a reinjected particle itself, so 
c its apparent prior position was outside grid.
c Compromise by simply using the end-of-step position which was inside.
c This may have unintended effects near the boundary. An alternative might
c be to reinject without a collision possibility.
                  if(rn.ge.r(NRFULL) .or. rn.le.1.)then
                     xp(1,i)=xpi1
                     xp(2,i)=xpi2
                     xp(3,i)=xpi3
c                  rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
c                  if(rn.le.1.)write(*,*)'Collide inside error',rn
                  endif
c Apply acceleration for the forward step duration. This is always
c starting from a newly collided particle. Fixes drift bias at high
c collisionality.
                  xp(4,i)=xp(4,i)+dtd*accel(1)
                  xp(5,i)=xp(5,i)+dtd*accel(2)
                  xp(6,i)=xp(6,i)+dtd*accel(3)
c Test for another collision
                  goto 1
               endif 
            endif
         endif
      enddo
      NCneutral=ncollide
      end
c******************************************************************
      subroutine mfpcollide(dt,colnwt,icolntype)
c Collision type icolntype=2: with specified mean free path.
c Input dt=timestep, colnwt=1/mfp (\propto target density)
c       ichoose=starting particle icycle=particle-number step.
      include 'piccom.f'

      icycle=10
      idum=1
c Decide which subset of particles to collide.
      ichoose=ran0(idum)*icycle+1
      if(ichoose.gt.icycle)ichoose=icycle
      do i=ichoose,iocprev,icycle
         if(ipf(i).gt.0)then
c     The velocity calculation cost about 5%.
c     It could be avoided if we stored it, perhaps. But actually padvan does
c     not take the square root, which is the dominant cost.
            v=sqrt(xp(4,i)**2+xp(6,i)**2+xp(6,i)**2)
c     Collision length. Should account for velocity dependence of sigma.
c     Enhance it by a factor of icycle to account for subsetting of collisions.
c     Account for the collision density, which perhaps should be set to give
c     collision length of one unit when it is one.
            xL=1./colnwt/icycle
c     Just a ran0 call increases the computational burden by about 5%.
            y=ran0(idum)
c     Here we need to invert a poisson distribution to tell if we have had
c     a collision. The cumulative poisson distribution is P(<x)=1-exp(-x/L)
c     where L=1/n\sigma is the collision length. Consequently,
c     if the random number is less than 1 - exp(v.dt/L) we have collided.
c     The exponentiation costs about 8.
            ytest=1.-exp(-v*dt/xL)
            if(y.lt.ytest)then
               write(*,*)'Collision:',y
            endif
         endif
      enddo
      end
c*******************************************************************
      subroutine fcalc_infdbl(dt)
      include 'piccom.f'
      real dt
      imin=1.

      call innerbc(imin,dt)

      do j=1,NTHUSED
         do k=1,NRUSED
            phi(k,j)=vprobe/rcc(k)
         enddo
      enddo
      do j=1,NTHUSED
         phi(0,j)=2.*phi(imin,j)-phi(imin+1,j)
      enddo
      do i=1,NRUSED
         phi(i,0)=phi(i,imin+1)
         phi(i,NTHUSED+1)=phi(i,NTHUSED-imin)
      enddo
      write(*,'($)')" "

      end

