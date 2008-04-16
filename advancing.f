c___________________________________________________________________________
c
c     This code is copyright (c)
c              Ian H Hutchinson    hutch@psfc.mit.edu.
c              Leonardo Patacchini patacchi@mit.edu
c
c     It may be used freely with the stipulation that any scientific or
c     scholarly publication concerning work that uses the code must give
c     an acknowledgement referring to the relevant papers
c
c     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
c     1953 (2002), vol 45, p 1477 (2003).
c
c     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
c     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
c
c     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
c     p013505 (2007)
c
c     The code may not be redistributed except in its original package.
c
c     No warranty, explicit or implied, is given. If you choose to build
c     or run the code, you do so at your own risk.
c___________________________________________________________________________

c     Version 3.5 Dec 07 a) Changed the integrator to account for the
c     external E-field in the drift part of the integrator. b) Include
c     storage of collected energy, and collected momentum at infinity
c     Version 3.0 April 07 Moving collisions into this advancing routine
c     Version 2.5; Jan 2005: subcycling of padvnc.  Version 2.5; Jan
c     2005: fixed reinjection flux option.  Advance the particles
      subroutine padvnc(dtin,icolntype,colnwt,step)

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
      logical lcollide,lcstep


c Choose the collision cycle here and set tau appropriately:
c icycle of 1 costs about 20% extra cf false. (Mostly alog, I'd guess).
      if(colnwt.gt.0.)then
c         icycle=1./(20.*colnwt*dtin)
c         icycle=1./(50.*colnwt*dtin)
         icycle=1
         if(.not.icycle.ge.1) icycle=1
         ichoose=0
c         ichoose=ran0(idum)*icycle
         tau=1./(colnwt*icycle)
         lcollide=.true.
      else
         lcollide=.false.
         icycle=1
         ichoose=1
         tau=1.e20
      endif

      idum=1

c Xp is the three x-coordinates followed by the 3 v coordinates.
c Use a leapfrog scheme, so interpret the v-coords as half a step
c behind the x-coords. 
      tisq=sqrt(Ti)
c If lsubcycle, use multiple fractional steps near inner boundary.
      dt=dtin
      rp2=r(1)**2
c Zero the sums.
      ncollide=0
      nrein=0
      nreintry=0
      spotrein=0.
      ninner=0
      fluxrein=0.
      ntrapre=0
      zmomprobe=0.
      enerprobe=0.
      collmom=0.
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

      ido=npart

c      write(*,*)'colnwt,tau,Eneutral,icycle',colnwt,tau,Eneutral,icycle
c End of setup
c------------------ Iterate over particles --------------------------
c No-subcycle default. Never gets changed w/o subcycling.
      dts=dtin
      isubcycle=1
      do i=1,ido

         if(ipf(i).gt.0) then
c ````````````````````````````````````````` Treatment of active slot.
c     Find the mesh position and the trigonometry.
c     Here we do need half quantities.
            ih=1
            hf=88.
            call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $           ,zetap,ih,hf)

c .................... Subcycle Loop .................
            remdt=dtin
            ic=0
            lcstep=.false.
c            do 81 ic=1,isubcycle    Obsolete.
c Here is the start of the modified loop, now explicit.
c We iterate till we have used up the whole time step dtin (remdt=0).
c Steps may be shortened by subcycling and collisions.
 80         ic=ic+1

c  Now we know where we are in radius rp. 
c  We decide the level of subcycling.
            if(lsubcycle) then
               isubcycle=r(nrfull)/rp
c          if(mod(i,1000).eq.0) write(*,'(i1,$)')isubcycle
               dts=dtin/isubcycle*1.00001

            endif

c If prior step was ended by a collision, restart the particle velocity.

               if(lcstep)then
                  call postcollide(i,tisq)
                  lcstep=.false.
c     Because postcollide selects the velocity and the position at the
c     same time, we need to set dtprec to zero, in order to offset v and
c     x by half a time step properly.
                  dtprec(i)=0
               endif
               dt=min(dts,remdt)
               if(lcollide .and. mod(i,icycle).eq.ichoose)then
c Here we calculate the time to next collision: cdt
c Based on random number draw and poisson distribution.
                  cdt= -alog(ran0(idum))*tau
c Using this approximation instead saves negligible time:
c                  cdt= ran0(idum)*tau
c So I conclude that the only loss of time is initialization.
                  if(cdt.lt.dt)then
c Collision at the end of cdt step.
                     dt=cdt
                     lcstep=.true.
                     ncollide=ncollide+1
                  endif
               endif
               if(.not.dt.lt.1000.)then
c Error trap 
                  write(*,*)'dt error: dt, cdt, dts, remdt',
     $                 dt, cdt, dts, remdt
               endif
               remdt=remdt-dt

c Except for the first time, find new position.
               if(ic.ne.1)then 
                  ih=1
                  hf=77.
                  call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)
               endif            
               call getaccel(i,accel,il,rf,ith,tf,ipl,pf,st,ct,
     $              sp,cp,rp,zetap,ih,hf)


c For acceleration, when dt is changing, use the average of prior and
c present values: dtnow.
c               if(dtprec(i).eq.0.)dtprec(i)=dt
               dtnow=0.5*(dt+dtprec(i))

               if(.not.verlet) then  

     
c     New integrators accounting for the fact that we know what an orbit
c     is in a uniform E-field parallel to a uniform B.

c     Kick (Because of the E-field from the probe)
                  do j=4,6
                     xp(j,i)=xp(j,i)+accel(j-3)*dtnow
                  enddo

c     Drift with Bz.ne.0 (Cyclotronic integrator) Perpendicular direction
                  if(Bz.ne.0) then
                     cosomdt=cos(Bz*dt)
                     sinomdt=sin(Bz*dt)
                     xp(1,i)=xp(1,i)+
     $                    (xp(5,i)*(1-cosomdt)+xp(4,i)*sinomdt)/Bz
                     xp(2,i)=xp(2,i)+
     $                    (xp(4,i)*(cosomdt-1)+xp(5,i)*sinomdt)/Bz
                     
                     temp=xp(4,i)
                     xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
                     xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt
                  else
                     do j=1,2
                        xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                     enddo          
                  endif

c     Additional drift of the z position and velocity due to the
c     external E-field (if neutral collisions)
                  xp(3,i)=xp(3,i)+xp(6,i)*dt+0.5*Eneutral*dt**2
                  xp(6,i)=xp(6,i)+Eneutral*dt

               else
c     Old Verlet integrator


c     Getaccel returns the accel based on the charge-field calculation.
                  accel(3)=accel(3)+Eneutral

                  if(Bz.eq.0.)then
c     Don't use split steps if Bz=0, for speed gain of 9%.
                     do j=4,6
                        xp(j,i)=xp(j,i)+accel(j-3)*dtnow
                     enddo
                     do j=1,3
                        xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                     enddo

                  else
c Old Boris integrator

c First half of velocity advance:    AccelPhi/2+AccelBz+AccelPhi/2
                     do j=4,6
                        xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
                     enddo
c B-field rotation
                     cosomdt=cos(Bz*dtnow)
                     sinomdt=sin(Bz*dtnow)         
                     temp=xp(4,i)
                     xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
                     xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt
c Second half of velocity advance
                     do j=4,6
                        xp(j,i)=xp(j,i)+accel(j-3)*dtnow/2
                     enddo
                     
                     do j=1,3
                        xp(j,i)=xp(j,i)+xp(j+3,i)*dt
                     enddo
                  endif

               endif

               dtprec(i)=dt
               rn2=0.
               xdv=0.
               v2=0.

c Position advance
               do j=1,3
                  rn2=rn2+xp(j,i)**2
                  xdv=xdv+xp(j,i)*xp(j+3,i)
                  v2=v2+xp(j+3,i)**2
               enddo
           
c The time prior to step end of closest approach
               tm=xdv/v2
               rn=sqrt(rn2)
c  Test if we went through the probe and came back out.
               if((0..lt.tm .and. tm.lt.dt .and.
     $              (rn2 - tm**2*v2).lt.rp2))then
c For a long time this had an error: used  tm**2/v2 erroneously. 
c Corrected 9 Apr 07.
                  if(rn.gt.r(1))then
c     write(*,*)'Through probe',tm,(rn2 - tm**2*v2)
                     rn=0.
                  endif
               endif

c-----------------------------------------------------------------               
c Handling boundaries :
               if(rn.le.r(1)) then
                  ninner=ninner+1

c     Solve for sphere crossing step fraction, s.
c It ought to be possible to do this with the tm-related information.
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
                  if(LCIC)then
                     icell=nint(ithc+tfc)
                  else
                     icell=ithc
                  endif
                  ninth(icell)=ninth(icell)+1
c     Collected momentum and energy
                  zmomprobe=zmomprobe+xp(6,i)
                  enerprobe=enerprobe+0.5*v2
                  collmom=collmom+vzinit(i)
               elseif(rn.ge.r(nr))then
c     Left the grid outer boundary.
                  zmout=zmout-xp(6,i)
c     Did not leave the grid. Jump to subcycle end.
               else
                  goto 81
               endif
c We left. If we haven't exhausted complement, restart particle i.
               if(nrein.lt.ninjcomp) then
                  call reinject(i,dtin,icolntype,bcr)
                  dtprec(i)=dtin
                  ipf(i)=1
                  zmout=zmout+xp(6,i)
                  vzinit(i)=xp(6,i)
                  if(i.le.norbits) then
                     if (.not.(orbinit))
     $                    iorbitlen(i)=0
                  endif
c If an ion is reinjected but was to collide outside, it should not collide
c after reinjection !
                  lcstep=.false.


c     New reinjection handling. Simply use the rest of the time step with
c     the new particle starting just at the edge. Get new position:
                  ih=1
                  hf=77.
                  call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,
     $                 sp,cp,rp,zetap,ih,hf)
c     Set the external step length, (isubcycle=1).
                  dts=dtin
c     The ion is reinjected with v and x synchronized. We set dtprec to
c     zero to offset v and x by half a timestep
                  dtprec(i)=0
c     Call the timestep fraction-remaining random.
                  remdt=dtin*ran0(idum)


c     Jump to subcycle end.
                  goto 81
               else
                  ipf(i)=0
               endif
c Break from subcycles after dealing with a particle that left.
               goto 82  
           
 81            continue
c     Explicit cycle controlled by remaining time in step:
               if(remdt.gt.0.) goto 80
c     .................... End of Subcycle Loop .................
c     Break jump point:
 82            continue
c -----------------------------------------------------------
               
               if(ldist) then
c Start of Various distribution diagnostics.
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
c------------------------End distribution diagnostics ---------------
            if(ipf(i).gt.0)iocthis=i
            
         elseif(nrein.lt.ninjcomp)then

c ```````````````````````````````````````` Treatment of INactive slot.
c Case for ipf(i) le 0 (empty slot) but still wanting to inject. 
c We should not come here unless .not.lfixedn.
c            write(*,*)'Reinjecting empty slot',i

            call reinject(i,dtin,icolntype,bcr)
            dtprec(i)=dtin
            ipf(i)=1
            iocthis=i
         elseif(i.ge.iocprev)then
c     Break if inactive slot and we have exhausted the complement of
c     injections.  And we have reached the maximum occupied slot of
c     previous run.
            goto 401
         endif
c---------------- End of padvnc particle iteration ------------------
      enddo
 401  continue

      NCneutral=ncollide
c      write(*,*)'ncollide=',ncollide,' icycle=',icycle
c We just want the diagnostics with the true particles for now
c      iocthis=min(iocthis,npartmax)

      iocprev=iocthis
c     if(.not.lfixedn)write(*,504)ninjcomp,nrein,i,iocprev
 504  format('  ninjcomp=',i6,'  nrein=',i6,'  i=',i6,
     $     '  iocprev=',i6)
 503  format('Orbit',i3,' length=',i5,' position=',4f7.3)
 502  format('Distrib. rn=',f6.3,' ithc=',i4,' vr=',f6.3)
 501  format('accel=',3f11.4,' xp=',3f11.4)

      end
c***********************************************************************
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

c parameters to find the floating potentail in presence of Bz
      integer nPhi
      data nPhi/300/
      real Irep(1:300)
      real ncs
      real LS(nthused)
      real z,iota,dpdr,Tau,eta,beta_e,beta_i
      data ncs/50./
c phispan is for floating potential if Bz.ne.0
      phispan=8
      beta_e=0
      iota=0


      if(linsulate.or.lfloat) then

         if(Bz.ne.0) then
            beta_i=Bz*sqrt(2/(Ti*pi))
            beta_e=beta_i*sqrt(Ti)*sqrt(rmtoz*1837.)
            z=beta_e/(1+beta_e)
            iota=1-0.0946*z-0.305*z**2+0.95*z**3-2.2*z**4+1.15*z**5
            
            do k=1,nthused
c Calculate the shielding length
               dpdr=1/phi(1,k)*(-phi(3,k)+4*phi(2,k)-3*phi(1,k))
     $              /((rcc(3)-rcc(1)))
               LS(k)=-1/(min(dpdr,-1.01)+1)

c               LS(k)=debyelen/sqrt(1+1/(Ti+rmtoz*vd**2))
c               LS(k)=debyelen/sqrt(1+1/Ti)
c               LS(k)=LS(k)+debyelen*log(1+1/debyelen)
            enddo
         endif

         flogfac=0.5*alog(2.*pi/(rmtoz*1837.))
         totflux=0.
         do j=1,nthused
c     Calculate the flux to each angular cell
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

               if(Bz.ne.0) then
                  do k=1,nPhi
                     Qth=1.-2*(j-1.)/(nthused-1.)
                     eta=phispan*k/nPhi/beta_e*
     $                    (1+beta_e/4*(1-exp(-4/(LS(j)*beta_e))))
                     Tau=eta/(1+eta)
                     A=0.678*Tau+1.543*Tau**2-1.212*Tau**3
                     Irep(k)=exp(-phispan*k/nPhi)
     $             *((2*(A+(1-A)*iota)-1)+2*(1-(A+(1-A)*iota))*abs(Qth))
                  enddo
               endif
               if(fluxofangle(j).gt.0.)then
                  if(Bz.ne.0) then
                     call invtfunc(Irep,nPhi,exp(flogfac)
     $                    *fluxofangle(j),x)
                     phi(imin,j)=
     $                    (-x*phispan/nPhi+(ncs-1.)*phi(imin,j))/ncs
                  else
                     phi(imin,j)=(alog(fluxofangle(j))+flogfac
     $                    +(ncs-1.)*phi(imin,j))/ncs
                  endif
               else
                  phi(imin,j)=phi(imin,j)
               endif
            endif
         enddo

         if(totflux.gt.0.)then
            if(Bz.eq.0) then
               vprobe=(alog(totflux/nthused)+
     $              flogfac+(ncs-1.)*vprobe)/ncs
c     comparison with old version
c               vprobe=alog(totflux/nthused)+flogfac
            else
c     calculate the total e- current to the sphere -> Irep
               do k=1,nPhi
                  Irep(k)=0
                  do j=1,nthused
                     axis=1
                     if (j.eq.1.or.j.eq.nthused) axis=0.5
                     Qth=1.-2*(j-1.)/(nthused-1.)
                     eta=phispan*k/nPhi/beta_e*
     $                    (1+beta_e/4*(1-exp(-4/(LS(j)*beta_e))))
                     Tau=eta/(1+eta)
                     A=0.678*Tau+1.543*Tau**2-1.212*Tau**3
                     Irep(k)=Irep(k)+axis*exp(-phispan*k/nPhi)
     $             *((2*(A+(1-A)*iota)-1)+2*(1-(A+(1-A)*iota))*abs(Qth))
                  enddo
                  Irep(k)=Irep(k)/(nthused-1)
               enddo
               call invtfunc(Irep,nPhi,exp(flogfac)
     $              *totflux/nthused,x)
               vprobe=(-x*phispan/nPhi+(ncs-1.)*vprobe)/ncs
            endif
         endif
         
c         write(*,*) vprobe,alog(totflux/nthused)+flogfac
         if(lfloat)then
            do j=1,nthused
               phi(imin,j)=vprobe
            enddo
         endif
c     write(*,*)
c     write(*,*)'fluxofangle=',(fluxofangle(j),j=1,NTHUSED)
c     write(*,*)'phi=',(phi(imin,j),j=1,NTHUSED)

c If prespecified probe potential
      else
         do j=1,NTHUSED
            phi(imin,j)=vprobe+Ezext*tcc(j)
         enddo     
      endif
      end

c****************************************************************** Set
c     the finite volumes coefficients for the outer boundary, as well as
c     the probe potential.
      subroutine fcalc_bc(dt,rshield,icolntype,colnwt)
      
      include 'piccom.f'
      real dt
c     If rshield is not an integer, problem with the MPI routines
      integer rshield
      real phislopeconst(0:nth+1),phislopefac(0:nth+1)
c Correct the type of temporary variable (IHH). Better use initial letter.
      integer c


c bpc=0 -> Use the spherical symmetry approximation (Hutch paper2)
c bpc=1 -> Quasineutrality on the 15% outer crone
c bpc=2 -> Phiout=0
c bpc=3 -> dPhi/dzout=0
c bcp=4 -> dPhi/drout=-Phiout/r (Use at high collisionality)

c Potential calculation in the shielding region
      if(bcphi.eq.1) then
         rshield=nint(NRUSED*.85)
         do j=1,nth
            do i=rshield,nr
               phi(i,j)=log(rho(i,j))
            enddo
         enddo
      else
         rshield=NRUSED
         if(bcphi.eq.2) then
            do j=1,nth
               phi(rshield,j)=0
            enddo
         endif
      endif
      
      if (bcphi.eq.0) then

c     Setting of the gpc array 
         redge= (rcc(NRFULL)+rcc(NRFULL-1))*0.5
         delredge=rcc(NRFULL)-rcc(NRFULL-1)
c     Screening k-number combines electrons and ions.
         if(debyelen.gt. 1.e-10) then
            el2=(1.+1./Ti)/debyelen**2
         else
            el2=2.e20
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
c         if(icolntype.eq.2)then
c Remove the deficit term when using simplistic rindex, otherwise 
c instability tends to result.
c         expE1=0.
c Simplistic trials.
c            rindex=(redge*el+1.)
c            rindex=1.
c            rindex=4.
c            rindex=10.
c         endif
         adeficit=0
c     Boundary slope factor calculations:
         do j=0,NTHUSED+1
c     Current fractional ion deficit due to collection.
c     Coefficient of 1/r^2 in modified shielding equation is
c     a = deficitj * r_edge^2 / \lambda_De^2
            deficitj=1-phi(NRUSED,j)/Ti -rho(NRUSED,j)
c Reduce the deficit term when collisionality is significant.
c Because it no longer applies. (Perhaps ought to account for vd).
            deficitj=deficitj/(1.+(colnwt*redge)**2/Ti)
            blfac1=(deficitj/debyelen**2) * redge
            adeficit=adeficit+blfac1
c     BC modification is (a/r_edge)[exp(EL*r) E_1(El*r)] given by approx.
            blfac=blfac1*expE1
            blfac=alpha*blfac
            phislopeconst(j)=blfac*redge*delredge/
     $           (redge+delredge*rindex*0.5)
            phislopefac(j)=(redge-delredge*rindex*0.5)/
     $           (redge+delredge*rindex*0.5)
c     Set gpc array
            gpc(j,1)=phislopefac(j)
            gpc(j,2)=0
            gpc(j,3)=0
            gpc(j,4)=-phislopeconst(j)
            gpc(j,5)=-1
         enddo
c     Actual a factor averaged over angles:
         adeficit=adeficit*redge/NTHUSED
         if(adeficit.lt.0.)then
c     write(*,*)'Negative adeficit',adeficit,' set to zero'
            adeficit=0.
         endif

      elseif(bcphi.eq.1.or.bcphi.eq.2) then
c     In this case, the potential outside the crone is prespecified
         do j=0,NTHUSED+1
            gpc(j,1)=0
            gpc(j,2)=0
            gpc(j,3)=0
            gpc(j,4)=0
            gpc(j,5)=0
         enddo

      elseif(bcphi.eq.3) then
c     This is the case where "dphi/dr=0" on the boundary
         delredge=rcc(NRFULL)-rcc(NRFULL-1)
         delcosth=2./(NTHUSED-1.)

         do j=0,NTHUSED+1
               Qth=1.-2*(j-1)/(NTHUSED-1.)
            gpc(j,1)=1
            gpc(j,2)=-delredge/(2*delcosth*rcc(NRUSED))*Qth*(1-Qth**2)
            gpc(j,3)=-gpc(j,2)
            gpc(j,4)=0
            gpc(j,5)=-1-delredge*((1-Qth**2)/(2.*rcc(NRUSED))
     $           +(1-Qth**2)**(1.5)/debyelen)
c            write(*,*) j,gpc(j,1),gpc(j,2),gpc(j,3),gpc(j,4),gpc(j,5)
         enddo
      elseif(bcphi.eq.4) then
c     This is to impose a 1/r potential at the outer edge, valid for
c     high collisionality.
         delredge=rcc(NRFULL)-rcc(NRFULL-1)
         do j=0,NTHUSED+1
            gpc(j,1)=1
            gpc(j,2)=0
            gpc(j,3)=0
            gpc(j,4)=0
            gpc(j,5)=-1-delredge/rcc(NRUSED)
c     write(*,*) j,gpc(j,1),gpc(j,2),gpc(j,3),gpc(j,4),gpc(j,5)
         enddo
         
      endif


c     Set inner Boundary conditions
      imin=1
      call innerbc(imin,dt)

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
c*******************************************************************
      subroutine fcalc_infdbl(dt)
      include 'piccom.f'
      real dt
      imin=1.

      call innerbc(imin,dt)
      do j=1,NTHUSED
         do k=1,NRUSED
            phi(k,j)=vprobe/rcc(k)*exp(-(rcc(k)-rcc(1))/debyelen)
         enddo
      enddo
      do j=1,NTHUSED
         phi(0,j)=2.*phi(imin,j)-phi(imin+1,j)
      enddo
      do i=1,NRUSED
         phi(i,0)=phi(i,imin+1)
         phi(i,NTHUSED+1)=phi(i,NTHUSED-imin)
      enddo

      end
c********************************************************************
      subroutine postcollide(i,tisq)
      include 'piccom.f'
      include 'colncom.f'
c Get new velocity; reflects neutral maxwellian shifted by vneutral.
      xp(4,i)=tisq*gasdev(idum)
      xp(5,i)=tisq*gasdev(idum)
      xp(6,i)=tisq*gasdev(idum)+ vneutral
      end
