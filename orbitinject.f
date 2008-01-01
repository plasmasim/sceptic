c***********************************************************************
c General version allows choice of reinjection scheme.
c***********************************************************************
      subroutine reinject(i,dt,icolntype,bcr)

      integer bcr 

      if(bcr.ne.0) then
         call maxreinject(i,dt,bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
         call fvreinject(i,dt,icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
         call ogenreinject(i,dt)
      else
c Injection from a shifted maxwellian at infinity
         call oreinject(i,dt)
      endif
      end
c***********************************************************************
      subroutine injinit(icolntype,bcr)

      integer bcr

      if(bcr.ne.0) then
c Injection from a maxwellian at boundary?
         call maxinjinit(bcr)
      elseif(icolntype.eq.1.or.icolntype.eq.5) then
c Injection from fv distribution at the boundary.
         call fvinjinit(icolntype)
      elseif(icolntype.eq.2.or.icolntype.eq.6)then
c Injection from a general gyrotropic distribution at infinity
         call ogeninjinit(icolntype)
      else
c Injection from a shifted maxwellian at infinity
         call oinjinit()
      endif
      end
c***********************************************************************
c***********************************************************************
c Other versions are in other source files.
      subroutine oreinject(i,dt)

      include 'piccom.f'
      integer i
      real dt
c Common data:
      parameter (eup=1.e-10)
      external pu
      logical istrapped
c Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist


c In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
      ilaunch=0
 1    continue
      ilaunch=ilaunch+1
      if(ilaunch.gt.1000)then
         write(*,*)'ilaunch excessive. averein=',averein,' brcsq=',
     $        brcsq,' ierr=',ierr,' rp=',rp
         stop
      endif
c Pick normal velocity from cumulative Pu
      y1=ran0(idum)
      call finvtfunc(pu,nvel,y1,u)
      iv=u
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(.not.dv.le.1)write(*,*)'Error in u calculation',y1,u,iv,dv
      vdist(iv)=vdist(iv)+1.
c Pick angle from cumulative Pc.
      y=ran0(idum)
c Here the drift velocity is scaled to the ion temperature.
      Uc=vd/vscale
      uu2=2.*Uc*u
      if(uu2.gt.50.) then
         crt=1.+alog(y)/uu2
      elseif(uu2.lt.-50.) then
         crt=-1.+alog(1-y)/uu2
      elseif(abs(uu2).lt.1.e-5)then
         crt=2.*y -1.
      else
         expuu2=exp(uu2)
c This expression is evaluated very inaccurately if expuu2 is nearly 1.
c That is why such a large limit on abs(uu2) is adopted.
         crt=alog(y*expuu2 + (1-y)/expuu2)/uu2
c The following do not do any better at solving this problem.
c         crt=alog( (y*expuu2 + (1-y)/expuu2)**(1./uu2))
c         crt=-1.+alog(1.+(expuu2**2-1.)*y)/uu2
      endif
      if(.not. abs(crt).le.1)then
c         write(*,*)'Strange crt:',crt,y,expuu2,uu2
c It seems impossible to avoid occasional strange results when uu2 is small.
         crt=2.*y-1.
      endif
c Testing angular distribution.
      if(LCIC)then
         icr=(1.+crt)*0.5*(NTHUSED-1) + 1.5
      else
         icr=(1.+crt)*0.5*(nth-1) + 1
      endif
      crdist(icr)=crdist(icr)+1.
c End of testing distribution monitor.
      srt=sqrt(1.- crt**2)
c Pick angle zt of poloidal impact and angle eta of impact parameter
      zt=ran0(idum)*2.*pi
      czt=cos(zt)
      szt=sin(zt)
      eta=ran0(idum)*2.*pi
      ceta=cos(eta)
      seta=sin(eta)
c Choose impact parameter, preventing overflow.
      chium2=-averein/Ti/(u+eup)**2
      if(chium2.le.-1.) then
c         write(*,*)'Impossible chium2=',chium2,' averein=', averein,
c     $        ' u=',u,' iv=',iv
         goto 1
c         stop
      endif
c      if(.not.lfixedn)chium2=0.
      brcsq=ran0(idum)*(1.+chium2)
c Reject a particle that will not reach boundary.
      if(brcsq.lt.0.) then
         goto 1
      endif
      brc=sqrt(brcsq)
c Get cosine and sine of impact angle relative to distant position.
c Based on integration.
      p2=brcsq*2.*Ti*u**2
      ierr=0
      if(debyelen.gt..001)then
c Orbit integration angle calculation.
c There is an overflow with this at zero debyelen. Ought to be properly fixed.
         call alphaint(p2,brcsq,cosal,sinal,ierr)
         if(ierr.ne.0)goto 1
c      write(*,'(4f9.4)')cosal-alcos(brc,chium2),sinal-alsin(brc,chium2)
c     Now ilaunch is the number of launches at infinity it took to get
c     one that reached the boundary.
      else
c Alternative based on analytic orbit calculation.
c Used for low debyelen, but really assumes negligible boundary potential.
c         call alcossin(brc,chium2,cosal,sinal)
         cosal=alcos(brc,chium2)
         sinal=alsin(brc,chium2)
      endif
c Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
      rs=r(nr)*0.99999
      xp(1,i)=rs*(czt*a1+szt*seta*sinal)
      xp(2,i)=rs*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rs*(-srt*ceta*sinal + crt*cosal)

c Obtain angle coordinate and map back to th for phihere.
      ct=xp(3,i)/rs
      call invtfunc(th(1),nth,ct,x)
      ic1=x
      ic2=ic1+1
      dc=x-ic1
c This expression should work for CIC And NGP.
      phihere=(phi(NRUSED,ic1)+phi(NRFULL,ic1))*0.5*(1.-dc)
     $        +(phi(NRUSED,ic2)+phi(NRFULL,ic2))*0.5*dc
c Section to correct the injection velocity and direction (but not the
c position) to account for local potential. 26 July 2004.
      if(localinj)then
         brcsq=(brcsq*(1.-phihere/Ti/(u+eup)**2)/(1.+chium2))
         if(brcsq.lt. 0.) then
c     This launch cannot penetrate at this angle. But it would have done
c     if the potential were equal here to averein. Thus it probably
c     should not be counted as part of the launch effort. So
            ilaunch=ilaunch-1
            goto 1
         endif
         chium2=-phihere/Ti/(u+eup)**2
         brc=sqrt(brcsq)
      endif
c Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
c Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
c Remove the following when using new advancing code.
c Increment the position by a random amount of the velocity.
c This is equivalent to the particle having started at an appropriately
c random position prior to reentering the domain.

c      xinc=ran0(idum)*dt
c      xinc=0.
c      vdx=0.

c Add magnetic field in the reinjection.
c ndivinj to add precision in the reinjection. rmoved for now
c      ndivinj=1
c      xinc=xinc/ndivinj

c      do k=1,ndivinj
c         do j=1,3
c            vdx=vdx+xp(j,i)*xp(j+3,i)
c         enddo
c         xp(3,i)=xp(3,i)+xp(6,i)*xinc
c         if(Bz.eq.0) then
c            do j=1,2
c               xp(j,i)=xp(j,i)+xp(j+3,i)*xinc
c            enddo
c         else
c            cosomdt=cos(Bz*xinc)
c            sinomdt=sin(Bz*xinc)
c            xp(1,i)=xp(1,i)+(xp(5,i)*(1-cosomdt)+xp(4,i)*sinomdt)/Bz
c            xp(2,i)=xp(2,i)+(xp(4,i)*(cosomdt-1)+xp(5,i)*sinomdt)/Bz
c            temp=xp(4,i)
c            xp(4,i)=temp*cosomdt+xp(5,i)*sinomdt
c            xp(5,i)=xp(5,i)*cosomdt-temp*sinomdt         
c         endif
c         
c         if(vdx.gt.0.)then
c            write(*,*)'Positive projection. u,phi=',u,phihere
c 601        format(a,5G10.5)
c      endif
c If we don't recalculate rp, then we don't trap NANs in the random choices.
      rcyl=xp(1,i)**2+xp(2,i)**2
      rp=rcyl+xp(3,i)**2
      rp=rs
c      write(*,*)'oreinject',rp
c Reject particles that are already outside the mesh.
      if(.not.rp.lt.r(nr)*r(nr))then
c      if(.not.rp.le.r(nr)*r(nr))then
         write(*,*)'Relaunch',rp,xp(1,i),xp(2,i),xp(3,i)
         goto 1
      else
c Do the outer flux accumulation.
c In order to accumulate the number of launches at infinity, rather than
c just the number of reinjections, we weight this by ilaunch
         spotrein=spotrein+phihere*ilaunch
         nrein=nrein+ilaunch
         fluxrein=fluxrein+1.
c Diagnostics of erroneous injects. Should not be needed:
c         if(.not. xp(1,i).le.400.)then
c            write(*,*)'Reinject overflow',i,xp(1,i),cosal,sinal
c     $        ,czt,szt,crt,srt,ceta,seta,brc,chium2,u,Ti,averein
c     $           ,y1,nvel
c         endif
         if(istrapped(i))then
            ntrapre=ntrapre+ilaunch
c            v=sqrt(xp(4,i)**2+xp(5,i)**2+xp(6,i)**2)
c            write(*,*)'Trapped',vdx/rp,u,v,sqrt(u**2-2.*averein)
c crt,czt,ceta,cosal
         endif
      endif
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine oinjinit()
c Common data:
      include 'piccom.f'

c Here the drift velocity is scaled to the ion temperature.
c And U's are in units of sqrt(2T/m), unlike vd.
      Uc=abs(vd)/sqrt(2.*Ti)
c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(Uc)

c Can't use these formulas for Uc exactly equal to zero.
      if(abs(Uc).lt.1.e-4)then
         if(Uc.lt.1.e-20) Uc=1.e-20
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            expu0=exp(-u0**2)
            pu2(i)=2.*Uc*expu0
            pu1(i)=0.5*(4.*u0**2*Uc + 2.*Uc)*expu0
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      else
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            uplus=u0+Uc
            uminus=u0-Uc
            pu2(i)=0.5*sqrt(pi)*(erfcc(uminus)-erfcc(uplus))
            pu1(i)=0.5*(-uminus*exp(-uplus**2)+uplus*exp(-uminus**2))
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      endif
      call srand(myid)
      end
c***********************************************************************
c***********************************************************************
c Calculate the cumulative probability for velocity index iu such that
c         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pu(iu)
      integer iu
c     averein is the average potential of reinjected particles, which is
c     used as an estimate of the potential at the reinjection boundary.
c     It is expressed in units of Te so needs to be scaled to Ti.
      include 'piccom.f'
      pudenom=pu1(1)-pu2(1)*averein/Ti
      pu=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
c********************************************************************
c Given a monotonic (increasing?) 
c function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine finvtfunc(Q,nq,y,x)
c Somehow this breaks the passing of a function reference.
c      implicit none
c      real external Q
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql
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
c Formerly .lt. which is an error.
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
c      x=(y-Ql)/(Qr-Ql)+iql
c Trap errors caused by flat sections.
      Qd=Qr-Ql
      if(Qd.eq.0.)then
         x=(iql+iqr)/2.
      else
         x=(y-Ql)/Qd+iql
      endif
      end
c**********************************************************************
C Inverse square law (phi\propto 1/r) injection functions:
c**********************************************************************
      subroutine alcossin(s,c,cosal,sinal)
      real s,c,cosal,sinal
      cosal=alcos(s,c)
      sinal=alsin(s,c)
      end
c**********************************************************************
      real function alcos(s,c)
      real s,c
      if(abs(s).le.1.e-12*abs(c))then
         alcos=-1.
         return
      else
         r=c/(2.*s)
         alcos=-(sqrt(1.+c-s**2)-(s-r)*r)/(1+r**2)
      endif
      end
c**********************************************************************
      real function alsin(s,c)
      real s,c
      if(s.le.1.e-12*c)then
         alsin=0.
         return
      else
         r=c/(2.*s)
         alsin=(sqrt(1.+c-s**2)*r+(s-r))/(1+r**2)
      endif
      end
c**********************************************************************
c**********************************************************************
c Return the angle cosine and sine for impact of a particle calculated
c by integrating the angle formula for a general central force.
c We integrate in the variable q=1/r from 0 to 1, provided that we don't
c encounter a barrier. If we do encounter one, we return an error.
c
c This version designed to work with uneven spacing if
c necessary, and trapezoidal integration. And uses init2ext
      subroutine alphaint(p2,b2,cosal,sinal,ierr)
c     The angular momentum and impact parameter squared
c     in units such that injection radius is 1.
      real p2, b2
c     The angle values returned.
      real cosal, sinal
c     The error return signal if non-zero.
      integer ierr
c Common data:
      include 'piccom.f'
      integer iqsteps
c This choice ensures iqsteps is large enough to accommodate a profile
c read in for processing with orbitint, but might not be the best choice for
c regular use.
c      parameter (iqsteps=nrsize+1)
      parameter (iqsteps=100+1)
      real phibye(iqsteps),phiei(iqsteps)
      real qp(iqsteps),pp(iqsteps)
      logical uninitialized
      data uninitialized/.true./
      save
c Statement function
c Inverse square law potential.
      extpot(q)=averein*q
c Not currently in use.

c First time through, do the initialization.
      if(uninitialized)then
         iqs=iqsteps
         uninitialized=.false.
c Screening length accounting for both ions and electrons.
         xlambda=debyelen/sqrt(1.+1./Ti)
         if(xlambda.lt.1.e-6)xlambda=1.e-6
         if(diagrho(1).eq.999.)then
c Special case to use the read in data.
c Some common data used in abnormal ways.
            qp(1)=0.
            phibye(1)=0.
            if(ninner.gt.iqsteps)stop 'Alphaint Too many r-cells read.'
            do i=1,ninner
               qp(i+1)=1./rcc(ninner+1-i)
               phibye(i+1)=diagphi(ninner+1-i)
            enddo
            iqs=ninner+1
            if(qp(iqs).ne.1.)then
               iqs=iqs+1
               qp(iqs)=1.
c Extrapolate linearly in q.
               phibye(iqs)=phibye(iqs-1)+
     $              (phibye(iqs-1)-phibye(iqs-2))*
     $              (qp(iqs)-qp(iqs-1))/(qp(iqs-1)-qp(iqs-2))
            endif
c The read-in case just uses phibye not phiei. And scales phibye(iqs) to 1
c putting the absolute edge value into averein.
            averein=phibye(iqs)
            adeficit=0.
            do i=1,iqs
               phibye(i)=phibye(i)/phibye(iqs)
            enddo
c End of read-in potential special case.
         else
c Specify the q-array. It goes from 0 to 1.
            do i=1,iqs
               qp(i)=((i-1.)/(iqs-1.))
            enddo
c            write(*,*)iqs,r(nr),xlambda,averein,adeficit
c     Since phi is specified in real space units, we tell the initialization
c     function what the rmax really is, and it does the transformation.
            call initext(iqs,qp,phibye,phiei,r(nr),xlambda)
         endif
c         write(*,'(2f8.4)')(qp(j),phibye(j),j=1,iqs)
c Hack to ignore outside:
c         do i=1,iqs
c            phibye(i)=0.
c            phiei(i)=phiei(i)*0.5
c         enddo
c Diagnostics
c         write(*,'(3f8.4)')(qp(j),phibye(j),phiei(j),j=1,iqs)
         if(averein.ne.0)then
c     When used by SCEPTIC averein is zero the first time, so this will
c     not be called.
            do i=1,iqs
               pp(i)=averein*phibye(i)-adeficit*phiei(i)
            enddo
            write(*,*)'adeficit=',adeficit
            call autoplot(qp,pp,iqs)
            call axlabels('q','potential')
            call pltend()
         endif
      endif
c End of initialization section.
c Do the integration for the orbit.
      ierr=0
      b2i=1./b2
      p2i2=2./p2
      sa=b2i - qp(1)**2 - p2i2*(averein*phibye(1)-adeficit*phiei(1))
      d1=(1./sqrt(sa))
c Inverse square case.
c      d1=1./sqrt(b2i)
      alpha=0.
c Trapezoidal rule.
      do i=2,iqs
         d2=d1
         sa=b2i - qp(i)**2 - p2i2*(averein*phibye(i)-adeficit*phiei(i))
c Inverse square case.
c         sa=b2i - qp(i)**2 - p2i2*extpot(qp(i))
         if(sa .le. 0.) goto 2
         d1=(1./sqrt(sa))
         alpha=alpha+(qp(i)-qp(i-1))*(d1+d2)*.5
      enddo
c      write(*,*)'alpha=',alpha
c Negative sign for definition of alpha relative to the forward direction.
      cosal=-cos(alpha)
      sinal=sin(alpha)
      return
 2    ierr=i
c      write(*,101)i,b2,p2,averein,adeficit
 101  format('Barrier: i=',i3,' b2=',f8.1,
     $     ' p2=',f8.4,' averein=',f8.4,' adeficit=',f8.4)
      end
c********************************************************************


