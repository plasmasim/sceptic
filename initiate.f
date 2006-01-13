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
c***********************************************************************
c r and th are the radius and cosine theta meshes.
c rcc and tcc are the center-cell values, where the charge assigned to
c that cell should be considered to be centered. For CIC this is the
c same as the r and th, but for NGP not.
c***********************************************************************
      subroutine meshinitngp(rmax)
      real rmax
c Common data:
      include 'piccom.f'

      r(0)=1.-(rmax-1.)/(NRFULL-1)
      do i=1,NRFULL
c linear r mesh   
         r(i)=1.+(i-1)*(rmax-1.)/(NRFULL-1)
         rcc(i-1)=0.5*(r(i-1)+r(i))
c distance from the probe surface, called \rho in notes.
         hr(i-1)=rcc(i-1)-r(1)
         zeta(i-1)=sqrt(abs(2.*hr(i-1)))
         if(hr(i-1).lt.0.)zeta(i-1)=-zeta(i-1)
      enddo
c Uniform r-mesh extrapolation.
      rcc(nr)=2.*rcc(nr-1)-rcc(nr-2)
      hr(nr)=rcc(nr)-r(1)
      zeta(nr)=sqrt(2.*hr(nr))
      zeta(nr+1)=sqrt(2.*(2.*hr(nr)-hr(nr-1)))
      do i=1,NRFULL
c Half-mesh quantities
         if(i.eq.1)then
            rat=1.
         elseif(i.eq.NRFULL)then
            rat=(sqrt(sqrt(2.*(2.*r(NRFULL)-r(NRFULL-1))))
     $           -sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         else
            rat=(sqrt(zeta(i+1))-sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         endif
         zetahalf(i)=0.5*(zeta(i)+zeta(i-1))
         cminus(i)=(rat-2./rat +1.)/6.
         cmid(i)=(rat+1./rat -2.)/6.
         cplus(i)=(1./rat -2.*rat +1)/6.
      enddo
c     The edge zetahalf must never permit ih=nr; so we kludge slightly.
      zetahalf(nr)=sqrt(2.*(rmax-r(1)))+1.e-6
c We should never need the following.
c      zetahalf(nr+1)=0.5*(zeta(nr)+zeta(nr+1))
      zetahalf(0)=-zetahalf(2)
c Avoid rounding errors 
      zetahalf(1)=0.
      do i=1,nth
c theta array including poles
c Uniform in cos theta
         th(i)=1.-2.*(i-1)/(nth-1)
      enddo
c Additional angle positions are given past the ends for the purposes of the
c boundary conditions. They are a distance beyond the ends equal to the
c last step.
      th(0)=2.*th(1)-th(2)
      th(NTHUSED+1)=2.*th(NTHUSED)-th(NTHUSED-1)
      do i=0,nth-1
         tcc(i)=0.5*(th(i)+th(i+1))
         thang(i)=acos(tcc(i))
      enddo
      tcc(NTHUSED+1)=2.*tcc(NTHUSED)-tcc(NTHUSED-1)
      thang(NTHUSED+1)=2.*thang(NTHUSED)-thang(NTHUSED-1)
      thang(0)=2.*thang(1)-thang(2)
      if(NRFULL.le.10 .and. nth.le.10) then
         write(*,*)'r,rcc,th,tcc,thang'
         write(*,*)(r(j),j=0,nrfull)
         write(*,*)(rcc(j),j=0,nrfull)
         write(*,*)(th(j),j=0,nthfull)
         write(*,*)(tcc(j),j=0,nthfull)
         write(*,*)(thang(j),j=0,nthfull)
      endif
c      write(*,*)'zeta=',zeta
c
      do i=1,NRUSED
         vol= r(i+1)**3-r(i)**3
         volinv(i)=3./(4.*pi*vol)
         voltot=voltot+vol
      enddo
c      if(myid.eq.0) write(*,*)'Voltot',voltot,'   Ratio to actual',
c          voltot/(r(NRFULL)**3-r(1)**3)
c     Zero the ninth storage.
      do k=1,nstepmax
         do j=1,nth
            ninthstep(j,k)=0
         enddo
      enddo
      call precalc()
      end
c***********************************************************************
c Interpolate onto the theta mesh. Return nearest index, fraction in thf.
      integer function interpth(ct,thf)
      include 'piccom.f'
      ithl=itpre(1+int((ct-th(1))*tfac))
      thf=(ct-th(ithl))/(th(ithl+1)-th(ithl))
      if(thf.gt.1.)then
         if(ithl+2.le.NTHFULL)then
            ithl=ithl+1
            thf=(ct-th(ithl))/(th(ithl+1)-th(ithl))
         else
            write(*,*)'INTERPTH error. ithl, thf incorrect'
            write(*,*)ithl,thf,ct
         endif
      endif
      interpth=ithl
      end
c***********************************************************************
      subroutine meshinitcic(rmax)
      real rmax
c Common data:
      include 'piccom.f'

      r(0)=1.-(rmax-1.)/(NRFULL-1)
      do i=1,NRFULL
c linear r mesh   
         r(i)=1.+(i-1)*(rmax-1.)/(NRFULL-1)
         rcc(i)=r(i)
c distance from the probe surface, called \rho in notes.
         hr(i)=r(i)-r(1)
         zeta(i)=sqrt(2.*hr(i))
      enddo
c Uniform r-mesh extrapolation.
      zeta(0)=-zeta(2)
      zetahalf(0)=-0.5*(zeta(2)+zeta(3))
      zeta(nr+1)=sqrt(2.*(2.*r(nr)-r(nr-1)-r(1)))
      do i=1,NRFULL
c Half-mesh quantities
         if(i.eq.1)then
            rat=1.
         elseif(i.eq.NRFULL)then
            rat=(sqrt(sqrt(2.*(2.*r(NRFULL)-r(NRFULL-1))))
     $           -sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         else
            rat=(sqrt(zeta(i+1))-sqrt(zeta(i)))/
     $        (sqrt(zeta(i))-sqrt(zeta(i-1)))
         endif
         zetahalf(i)=0.5*(zeta(i)+zeta(i-1))
         cminus(i)=(rat-2./rat +1.)/6.
         cmid(i)=(rat+1./rat -2.)/6.
         cplus(i)=(1./rat -2.*rat +1)/6.
      enddo
      zetahalf(nr+1)=0.5*(zeta(nr)+zeta(nr+1))
      do i=1,nth
c theta array including poles
c Uniform in cos theta
         th(i)=1.-2.*(i-1)/(nth-1)
      enddo
c Additional angle positions are given past the ends for the purposes of the
c boundary conditions. They are a distance beyond the ends equal to the
c last step.
      th(0)=2.*th(1)-th(2)
      th(NTHUSED+1)=2.*th(NTHUSED)-th(NTHUSED-1)
      do i=1,nth
c     Cic version
         tcc(i)=th(i)
         thang(i)=acos(th(i))
      enddo
      tcc(NTHUSED+1)=2.*tcc(NTHUSED)-tcc(NTHUSED-1)
      thang(NTHUSED+1)=2.*thang(NTHUSED)-thang(NTHUSED-1)
      thang(0)=2.*thang(1)-thang(2)
      if(NRFULL.le.10 .and. nth.le.10) then
         write(*,*)'r,rcc,th,tcc,thang'
         write(*,*)(r(j),j=0,nrfull)
         write(*,*)(rcc(j),j=0,nrfull)
         write(*,*)(th(j),j=0,nthfull)
         write(*,*)(tcc(j),j=0,nthfull)
         write(*,*)(thang(j),j=0,nthfull)
      endif
c      write(*,*)'th=',th

c Calculate the mesh volumes
      rim=0.
      rm2=0.
      rm3=0.
      voltot=0.
c Silence warnings. Not otherwise necessary.
      ri=r(1)
      ri1=r(2)
c
      do i=1,NRUSED
         if(i.lt.NRUSED)then
            ri=r(i)
            ri1=r(i+1)
            rs3=ri1**3+ri1**2*ri+ri1*ri**2+ri**3
            rs2=ri1**2+ri1*ri+ri**2
         else
            rs3=0.
            rs2=0.
         endif
         vol= ri1*rs2-0.75*rs3 + 0.75*rm3-rim*rm2
         volinv(i)=3./(4.*pi*vol)
         rim=ri
         rm2=rs2
         rm3=rs3
         voltot=voltot+vol
      enddo
c      if(myid.eq.0) write(*,*)'Voltot',voltot,'   Ratio to actual',
c     $     voltot/(r(NRFULL)**3-r(1)**3)
c      write(*,*)'Volinv',volinv
c Zero the ninth storage.
      do k=1,nstepmax
         do j=1,nth
            ninthstep(j,k)=0
         enddo
      enddo
      call precalc()
      end
c***********************************************************************
c Initializing particles.
      subroutine pinit()
c Common data:
      include 'piccom.f'
      logical istrapped

c For now use the whole array.
      nrealin=0
      ntries=0
      ntrapped=0
      rmax=r(NRFULL)
      rsp=r(rsplit)
      rmax2=rmax*rmax
      rsp2=rsp*rsp
      idum=1
      if(rmax2.le.1.) stop 'Error: rmax is less than 1.'

c     We initialize the 'true' particles'
      do i=1,npart
         ipf(i)=1
 1       continue
         ntries=ntries+1
         xp(1,i)=rmax*(2.*ran0(idum)-1.)
         xp(2,i)=rmax*(2.*ran0(idum)-1.)
         xp(3,i)=rmax*(2.*ran0(idum)-1.)
         rc=0.
         do j=1,3
            rc=rc+xp(j,i)**2
         enddo
c If we are not in the plasma region, try again.
         if(rc.ge.rmax2 .or. rc.le.1.) goto 1
         Ti0=Ti
         tisq=sqrt(Ti0)
         xp(4,i)=tisq*gasdev(idum)
         xp(5,i)=tisq*gasdev(idum)
         xp(6,i)=tisq*gasdev(idum) + vd
         if(istrapped(i))then
            ntrapped=ntrapped+1
c If this goto is included then trapped particles are rejected.
c But that tends to deplete the region close to the probe.
c            goto 1
         endif

c Counting the number of real particles in the inner domain
         if (sqrt(rc).le.r(rsplit)) then
            nrealin=nrealin+1
         endif

      enddo

c     We now initialize the add particles
      if (dsub.and.(npartadd.ge.1)) then
         do i=npartmax+1,npartmax+npartadd
            ipf(i)=1
 2          continue
            ntries=ntries+1
            xp(1,i)=rsp*(2.*ran0(idum)-1.)
            xp(2,i)=rsp*(2.*ran0(idum)-1.)
            xp(3,i)=rsp*(2.*ran0(idum)-1.)
            rc=0.
            do j=1,3
               rc=rc+xp(j,i)**2
            enddo
c     If we are not in the inner region, try again.
            if(rc.ge.rsp2 .or. rc.le.1.) goto 2
            Ti0=Ti
            tisq=sqrt(Ti0)
            xp(4,i)=tisq*gasdev(idum)
            xp(5,i)=tisq*gasdev(idum)
            xp(6,i)=tisq*gasdev(idum) + vd
            if(istrapped(i))then
               ntrapped=ntrapped+1
c     If this goto is included then trapped particles are rejected.
c     But that tends to deplete the region close to the probe.
c     goto 2
            endif
               
c     Initialize the storage arrays
            if (i-npartmax.le.10) then
               do k=1,addhist
                  do l=1,6
                     xpstorage(l,i-npartmax,k)=xp(l,i)
                  enddo
                  xpstonum(k)=10
               enddo
            endif
         enddo
      endif

c Set flag of unused slots to 0
      do i=npart+1,npartmax
         ipf(i)=0
      enddo
      
c      write(*,*)'Initialized ','id=',myid,
c     $     '  n=',npart,'  ntries=',ntries,'  ntrapped=',ntrapped
c Initialize rhoinf:
      rhoinf=numprocs*npart/(4.*pi*r(NRFULL)**3/3.)
c Initialize orbit tracking
      do ko=1,norbits
         iorbitlen(ko)=0
      enddo

      end
c***********************************************************************
c Initializing the fields
      subroutine finit()
      include 'piccom.f'
      if(lfext)then
c Fields read from external file stored in diagphi
c cic boundary is at i=1, ngp at 0+(1/2) (sort of).
         if(LCIC)then
            imin=1
         else
            imin=0
         endif
         do j=1,NTHUSED
            phi(imin,j)=vprobe
            do i=imin+1,NRUSED
               phi(i,j)=diagphi(i)
            enddo
            phi(0,j)=2.*phi(imin,j)-phi(imin+1,j)
c Kludge the outer boundary for NGP
            if(NRFULL.gt.NRUSED)
     $           phi(NRFULL,j)=2.*phi(NRUSED,j)-phi(NRUSED-1,j)
         enddo
         do i=1,NRUSED
            phi(i,0)=phi(i,imin+1)
            phi(i,NTHUSED+1)=phi(i,NTHUSED-imin)
         enddo
      else
c Usual case.
         do j=0,NTHFULL
            do i=0,NRFULL
c     Free-space initialization.
               phi(i,j)=vprobe*r(1)/r(i)
c     Trivial initialization was used for a long time. Hardly different.
c     phi(i,j)=0.
               if(i.eq.0)phi(i,j)=vprobe
            enddo
         enddo
      endif
      end
c************************************************************************
      subroutine precalc()
c     Precalculation functions
      include 'piccom.f'

      rfac=(nrpre-1.)/(r(NRFULL)-r(1))
      tfac=(ntpre-1.)/(th(nth)-th(1))
      do j=1,ntpre
c     finding the theta precalculated mesh.
         thp=(j-1.)/tfac+th(1)
         thl=th(1)
         thr=th(nth)
         itl=1
         itr=nth
 200     if(itr-itl.le.1)goto 210
         itx=(itr+itl)/2
         thx=th(itx)
         if(thx.ge.thp) then
            thl=thx
            itl=itx
         else
            thr=thx
            itr=itx
         endif
         goto 200
 210     continue
         itpre(j)=itl
      enddo
c     r grid may be nonlinear find grid position by bisection.
      do j=1,nrpre
         rp=(j-1.)/rfac+r(1)
         rl=r(1)
         rr=r(NRFULL)
         il=1
         ir=NRFULL
 201     ix=(ir+il)/2
         rx=r(ix)
         if(rx.le.rp) then
            rl=rx
            il=ix
         else
            rr=rx
            ir=ix
         endif
         if(ir-il.gt.1)goto 201
c     Now il and ir, rl and rr bracket the radius.
         irpre(j)=il
      enddo
c      write(*,*)'Precalculated the r and theta mesh lookups.'
c      write(*,*)itpre,tfac


c Now irpre(1+int(rp-r(1))*rfac)) is the irl except for rounding etc.
c The irpre spacing must be small enough that the maximum increment of
c irpre from j to j+1 is 1. Then it is possible that the downward rounding
c causes irpre to be at most 1 too small.
c The same applies to itpre

      end
c***********************************************************************
c Initialize the poison iteration coefficients. Must be done after
c mesh initiation.
      subroutine poisinitngp()

c Common data:
      include 'piccom.f'

      do i=1,NRUSED
         ri=r(i)
         rip1=r(i+1)
         rave=(rip1**2+2.*rip1*ri+ri**2)/4.
         apc(i)=debyelen**2 *rip1**2/rave/(rcc(i+1)-rcc(i))/(rip1-ri)
         bpc(i)=debyelen**2 *ri**2/rave/(rcc(i)-rcc(i-1))/(rip1-ri)
         if(i.eq.1)then
c During iteration, phi(0) is set to phiprobe. So we have to have a different
c form for bpc. The following is not correct to second order:
c            bpc(i)=2.*bpc(i)
c I hope this is:
            raj=(rip1-ri)*.25
            bpc(1)=bpc(1)*(1.+((ri+raj)/(ri-raj))**2)
         endif

         do j=1,NTHUSED
            cpc(i,j)=debyelen**2/rave/(thang(j+1)-thang(j))**2
            dpc(i,j)=debyelen**2/rave/(thang(j)-thang(j-1))**2
            if(j.eq.1) dpc(i,j)=0.
            if(j.eq.NTHUSED) cpc(i,j)=0.
            fpc(i,j)=apc(i)+bpc(i)+cpc(i,j)+dpc(i,j)
         enddo

      enddo

      end

c***********************************************************************
c Initialize the poison iteration coefficients. Must be done after
c mesh initiation.
      subroutine poisinitcic()

c Common data:
      include 'piccom.f'

      do i=1,NRUSED-1
c These two statements are different for CIC.
         ri=  (r(i)+r(i-1))/2.
         rip1=(r(i+1)+r(i))/2.
c
         rave=(rip1**2+2.*rip1*ri+ri**2)/4.
         apc(i)=debyelen**2 *rip1**2/rave/(rcc(i+1)-rcc(i))/(rip1-ri)
         bpc(i)=debyelen**2 *ri**2/rave/(rcc(i)-rcc(i-1))/(rip1-ri)
c         if(i.eq.NRUSED-1)then
c CIC boundary condition for inverse square behaviour.
c            apc(i)=apc(i)*4.*(rcc(i+1)-rcc(i))/(3.*r(i+1)-r(i))
c         endif
         do j=1,NTHUSED
            cpc(i,j)=debyelen**2/rave
     $           *(0.5*(nth-1)*sin(acos((th(j+1)+th(j))/2.)))**2
            dpc(i,j)=debyelen**2/rave
     $           *(0.5*(nth-1)*sin(acos((th(j)+th(j-1))/2.)))**2
c old
c            cpc(i,j)=debyelen**2/rave/(thang(j+1)-thang(j))**2
c            dpc(i,j)=debyelen**2/rave/(thang(j)-thang(j-1))**2
c            fpc(i,j)=apc(i)+bpc(i)+cpc(i,j)+dpc(i,j)
            if(j.eq.1)then
               cpc(i,j)=2.*cpc(i,j)
               dpc(i,j)=0.
            elseif(j.eq.NTHUSED)then
               dpc(i,j)=2.*dpc(i,j)
               cpc(i,j)=0.
            endif
            fpc(i,j)=apc(i)+bpc(i)+cpc(i,j)+dpc(i,j)
         enddo
c Now set the boundary coefficient to zero once fpc is correctly calculated.
c That way we can use whatever boundary phi is most convenient.
c         if(i.eq.NRUSED-1)apc(i)=0.
      enddo
c Test section exploring angle approximations.
c      write(*,*)' sin(thang+1/2), -2/(nth-1)(th(+)-th)'
c      do j=1,NTHUSED-1
c         write(*,*)sin(acos((th(j+1)+th(j))/2.)),
c     $        2./(nth-1)/(thang(j+1)-thang(j))
c      enddo
c      read(*,*)iin
      end
c***********************************************************************
      logical function istrapped(i)

c     Return as logical whether the particle i is trapped or not.  It is
c     considered trapped if the energy available for radial velocity at
c     the outer boundary and at the probe, conserving angular momentum
c     and energy is negative. The assumption of angular momentum
c     conservation is false for non-symmetric situations.

      include 'piccom.f'

      ih=0
      hf=66.
      call ptomesh(i,il,rf,ith,tf,ipl,pf,st,ct,sp,cp,rp
     $     ,zetap,ih,hf)

      rn=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)

c The interpolation here might not be correct for both schemes.
      phin=(phi(il,ith)*(1.-tf)+phi(il,ith+1)*tf)*(1.-rf) +
     $     (phi(il+1,ith)*(1.-tf)+phi(il+1,ith+1)*tf)*rf

c Definition as being that the particle does not leave the domain
c      phie= phi(NRUSED,ith)
c Definition that the particle does not reach infinity.
      phie=0.
      phip= vprobe

      vr2=(xp(4,i)*xp(1,i)+xp(5,i)*xp(2,i)+xp(6,i)*xp(3,i))**2
     $     /rn
      v2=xp(4,i)**2+xp(5,i)**2+xp(6,i)**2
      vt2=v2-vr2

c Domain definition.
c      vte2=(vt2*(rn/rcc(NRUSED))**2)
      vte2=0. 
      vtp2=vt2*rn**2

      if( (vte2 .gt. v2 + 2.*(phin-phie) ) .and.
     $    (vtp2 .gt. v2 + 2.*(phin-phip) ) ) then
         istrapped=.true.
      else
         istrapped=.false.
      endif
c      write(*,*)'phin=',phin,'  phie=',phie,'  phip=',phip
c      write(*,*)'vte=',vte,' vtp=',vtp

      end
c********************************************************************
      function mtrapped()
      include 'piccom.f'
      logical istrapped

      mtrapped=0
      do i=1,iocprev
         if(ipf(i).gt.0) then
            if(istrapped(i)) mtrapped=mtrapped+1
         endif
      enddo

      end
