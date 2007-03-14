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
      subroutine chargetomesh()
c Common data:
      include 'piccom.f'

c Corretion factor to calculate the density in the inner domain
      real corr

c      ninner=0
      do i=0,nr
         do j=0,nth+1
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
c Perhaps this needs to be larger than npart for .not.lfixed.
c      write(*,*)'Starting chargetomesh',npart
      do i=1,iocprev
         if(ipf(i).gt.0)then
c         if(i.lt.10000)write(*,'(i6,$)')i
c Use fast ptomesh, half-quantities not needed.
            ih=0
            hf=99.
            call ptomesh(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp
     $           ,zetap,ih,hf)
            if(rf.lt.0..or.rf.gt.1.)then
               rp=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
               write(*,*)'Outside mesh, rf error in chargetomesh',
     $              rf,irl,i,rp
            else
               call chargeassign(i,irl,rf,ithl,thf,
     $              ipl,pf,st,ct,sp,cp,rp)
            endif
         endif
      enddo

c Correction for --ds
      
      if(dsub) then
         do ir=0,nr
            corr=1.
            if (ir.lt.rsplit) then 
               corr=1.*nrealin/(nrealin+npartadd)
            elseif (ir.eq.rsplit) then
               dr=0.5*(r(rsplit)-r(rsplit-1))
               cvol=( r(rsplit)**3 - (r(rsplit)-dr)**3 )/
     $              ( (r(rsplit)+dr)**3 - (r(rsplit)-dr)**3 )
               corr=1.*nrealin/(nrealin+npartadd*cvol)
            endif

            do ith=0,nth+1
               psum(ir,ith)=psum(ir,ith)*corr
               vrsum(ir,ith)=vrsum(ir,ith)*corr
               vtsum(ir,ith)=vtsum(ir,ith)*corr
               vpsum(ir,ith)=vpsum(ir,ith)*corr
               vr2sum(ir,ith)=vr2sum(ir,ith)*corr
               v2sum(ir,ith)=v2sum(ir,ith)*corr
               vtp2sum(ir,ith)=vtp2sum(ir,ith)*corr
               vzsum(ir,ith)=vzsum(ir,ith)*corr
            enddo
         enddo
      endif
      
      
      end
c***********************************************************************
c Accumulate particle charge into rho mesh and other diagnostics.
      subroutine chargeassign(i,irl,rf,ithl,thf,ipl,pf,st,ct,sp,cp,rp)
c      implicit none
      integer i
c Common data:
      include 'piccom.f'
c Assign as if square for now. Area weighting might be better.
c Charge summation.
      psum(irl,ithl)=psum(irl,ithl) + (1.-rf)*(1.-thf)
      psum(irl+1,ithl)=psum(irl+1,ithl) + rf*(1.-thf)
      psum(irl,ithl+1)=psum(irl,ithl+1) + (1.-rf)*thf
      psum(irl+1,ithl+1)=psum(irl+1,ithl+1) + rf*thf

      if(diags .or. irl.le.2) then
      vxy=xp(4,i)*cp + xp(5,i)*sp
      vr=vxy*st + xp(6,i)*ct
      vrsum(irl,ithl)=vrsum(irl,ithl) + (1.-rf)*(1.-thf)*vr
      vrsum(irl+1,ithl)=vrsum(irl+1,ithl) + rf*(1.-thf)*vr
      vrsum(irl,ithl+1)=vrsum(irl,ithl+1) + (1.-rf)*thf*vr
      vrsum(irl+1,ithl+1)=vrsum(irl+1,ithl+1) + rf*thf*vr

      vt= vxy*ct - xp(6,i)*st
      vtsum(irl,ithl)=vtsum(irl,ithl) + (1.-rf)*(1.-thf)*vt
      vtsum(irl+1,ithl)=vtsum(irl+1,ithl) + rf*(1.-thf)*vt
      vtsum(irl,ithl+1)=vtsum(irl,ithl+1) + (1.-rf)*thf*vt
      vtsum(irl+1,ithl+1)=vtsum(irl+1,ithl+1) + rf*thf*vt

      vp=-xp(4,i)*sp  + xp(5,i)*cp
      vpsum(irl,ithl)=vpsum(irl,ithl) + (1.-rf)*(1.-thf)*vp
      vpsum(irl+1,ithl)=vpsum(irl+1,ithl) + rf*(1.-thf)*vp
      vpsum(irl,ithl+1)=vpsum(irl,ithl+1) + (1.-rf)*thf*vp
      vpsum(irl+1,ithl+1)=vpsum(irl+1,ithl+1) + rf*thf*vp

      vr2=vr*vr
      vr2sum(irl,ithl)=vr2sum(irl,ithl) + (1.-rf)*(1.-thf)*vr2
      vr2sum(irl+1,ithl)=vr2sum(irl+1,ithl) + rf*(1.-thf)*vr2
      vr2sum(irl,ithl+1)=vr2sum(irl,ithl+1) + (1.-rf)*thf*vr2
      vr2sum(irl+1,ithl+1)=vr2sum(irl+1,ithl+1) + rf*thf*vr2

      v2=(xp(4,i)*xp(4,i) +xp(5,i)*xp(5,i) +xp(6,i)*xp(6,i))
      v2sum(irl,ithl)=v2sum(irl,ithl) + (1.-rf)*(1.-thf)*v2
      v2sum(irl+1,ithl)=v2sum(irl+1,ithl) + rf*(1.-thf)*v2
      v2sum(irl,ithl+1)=v2sum(irl,ithl+1) + (1.-rf)*thf*v2
      v2sum(irl+1,ithl+1)=v2sum(irl+1,ithl+1) + rf*thf*v2

      vtp2=v2-vr2
      vtp2sum(irl,ithl)=vtp2sum(irl,ithl) + (1.-rf)*(1.-thf)*vtp2
      vtp2sum(irl+1,ithl)=vtp2sum(irl+1,ithl) + rf*(1.-thf)*vtp2
      vtp2sum(irl,ithl+1)=vtp2sum(irl,ithl+1) + (1.-rf)*thf*vtp2
      vtp2sum(irl+1,ithl+1)=vtp2sum(irl+1,ithl+1) + rf*thf*vtp2

      vz=xp(6,i)
      vzsum(irl,ithl)=vzsum(irl,ithl) + (1.-rf)*(1.-thf)*vz
      vzsum(irl+1,ithl)=vzsum(irl+1,ithl) + rf*(1.-thf)*vz
      vzsum(irl,ithl+1)=vzsum(irl,ithl+1) + (1.-rf)*thf*vz
      vzsum(irl+1,ithl+1)=vzsum(irl+1,ithl+1) + rf*thf*vz
      endif
      end
c***********************************************************************
c Calculate potential phi from rho.
      subroutine fcalc(dt)
c Common data:
      include 'piccom.f'
      real phi0mphi1(0:nthsize),delphi0(0:nthsize)
      parameter (nthp1=nthsize+1)
c      real phi1ave
      real bcifac,bcpfac,bci,bcp,bvf
      real relax
      real cs(nthsize),csd(nthsize)
      real ncs
      logical first
      data relax/1./
      data bcifac/.2/bcpfac/.1/
      data bvf/1.2071/
      data first/.true./
      data phi0mphi1/nthp1*0./
      data delphi0/nthp1*0./
      data ncs/50./
      save
      cerr=0.

      do j=1,nth
         do i=1,nr
            if(rho(i,j).le.0.)then
               write(*,*)'rho=0',i,j
               stop
            endif
c Simplistic Boltzmann scheme. May need relaxation.
            delta=phi(i,j)-log(rho(i,j))
            if(abs(delta).le.cerr)cerr=abs(delta)
            phi(i,j)=phi(i,j)-relax*delta
         enddo
      enddo
c Probe boundary condition.
      if(first)then
         do j=1,nth
            cs(j)=-sqrt(1.+Ti)
         enddo
         first=.false.
      endif
c      write(*,*)'p1,p2,v1,v2,csd,cs,vs,phi0,phi1,delphinew'
c      write(*,*)'cs=',(cs(kk),kk=1,nth)
      do j=1,nth
         p1=vr2sum(1,j)*psum(1,j)-vrsum(1,j)**2
         v1=vrsum(1,j)
         if(p1.ne.0)then
            if(psum(1,j).gt.0)then
               p1=p1/psum(1,j)
               v1=v1/psum(1,j)
            else
               write(*,*)'psum(1,j)=0'
               stop
            endif
         endif
         p2=vr2sum(2,j)*psum(2,j)-vrsum(2,j)**2
         v2=vrsum(2,j)
         if(p2.ne.0)then
            if(psum(2,j).gt.0)then
               p2=p2/psum(2,j)
               v2=v2/psum(2,j)
            else
               write(*,*)'psum(2,j)=0'
               stop
            endif
         endif
         vs=(1.+bvf)*v1-bvf*v2
c Fix the psum difference so it can't be zero.
c Original derivative form
         csd(j)=(p2 - p1)/(psum(2,j)-psum(1,j)+.5)+1.
c Gamma Ti form with gamma=3.
c         csd(j)=3.*p1/(psum(1,j)+.5) + 1.
         if(csd(j).lt.0.)csd(j)=0.
         csd(j)=-sqrt(csd(j))
c Clip the excursion symmetrically.
         if(csd(j).lt.2.*cs(j))  csd(j)=2.*cs(j)
c Average the sound-speed value over ncs steps.
         cs(j)=((ncs-1.)*cs(j)+csd(j))/ncs
         if(cs(j).gt.0)then
            write(*,*)'cs positive',j,cs(j),csd(j)
     $           ,p1,p2,psum(1,j),psum(2,j),ncs
            stop
         endif
         bci=-bcifac*cs(j)**2*dt/(r(2)-r(1))
         bcp=-bcpfac*(r(2)-r(1))/(dt*cs(j))
         delphinew=(vs - cs(j))*bci
         phi(0,j)=phi(1,j)+phi0mphi1(j)+delphinew+
     $        (delphinew-delphi0(j))*bcp
         delphi0(j)=delphinew
         if(phi(0,j).gt.phi(1,j))phi(0,j)=phi(1,j)
         phi0mphi1(j)=phi(0,j)-phi(1,j)
c         write(*,'(10f8.3)')p1,p2,v1,v2,csd(j),cs(j),vs,
c     $        phi(0,j),phi(1,j),delphinew
c Adjusting the potential of the first cell.



         phi(1,j)=phi(0,j)



         if(.not.abs(phi(1,j)).lt.1.e20)then
            write(*,*)'phi1 overflow',phi(1,j),bcp,cs(j),ncs,delphinew
     $           ,p1,p2,psum(1,j),psum(2,j),csd(j)
            stop
         endif
      enddo

      do i=1,nr
         phi(i,0)=phi(i,1)
         phi(i,nth+1)=phi(i,nth)
      enddo
c      write(*,*)'At end cs=',(cs(kk),kk=1,nth)
c      write(*,'(10f7.4)')((phi(i,ih),i=0,9),ih=0,9)
      cerr=1.
      end
c***********************************************************************
c cic version.
      subroutine getaccel(i,accel,il,rf,ith,tf,ipl,pf,
     $     st,ct,sp,cp,rp,zetap,ih,hf)
c Evaluate the cartesian acceleration into accel. Using half-mesh
c parameters.
c accel is minus the gradient of phi for the ith particle.
c Be careful with variables in this routine.

      implicit none
      integer i
      real accel(3)
c Common data:
      include 'piccom.f'
      real ar,at
      real ct,st,cp,sp,rp
      real zetap,hf
      integer ih
      real dp,dth
c      parameter (dp=2.*pi/np,dth=pi/(nth-1))
      real dpinv,dthinv
c      parameter (dpinv=np/2./pi,dthinv=(nth-1)/pi)
      integer ipl,ith
      integer il,ir,ithp1,ithp2,ithm1,ilm1
      real rlm1,rf,tf,pf,rr,rl
      real philm1t,philm1p,phihp1t,phihp1p
      real phihp1m,phihp12
      data dp/0./

c Don't use parameter statements. Just set values first time
      if(dp.eq.0.)then
         dp=2.*pi/np
         dth=pi/(nth-1)
         dpinv=np/2./pi
         dthinv=(nth-1)/pi
      endif

c Set up indexes. Here we are using the half-mesh for r.

      rl=r(ih)
      ir=ih+1
c      if(hf.gt.1. .or. hf.lt.0.)write(*,*)'hf=',hf

      ithp1=ith+1
      ithp2=ith+2
      ithm1=ith-1
c Deal with r-boundary conditions.
c Constant slope at r-boundary. (Zero second derivative).
      if(ih.eq.nr)then
         phihp1t=phi(ih,ith)+(phi(ih,ith)-phi(ih-1,ith))
         phihp1p=phi(ih,ithp1)+(phi(ih,ithp1)-phi(ih-1,ithp1))
         phihp1m=phi(ih,ithm1)+(phi(ih,ithm1)-phi(ih-1,ithm1))
         phihp12=phi(ih,ithp2)+(phi(ih,ithp2)-phi(ih-1,ithp2))
         rr=rl+(rl-r(ih-1))
      else
         phihp1t=phi(ir,ith)
         phihp1p=phi(ir,ithp1)
         phihp1m=phi(ir,ithm1)
         phihp12=phi(ir,ithp2)
         rr=r(ir)
      endif
      ilm1=ih-1
c Here we control whether we use the zeta or r interpolation.
      if(debyelen.lt.1.e-2)then
c      if(.true.)then
c Linear approx to sqrt form at boundary.
         if(ih.eq.1)then
c            rr=r(ir)
            philm1t=phi(il,ith)-bdyfc*sqrt(2.*(rr-rl))*0.25
            philm1p=phi(il,ithp1)-bdyfc*sqrt(2.*(rr-rl))*0.25
            rlm1=2.*rl - rr
c     Constant slope
         else
            philm1t=phi(ilm1,ith)
            philm1p=phi(ilm1,ithp1)
            rlm1=r(ilm1)
         endif
c Uniform interpolation of d\phi/d\zeta times 1/\zeta
         if(zetap.le.1.e-2)zetap=1.e-2
         ar=(  ( (phihp1t-phi(ih,ith))/(zeta(ir)-zeta(ih))*hf +
     $        (phi(ih,ith)-philm1t)/(zeta(ih)-zeta(ilm1))*(1.-hf)
     $           )*(1.-tf)
     $        +( (phihp1p-phi(ih,ithp1))/(zeta(ir)-zeta(ih))*hf +
     $        (phi(ih,ithp1)-philm1p)/(zeta(ih)-zeta(ilm1))*(1.-hf)
     $         )*tf )/zetap
      else
c Interpolate in r.
         philm1t=phi(ilm1,ith)
         philm1p=phi(ilm1,ithp1)
         rlm1=r(ilm1)
         ar=(  ( (phihp1t-phi(ih,ith))/(rr-rl)*hf +
     $        (phi(ih,ith)-philm1t)/(rl-rlm1)*(1.-hf)
     $           )*(1.-tf)
     $        +( (phihp1p-phi(ih,ithp1))/(rr-rl)*hf +
     $        (phi(ih,ithp1)-philm1p)/(rl-rlm1)*(1.-hf)
     $         )*tf )
c
      endif
         
      ar=-ar

      if(tf.le.0.5)then
         at= ( (phi(ih,ithp1)-phi(ih,ith))*(tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith)))
     $        +(phi(ih,ithp1)-phi(ih,ithm1))*(0.5-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ithm1))) ) * (1.-rf)
     $        + ( (phihp1p-phihp1t)*(tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith)))
     $        +(phihp1p-phihp1m)*(0.5-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ithm1))) ) * rf
      else
         at= ( (phi(ih,ithp2)-phi(ih,ith))*(tf-0.5)*2.
     $          /(rl*(thang(ithp2)-thang(ith)))
     $        +(phi(ih,ithp1)-phi(ih,ith))*(1.-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith))) ) * (1.-rf)
     $        + ( (phihp12-phihp1t)*(tf-0.5)*2.
     $          /(rl*(thang(ithp2)-thang(ith)))
     $        +(phihp1p-phihp1t)*(1.-tf)*2.
     $          /(rl*(thang(ithp1)-thang(ith))) ) * rf
      endif
      
      at=-at
c      ap=0.
      if(lat0)at=0.

 501  format(a,6f10.4)
      accel(3)=ar*ct - at*st
      accel(2)=(ar*st+ at*ct)*sp
      accel(1)=(ar*st+ at*ct)*cp
c Trap errors.

      if(.not.accel(1).lt.1.e5)then
         write(*,*) 'i: ',i,' x: ',xp(1,i),' y: ',xp(2,i),' z: ',xp(3,i)
         write(*,*)'Accel Excessive: ar,at,st,ct,ih,hf,rf,ith,tf'
         write(*,*) ar,at,st,ct,ih,hf,rf,ith,tf
         write(*,*) 'phi at ith and ithp1:'
         write(*,*) philm1t,phi(ih,ith),phihp1t
         write(*,*) philm1p,phi(ih,ithp1),phihp1p
         write(*,*) 'zetap=',zetap,'  bdyfc=',bdyfc
         write(*,*) 'zeta '
         write(*,*) zeta(ilm1),zeta(ih),zeta(ir)
         write(*,'(10f7.4)')((phi(i,ih),i=1,10),ih=1,10)
         stop
      endif
      end
c***********************************************************************

c**********************************************************************
      subroutine esforce(ir,qp,fz,epz)
      include 'piccom.f'
c Version to specify radius node ir at which to calculate force.
c Return the charge qp, esforce fz, and electron pressure force epz.
      
      real ercoef(nthsize),etcoef(nthsize),ertcoef(nthsize)
      real qpcoef(nthsize)

      logical lnotinit
      data lnotinit/.true./
      save

      if(lnotinit)then
c Initialize coefficient arrays
         do j=1,nthused
            if(j.eq.1)then
               qpcoef(j)=0.5*(th(j+1)-th(j))
               ercoef(1)=th(2)*(th(2)+th(1))/2.
     $              - (th(2)*(th(2)+th(1))+th(1)*th(1))/3.
               etcoef(j)=(th(j+1)**2/2.-th(j+1)**4/4.
     $              -th(j)**2/2.+th(j)**4/4.)/(th(j+1)-th(j))**2
               ertcoef(j)=1.-(th(j+1)*(th(j+1)+th(j))+th(j)**2)/3.
            elseif(j.eq.nthused)then
               qpcoef(j)=0.5*(th(j)-th(j-1))
               etcoef(j)=0.
               ertcoef(j)=0.
               ercoef(j)=(th(j)*(th(j)+th(j-1))+th(j-1)**2)/3.
     $              - th(j-1)*(th(j)+th(j-1))/2.
            else
               qpcoef(j)=(th(j+1)-th(j))
               ercoef(j)=(th(j+1)*(th(j+1)+th(j))
     $              -th(j-1)*(th(j)+th(j-1)))/6.
               etcoef(j)=(th(j+1)**2/2.-th(j+1)**4/4.
     $              -th(j)**2/2.+th(j)**4/4.)/(th(j+1)-th(j))**2
               ertcoef(j)=1.-(th(j+1)*(th(j+1)+th(j))+th(j)**2)/3.
            endif
c     Multiply by calibrations, but not radius factors.
            qpcoef(j)=-qpcoef(j)*2.*pi
            ercoef(j)=-0.5*ercoef(j)*2.*pi
            etcoef(j)=-0.5*etcoef(j)*2.*pi
            ertcoef(j)=-ertcoef(j)*2.*pi
         enddo
c         write(*,*)'qpcoef,ercoef, etcoef,ertcoef'
c         write(*,'(i3,4f10.5)')(j,qpcoef(j),ercoef(j),etcoef(j),
c     $        ertcoef(j),j=1,nthused)

         lnotinit=.false.
      endif

c Calculate charge and electrostatic force on probe surface
c We do this at a specified radius node.
      k=ir
      if(k.gt.nrused)then
         write(*,*)'esforce radial node number too large. Reset.'
         k=nrused
      endif
      delr=(r(2)-r(1))
      j=0
      if(k.eq.1)then
         rkp=0.5*(rcc(k)+rcc(k+1))
         rkp2=0.5*(rcc(k+1)+rcc(k+2))
            erp=-(rkp*rkp*(phi(k+1,j+1)-phi(k,j+1))*(1+.5)-
     $         .5*rkp2*rkp2*(phi(k+2,j+1)-phi(k+1,j+1)))/delr
c         erp=-(-.5*phi(k+2,1)+2.*phi(k+1,1)-1.5*phi(k,1))/delr
      elseif(k.eq.nrused)then
         rkm=0.5*(rcc(k)+rcc(k-1))
         rkm2=0.5*(rcc(k-1)+rcc(k-2))
         erp=-(rkm*rkm*(phi(k,j+1)-phi(k-1,j+1))*(1+.5)-
     $        .5*rkm2*rkm2*(phi(k-1,j+1)-phi(k-2,j+1)))/delr
c         erp=-( .5*phi(k-2,1)-2.*phi(k-1,1)+1.5*phi(k,1))/delr
      else
         rkp=0.5*(rcc(k)+rcc(k+1))
         rkm=0.5*(rcc(k)+rcc(k-1))
            erp=-0.5*(rkp*rkp*(phi(k+1,j+1)-phi(k,j+1))+
     $           rkm*rkm*(phi(k,j+1)-phi(k-1,j+1)))/delr
c         erp=-0.5*(phi(k+1,1)-phi(k-1,1))/delr
      endif
      qp=qpcoef(1)*erp
      fz=ercoef(1)*erp*erp/rcc(k)**2
      epz=ercoef(1)*2.*exp(phi(k,1))
      do j=1,nthused-1
c     Radial field extrapolated to index position if necessary.
         er=erp
         if(k.eq.1)then
c            erp=-(-.5*phi(k+2,j+1)+2.*phi(k+1,j+1)-1.5*phi(k,j+1))/delr
            erp=-(rkp*rkp*(phi(k+1,j+1)-phi(k,j+1))*(1+.5)-
     $         .5*rkp2*rkp2*(phi(k+2,j+1)-phi(k+1,j+1)))/delr
         elseif(k.eq.nrused)then
c            erp=-( .5*phi(k-2,j+1)-2.*phi(k-1,j+1)+1.5*phi(k,j+1))/delr
            erp=-(rkm*rkm*(phi(k,j+1)-phi(k-1,j+1))*(1+.5)-
     $         .5*rkm2*rkm2*(phi(k-1,j+1)-phi(k-2,j+1)))/delr
         else
            erp=-0.5*(rkp*rkp*(phi(k+1,j+1)-phi(k,j+1))+
     $           rkm*rkm*(phi(k,j+1)-phi(k-1,j+1)))/delr
         endif
         qp=qp+qpcoef(j+1)*erp
         fz=fz+ercoef(j+1)*erp*erp/rcc(k)**2
         eth=phi(k,j+1)-phi(k,j)
         fz=fz+etcoef(j)*eth*eth
         fz=fz - ertcoef(j)*0.5*(erp+er)*eth/rcc(k)
         epz=epz + ercoef(j+1)*2.*exp(phi(k,j+1))
c      write(*,*)j,erp,eth,ercoef(j+1)*erp*erp,etcoef(j)*eth*eth,
c     $           - ertcoef(j)*0.5*(erp+er)*eth
c      write(*,*)fz
      enddo
      fz=fz
      qp=qp
      epz=-epz*r(k)**2
      end
