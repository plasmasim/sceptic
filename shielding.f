c******************************************************************
      subroutine fcalc_shielding(dt,icolntype,colnwt,sor_comm,myid2)

      include 'piccom.f'
      include 'mpif.h'
      real dt
      real relax
      real delta
c     If rshield is not an integer, problem with the MPI routines
      integer rshield
c     sor_comm is the subset of MPI_COMM_WORLD communicator used for
c     the bloc sor
      integer sor_comm
      real phislopeconst(0:nth+1),phislopefac(0:nth+1)
c Correct the type of temporary variable (IHH). Better use initial letter.
      integer c

c Chebychev acceleration. Wild guess at the Jacoby convergence radius.

      rjac=1.-4./max(10,NRUSED)**2
      omega=1.
      maxits=2.5*NRUSED
      dconverge=1.e-5
      imin=1
      relax=0.5


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


      call innerbc(imin,dt)
      
c     Do the iterations in serial
      if (.not.sorparallel.and.myid2.eq.0) then
c      if (.true.) then
c SOR iteration
         do k=1,maxits
c Use over-relaxation if debyelen is large, or straight newton otherwise.
            relax=(omega*debyelen**2+1.)/(debyelen**2+1.)
            deltamax=0.
c Changed to chessboard calculation. Twice as fast
            c=mod(k,2)


            do j=1,NTHUSED
c     We only go to rshield-1
               do i=imin+1+mod(j+c,2),rshield,2

c     Boundary at rshield
                  if (i.eq.rshield) then
                     if(Ezext.eq.0)then
                        delta=gpc(j,1)*phi(i-1,j)+gpc(j,2)*phi(i,j-1)+
     $                   gpc(j,3)*phi(i,j+1)+gpc(j,4)+gpc(j,5)*phi(i,j)
                        if(abs(delta).gt.abs(deltamax))deltamax=delta
                        phi(i,j)=phi(i,j)+relax*delta
                     else
c Implement edge-potential-control if set.
                        phi(i,j)=Ezext*tcc(j)*r(i)
                     endif
                  else
                     expphi=exp(phi(i,j))
                    dnum= apc(i)*phi(i+1,j)+bpc(i)*phi(i-1,j) + cpc(i,j)
     $               *phi(i,j+1)+dpc(i,j)*phi(i,j-1) -fpc(i,j)*phi(i,j)
     $                    + rho(i,j) - expphi
                     dden=fpc(i,j) + expphi
                     delta=relax*dnum/dden
                     if(abs(delta).gt.abs(deltamax))deltamax=delta
                     phi(i,j)=phi(i,j)+delta
                  endif

               enddo           
            enddo

            if(abs(deltamax).lt.dconverge.and.k.ge.2) goto 11
            if(k.eq.1)then
               omega=1./(1.-0.5*rjac**2)
            else
               omega=1./(1.-0.25*rjac**2*omega)
            endif
         enddo

      elseif(sorparallel) then
c fcalc_sor call : Do the iterations in parallel
c nrsize+1 is the 1st dimension of the mesh (phi, rho, ...)
c rshield+1 is 1+ the r dimension really used. From 1 to rshield, but we add
c 1 to ensure bloc communications of the BC at each iteration
c NTHUSED+2 is the th dimension really used
         call fcalc_sor(nrsize+1,nthsize+1,rshield+1,NTHUSED+2,rjac,
     $        maxits,dconverge,k,sor_comm)
     
      endif
 11   continue

      if (myid2.eq.0) then
         write(*,'('':'',i3,$)')k
c     write(*,201)k,deltamax,relax
 201     format(' SOR iteration',I4,' delta:',f10.6,' relax=',f8.4)
c     Calculate electric force on probe. Moved to main.
      endif

c     Inner Boundary values
      do j=1,NTHUSED
         phi(0,j)=2.*phi(imin,j)-phi(imin+1,j)
      enddo
      do i=1,NRUSED
         phi(i,0)=phi(i,imin+1)
         phi(i,NTHUSED+1)=phi(i,NTHUSED-imin)
      enddo
c     write(*,*)'phi(rmax)=',phi(NRFULL,NTHUSED/2)

 123  continue

      end 

