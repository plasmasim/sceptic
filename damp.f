
c***********************************************************************
c Calculate auxiliary damping potential from time-derivative of density.
c Adjust the total potential accordingly. 
      subroutine damp(dt,damplen)
c Time-step, dt; damping length, damplen, in normalized units.
c Algorithm is designed to give a damping time equal to damplen divided
c by the sound speed, at low frequency.
      real dt,damplen
c Common data:
      include 'piccom.f'
      real relax
c Local data:
c Prior density, and working density derivative.
      real rhoprevious(0:nrsize,0:nthsize)
c Adjustment to potential
      real potadj(0:nrsize,0:nthsize)
      logical lfirstcall
      data lfirstcall/.true./
      save

      if(damplen.eq.0. .or. dt.eq.0.) stop 'damp: zero argument(s)'
      wmu=debyelen**2/dt

c Set up arrays
      do j=1,NTHUSED
         do i=imin+1,NRFULL-1
            if(lfirstcall)then
               rhoprevious(i,j)=0.
               potadj(i,j)=0.
            else
               rhoprevious(i,j)=(rho(i,j)-rhoprevious(i,j))*wmu
            endif
         enddo
      enddo
      lfirstcall=.false.

c cic boundary is at i=1, ngp at 0+(1/2) (sort of).
      if(LCIC)then
         imin=1
      else
         imin=0
      endif

c   Boundary conditions.
      do j=1,NTHUSED
         potadj(imin,j)=0.
         potadj(NRFULL,j)=0.
      enddo     

c Chebychev acceleration. Wild guess at the Jacoby convergence radius.
c adjustable fudge factor to optimize convergence.
      rjfac=3.
      rjac=1.-rjfac*4./max(10,NRUSED)**2
      omega=1.
      maxits=2.5*NRUSED
c No need for very fine convergence here.
      dconverge=1.e-3

c SOR iteration.
      do k=1,maxits
         relax=omega
         deltamax=0.
c Alternate iteration directions for better convergence.
         if(mod(k/2,2).eq.0)then
         do j=1,NTHUSED
            do i=imin+1,NRFULL-1
               dnum= apc(i)*potadj(i+1,j)+bpc(i)*potadj(i-1,j)+cpc(i,j)
     $              *potadj(i,j+1)+dpc(i,j)*potadj(i,j-1)
     $              -fpc(i,j)*potadj(i,j)
     $              + rhoprevious(i,j)
               dden=fpc(i,j)
               delta=relax*dnum/dden
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               potadj(i,j)=potadj(i,j)+delta
            enddo
         enddo
         else
         do j=NTHUSED,1,-1
            do i=NRFULL-1,imin+1,-1
               dnum= apc(i)*potadj(i+1,j)+bpc(i)*potadj(i-1,j)+cpc(i,j)
     $              *potadj(i,j+1)+dpc(i,j)*potadj(i,j-1)
     $              -fpc(i,j)*potadj(i,j)
     $              + rhoprevious(i,j)
               dden=fpc(i,j)
               delta=relax*dnum/dden
               if(abs(delta).gt.abs(deltamax))deltamax=delta
               potadj(i,j)=potadj(i,j)+delta
            enddo
         enddo
         endif
         if(abs(deltamax).lt.dconverge.and.k.ge.2)goto 11
         if(k.eq.1)then
            omega=1./(1.-0.5*rjac**2)
         else
            omega=1./(1.-0.25*rjac**2*omega)
        endif
c        write(*,201)k,deltamax,relax
      enddo
c      write(*,*)'SOR not converged. deltamax=',deltamax
 11   continue
      write(*,'('':'',i3,$)')k
 201  format(' Damp SOR iteration',I4,' delta:',f10.6,' relax=',f8.4)
      do i=1,NRUSED
         potadj(i,0)=potadj(i,imin+1)
         potadj(i,NTHUSED+1)=potadj(i,NTHUSED-imin)
      enddo
      pamax=0.
      paave=0.
      ipa=0
      do j=1,NTHUSED
         do i=imin+1,NRFULL-1
            rhoprevious(i,j)=rho(i,j)
            phi(i,j)=phi(i,j)+potadj(i,j)/damplen
            ipa=ipa+1
            paave=paave+potadj(i,j)
            if(abs(potadj(i,j)).gt.pamax) pamax=abs(potadj(i,j))
         enddo 
      enddo 
      paave=paave/ipa
c      write(*,'('' pamax='',f8.4,''  paave='',f8.5)')pamax,paave
      end
c*******************************************************************
