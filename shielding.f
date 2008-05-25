
      subroutine fcalc_shielding(dt,rshield)

c Serial chess-board SOR solver

      include 'piccom.f'
      real dt
      real relax
      real delta
c     If rshield is not an integer, problem with the MPI routines
      integer rshield
c Correct the type of temporary variable (IHH). Better use initial letter.
      integer c

c Chebychev acceleration. Wild guess at the Jacoby convergence radius.

      rjac=1.-4./max(10,NRUSED)**2
      omega=1.
      maxits=2.5*NRUSED
      dconverge=1.e-5
      imin=1
      relax=0.5


c     Start SOR iteration
      do k=1,maxits
c     Use over-relaxation if debyelen is large, or straight newton otherwise.
         relax=(omega*debyelen**2+1.)/(debyelen**2+1.)
         deltamax=0.
c     Changed to chessboard calculation. Twice as fast
         c=mod(k,2)
         

         do j=1,NTHUSED
c     We only go to rshield-1
            do i=imin+1+mod(j+c,2),rshield,2
               
c     Boundary at rshield
               if (i.eq.rshield) then
                  if(Ezext.eq.0)then
                     delta=gpc(j,1)*phi(i-1,j)+gpc(j,2)*phi(i,j-1)+
     $                    gpc(j,3)*phi(i,j+1)+gpc(j,4)+gpc(j,5)*phi(i,j)
                     if(abs(delta).gt.abs(deltamax))deltamax=delta
                     phi(i,j)=phi(i,j)+relax*delta
                  else
c     Implement edge-potential-control if set.
                     phi(i,j)=Ezext*tcc(j)*r(i)
                  endif
               else
                  expphi=exp(phi(i,j))
                  dnum= apc(i)*phi(i+1,j)+bpc(i)*phi(i-1,j) + cpc(i,j)
     $                 *phi(i,j+1)+dpc(i,j)*phi(i,j-1) -fpc(i,j) *phi(i
     $                 ,j)+ rho(i,j) - expphi
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
 11   continue
      

      write(*,'('':'',i3,$)')k
c     write(*,201)k,deltamax,relax
 201  format(' SOR iteration',I4,' delta:',f10.6,' relax=',f8.4)
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

