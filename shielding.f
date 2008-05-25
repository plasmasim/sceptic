c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
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
c
c     Version:  1.25   Sun May 25 19:10:02 EDT 2008
c
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___

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

