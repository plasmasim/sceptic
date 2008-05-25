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

c******************************************************************

c Do the SOR iteration in parallel

      subroutine fcalc_shielding_par(dt,rshield,sor_comm,myid2)

      include 'piccom.f'
      real dt
c     If rshield is not an integer, problem with the MPI routines
      integer rshield
c     sor_comm is the subset of MPI_COMM_WORLD communicator used for
c     the bloc sor
      integer sor_comm
c Correct the type of temporary variable (IHH). Better use initial letter.
      integer c

c Chebychev acceleration. Wild guess at the Jacoby convergence radius.

      rjac=1.-4./max(10,NRUSED)**2
      omega=1.
      maxits=2.5*NRUSED
      dconverge=1.e-5
      imin=1
      relax=0.5

c     fcalc_sor call : Do the iterations in parallel nrsize+1 is the 1st
c     dimension of the mesh (phi, rho, ...)  rshield+1 is 1+ the r
c     dimension really used. From 1 to rshield, but we add 1 to ensure
c     bloc communications of the BC at each iteration NTHUSED+2 is the
c     th dimension really used


      call fcalc_sor(nrsize+1,nthsize+1,rshield+1,NTHUSED+2,rjac,
     $     maxits,dconverge,k,sor_comm,myid2)

      
      
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

