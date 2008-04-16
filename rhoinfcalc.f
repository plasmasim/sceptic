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

c***********************************************************************
c Calculate the density at infinity based on the number of injections
c per step. 
      subroutine rhoinfcalc(dt,icolntype,colnwt)
      real dt
      integer icolntype
      real colnwt,driftout,densred
      include 'piccom.f'
      include 'fvcom.f'
      real riave
      real fluxofangle(nthsize)
      real kappa
      save

c This allows us to restart with nstepsave .ne. 1 if rhoinf is set.
      if(riave.eq.0)riave=rhoinf
      if(finnerave.eq.0)finnerave=ninner

      if(nrein .gt. 0) then
c Combination equivalent to phihere usage.
c The following statement requires us to call this routine after 
c (part at least of) chargediag has been done. 
         averein=(diagphi(NRFULL)+diagphi(NRUSED))*.5

         if(averein.gt.0.5*Ti)then
c This is necessary to prevent smaxflux errors. smaxflux is not correct
c for repulsive potentials.
c            write(*,*)'Excessive averein',averein,' capped'
            averein=0.5*Ti
         endif

c     We have to calculate rhoinf consistently with the reinjection

c Begin with the standard case
         if(bcr.eq.0) then
            
         if(icolntype.eq.1.or.icolntype.eq.5) then
c Using general fvinject, injecting at the computational boundary.
c  qthfv(nthfvsize) contains the one-way flux density
c integrated dcos(theta) in 2Ti-normalized units.

            riest=(nrein/dt) /
     $           ((sqrt(2.*Ti)/2*qthfv(nthfvsize))
     $           *4.*3.141593*r(NRFULL)**2)

            if(bcphi.eq.4) then

c     bcphi4 is for the stationary continuum regime. Need the total flux
c     to use fluid boundary conditions on the density at the outer boundary
               totflux=0.
               do j=1,nthused
c     Calculate the flux to each angular cell
                  if(lcic)then
                     fluxofangle(j)=finthave(j)*(nthused-1.)/
     $                    (4.*pi*rhoinf*dt*r(1)**2)
                     if(j.eq.1 .or. j.eq.nthused)
     $                    fluxofangle(j)=fluxofangle(j)*2.
                  else
                     fluxofangle(j)=finthave(j)*(nthused)/
     $                    (4.*pi*rhoinf*dt*r(1)**2)
                  endif
                  totflux=totflux+fluxofangle(j)
               enddo
               totflux=totflux/nthused

c     For stationary strongly collisional plasmas, the outer
c     distribution function is, for rb>>rp, a radially shifted
c     Maxwellian with velocity driftout.
               driftout=totflux/rcc(nrused)**2
c     But riest is calculated assuming a stationary Maxwellian (valid at
c     nu=0). Hence rescale riest to account for Gamma\sim
c     Gamma(vd=0)+vd/2, with vd=totflux/rcc(nrused)**2
               densred=1+driftout*sqrt(pi)/sqrt(2*Ti)   
               kappa=sqrt(pi)*sqrt(2*Ti)/colnwt
c     However in the strongly collisional case the outer density is not 1:
c               nout=ninf*(1-Ti/(1+Ti)*Gamma/(kappa*Gamma0)/rcc(nrused))
               densred=densred*(1-Ti/(1+Ti)*totflux/(kappa*sqrt(2*Ti)/(2
     $              *sqrt(pi)))/rcc(nrused))
               if(densred<0.5) densred=0.5
               
               riest=riest/densred
            endif

         elseif(icolntype.eq.2 .or. icolntype.eq.6)then
c ogeninject from infinity with a general distribution numerically
c Flux/(2\pi riest v_n rmax^2) = pu1(1) - averein*pu2(1)
            riest=(nrein/dt) /
     $           (sqrt(2.*Ti)*
     $           2.*3.1415926**2*(pu1(1)-pu2(1)*averein/Ti)
     $           *r(NRFULL)**2 )
c Correction for collisional drag solution in outer region.
c            riest=riest
c     $           +finnerave*colnwt/(4.*3.1415926*dt*(1.+Ti)*r(NRFULL))
c last is correction for diffusion.
         else
c smaxflux returns total flux in units of Ti (not 2Ti)
            riest=(nrein/dt) /
     $           (sqrt(Ti)*
     $           smaxflux(vd/sqrt(2.*Ti),(-averein/Ti))
     $           *r(NRFULL)**2 )
         endif

c Now deal with custom reinject. For example if the neutral velocity is vd,
c we need to reinject a Maxwellian, not the collisional distribution
         else
            if (bcr.eq.1) then
c Set averein =0, because we should not consider the potential at the outer
c boundary for the reinjection
               averein=0
               riest=(nrein/dt) /
     $              (sqrt(Ti)*
     $              smaxflux(vd/sqrt(2.*Ti),(-averein/Ti))
     $              *r(NRFULL)**2 )
            elseif (bcr.eq.2) then
               riest=(nreintry/dt) /
     $              (sqrt(Ti)*
     $              smaxflux(vd/sqrt(2.*Ti),(-averein/Ti))
     $              *r(NRFULL)**2 )
c     write(*,*) riest
            endif
            
         endif

c         write(*,*)'nrein=',nrein
c     $        ,'  averein=',averein
c     $        ,' riest=',riest,' rhoinf=',rhoinf
      else
c If no valid estimate, just keep it the same.
c         write(*,*) 'nrein=0'
         if(.not. rhoinf.gt.1.e-4)then
            write(*,*)'rhoinf error in chargediag:',rhoinf,
     $           ' approximated.'
c estimate rhoinf approximately:
            rhoinf=numprocs*npart/(4.*pi*r(NRFULL)**3/3.)
         endif
         averein=0.
         riest=rhoinf
      endif
c Average the flux of particles to the inner, because it's a small number.
      finnerave=(finnerave*(nstepsave-1) + float(ninner))/nstepsave 
c Use an averaging process to calculate rhoinf (=riave)
      riave=(riave*(nstepsave-1) + riest)/nstepsave 
c      if(myid.eq.0) write(*,'(a,2f9.2,a,i6,f8.1,f8.1)')
c     $     'riest,riave',riest,riave,
c     $     '  ninner,finave,adj',ninner,finnerave,
c     $     finnerave*colnwt/(4.*3.1415926*dt*(1.+Ti)*r(NRFULL))
c*******
      rhoinf=riave
      end
