c___________________________________________________________________________
c
c This code is copyright (c) Ian H Hutchinson (2004) hutch@psfc.mit.edu.
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build or run
c the code, you do so at your own risk.
c___________________________________________________________________________
c Version 2.5: Jan 2005

C Main program for pic code
      program sceptic

c     For now, put maxsteps in piccom
      integer ipstep
      logical success,finaldiags
      real rmax
      character*100 string
      logical readpart,writepart
      logical lsmooth,lcolcont,lpstore
      logical rmaxauto,lfcalcl

c Communicator and id for the sor communicator
c In the case of sorparallel=false, myid2=0 only
c does the potential calculation
      integer sor_comm,myid2
      
c Common storage
      include 'piccom.f'
      include 'colncom.f'
      data maxsteps/500/
      data rmax/5./
      data dtf/0.025/bdt/1./
      data success/.false./
      data readpart/.false./
      data writepart/.false./
      data sorparallel/.false./
      data lfcalcl/.false./
      data rmaxauto/.false./
c Parallel processing MPI options.
#ifdef MPI
      include 'mpif.h'
      double precision sortime
      integer rc
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )

#else
      myid=0
      numprocs=1
#endif

c Full size arrays by default. Can be changed later by switches.
      nr=nrsize
      nth=nthsize
      if(LCIC)nth=nthsize-1
c Common defaults. Avoid block data.
      vd=.000
      debyelen=.1
      vprobe=-4.
      Ezext=0.
      bdyfc=4.
      Ti=1.
      diags=.true.
      finaldiags=.true.
      ldist=.false.
      avelim=0.6
      nstepsave=0
      nsamax=40
      lsmooth=.false.
      lcolcont=.true.
      ircell=1
      itcell=1
      norbits=0
      npart=npartmax
      ninjcomp0=npartmax
      ninjcomp=ninjcomp0
      linsulate=.false.
      lfloat=.false.
      rmtoz=1.
      ierad=1
      rhomin=0.
      rhomax=0.
      lat0=.false.  

      lfext=.false.
      colnwt=0.
      icolntype=0
      localinj=.false.
      damplen=0.
      ipfsw=3
      call pfset(0)
      lfixedn=.true.
      ipstep=1
      vneutral=0.
c     bcphi : BC on potential
      bcphi=0
c     bcr : BC on reinjection
      bcr=0
      dsub=.false.
      Bz=0.
      infdbl=.false.
      orbinit=.false.

c Deal with arguments
      if(iargc().eq.0) goto 51
      do 1 i=1,iargc()
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:8) .eq. '--rhoinf')then
            read(string(9:),*)rhoinf
            lfixedn=.false.
         endif
         if(string(1:2) .eq. '-r') readpart=.true.
         if(string(1:2) .eq. '-w') writepart=.true.
         if(string(1:3) .eq. '-gc') then
            lcolcont=.false.
         elseif(string(1:3) .eq. '-ge') then
         elseif(string(1:2) .eq. '-g')then
            if(string(3:3).ne.' ')then
               read(string(3:),*,err=259,end=259)ipstep
            else
               diags=.false.
            endif
            goto 260
 259        ipstep=1
 260        continue
         endif
         if(string(1:7) .eq. '--fcalc') lfcalcl=.true.
         if(string(1:2) .eq. '-f') finaldiags=.false.
         if(string(1:3) .eq. '-er') then
            read(string(4:),*)ierad
         elseif(string(1:2) .eq. '-e') then
            ldist=.true.
c try to read the cell from the rest of the string.
            read(string(3:),*,err=261,end=261)ircell,itcell
            goto 262
c restore defaults on read error.
 261        ircell=1
            itcell=1
 262        continue
         endif
         if(string(1:4) .eq. '--bz') read(string(5:),*)Bz
         if(string(1:4) .eq. '--ds') dsub=.true.
         if(string(1:4) .eq. '--sp') sorparallel=.true.
         if(string(1:3) .eq. '-sa') then
            rmaxauto=.true.
         elseif(string(1:2) .eq. '-s') then
            read(string(3:),*)maxsteps
         endif
         if(string(1:2) .eq. '-b') read(string(3:),*)bdyfc
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
         if(string(1:6) .eq. '--damp')read(string(7:),*)damplen
         if(string(1:2) .eq. '-d') read(string(3:),*)dtf
         if(string(1:5) .eq. '-cmin')then
            read(string(6:),*)rhomin
         elseif(string(1:5) .eq. '-cmax')then
            read(string(6:),*)rhomax
         elseif(string(1:2) .eq. '-c') then
            read(string(3:),*)bdt
         endif
         if(string(1:2) .eq. '-x') read(string(3:),*)rmax
         if(string(1:3) .eq. '-vn') then
            read(string(4:),*)vneutral
         elseif(string(1:2) .eq. '-v')then
            read(string(3:),*)vd
         endif
         if(string(1:3) .eq. '-pf') then
            lfloat=.true.
         elseif(string(1:3) .eq. '-pi') then
            linsulate=.true.
         elseif(string(1:2) .eq. '-p') then
            read(string(3:),*)vprobe
         endif
         if(string(1:2) .eq. '-z') read(string(3:),*)Ezext
         if(string(1:2) .eq. '-l') read(string(3:),*)debyelen
         if(string(1:3) .eq. '-nr') read(string(4:),*)nr
         if(string(1:3) .eq. '-nt') read(string(4:),*)nth
         if(string(1:3) .eq. '-ni') read(string(4:),*)npart
         if(string(1:2) .eq. '-m') read(string(3:),*)rmtoz
         if(string(1:2) .eq. '-a') then
            read(string(3:),*,err=251,end=251)ipfsw
            call pfset(ipfsw)
            goto 252
 251        ipfsw=0
 252        continue
         endif
         if(string(1:5) .eq. '--bcr') read(string(6:),*)bcr
         if(string(1:7) .eq. '--bcphi') read(string(8:),*)bcphi
         if(string(1:7) .eq. '-infdbl') infdbl=.true.
         
         if(string(1:2) .eq. '-o') read(string(3:),*)norbits
         if(string(1:3) .eq. '-kt')then
            read(string(4:),*)icolntype
         elseif(string(1:2) .eq. '-k') then
            read(string(3:),*)colnwt
         endif
         if(string(1:2) .eq. '-?') goto 51
         if(string(1:8) .eq. '--smooth') lsmooth=.true.
         if(string(1:5) .eq. '--at0') lat0=.true.
         if(string(1:10) .eq. '--localinj') localinj=.true.
         if(string(1:6) .eq. '--fext')then
c Read potential from a file
            call potread(string(7:),nrTi,rccmax)
c            call exit(0)
         endif
 1    continue
 3    continue
c PROVIS
      if (rmaxauto) then
         maxsteps=nint(4*rmax/(sqrt(2*Ti)*dtf))
      endif




c      write(*,*)linsulate,diags
c Set Array sizes, allowed variable.
      if(nr.gt.nrsize)then
         write(*,*)'Too many radial points:',nr,'  Set to',nrsize
         nr=nrsize
      endif
      if(npart.gt.npartmax)then
         write(*,*)'Too many ions:',npart,'  Set to',npartmax
         npart=npartmax
      endif
      if(LCIC)then
         if(nth.gt.nthsize-1)then
            write(*,*)'Too many theta points:',nth,'  Set to',nthsize-1
            nth=nthsize-1
         endif
         NRUSED=nr
         NTHUSED=nth
         NRFULL=nr
         NTHFULL=nth+1
      else
         if(nth.gt.nthsize)then
            write(*,*)'Too many theta points:',nth,'  Set to',nthsize
            nth=nthsize
         endif
         NRUSED=nr-1
         NTHUSED=nth-1
         NRFULL=nr
         NTHFULL=nth
      endif

c Inner domain parameter
      rsplit=nint(NRUSED/2.)

c Orbit tracking setting
      trackinit=maxsteps/2
      
      if(myid.eq.0)then
c         write(*,*)'Command line and defaults:'
         write(*,505)nr,nth,npart,Bz
         write(*,502)debyelen,Ti,maxsteps,dtf,rmax,vd,vprobe
 505     format(' Mesh: ',i3,' x',i3,'  Particles:',i7,' Bz:',f4.2)
 502     format('l_d=',f10.3,' Ti=',f7.3,' Steps=',i4,
     $        ' dt=',f6.4,' rmax=',f4.1,' vd=',f6.3,' vp=',f8.4)
      endif
      if(maxsteps.gt.nstepmax-1) Stop 'Too many steps requested.'
      lplot=diags.and.(myid.eq.0)
      k=0

c     Parallel solver stuff

c Don't need the parallel solver for 0 Debye length
      if ((debyelen.eq.0).or.infdbl) sorparallel=.false.

#ifdef MPI
      if (sorparallel) then
         call sorparinit(myid2,sor_comm)
      elseif (myid.eq.0) then
         myid2=0
      else
         myid2=-1
      endif
#else
      myid2=0
#endif

c Initialize the mesh and poisson coefficients
      if(LCIC)then
         call meshinitcic(rmax)
         call poisinitcic()
      else
         call meshinitngp(rmax)
         call poisinitngp()
      endif

c Checks when using an external file for potential
      if(lfext)then
         if(nrTi.ne.NRUSED)then
            write(*,*)'Incompatible r-array lengths',nrTi,NRUSED
            call exit()
         endif
         if(rccmax.ne.rcc(NRUSED))then
            write(*,*)'Incompatible r-domains',rccmax,rcc(NRUSED)
            call exit()
         endif
      endif

      if(nr.eq.10)then            
         write(*,*)'apc'
         write(*,521)apc
         write(*,*)'bpc'
         write(*,521)bpc
         write(*,*)'cpc'
         write(*,521)((cpc(i,j),i=1,10),j=0,NTHFULL)
         write(*,*)'dpc'
         write(*,521)((dpc(i,j),i=1,10),j=0,NTHFULL)
         write(*,*)'fpc'
         write(*,521)((fpc(i,j),i=1,10),j=0,NTHFULL)
 521     format(10g8.2)
      endif

c Initialize the fields. (Changed the order of finit and injinit)
      call finit()
c Initialize the random functions for reinjection.
      call injinit(icolntype,bcr)
c Initialize velocity diagnostics
      vrange=8.*sqrt(Ti)+1.4*abs(vd)+1.4*sqrt(abs(vprobe))
      do kk=1,nvmax
         nvdiag(kk)=0
         vrdiagin(kk)=0
         vtdiagin(kk)=0
         vdiag(kk)=vrange*(float(kk-1)/nvmax-0.499)
      enddo

      if(.not.lfixedn)then
         ninjcomp0=dtf*rhoinf*sqrt(Ti)*
     $        smaxflux(vd/sqrt(2.*Ti),0.)*r(NRFULL)**2
         ninjcomp=ninjcomp0
         if(myid.eq.0) write(*,*)'Initial:  rhoinf=',rhoinf
     $        ,'  ninjcomp=',ninjcomp
      endif
      iocprev=npart

c Initialize collisions and related fields.
      call colninit(colnwt,icolntype)
c Read in the previous particle distribution and averages.
      if(readpart) call partrd(success)
c Or Initialize (load) particles
      if(.not.success) call pinit()
c Assign charge to mesh
      call chargetomesh()
c Collect the partial sums of moments of distribution.
      call sumreduce()
      call rhocalc(lsmooth)
c Print out diagnostics for tiny meshes.
      if(nr.eq.10) then
         write(*,*)'rho,phi,psum'
         write(*,504)((rho(iw,jw),iw=1,nr),jw=1,NTHUSED)
         write(*,504)((phi(iw,jw),iw=1,nr),jw=0,NTHFULL)
         write(*,503)((psum(iw,jw),iw=1,nr),jw=0,NTHFULL)
      endif
      if(lplot)then
         if(ldist) then
            call multiframe(2,3,3)
         else
            call multiframe(2,2,3)
         endif
      endif

      time=0.
c Save the permanent plot switch.
      lpstore=lplot
      maccel=maxsteps/3
      nsamax=min(39,maxsteps/10)+1
      if(success)nstepsave=nsamax

c     Timing
#ifdef MPI
      sortime=MPI_WTIME()
#endif
      
      if (myid.eq.0) then
         write(*,*) "Maxsteps : ",maxsteps
      endif
c Main Stepping loop.
      do i=1,maxsteps

c Plot at some subset of steps.
c         write(*,*)mod(i-1,ipstep),i,ipstep
         lplot=lpstore.and.(mod(i-1,ipstep).eq.0)
c Acceleration for first 1/3rd of steps.
         bdtnow=max(1.,bdt*(maccel-i)/(maccel-1.))
         dt=bdtnow*dtf
         ninjcomp=bdtnow*ninjcomp0
c     Here there is a problem if we have ipstep!=1 because it is
c     essential that right after changing pfset one should call
c     pltinit. Otherwise initialization may be incorrect. This problem
c     is avoided by ensuring that the lplot is not changed back to true,
c     on leaving the loop, until we are about to reinitialize plotting.
c     As a result the last plot is only saved if it is on the maxsteps step.

c     If this flag is on, then we force the tracked particles to have the
c     following properties.

         if (myid.eq.0.and.i.eq.trackinit.and.orbinit) then 
            h=1./2*rmax
            do k=1,norbits
               iorbitlen(k)=0
               xp(3,k)=-2./3*rmax
c               xp(1,k)=h*k/norbits*cos(2*pi*k/norbits)
c               xp(2,k)=h*k/norbits*sin(2*pi*k/norbits)
               xp(1,k)=h*k/norbits
               xp(2,k)=0
               xp(6,k)=vd*1.
               xp(4,k)=0.
               xp(5,k)=0.
            enddo
         endif
               



         if(i.eq.maxsteps)call pfset(ipfsw)
c Assign charge to mesh

         call chargetomesh()
c Collect the partial sums of moments of distribution.
         call sumreduce()
         call aveupstep()

c myid.eq.0 for the master node.
         if(myid.eq.0)then
            
c     Refresh precalculated fuctions
            
c     if ((mod(i,20).eq.0)) then
c               call injinit(icolntype,bcr)
c     endif
            
c     Plot density contours.
            if(lplot)then
               if(lcolcont)then
                  call rhodisplay(10,rhomin,rhomax)
                  if(norbits.gt.0) call plotorbits()
               else
                  call pltinit(0.,1.,0.,1.)
               endif
            endif
            call rhocalc(lsmooth)
         endif

#ifdef MPI
c     Broadcast back to the slaves.
         if (sorparallel) then
            call MPI_BCAST(rho,(nrsize+1)*(nthsize+1),
     $           MPI_REAL,0,MPI_COMM_WORLD,ierr)
         endif
#endif

         if(myid.eq.0) then          
c Plot the slices through the probe.


            if(lplot) call slices(jstepth,rhomin,rhomax,time)

c Document charge accumulated; calculate diagrho.
c This is where rhoinf and finthave get calculated.
            call chargediag(dt,i)
         endif

c If we use floating potential, the slaves may need rhoinf and fluxes
#ifdef MPI
         if (sorparallel.and.(linsulate.or.lfloat)) then
            call MPI_BCAST(rhoinf,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(finthave,nthused,MPI_REAL,0,
     $           MPI_COMM_WORLD,ierr)
         endif
#endif
         if(myid.eq.0) then
c  Write step information.
            if(mod(i,5).eq.1) then
               if(i.gt.1)then
                  mtrap=mtrapped()
                  if(linsulate.or.lfloat)then
                     yvpre=vprobe
                  else
                     if(lfixedn)then
c                        yvpre=fluxrein/nrein
                        yvpre=float(nrein)/npart
                     else
                        yvpre=float(nrein)/ninjcomp
                     endif
                  endif
                  write(*,'(2f7.3,i7,f7.3)')
     $              fluxprobe(i-1)/(4.*pi*r(1)**2)/rhoinf/dtf,
     $              (phi(NRFULL,NTHUSED/2)+phi(NRUSED,NTHUSED/2))*.5,
     $                 mtrap,yvpre
               endif
               if(.not.lfixedn) write(*,601)ninjcomp,nrein,rhoinf,
     $              bdtnow,iocprev
 601           format('ninjcomp=',i6,'  nrein=',i6,'  rhoinf=',f8.2,
     $              '  bdtnow',f7.4
     $              ,'  iocprev=',i6)
               if(debyelen.ne.0. .and. mod(i,100).eq.1)then
                  write(*,*)
     $              'Step:SOR Step:SOR Step:SOR',
     $              ' Step:SOR Step:SOR',
     $              '  Flux, Phiedge, Trap, Vp or Rein'
               endif
            endif
            write(*,'(i5,$)')i
         endif

c Avoid uggly printing on the screen :
#ifdef MPI
         if (NRUSED.lt.150) then
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         endif
#endif


c     Calculate the potential field, unless it is fixed.
         if(.not.lfext)then
            if (myid2.eq.0.and.infdbl) then
               call fcalc_infdbl(dt)
            elseif (myid2.eq.0.and.debyelen.eq.0) then
               call fcalc(dt)
            elseif (myid2.ge.0)then
               if(lfcalcl.and.(.not.sorparallel)) then
                  call fcalc_lambda(dt)
               else
                  call fcalc_shielding(dt,sor_comm,myid2)
               endif
            endif
         endif


c     Broadcast back to the slaves.
#ifdef MPI
         call MPI_BCAST(phi,(nrsize+1)*(nthsize+1),
     $        MPI_REAL,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(adeficit,1,
     $        MPI_REAL,0,MPI_COMM_WORLD,ierr)
c     Averein is always broadcasted from 0 since it is calculated
c     in diags, called only by myid=0, not myid2=0
         call MPI_BCAST(averein,1,
     $        MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif

         if (myid.eq.0) then
c     Damping on ion motion by artificial potential modification.
            if(.not.lfext.and.damplen.gt.0.) call damp(dt,damplen)
            
c     Evaluate total charge and z-forces. Not yet implemented for NGP.
            if(LCIC)then
               call esforce(ierad,qp,fz,epz)
               zmom(i,enccharge,1)=qp
               zmom(i,fieldz,1)=fz
               zmom(i,epressz,1)=epz
               call esforce(nrused,qp,fz,epz)
               zmom(i,enccharge,2)=qp
               zmom(i,fieldz,2)=fz
               zmom(i,epressz,2)=epz
            endif
            
         endif

         
c     Plot radial average velocity and temperatures.
         if(lplot)then
            call vrdiag()
            if(ldist)then
               call vdiagsout(.false.)
               call vdiagsin(.false.)
            endif
         endif
c     Main particle advance.
         call padvnc(dt,icolntype,i)
c         call padvnc(dt,icolntype)
c     Collide the particles. See routine in advancing.f.
         if(icolntype.ne.0) call collide(dt,colnwt,icolntype)
c     Reduce back the flux and distribution data from the particle advance.
         call partreduce(i)

c     Adjust to the flux that would have occurred for standard step size.
         fluxprobe(i)=fluxprobe(i)/bdtnow
c     
         zmom(i,partz,1)=zmom(i,partz,1)/dt
         zmom(i,partz,2)=zmom(i,partz,2)/dt
         
c Small mesh written diagnostics.   
    
         if(nr.eq.10) then
            write(*,*)'Rho,phi'
            write(*,504)((rho(iw,jw),iw=1,nr),jw=0,NTHFULL)
            write(*,504)((phi(iw,jw),iw=1,nr),jw=0,NTHFULL)
         endif
         time=time+dt
 503     format(10f8.1)
 504     format(10f8.3) 
      enddo
c End of particle stepping section.

      itotsteps=maxsteps
      if(myid.eq.0)then
c         write(*,'(/,a)')'Exhausted maxsteps'
         write(*,*)
         write(*,505)nr,nth,npart,Bz
         write(*,502)debyelen,Ti,maxsteps,dtf,rmax,vd,vprobe
         if(iocprev.ne.npart)write(*,*)'iocprev=',iocprev
         if(rmtoz.ne.1.) write(*,'(''rmtoz='',f10.4)')rmtoz
         if(colnwt.ne.0) write(*,
     $        '(''Collisions: type='',i4,'' density='',f10.5)')
     $        icolntype,colnwt
c Average various fluxes and fields:
         call avefluxes(itotsteps,dt,fave,
     $        zmomave,fezave,zmoutave,qprobeave,epzave)
         if(lplot.and.finaldiags) then
            write(string,'(''Probe Flux Density='',f8.4)')fave
            call jdrwstr(0.05,0.7,string,1.)
         endif
      endif 
      if(lplot) call pltend()
c Restore the permanent plotting switch.
      lplot=lpstore

      if(writepart) call partwrt()

      diags=.true.
      call chargetomesh()
      call sumreduce()
c Draw some final diagnostic figures, if flag set.    
      if(finaldiags.and.myid.eq.0)then
         call multiframe(0,0,0)
         call yautoplot(fluxprobe,itotsteps)
         call axlabels('step','Particles to probe in step')
         call pltend()
         call multiframe(1,1,0)
         call rhodisplay(15,rhomin,rhomax)
         if(norbits.gt.0) call plotorbits()
         call pltend()
         if(ldist)then
            call charsize(0.025,0.025)
            call vdiagsout(.true.)
            call pltend()
            call vdiagsin(.true.)
            call pltend()
            call charsize(0.0,0.0)
         endif
      endif
      if(myid.eq.0)then
c Write the output file.
         call output(dt,damplen,i-1,fave)
         if (norbits.ge.1) call orbitoutput()
      endif

#ifdef MPI
      if (myid2.eq.0) then
         write(*,*) "time : ",MPI_WTIME()-sortime
      endif

      call MPI_FINALIZE(rc)
#endif
c Use call exit to avoid silly fortran message.
      call exit(0)

 51   continue
c Help section
      write(*,*)
     $ 'Usage: sceptic [-t.. -s.. -d.. -c..',
     $     ' -x.. -v.. -p.. -l.. -e[..] -r -w -g -?]'
      write(*,*)' e.g. ./sceptic -v0.2 -nr50 -nt30 -d.05'
      write(*,*)' Switch arguments (defaults):[.ff=float, nnn=integer]'
      write(*,*)' -t.ff ion temp Ti (1);',
     $     ' -snnn max time steps (500);',
     $     ' -sa Autodetermination of max steps'
      write(*,*)' -d.ff time step dt (.025);',
     $     ' -c.ff initial accel factor (1)'
      write(*,*)' -v.ff drift velocity (0.0);',
     $     ' -x.ff max r of domain (5)'
      write(*,*)' -p.ff probe potential (-4.) -pf float -pi insulate;',
     $     ' -b.ff Boundary factor (4)'
      write(*,*)' -l.ff Debye length (.1); -onnn Orbits to plot (0)',
     $     '; -z.ff External Ez (0)'
      write(*,*)' -e[ir,it] distribution function diagnostics',
     $     ' [in cell ir,it](1,1);'
      write(*,*)' -r read particle data (no), -w write ptcl. data (no)'
      write(*,*)' -g<nnn> diag plots only on nnn th step;',
     $     ' -a<n> save plots to disk (pfset n)'
      write(*,*)' -g no diags, -f no final diags -? Print this help.'
      write(*,*)' -nrnnn, -ntnnn, -ninnn, set radial, angle mesh-size',
     $     ', particle number.'
      write(*,*)' -ernnn radius at which to calculate q,E-force.',
     $     ' -mfff ratio of mass to Z (1.)'
      write(*,*)' -ktnnn collision type (0: none, 1 nu, 2 mfp const).'
     $     ,' -k.ff collision density.'
      write(*,*)' -vn.ff set neutral drift velocity for collisions.'
      write(*,*)' --at0 set theta acceleration artificially to zero'
      write(*,*)' --localinj use local potential for injection.',
     $     ' --smooth symmetrize in angle.'
      write(*,*)' --fextccccc.. use symmetric potential read from file.'
      write(*,*)' --fcalc use old potential solver (for validation etc)'
      write(*,*)' --damp.ff specify damping length (typically ~rmax).'
      write(*,*)' --rhoinf.ff specify constant injection rate, rather',
     $     ' than particle number.'
      write(*,*)' --bz.ff Bz (0),--ds subdivide the domain, --sp use',
     $     'parallel bloc solver'
      write(*,*)' --bcphi(0) bc potential (0:spherical sym, 1: Quasi',
     $     ' neutrality on outer crone, 2: Phiout=0, 3: dPhiout/dr=0)'
      write(*,*)' --bcr(0) bc reinject (0:spherical sym, 1: Simple',
     $     'Maxwellian, 2: Adiabatic)'

      end

c***********************************************************************
c End of Main Program.
c***********************************************************************

c Calculate rho from the psum.
      subroutine rhocalc(lsmooth)

c If lsmooth, then symmetrize the density, by averaging.
      logical lsmooth

c Common data:
      include 'piccom.f'
 501  format(10f8.1)
c Now we have added up fractional charges assigned to each mesh point. 
c We need to divide by the volume corresponding to each.
      if(.not.rhoinf.gt.1.e-4) then
c Initialize rhoinf approximately:
         write(*,*)'Rhoinf in rhocalc too small',rhoinf
         rhoinf=numprocs*npart/(4.*pi*r(NRFULL)**3/3.)
         write(*,*)'Rhoinf in rhocalc approximated as',rhoinf
      endif

c Correction for --ds
      do ir=1,NRUSED

c     Volumes are now calculated in meshinit.
         do ith=1,NTHUSED
            if(psum(ir,ith).eq.0. .and. debyelen.eq.0.)then
c               write(*,*)'Rhocalc fixup: zero particles in cell',ir,ith
               rho(ir,ith)=0.01*volinv(ir)*(nth-1.)*np/rhoinf
            else
               rho(ir,ith)=
     $              psum(ir,ith)*volinv(ir)*(nth-1.)*np/rhoinf
            endif
         enddo
         
         if(LCIC)then
c Fix the theta boundaries. The extra small quantity prevents zeros, I hope.
            rho(ir,1)=2.*rho(ir,1)+1.e-4
            rho(ir,nth)=2.*rho(ir,nth)+1.e-4
         endif

      
c     Smoothing in angle
c     Total symmetry case:
         if(lsmooth) then
            rhotot=0.
            do ith=1,NTHUSED
               rhotot=rhotot+rho(ir,ith)
            enddo
            rhotot=rhotot/nth
            do ith=1,NTHUSED
               rho(ir,ith)=rhotot
            enddo
         endif
         
      enddo

c     End of smoothing
      if(.not.LCIC)then
c     Fix theta ends for plotting.
         rho(ir,0)=rho(ir,1)
         rho(ir,nth)=rho(ir,nth-1)
      endif
      if(LCIC.and.debyelen.eq.0)then
c     Section to fix probe boundary bias arising from use of half a cell
c     in a region where there is strong gradient.
c     Not correct for large lambda. Because the
c     density can be increasing, not decreasing. 30 Aug 2003.
         do ith=1,nth
c     rho(1,ith)=rho(1,ith)/(1.+(2./3.)*sqrt(r(2)-r(1)))
            rho(1,ith)=rho(1,ith)/(1.+0.75*sqrt(r(2)-r(1)))
         enddo
      endif
      
      end
 
c*********************************************************************

      subroutine sumreduce()
#ifdef MPI
c The particle number total
      include 'piccom.f'
      include 'mpif.h'
      real ptot(0:nrsize,0:nthsize),vrtot(0:nrsize,0:nthsize)
      real vttot(0:nrsize,0:nthsize)
      real vr2tot(0:nrsize,0:nthsize),v2tot(0:nrsize,0:nthsize)
      real vtp2tot(0:nrsize,0:nthsize)
      call MPI_REDUCE(psum,ptot,(nrsize+1)*(nthsize+1),
     $     MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(vrsum,vrtot,(nrsize+1)*(nthsize+1),
     $     MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(vr2sum,vr2tot,(nrsize+1)*(nthsize+1),
     $     MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      if(diags)then
         call MPI_REDUCE(vtsum,vttot,(nrsize+1)*(nthsize+1),
     $        MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(v2sum,v2tot,(nrsize+1)*(nthsize+1),
     $        MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(vtp2sum,vtp2tot,(nrsize+1)*(nthsize+1),
     $        MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
      endif
      if(myid.eq.0)then
         do i1=0,nr
            do i2=0,nth
               psum(i1,i2)=ptot(i1,i2)
               vrsum(i1,i2)=vrtot(i1,i2)
               vr2sum(i1,i2)=vr2tot(i1,i2)
               if(diags)then
                  vtsum(i1,i2)=vttot(i1,i2)
                  v2sum(i1,i2)=v2tot(i1,i2)
                  vtp2sum(i1,i2)=vtp2tot(i1,i2)
               endif
            enddo
         enddo
      endif
#endif
      end
c***********************************************************************
      subroutine partreduce(i)
      integer i
      include 'piccom.f'
#ifdef MPI
      include 'mpif.h'
      real nvdiagtot(nvmax)
         call MPI_REDUCE(nrein,nreintot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(nreintry,nreintrytot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(spotrein,spotreintot,1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(fluxrein,fluxreintot,1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(ninner,nintot,1,MPI_INTEGER,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(ninth,ninthstep(1,i),nthsize,MPI_INTEGER,
     $        MPI_SUM,0,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(zmomprobe,zmom(i,partz,1),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(zmout,zmom(i,partz,2),1,MPI_REAL,MPI_SUM,0,
     $        MPI_COMM_WORLD,ierr)
         if(diags)then
            call MPI_REDUCE(nvdiag,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  nvdiag(kk)=nvdiagtot(kk)
               enddo
            endif
         endif
         if(ldist)then
            call MPI_REDUCE(vrdiagin,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  vrdiagin(kk)=nvdiagtot(kk)
               enddo
            endif
            call MPI_REDUCE(vtdiagin,nvdiagtot,nvmax,MPI_REAL,
     $           MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
               do kk=1,nvmax
                  vtdiagin(kk)=nvdiagtot(kk)
               enddo
            endif
         endif
#else
c This shuffle is necessary to accommodate the reduce.
         nreintot=nrein
         nreintrytot=nreintry
         spotreintot=spotrein
         fluxreintot=fluxrein
         nintot=ninner
         do j=1,nth
            ninthstep(j,i)=ninth(j)
         enddo
         zmom(i,partz,1)=zmomprobe
         zmom(i,partz,2)=zmout
#endif
c         write(*,*)'zmomstep=',zmomstep(i),'  zoutstep=',zmoutstep(i)
         nrein=nreintot
         nreintry=nreintrytot
         spotrein=spotreintot
         fluxrein=fluxreintot
         fluxprobe(i)=nintot
         end
c***************************************************************
      subroutine aveupstep()
      include 'piccom.f'
c Increment the averaging processes only here.
      if(nstepsave.lt.nsamax)nstepsave=nstepsave+1
      do j=1,nvmax
         nvdiag(j)=nvdiag(j)*(nstepsave-1.)/nstepsave
      enddo
      do j=1,nvmax
         vtdiagin(j)=vtdiagin(j)*(nstepsave-1.)/nstepsave
         vrdiagin(j)=vrdiagin(j)*(nstepsave-1.)/nstepsave
      enddo
      end
c*****************************************************************
      subroutine avefluxes(itotsteps,dt,fave,
     $     zmomave,fezave,zmoutave,qprobeave,epzave)
      include 'piccom.f'
      real phiout
c     Average the flux to the probe from the second half of the run.
c Also average the z-momentum components for the second half and 
c place them in zmom(nstepmax,*,*). maxsteps is limited to be nstepmax-1.
c The average is what is output. 
      fave=0.
      zmomave=0.
      zmoutave=0.
      zmoutvar=0.
      fezave=0.
      qprobeave=0.
      epzave=0.
      do j=1,4
         zmom(nstepmax,j,1)=0.
         zmom(nstepmax,j,2)=0.
      enddo
      if(itotsteps.ge.2)then
         m2=itotsteps/2
         do k=m2,itotsteps
            fave=fave+fluxprobe(k)
            do j=1,4
               zmom(nstepmax,j,1)=zmom(nstepmax,j,1)+zmom(k,j,1)
               zmom(nstepmax,j,2)=zmom(nstepmax,j,2)+zmom(k,j,2)
            enddo
            zmomave=zmomave+zmom(k,partz,1)
            fezave=fezave+zmom(k,fieldz,1)
            zmoutave=zmoutave+zmom(k,partz,2)
            zmoutvar=zmoutvar+zmom(k,partz,2)**2
            qprobeave=qprobeave+zmom(k,enccharge,1)
            epzave=epzave+zmom(k,epressz,2)
         enddo
         fave=fave/(itotsteps-m2+1)
c     normalize to rhoinf and the probe area to give normalized flux density
         fave=fave/(4.*pi*r(1)**2)/rhoinf/dt
         do j=1,4
            zmom(nstepmax,j,1)=zmom(nstepmax,j,1)/(itotsteps-m2+1)
            zmom(nstepmax,j,2)=zmom(nstepmax,j,2)/(itotsteps-m2+1)
         enddo
         zmom(nstepmax,partz,1)=zmom(nstepmax,partz,1)/rhoinf
         zmom(nstepmax,partz,2)=zmom(nstepmax,partz,2)/rhoinf
         zmomave=zmomave/(itotsteps-m2+1)/(rhoinf)
         fezave=fezave/(itotsteps-m2+1)
         qprobeave=qprobeave/(itotsteps-m2+1)
         zmoutave=zmoutave/(itotsteps-m2+1)/(rhoinf)
         zmoutvar=zmoutvar/(itotsteps-m2+1)/(rhoinf**2)-zmoutave**2
         epzave=epzave/(itotsteps-m2+1)
         write(*,501)fave
 501     format('Probe flux density=',f8.4)
         total1=zmom(nstepmax,fieldz,1)*debyelen**2
     $        +zmom(nstepmax,epressz,1)+zmom(nstepmax,partz,1)
         total2=zmom(nstepmax,fieldz,2)*debyelen**2
     $        +zmom(nstepmax,epressz,2)+zmom(nstepmax,partz,2)
         write(*,*)'Charge      E-field       Electrons',
     $        '      Ions     Total'
         write(*,*)(zmom(nstepmax,j,1),j=1,4),total1
         write(*,*)(zmom(nstepmax,j,2),j=1,4),total2

         phiout=0

         do j=1,NTHUSED
            phiout=phi(NRUSED,j)+phiout
         enddo
         phiout=phiout/NTHUSED
 
         write(*,*) 'rhoinf', rhoinf
         write(*,*) 'Outer potential',phiout


c tests
c         write(*,*)' Steps zmominner,outer'
c         do k=1,itotsteps
c            write(*,*)k,zmom(k,partz,1)/rhoinf,zmom(k,partz,2)/rhoinf
c         enddo
         write(*,'(''Ion Forces: Inner, outer, uncertainty,'',2f10.4,
     $        ''+-'',f7.4)')
     $        zmomave,zmoutave,sqrt(zmoutvar/(itotsteps-m2+1))
      endif
      end
c*******************************************************************
      subroutine potread(string,nrTi,rccmax)
      character*(*) string
      character*100 charin
      include 'piccom.f'
c Try to read file:
      write(*,*)'Trying to read file',string
      open(10,file=string,status='old',err=101)
c     Line for nothing.
      read(10,*)charin
      read(10,'(2f8.5,f8.4,i6,f12.4,f12.6,f7.4,2f14.5)',err=201)
     $     dtf,vd,Ti,istepsf,rhoinff,phiinff,
     $     favef,debyelen,vprobe
 201  continue
      write(*,*)'dtf,vd,Ti,istepsf,rhoinff,phiinff,favef,debyelen,Vp'
      write(*,*)dtf,vd,Ti,istepsf,rhoinff,phiinff,favef,debyelen,vprobe
      read(10,*,err=101)nrTi
      write(*,*)nrTi
      if(nrTi.gt.nrsize)then
         write(*,*)'r-array too long for storage'
         call exit()
      endif
c     write(*,*)'nrTi=',nrTi
      do ir=1,nrTi
         read(10,*,err=101)rcc(ir),diagphi(ir)
c         write(*,*)rcc(ir),diagphi(ir)
      enddo
      rccmax=rcc(nrTi)
      lfext=.true.
      close(10)
      return
 101  write(*,*)'Failed to find external file.'
      call exit()
      end

c***************************************************

      subroutine sorparinit(myid2,sor_comm)

c     Creates a new communicator, sor_comm, who contains a subset of
c     MPI_COMM_WORLD, because we don't need all the processes for the
c     bloc SOR

c     Not to be confused with mpicommcart, which is a reorganisation of
c     sor_comm

c     myid2 is the process id in the new communicator. If myid2 is negative,
c     the process does not belong to the communicator

      integer myid2,nprocsor
      integer members(0:100)
      integer group_world,sor_group,sor_comm

      include 'piccom.f'
      include 'mpif.h'

      idim1=nint(sqrt(numprocs*NRUSED/NTHUSED*1.))
      idim2=nint(numprocs/idim1*1.)

      nprocsor=idim1*idim2

      do i=0,nprocsor-1
         members(i)=i
      enddo

      call MPI_COMM_GROUP(MPI_COMM_WORLD,group_world,ierr)
      call MPI_GROUP_INCL(group_world,nprocsor,members,sor_group,ierr)
      call MPI_GROUP_RANK(sor_group,myid2,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,sor_group,sor_comm,ierr)
      if (myid2.ge.0) then
         call MPI_COMM_RANK(sor_comm,myid2,ierr)
      endif

      if (myid2.eq.0) then
c         write(*,*) "sorp : ",idim1,idim2,nprocsor,numprocs
         write(*,510) idim1,idim2,nprocsor
 510     format('BlocX=',i2,'  BlocY=',i2,'  SorProc=',i2)
      endif

      end
