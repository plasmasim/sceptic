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
c Version 2.6 outputs two files. T...frc traces the force evolution.


c     Writes a txt file with the orbits of the traced particules
      subroutine orbitoutput()

c     Common data
      include 'piccom.f'


      character*30 filename
      integer iti,it2

c Construct a filename that contains many parameters

      write(filename,'(a)')'T'
      iti=nint(alog10(Ti)-0.49)
      it2=nint(Ti/10.**iti)
      write(filename(2:2),'(i1.1)')it2
      if(iti.lt.0) then
         filename(3:3)='m'
         iti=-iti
      else
         filename(3:3)='e'
      endif
      write(filename(4:4),'(i1.1)')iti

      filename(5:5)='v'
      write(filename(6:8),'(i3.3)')nint(100*vd)
      filename(9:9)='r'
      write(filename(10:11),'(i2.2)')ifix(r(nr))
      filename(12:12)='P'
      write(filename(13:14),'(i2.2)')ifix(abs(Vprobe))

      filename(15:15)='L'
      if(debyelen.gt.1.e-10)then
         iti=nint(alog10(debyelen)-0.49)
         it2=nint(debyelen/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(16:16),'(i1.1)')it2
      if(iti.lt.0) then
         filename(17:17)='m'
         iti=-iti
      else
         filename(17:17)='e'
      endif
      write(filename(18:18),'(i1.1)')iti

      filename(19:19)='B'
      if(Bz.gt.1.e-10)then
         iti=nint(alog10(Bz)-0.49)-1
         it2=nint(Bz/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(20:21),'(i2.2)')it2
      if(iti.lt.0) then
         filename(22:22)='m'
         iti=-iti
      else
         filename(22:22)='e'
      endif
      write(filename(23:23),'(i1.1)')iti

      filename(24:27)='.orb'



      open(15,file=filename)
      write(15,*) 'Number of orbits'
      write(15,*) norbits
      do k=1,norbits
         write(15,*) k,'th orbit'
         write(15,*) iorbitlen(k)
         do i=1,iorbitlen(k)
            write(15,590) xorbit(i,k),yorbit(i,k),zorbit(i,k),
     $      vxorbit(i,k),vyorbit(i,k),vzorbit(i,k)
         enddo
      enddo
 590  format(6f9.4)
      close(15)
      end


c     Writes the main output file
      subroutine output(dt,damplen,i,fave)
c Common data:
      include 'piccom.f'
      character*30 filename
      integer iti,it2
c Construct a filename that contains many parameters
      write(filename,'(a)')'T'
      iti=nint(alog10(Ti)-0.49)
      it2=nint(Ti/10.**iti)
      write(filename(2:2),'(i1.1)')it2
      if(iti.lt.0) then
         filename(3:3)='m'
         iti=-iti
      else
         filename(3:3)='e'
      endif
      write(filename(4:4),'(i1.1)')iti

      filename(5:5)='v'
      write(filename(6:8),'(i3.3)')nint(100*vd)
      filename(9:9)='r'
      write(filename(10:11),'(i2.2)')ifix(r(nr))
      filename(12:12)='P'
      write(filename(13:14),'(i2.2)')ifix(abs(Vprobe))

      filename(15:15)='L'
      if(debyelen.gt.1.e-10)then
         iti=nint(alog10(debyelen)-0.49)
         it2=nint(debyelen/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(16:16),'(i1.1)')it2
      if(iti.lt.0) then
         filename(17:17)='m'
         iti=-iti
      else
         filename(17:17)='e'
      endif
      write(filename(18:18),'(i1.1)')iti

      filename(19:19)='B'
      if(Bz.gt.1.e-10)then
         iti=nint(alog10(Bz)-0.49)-1
         it2=nint(Bz/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(20:21),'(i2.2)')it2
      if(iti.lt.0) then
         filename(22:22)='m'
         iti=-iti
      else
         filename(22:22)='e'
      endif
      write(filename(23:23),'(i1.1)')iti

      filename(24:27)='.dat'

c Write out averaged results.
      open(10,file=filename)
      write(10,'(a,a)')'  dt       vd       Ti  steps    rhoinf  
     $          phiinf   fave  debyelen Vp damplen / nr/ phi  Bz'
      write(10,'(2f8.5,f8.4,i6,f12.4,f11.5,f8.4,2f14.5,f8.3,f8.4)')
     $     dt,vd,Ti,i,rhoinf,log(rhoinf),fave,debyelen,vprobe,damplen,Bz
      write(10,*)NRUSED
      do j=1,NRUSED
         write(10,*)rcc(j),diagphi(j),diagrho(j)
c -log(rhoinf)
      enddo
      write(10,'(a)')'Number of steps, Particles to probe each step'
      write(10,*)i
      write(10,*)(fluxprobe(j),j=1,i)
      write(10,'(a)')'Number of theta cells, Number of steps, Particles'
      write(10,*)NTHUSED,i
      do j=1,NTHUSED
         ninth(j)=0
      enddo
      nastep=0
      do k=1,i
         write(10,*)(ninthstep(j,k),j=1,NTHUSED)
         if(k.gt.3*i/4)then
            nastep=nastep+1
            do j=1,NTHUSED
               ninth(j)=ninth(j)+ninthstep(j,k)
            enddo
         endif
      enddo
      write(10,'(a,a)')'Particle angular distrib summed over last'
     $     ,' half of steps, numbering:'
      write(10,*)nastep
      write(10,*)(ninth(j),j=1,NTHUSED)
      write(10,'(a,i4,i4)')'Mesh potential. Grid',NRUSED,NTHUSED
      do j=1,NRUSED
         write(10,'(10f8.3)')(phi(j,k),k=1,NTHUSED)
      enddo
      write(10,'(a,i4,i4)')'Mesh density/infinity. Grid',NRUSED,NTHUSED
      do j=1,NRUSED
         write(10,'(10f8.3)')(rho(j,k),k=1,NTHUSED)
      enddo
      write(10,'(a,i4,i4)')'Volinv. Grid',NRUSED
      write(10,'(10f8.3)')(volinv(k),k=1,NRUSED)

      call outsums(dt,i+1)

c Output time-averages of z-force components stored in zmom(nstepmax,*,*).
c Particle units nTr^2, Electric nT lambda_D^2. 
      total1=zmom(nstepmax,fieldz,1)*debyelen**2
     $     +zmom(nstepmax,epressz,1)+zmom(nstepmax,partz,1)
      total2=zmom(nstepmax,fieldz,2)*debyelen**2
     $     +zmom(nstepmax,epressz,2)+zmom(nstepmax,partz,2)
      write(10,*)'Charge      E-field       Electrons',
     $     '      Ions     Total'
      write(10,*)(zmom(nstepmax,j,1),j=1,4),total1
      write(10,*)(zmom(nstepmax,j,2),j=1,4),total2
      if(rmtoz.ne.1.) write(10,'(''rmtoz='',f10.4)')rmtoz
      if(icolntype.ne.0) write(10,
     $     '(''Collisions: type='',i4,'' weight='',f10.5)')
     $     icolntype,colnwt

c End of output file.
      close(10)
      filename(24:27)='.frc'
      call outforce(filename,i)
      end
c*********************************************************************
c Write a second file with the force data as a function of step
      subroutine outforce(filename,istepmax)
      character*(*) filename
      include 'piccom.f'
      real zmn(nstepmax,4,2)
c Apply normalization factors but don't change the zmom.
c Perhaps this extra storage is unnecessary.
      do i=1,istepmax
         do k=1,2
            zmn(i,enccharge,k)=zmom(i,enccharge,k)
            zmn(i,partz,k)=zmom(i,partz,k)/rhoinf
            zmn(i,epressz,k)=zmom(i,epressz,k)
            zmn(i,fieldz,k)=zmom(i,fieldz,k)*debyelen**2
         enddo
      enddo
      open(9,file=filename)
      write(9,*)istepmax
      write(9,*)'Step    Charge     E-field      Electrons',
     $        '      Ions   Total Force'
      write(9,'(i5,5f12.5)')((i,(zmn(i,j,k),j=1,4)
     $     ,zmn(i,partz,k)+zmn(i,fieldz,k)+zmn(i,epressz,k)
     $     ,k=1,2),i=1,istepmax)
      close(9)
      end
c**********************************************************************
c Write out the particle data.
      subroutine partwrt()
c Common data:
      include 'piccom.f'
      character*11 filename

      write(filename,'(''part'',i3.3,''.dat'')')myid
c Delete the file first to help with nfs problems.
      open(11,file=filename,status='unknown')
      close(11,status='delete')
c
      open(11,file=filename,status='unknown')
      write(11,*)npartmax,npart,nr,nth,ndim,np
      write(11,*)xp
      write(11,*)rhoinf,spotrein,averein
      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
      close(11)
      end

c**********************************************************************
c Read in the particle data.
      subroutine partrd(success)
      logical success
c Common data:
      include 'piccom.f'
      character*11 filename

      write(filename,'(''part'',i3.3,''.dat'')')myid
      success=.false.
      open(11,file=filename,status='old',err=101)
      read(11,*,err=100,end=100)ipartmax,ipart,ir,ith,idim,ip
      if(ipartmax.eq.npartmax .and. ipart.eq.npart
     $     )then
c     $ .and. ir.eq.nr .and. ith.eq.nth )then
         write(*,*)'Using saved particle data.'
         read(11,*,err=100,end=100)xp
         read(11,*,err=100,end=100)rhoinf,spotrein,averein
      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
         success=.true.
      else
         write(*,*)'Particle data mismatch',ipartmax,npartmax,ipart,
     $        npart,ir,nr,ith,nth
      endif
      close(11)
      return
 100  close(11)
      write(*,*) 'Error reading pardata.dat'
      return
 101  write(*,*) 'No particle file to read.'
      end
c**********************************************************************
c Get the average and slope over the rmesh range i1,i2.
      subroutine slopegen(phi,r,nr,i1,i2,slope,average)
      integer nr
      real phi(nr),r(nr)

c Assume r-mesh is linear
      rmom0=0.
      rmom1=0.
      rmid=(r(i2)+r(i1))/2.
      do i=i1,i2
         rmom0=rmom0+phi(i)
         rmom1=rmom1+(r(i)-rmid)*phi(i)
      enddo
      average=rmom0/(i2-i1+1)
c      rave=rmom1/(i2-i1+1)
      rlen=r(i2)-r(i1)
      slope=12.*(rmom1)/(rlen*rlen)/(i2-i1+1)
c      write(*,*)rmom0,rmom1,r(i1),r(i2),rmid,rlen
      end

c**********************************************************************
      subroutine outsums(dt,i)
c Common data:
      include 'piccom.f'
c Write out summed results.
      nrhere=NRUSED
      nthhere=NTHUSED
      nphere=np
c Combined files. Don't have to open.
c      open(10,file=filename)
      write(10,'(a,a)')
     $     '    dt       vd       Ti  steps  rmax',
     $     ' rhoinf debyelen Vp /nr,nth,np; sums'
      write(10,'(2f8.5,f8.4,i6,f8.3,f12.3,2f14.5)')
     $     dt,vd,Ti,i,r(nr),rhoinf,debyelen,vprobe
      write(10,*)nrhere,nthhere,nphere
      write(10,*)'psum'
      write(10,*)((psum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vrsum'
      write(10,*)((vrsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vtsum'
      write(10,*)((vtsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vpsum'
      write(10,*)((vpsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'v2sum'
      write(10,*)((v2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vr2sum'
      write(10,*)((vr2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vtp2sum'
      write(10,*)((vtp2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'diagvr'
      write(10,*)((diagvr(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'r[cc]'
      write(10,*)(rcc(k1),k1=1,nrhere)
      write(10,*)'volinv'
      write(10,*)(volinv(k1),k1=1,nrhere)
      write(10,*)'t[cc]'
      write(10,*)(tcc(k2),k2=1,nthhere)
      end




