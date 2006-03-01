c Tools for postprocessing a serie of log files with different Bz, but
c same other parameters.
c The launch syntax is postprocB [switches] T...B*.dat
c Must launch at least two files to work properly



c********************************************************************
      character*100 string,filename
      include 'piccompost.f'
      real rholocal(0:NRFULL,0:NTHFULL)
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2
      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      integer nti0
      parameter (nti0=100)
      character*100 charin
      logical lpcic,lreaddiag
      logical ledge
      logical output,sh

c     Data tables to be collected for the plotting
      real Blist(40),fluxlist(40),probeplist(40)
      integer narg
      real gamma(40)
      real temp1,temp2,randflux
      integer oaverage,temp
      real outdens(250,40),outphi(250,40),inphi1(250,40),thcos(250,40)
      real thflux(250,40),inphi2(250,40),indens1(250,40),indens2(250,40)
      real tempt1(250),tempt2(250),inphi3(250,40),indens3(250,40)


c     Tables for the sheath calculation. Doesn't make sense
c     for high B field
      real xsheath(250,40)
      real ysheath(250,40)
      real sheathrho(250,250)
      real sheathr(250)


      data lpcic/.false./
      data lreaddiag/.false./
      data ledge/.false./
      data output/.false./
      data oaverage/10/
      data sh/.false./
c     Temporary data, for the bubble sort or the plotting
      temp=1.
      temp1=1
      temp2=1

c Deal with arguments, and store all the values in the table for each B
      narg=iargc()
      if (iargc().eq.0) goto 51
      do i=1,narg
 11      call getarg(i,string)

         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:2).eq.'-f') then
            output=.true.
            goto 3
         elseif(string(1:3).eq.'-sh') then
            sh=.true.
            goto 3
         else
            filename=string
            do k=1,len(filename)-2
               if(filename(k:k+1).eq.'Sp')then
                  filename(k:k+1)='Ti'
                  write(*,*) 'Using Ti file ',filename
                  goto 110
               endif
            enddo
 110        continue
         endif
c     Read the outputfile
         call readoutput(lreaddiag,lpcic,ledge,
     $        filename,rholocal,nrhere,nthhere,nphere,
     $        phipic,rhopic,rpic,rpicleft,phicos,
     $        rhomax,rhomin,nrti,phiinf,nastep,nsteps,dt,rmax,fave,
     $        debyelen,vprobe,
     $        ierr)

c     Creating a list with all the Bz to analyse, and the corresponding
c     potential. Useful in the float potential case
         Blist(i)=Bz
         probeplist(i)=vprobe

c         open(15,file='flux.txt',access='append')
c         write(15,*) 'B: ',Bz,' Vp: ',vprobe
c         close(15)
c     Creating a list of the flux as a function of theta
         fluxlist(i)=0.
         write(*,*) "nastep",nastep
         do j=1,nthhere
            thflux(j,i)=ninth(j)*(nthhere-1.)/
     $           (4.*pi*rhoinf*dt*nastep)
            if(j.eq.1 .or. j.eq.nthhere) thflux(j,i)=thflux(j,i)*2.
            fluxlist(i)=fluxlist(i)+thflux(j,i)
         enddo
         randflux=2*sqrt(2*Ti)*sqrt(pi)*(1-0/Ti)
         fluxlist(i)=4*pi*fluxlist(i)/nthhere/randflux

c     Creating a list of the flux ratios for all the files
         gamma(i)=1.*ninth(1)/ninth(nthhere)

c     Creating a list of the outer density as a funtion of theta
         do k=1,nthhere
            outdens(k,i)=0.
            do l=0,oaverage-1
               outdens(k,i)=outdens(k,i)+rholocal(nrhere-l,k)
            enddo
            outdens(k,i)=outdens(k,i)/(oaverage)
            thcos(k,i)=tcc(k)
         enddo

c     Creating a list of the outer potential as a function of theta
         do k=1,nthhere
            outphi(k,i)=0.
            do l=0,oaverage-1
               outphi(k,i)=outphi(k,i)+phi(nrhere-l,k)
            enddo
            outphi(k,i)=outphi(k,i)/(oaverage)
         enddo

c     Creating a list of the potential and density as a function of r
         do k=1,nrhere
            inphi1(k,i)=phi(k,1)
            inphi2(k,i)=phi(k,int(nthhere/2))
            inphi3(k,i)=phi(k,nthhere)
            indens1(k,i)=rho(k,1)
            indens2(k,i)=rho(k,int(nthhere/2))
            indens3(k,i)=rho(k,nthhere)
         enddo


c     calculates the sheath size
         ptaverage=1.
         do j=1,nthhere
            do k=nrhere,1,-1
               ptaverage=1.
               sheathrho(k,j)=rholocal(k,j)
               if (j.ge.2) then
                  ptaverage=ptaverage+1
                  sheathrho(k,j)=sheathrho(k,j)+rholocal(k,j-1)
               endif
               if (j.le.(nthhere-1)) then
                  ptaverage=ptaverage+1
                  sheathrho(k,j)=sheathrho(k,j)+rholocal(k,j+1)
               endif
               if (k.ge.2) then
                  ptaverage=ptaverage+1
                  sheathrho(k,j)=sheathrho(k,j)+rholocal(k-1,j)
               endif
               if (k.le.(nrhere-1)) then
                  ptaverage=ptaverage+1
                  sheathrho(k,j)=sheathrho(k,j)+rholocal(k+1,j)
               endif
               sheathrho(k,j)=sheathrho(k,j)/ptaverage
               if ((-exp(phi(k,j))+sheathrho(k,j)).le.(0.1)) then
                  sheathr(j)=rcc(k)
               endif
            enddo
         enddo
         do j=1,nthhere
            xsheath(j,i)=sheathr(j)*tcc(j)
            ysheath(j,i)=sheathr(j)*sqrt(1-tcc(j)**2)
         enddo

3       continue
      enddo


c     We sort all the lists to have them in rising B order
c     We do a simple bubble sort

      temp=1
 900  if(Blist(temp).gt.Blist(temp+1)) then


         temp1=Blist(temp+1)
         Blist(temp+1)=Blist(temp)
         Blist(temp)=temp1

         temp1=probeplist(temp+1)
         probeplist(temp+1)=probeplist(temp)
         probeplist(temp)=temp1

         temp1=gamma(temp+1)
         gamma(temp+1)=gamma(temp)
         gamma(temp)=temp1

         temp1=fluxlist(temp+1)
         fluxlist(temp+1)=fluxlist(temp)
         fluxlist(temp)=temp1
         
         do k=1,nthhere
            temp1=outdens(k,temp+1)
            outdens(k,temp+1)=outdens(k,temp)
            outdens(k,temp)=temp1

            temp1=thcos(k,temp+1)
            thcos(k,temp+1)=thcos(k,temp)
            thcos(k,temp)=temp1

            temp1=outphi(k,temp+1)
            outphi(k,temp+1)=outphi(k,temp)
            outphi(k,temp)=temp1

            temp1=xsheath(k,temp+1)
            xsheath(k,temp+1)=xsheath(k,temp)
            xsheath(k,temp)=temp1

            temp1=ysheath(k,temp+1)
            ysheath(k,temp+1)=ysheath(k,temp)
            ysheath(k,temp)=temp1
            
            temp1=thflux(k,temp+1)
            thflux(k,temp+1)=thflux(k,temp)
            thflux(k,temp)=temp1
         enddo

         do k=1,nrhere
            temp1=inphi1(k,temp+1)
            inphi1(k,temp+1)=inphi1(k,temp)
            inphi1(k,temp)=temp1
            temp1=inphi2(k,temp+1)
            inphi2(k,temp+1)=inphi2(k,temp)
            inphi2(k,temp)=temp1
            temp1=inphi3(k,temp+1)
            inphi3(k,temp+1)=inphi3(k,temp)
            inphi3(k,temp)=temp1

            temp1=indens1(k,temp+1)
            indens1(k,temp+1)=indens1(k,temp)
            indens1(k,temp)=temp1
            temp1=indens2(k,temp+1)
            indens2(k,temp+1)=indens2(k,temp)
            indens2(k,temp)=temp1
            temp1=indens3(k,temp+1)
            indens3(k,temp+1)=indens3(k,temp)
            indens3(k,temp)=temp1
         enddo

         temp=0
      endif
      if (temp.eq.(narg-1)) goto 901
      temp=temp+1
      goto 900
 901  continue




c Writes on the output file if necessary
      if (output) then
         open(15,file='flux.txt',access='append')
         write(15,*) 'Ti:',Ti,'   Lambda:',debyelen,'   Vp:',vprobe
         write(15,*) 'B        flux       Vp'
         do k=1,narg
            write(15,'(3f8.4)') Blist(k),fluxlist(k),probeplist(k)
         enddo
         write(15,*) ""
         close(15)
      endif

c Start of Plotting:
c Set arrow scale
      v1=max(1.,vd)
c Change the B scale for the sake of plotting
      do k=1,narg
         Blist(k)=Blist(k)/(1+Blist(k))
      enddo

      call pfset(3)


      call multiframe(0,0,0)

c     Plot floating potential as a function of B
      call pltinit(0.,1.,-5.,0.)
      call axis()
      call axis2()
      call winset(.true.)
      call polyline(Blist,probeplist,narg)
      call winset(.false.)
      call fwrite(vd,iwdth,2,charin)
      call jdrwstr(0.05,.70,
     $     'v!dd!d='//charin(1:iwdth)//char(0),1.)
      call fwrite(Ti,iwdth,2,charin)
      call jdrwstr(0.05,0.60,
     $     'T!di!d='//charin(1:iwdth)//char(0),1.)
      call axlabels('B/(1+B)','Probe floating potential')
      call pltend()

c     Plot outer flux as a function of B
      call pltinit(0.,1.0,0.,-vprobe/Ti+1.)
      call axis()
      call axis2()
      call winset(.true.)
      call polyline(Blist,fluxlist,narg)
      call winset(.false.)
      call fwrite(vd,iwdth,2,charin)
      call jdrwstr(0.05,.70,
     $     'v!dd!d='//charin(1:iwdth)//char(0),1.)
      call fwrite(Ti,iwdth,2,charin)
      call jdrwstr(0.05,0.60,
     $     'T!di!d='//charin(1:iwdth)//char(0),1.)
      call axlabels('B/(1+B)','Flux to probe')
      call pltend()
      
c     Plot flux ratio as a function of B
      call autoplot(Blist,gamma,narg)
      call fwrite(vd,iwdth,2,charin)
      call jdrwstr(0.05,.70,
     $     'v!dd!d='//charin(1:iwdth)//char(0),1.)
      call fwrite(Ti,iwdth,2,charin)
      call jdrwstr(0.05,0.60,
     $     'T!di!d='//charin(1:iwdth)//char(0),1.)
      call axlabels('B','Flux ratio')
      call pltend()


c     Plot outer density as a function of cos(theta) for all B
      call axptset(0.,0.)

      call minmax(outdens(1,1),nthhere*narg,outmin,outmax)
      call pltinit(-1.0,1.0,outmin-0.1,outmax+0.1)
      call axis()
      call axis2()
      call axlabels('cos(!Aq!@)','rho boundary')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nthhere
            tempt1(k)=outdens(k,j)
            tempt2(k)=thcos(k,j)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nthhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo
      call pltend()

c     Plot flux as a funtion of theta
      call minmax(thflux(1,1),nthhere*narg,outmin,outmax)
      call pltinit(-1.0,1.0,outmin-0.1,outmax+0.1)
      call axis()
      call axis2()
      call axlabels('cos(!Aq!@)','flux to probe')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nthhere
            tempt1(k)=thflux(k,j)
            tempt2(k)=thcos(k,j)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nthhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
         if (output) then
            open(15,file='flux.txt',access='append')
            write(*,*) Blist(j)
            write(15,*) tempt1
            write(15,*) ""
            close(15)
         endif
      enddo
      
      call pltend()

c     Plot outer potential as a function of cos(theta)
      call minmax(outphi(1,1),nthhere*narg,outmin,outmax)
      call pltinit(-1.0,1.0,outmin-0.1,outmax+0.1)
      call axis()
      call axis2()
      call axlabels('cos(!Aq!@)','phi boundary')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nthhere
            tempt1(k)=outphi(k,j)
            tempt2(k)=thcos(k,j)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nthhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo
      call pltend()
      
      call multiframe(2,1,3)

c     Plot potential at theta=0 and Pi/2 as function of r
      call minmax(inphi1(1,1),nrhere*narg,outmin,outmax)
      call pltinit(0.,rcc(nrhere)+0.,outmin-0.1,0.2)
      call axis()
      call axis2()
      call axlabels('r','phi !Aq!@=Pi/2')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nrhere
            tempt1(k)=inphi1(k,j)
            tempt2(k)=rcc(k)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nrhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo

      call color(15)
      call pltinit(0.,rcc(nrhere)+0.,outmin,outmax)
      call axis()
      call axis2()
      call axlabels('r','phi !Aq!@=0')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nrhere
            tempt1(k)=inphi2(k,j)
            tempt2(k)=rcc(k)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nrhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo
      call pltend()

c     Plot density at theta=0 and Pi/2 as function of r
      call minmax(indens1(1,1),nrhere*narg,outmin,outmax)
      call pltinit(0.,rcc(nrhere)+0.,outmin-0.1,outmax+0.2)
      call axis()
      call axis2()
      call axlabels('r','rho !Aq!@=Pi/2')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nrhere
            tempt1(k)=indens1(k,j)
            tempt2(k)=rcc(k)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nrhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo
      call minmax(indens2(1,1),nrhere*narg,outmin,outmax)
      call color(15)
      call pltinit(0.,rcc(nrhere)+0.,outmin-0.1,outmax+0.2)
      call axis()
      call axis2()
      call axlabels('r','rho !Aq!@=0')
      do j=1,narg
         call color(mod(j,15)+1)
         do k=1,nrhere
            tempt1(k)=indens2(k,j)
            tempt2(k)=rcc(k)
         enddo
         call winset(.true.)
         call polyline(tempt2,tempt1,nrhere)
         call winset(.false.)
         write(charin,'(f4.2)') Blist(j)
         call legendline(-.48,0.05*(j-1),0,
     $        charin)
      enddo
      call pltend()
      call multiframe(0,0,0)

c     Plot sheaths
      if (sh) then
         call minmax(xsheath(1,1),nthhere*narg,outmin,outmax)
         call pltinit(-outmax*1.1,outmax*1.1,0.,outmax*2.2)
         call axis()
         do j=1,narg
            call color(mod(j,10)+1)
            do k=1,nthhere
               tempt1(k)=xsheath(k,j)
               tempt2(k)=ysheath(k,j)
            enddo
            call winset(.true.)
            call polyline(tempt1,tempt2,nthhere)
            call winset(.false.)
            write(charin,'(f4.2)') Blist(j)
            call legendline(-.48,0.05*(j-1),0,
     $           charin)
            
         enddo
         call pltend()
      endif

      goto 52
 51   continue
c Help section
      write(*,*)
     $     'Usage: postprocB [-sh..-f..] filenameB*.dat'
      write(*,*) 'e.g. ./postprocB -f T1e0v000r05P04L1m1B00e0.dat'
      write(*,*) 'Switch arguments (defaults)'
      write(*,*) '-sh(false) Plot sheaths (Quasineutrality break down)'
      write(*,*) 'f(false) Writes output file'

 52   end

c Data reading subroutine
      subroutine readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,nphere,
     $     phipic,rhopic,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,debyelen,vprobe,
     $     ierr)
      logical lreaddiag,lpcic,ledge
      character*100 string,filename
      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      include 'piccompost.f'
      real rholocal(0:NRFULL,0:NTHFULL)
      character*100 charin
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2


      ierr=0
c Read the data file:
c__________________________________________________________________
      open(10,file=filename,status='old',err=101)
c Line for nothing.
      read(10,*)charin

      read(10,'(2f8.5, f8.4, i6, f12.4, f11.5, f8.4, 2f14.5, f8.3,
     $     f8.4)',err=201) dt,vd,Ti,isteps,rhoinf,phiinf,fave,
     $     debyelen,vprobe,damplen,Bz
 201  continue
      write(*,*)'Bz,dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,Vp'
      write(*,*)Bz,dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
      read(10,*,err=202)nrTi
      nrhere=nrTi
c      write(*,*)'nrTi=',nrTi
      do i=1,nrTi
         read(10,*,err=203)rpic(i),phipic(i)
      enddo
      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*,err=204)nsteps
c      write(*,*)nsteps
      if(nsteps.gt.nstepmax) then
         write(*,*)'Number of steps',nsteps,
     $        ' exceeds allocation',nstepmax
         call exit
      endif
      read(10,*)(fluxprobe(j),j=1,nsteps)
c Read theta cells
      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*)nthhere,nsteps
c      write(*,*)nthhere,nsteps
      do i=1,nsteps
         read(10,*)(ninthstep(j,i),j=1,nthhere)
      enddo

      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*)nastep
c      write(*,*)'nastep',nastep
      nread=nthhere
c This is not necessary here and indeed breaks the combined read because
c lpcic has not yet been set.
c      if(.not.lpcic)nread=nread+1
      read(10,*)(ninth(j),j=1,nread)
      if(lreaddiag)then
         write(*,*)'nastep',nastep,' ninth:'
         write(*,*)( ninth(j),j=1,nthhere)
      endif

      read(10,*)charin
      do j=1,nrhere
         read(10,*)(phi(j,k),k=1,nthhere)
      enddo
      read(10,*)charin
      do j=1,nrhere
         read(10,'(10f8.3)')(rholocal(j,k),k=1,nthhere)
      enddo
      read(10,*)charin
      read(10,*)(volinv(k),k=1,nrhere)
c We don't use open and close for combined files.
      if(filename(1:2).eq.'Ti')then
c But this must be using the split version.
         close(10)
c Read in  summed results.
         filename(1:2)='Sp'
         open(10,file=filename,err=210,status='old')
      endif
      read(10,'(a)')string
      read(10,'(2f8.5,f8.4,i6,f8.3,f12.3,2f14.5)',err=200)
     $     dt,vd,Ti,i,rmax,rhoinf,debyelen,vprobe
 200  continue
      read(10,*)nrhere,nthhere,nphere
      if(nrhere.gt.NRUSED .or. nthhere.gt.NTHUSED)then
         write(*,*)'Required dimensions: nr',nrhere,' nth',nthhere
         write(*,*)'are too large for the allocated values:'
     $        ,NRUSED,NTHUSED
         stop
      endif
      read(10,*)string
      read(10,*)((psum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vrsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vtsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vpsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((v2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vr2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vtp2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((diagvr(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)(rcc(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(volinv(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(tcc(k2),k2=1,nthhere)      
      read(10,*,err=402,end=402)string
      read(10,*)charge1,ffield1,felec1,fion1,ftot1
      read(10,*)charge2,ffield2,felec2,fion2,ftot2
 402  close(10)
      write(*,*)'nrhere,nthhere'
      write(*,*)nrhere,nthhere
      if(lreaddiag)then
         write(*,*)'Finished reading'
         write(*,*)'vrsum(1)'
         write(*,501)(vrsum(1,k2), k2=1,nthhere)
         write(*,*)'psum(1)'
         write(*,501)(psum(1,k2), k2=1,nthhere)
         write(*,*)'vr(1)'
         write(*,501)(vrsum(1,k2)/psum(1,k2), k2=1,nthhere)
         write(*,*)'diagvr(1)'
         write(*,501)(diagvr(1,k2), k2=1,nthhere)
         write(*,*)'rcc'
         write(*,*)(rcc(k1),k1=1,nrhere)
         write(*,*)'volinv'
         write(*,*)(volinv(k1),k1=1,nrhere)
         write(*,*)'tcc'
         write(*,*)(tcc(k2),k2=1,nthhere)
      endif
 501  format(10f8.3)
c__________________________________________________________________
c End of reading the data section

c__________________________________________________________________
c Fix up data 
      if(tcc(1).eq.1)lpcic=.true.
c Correct the outside angle centers if necessary.
      if(lpcic)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.

         thcells=nthhere-1.
      else
         thcells=nthhere
      endif


c In postproc, th is used as the weighting of the cells, essentially 
c the delta of cos theta that corresponds to each cell.
      do i=1,nthhere
         th(i)=2./(thcells)
      enddo
c The ends are half-weighted for cic.
      if(lpcic)then
         th(1)=0.5*th(1)
         th(nthhere)=0.5*th(nthhere)
      endif

c      write(*,*)'      tcc      th      thang','  lpcic=',lpcic
c      write(*, '(3f10.4)')(tcc(kk),th(kk),thang(kk),kk=0,nthhere+1)

c      call autoplot(th,diagvr,nthhere)

c Calculate rho(i,j) from the psum, volinv. Normalized by rhoinf.
      rhomax=0.
      do k1=1,nrhere
         do k2=1,nthhere
            rho(k1,k2)=psum(k1,k2)*volinv(k1)*thcells*nphere/rhoinf
            if(rho(k1,k2).gt.rhomax)rhomax=rho(k1,k2)
         enddo
         if(lpcic)then
c fix up rho as double on boundary.
            rho(k1,1)=2.*rho(k1,1)
            rho(k1,nthhere)=2.*rho(k1,nthhere)
         endif
c fix angle ends of rho and phi
         rho(k1,0)=rho(k1,1)
         rho(k1,nthhere+1)=rho(k1,nthhere)
         phi(k1,0)=phi(k1,1)
         phi(k1,nthhere+1)=phi(k1,nthhere)
      enddo
      ir=10.
      if(ledge)then
         rhomax=min(rhomax,1.5)
         rhomin=.5
      else
         rhomin=0.
      endif

      if(lreaddiag)then
         write(*,*)'rho   ','rholocal',
     $     ' ;  rho is from psum, rholocal from Ti file'
         do i=1,nrhere
            write(*,*)rho(i,1),rholocal(i,1)
         enddo
         write(*,*)'End of rho comparison'
      endif

      phiinf=0.

      jmin=1
      jmax=nthhere
      if(lreaddiag)write(*,*)'jmin,jmax',jmin,jmax
      do i=1,nrhere
         rhopic(i)=0.
         phicos(i)=0.
c         write(*,*)'th   tcc   phi'
         do j=jmin,jmax
            rhopic(i)=rhopic(i)+rholocal(i,j)
c     rhopic(i)=rhopic(i)+rho(i,j)
c \int cos(\theta) \phi(\theta) d\cos(\theta)
            phicos(i)=phicos(i)+th(j)*tcc(j)*phi(i,j)
c            write(*,'(4f10.4)')th(j),tcc(j),phi(i,j),phicos(i)
         enddo
         rhopic(i)=rhopic(i)/float(jmax-jmin+1)
      enddo

c     rescale rho; but usually this is the identity transformation.
      do i=1,nrhere
         rpicleft(i)=-rpic(i)
      enddo

      return
c End of data fix-up section
c__________________________________________________________________
 202  write(*,*)"nr error"
      call exit
 203  write(*,*)"rpicphipic error"
      call exit
 204  write(*,*)"nsteps error"
      call exit
 210  write(*,*)'Error opening file: ',filename(1:50)
      call exit
 101  ierr=101
      end


      
