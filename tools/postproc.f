C Postprocessing of a pair of output files.
c Will still work without the Ti file in a pinch.
c
c********************************************************************
c Split out the file reading part from the rest. July 2004.
c Removed infinity recalculation, and related code.
c Made the angle average the full pi not just upstream.
c********************************************************************
      character*100 string,filename
      include 'piccompost.f'
c      include 'cic/piccompost.f'

      real rholocal(0:NRFULL,0:NTHFULL)
c      real thanglocal(0:NTHFULL)
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      real phiyukawa(1000)
      integer nti0
      parameter (nti0=100)
      real rti0(nti0)
      real phiti0(nti0)
      character*100 charin
      real fluxofangle(nth),cflux(nth)
      integer jstepth
      logical lpcic,ltempc,lphip,lreaddiag,lgraph,larrows,lconline
      logical lvfwrt,lseptp,lunlabel,ledge,ldens,langle,lnlog
      data lpcic/.false./
      data ltempc/.false./
      data lconline/.false./
      data lphip/.false./
      data lreaddiag/.false./
      data lgraph/.true./
      data larrows/.false./
      data lvfwrt/.true./
      data lseptp/.false./
      data lunlabel/.false./
      data ledge/.false./
      data ldens/.false./
      data langle/.false./
      data lnlog/.false./
      data jstepth/1/

c Deal with arguments
      do 1 i=1,iargc()
         call getarg(i,string)
C         write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
         if(string(1:2) .eq. '-r') lreaddiag=.true.
         if(string(1:2) .eq. '-x') lgraph=.false.
c Legacy usage for summarize:         
         if(string(1:2) .eq. '-p') lgraph=.false.
         if(string(1:2) .eq. '-c') lpcic=.true.
         if(string(1:2) .eq. '-f') lphip=.true.
         if(string(1:2) .eq. '-n') ldens=.true.
         if(string(1:2) .eq. '-t') ltempc=.true.
         if(string(1:2) .eq. '-l') lconline=.true.
         if(string(1:2) .eq. '-a') larrows=.true.
         if(string(1:2) .eq. '-v') lvfwrt=.false.
         if(string(1:2) .eq. '-s') lseptp=.true.
         if(string(1:2) .eq. '-u') lunlabel=.true.
         if(string(1:2) .eq. '-e') ledge=.true.
         if(string(1:2) .eq. '-b') langle=.true.
         if(string(1:3) .eq. '-o') lnlog=.true.
         if(string(1:2) .eq. '-j') read(string(3:),*)jstepth
         if(string(1:2) .eq. '-?') goto 51
         else
            filename=string
         endif
 1    continue
 3    continue
      if(i.eq.1)goto 51

      do i=1,len(filename)-2
         if(filename(i:i+1).eq.'Sp')then
            filename(i:i+1)='Ti'
            write(*,*) 'Using Ti file ',filename
            goto 110
         endif
      enddo
 110  continue

c Read the outputfile
      call readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,nphere,
     $     phipic,rhopic,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,debyelen,vprobe,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)
      if(ierr.eq.101) goto 101


c Set arrow scale
      v1=max(1.,vd)
c Start of Plotting:

      call pfset(3)
      if(lgraph)then
c Now we make the first plot a time trace.
         call yautoplot(fluxprobe,nsteps)
         call axlabels('step','Particles to probe')
         call pltend()

c Old start of plotting.
         call multiframe(2,1,3)
         call autoplot(rpic,rhopic,nrhere)
         call axlabels('r','angle averaged density')
         call winset(.true.)
         call vecw(rpic(1),1.,0)
         call vecw(rpic(nrhere),1.,1)
         call winset(.false.)
         call minmax(rho(1,1),nthhere*nrhere,rhmin,rhmax)
         if(ledge)rhmax=min(rhmax,1.5)
         call pltinit(-rpic(nrhere),rpic(nrhere),0.,max(rhmax,1.1))
         call axis()
         call axis2()
         call axlabels('radial position','n/n!A!d;!d!@')
         do j=1,nthhere/2,jstepth
            call color(mod(j,15)+1)
            call winset(.true.)
            call polyline(rpicleft,rho(1,nthhere-j+1),nrhere)
            call polyline(rpic,rho(1,j),nrhere)
            call winset(.false.)
            write(charin,'(i3,f5.2)')j,tcc(j)
            call legendline(1.05,0.08*(j-1)/jstepth,
     $           0,charin(1:8)//char(0))
         enddo
         call color(15)
         call legendline(1.05,0.08*(j-1)/jstepth,258,
     $        'angle:   j |cos!Aq!@|')
         call vecw(-rpic(nrhere),1.,0)
         call vecw(rpic(nrhere),1.,1)
         call fwrite(vd,iwdth,2,charin)
         if(lvfwrt) call jdrwstr(0.05,0.72,
     $        'v!dd!d='//charin(1:iwdth)//char(0),1.)
         call pltend()
         call multiframe(0,0,0)
      endif
      if(lnlog)then
c         call lautoplot(rpic,rhopic,nrhere,.true.,.true.)
         call pltinit(0.,1.,0.,1.)
         call scalewn(1.,200.,1.,100.,.true.,.true.)
         call polyline(rpic,rhopic,nrhere)
         call axis()
         call axlabels('r','angle averaged density')
         call pltend()
         open(15,status='unknown',file='rhopic.dat')
         write(15,*)'dt,      vd,      Ti,      rmax,',
     $        '   fave, debyelen,    Vp [icoln,colnwt]'
         write(15,'(2f8.5,f8.4,f8.3,f8.3,f12.5,f10.5,i4,e13.4)')
     $     dt,vd,Ti,rmax,fave,debyelen,vprobe,icolntype,colnwt
         write(15,*)nrhere
         write(15,'(2f12.5)')(rpic(jj),rhopic(jj),jj=1,nrhere)
         close(15)
      endif
c Ti=0 quasineutral case:
      do j=1,nti0
         phiti0(j)=-0.5*j/float(nti0)
         rti0(j)=sqrt(exp(-0.5-phiti0(j))/sqrt(-2.*phiti0(j)))
      enddo
      do j=1,nrTi
         phipic(j)=phipic(j)-phiinf
      enddo
      goto 102
 101  write(*,*) 'Does not seem to be a Ti... file here.'
 102  continue

C End of stuff dependent on Ti file reading.
      open(13,status='unknown',file='phiout.dat')
      write(13,*)'dt,      vd,      Ti,      rmax,',
     $     '   fave, debyelen,    Vp [icoln,colnwt]'
      write(13,'(2f8.5,f8.4,f8.3,f8.3,f12.5,f10.5,i4,e13.4)')
     $     dt,vd,Ti,rmax,fave,debyelen,vprobe,icolntype,colnwt
      write(13,*)nrhere
      write(13,'(2f12.5)')(rpic(j),phipic(j),j=1,nrhere)

      if(lphip)then
         call minmax2(phi(1,1),nr+1,nrhere,nthhere,pmin,pmax)
c         write(*,*)'minmax',nr,nrhere,nthhere,pmin,pmax
         cscale=0.02
         if(abs(phipic(1)) .lt. 5.)cscale=0.1
         ppmax=1.2
         call pltinit(rpic(1),rpic(nrhere),pmin,ppmax)
         call polyline(rpic,phipic,nrhere)
         call axis()
         call vecw(rpic(1),0.,0)
         call vecw(rpic(nrhere),0.,1)
         vt2=vd**2+3.*Ti
c Linearized shielding length corrected for finite size.
         slambda=sqrt(debyelen**2*vt2/(3.+vt2) + 1.)
c         slambda=8.5
c         slambda=debyelen
c         write(*,*)'vt2,slambda,vprobe,phipic(1)'
c     $        ,vt2,slambda,vprobe,phipic(1)
         if(slambda.eq.0.) slambda=1.e-4
         do kk=1,nrhere
            phiyukawa(kk)=(vprobe*rpic(1)/exp(-rpic(1)/slambda))
     $           *exp(-rpic(kk)/slambda)/rpic(kk)
c            write(*,*)kk,rpic(kk),phiyukawa(kk)
         enddo
         call dashset(1)
         call polyline(rpic,phiyukawa,nrhere)
         call dashset(0)
         call axlabels('r',
     $        '<!Af!@>!A=Jf!@ dcos!Aq!@/2')
         call scalewn(rpic(1),rpic(nrhere),pmin*cscale,ppmax*cscale,
     $        .false.,.false.)
         call color(iblue())
         call winset(.true.)
         call polyline(rpic,phicos,nrhere)
         call winset(.false.)
         call axptset(1.,0.)
         call ticrev()
         call altyaxis(1.,1.)
         call ticrev()
         call axlabels('','!AJf!@cos!Aq!@ dcos!Aq!@')
         call color(15)
         call pltend()
c Potential slices
c         call multiframe(0,0,0)
         call axptset(0.,0.)
         call minmax(phi(1,1),nthhere*nrhere,rhmin,rhmax)
         if(ledge)rhmin=max(rhmin,-.4)
         call pltinit(-rpic(nrhere),rpic(nrhere),rhmin,max(rhmax,0.1))
         call axis()
         call axis2()
         call axlabels('radial position','!Af!@')
         do j=1,nthhere/2,jstepth
            call color(mod(j,15)+1)
            call winset(.true.)
            call polyline(rpicleft,phi(1,nthhere-j+1),nrhere)
            call polyline(rpic,phi(1,j),nrhere)
            call winset(.false.)
            write(charin,'(i3,f5.2)')j,tcc(j)
            call legendline(-.48,0.08*(j-1)/jstepth,
     $           0,charin(1:8)//char(0))
         enddo
         call color(15)
         call legendline(-.48,0.08*(j-1)/jstepth,258,
     $        'angle:   j |cos!Aq!@|')
         call vecw(-rpic(nrhere),1.,0)
         call vecw(rpic(nrhere),1.,1)
         call fwrite(vd,iwdth,2,charin)
         if(lvfwrt) call jdrwstr(.05,0.72,
     $        'v!dd!d='//charin(1:iwdth)//char(0),1.)
         call winset(.true.)
         call dashset(4)
         call polyline(rpicleft,phipic,nrhere)
         call polyline(rpic,phipic,nrhere)
         call winset(.false.)
         call legendline(-.48,-0.08,
     $           0,'Upstream angle/time Average')
         call dashset(0)
         call pltend()
c Contouring
c         call condisphi(ir,jstepth,0.,vprobe,
c     $     nrhere,nthhere,v1,larrows,lconline)
         call condisphi(ir,jstepth,pmin,pmax,
     $     nrhere,nthhere,v1,larrows,lconline,lpcic,ledge)
         call pltend()
      endif

c Contouring:
      if(lseptp)then
         if(ltempc)then
            call condisplay2(ir,jstepth,rhomax,rhomin,
     $           nrhere,nthhere,v1,larrows,lconline)
         endif
      else
         if(ldens)then
            call conrho(ir,jstepth,rhomax,rhomin,
     $           nrhere,nthhere,v1,larrows,lconline,rholocal)
         endif
         if(ltempc)then
            call contemp(ir,jstepth,rhomax,rhomin,
     $           nrhere,nthhere,v1,larrows,lconline,ledge)
         endif
      endif

      if(lunlabel)then
         call condisunlabel(ir,jstepth,rhomax,rhomin,
     $           nrhere,nthhere,v1,larrows,lconline)
         call pltend()
      endif

      if(lreaddiag)write(*,*)'tcc, fluxofangle'
      if(.not.lgraph)write(*,*)'cos(theta), flux, coulflux'
c Calculate the flux as a function of angle.
      do j=1,nthhere
c Seems to be an error here for the ngp version.
         fluxofangle(j)=ninth(j)*(nthhere-1.)/
     $        (4.*pi*rhoinf*dt*nastep)
         if(lpcic)then
            if(j.eq.1 .or. j.eq.nthhere)fluxofangle(j)=fluxofangle(j)*2.
         else
c Correct the ngp error: 2 Sep 2002
            fluxofangle(j)=fluxofangle(j)*nthhere/(nthhere-1.)
         endif
c         write(*,*)rhoinf,dt,nastep,ninth(j)
         cthang=acos(tcc(j))
         cunnorm=coulflux(cthang,-vprobe/Ti,vd/sqrt(Ti))
         cflux(j)=cunnorm
     $        /sqrt(2.*3.14159/Ti)
         if(.not.lgraph) write(*,'(4f10.5)')tcc(j),fluxofangle(j)
     $        ,cflux(j),cunnorm
c     $        ,cthang
      enddo

c Plot flux as a function of angle.
      if(lgraph .and. langle)then
         call multiframe(0,0,0)
         call autoplot(tcc(1),fluxofangle,nthhere)
         call axlabels('cos!Aq!@','flux density (normalized)')
         call axis2()
c         if(vprobe.eq.0)
         call polymark(tcc(1),cflux,nthhere,1)
         call fwrite(vd,iwdth,2,charin)
         if(lvfwrt) call jdrwstr(0.05,0.65,
     $        'v!dd!d='//charin(1:iwdth)//char(0),1.)
         call pltend()
      endif

c Plot angle for differents steps
      if(lgraph .and. langle)then
         ifirst=nsteps/2
         ispace=20
         do k=ifirst,nsteps,ispace
            do kk=1,nthhere
               cflux(kk)=ninthstep(kk,k)
            enddo
            if(lpcic)then
c     fix up as double on boundary.
               cflux(1)=2.*cflux(1)
               cflux(nthhere)=2.*cflux(nthhere)
            endif
            if(k.eq.ifirst)then
               call autoplot(tcc(1),cflux(1),nthhere)
               call axlabels('cos!Aq!@','particles per step')
               call axis2()
               call fwrite(vd,iwdth,2,charin)
               call jdrwstr(0.05,0.65,
     $              'v!dd!d='//charin(1:iwdth)//char(0),1.)
               call winset(.true.)
            endif
            call color((k-ifirst)*14/(nsteps-ifirst))
            call polyline(tcc(1),cflux(1),nthhere)
         enddo
         call pltend()
      endif

      call exit
 51   write(*,*)"Usage: postproc [-f -x ... ] filename"
      write(*,*)"-f plot angle-averaged phi and phi contours.",
     $     "-n contour n.      -t contour T.     -l contour lines.",
     $     "-c CIC output. -j<n> theta plot step.",
     $     "-r print diagnostics on file reading.",
     $     "-i recalculate rho-infinity, -v do not write v_f value",
     $     "-a put velocity arrows on T plots. -s separate T plots",
     $     "-u plot unlabelled full sphere.  -x do no line graphs.",
     $     "-e use density contours spaced closer to 1, for edge.",
     $     "-b plot angle distributions.",
     $     "-o plot log(n) versus log(r), angle averaged."
      
      call exit
      end
c***************************************************************************
c Contouring of the charge density, rho, on distorted mesh.
      subroutine conrho(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline,rholocal)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline
c Common data:
      include 'piccompost.f'
      real rholocal(0:NRFULL,0:NTHFULL)
c      include 'cic/piccompost.f'
c      save
      character*20 cstring
      character*30 tstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Conrho error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif
c Correct the outside angle centers if necessary.
      if(tcc(1).eq.1)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
      enddo

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax


      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+rhomin
      enddo
      icl=-ncont

c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,25000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rholocal(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)
c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
            call fitrange(rhomin,rhomax,ncont,ipow,
     $           fac10,delta,first,xlast)
c               write(*,*)'rhomax,fac10,first,delta',
c     $              rhomax,fac10,first,delta
            do j=1,ncont
               zclv(j)=(first+j*delta)
            enddo
            ntype=2
            icl=(ncont-1)
            call contourl(rholocal(1,0),cworka,NRFULL+1,nrhere,
     $           nthhere+2,zclv,icl,zrho,xrho,ntype)
            write(*,'(a,30f5.2)')'Density Contours=',zclv
            call fwrite(delta,iwd,1,cstring)
            tstring=' contour spacing: '//cstring
c            call legendline(-.1,-.22,258,tstring)
         endif
c      endif
c Fit closer than boxtitle
      call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
c      call boxtitle('n/n!A!d;!d!@')
      call axis()
      call axlabels('z','r sin!Aq!@')
  
      
      if(larrows) then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere,it
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              '!A_!@',0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0))
         call charsize(0.,0.)
         write(cstring,'(''  v='',f4.1)') v1
         call color(15)
         call legendline(0.8,0.95,258,cstring(1:8)//char(0))
      endif
      
      call pltend()
      end

c***************************************************************************
c Contouring of the temperature.
      subroutine contemp(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline,ledge)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline,ledge
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
c      character*20 cstring
c      character*30 tstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(ncont),zclv2(ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      real Trave(NRFULL+1),Ttpave(NRFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif
c Correct the outside angle centers if necessary.
      if(tcc(1).eq.1)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
      enddo
      rpmax=rcc(nrhere)

c Temperature plots.
      call multiframe(2,1,3)
      Tmax=0.
      Tmin=1000.
      Trmax=0.
      Trmin=1000.
      do i=1,nrhere
         Trave(i)=0.
         Ttpave(i)=0.
         do j=1,nthhere
            Tr(i,j)=(vr2sum(i,j)*psum(i,j)-vrsum(i,j)**2)
     $           /(psum(i,j)**2+1.e-5)
            Ttp(i,j)=0.5*(vtp2sum(i,j)*psum(i,j)-
     $           (vpsum(i,j)**2+vtsum(i,j)**2) )
     $           /(psum(i,j)**2+1.e-5)
            if(Tr(i,j).lt.Trmin)Trmin=Tr(i,j)
            if(Tr(i,j).gt.Trmax)Trmax=Tr(i,j)
            if(Ttp(i,j).lt.Tmin)Tmin=Ttp(i,j)
            if(Ttp(i,j).gt.Tmax)Tmax=Ttp(i,j)
            Trave(i)=Trave(i)+Tr(i,j)
            Ttpave(i)=Ttpave(i)+Ttp(i,j)
         enddo
         Trave(i)=Trave(i)/nthhere
         Ttpave(i)=Ttpave(i)/nthhere
         Tr(i,0)=Tr(i,1)
         Tr(i,nthhere+1)=Tr(i,nthhere)
         Ttp(i,0)=Ttp(i,1)
         Ttp(i,nthhere+1)=Ttp(i,nthhere)
      enddo
      Tmax=0.8*Tmax
      
c Tr plot
c      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
c         zclv(j)=Trmin+(Tmax-Trmin)*(1.*(j-1)/float(ncont-1))
c      enddo
c      icl=-ncont
c
      icl2=0
      if(ledge) then
         Trmin=0.5*Ti
         Tmin=Trmin
         Trmax=5*Ti
         Tmax=Trmax
         icl2=20
         do ii=1,icl2
            zclv2(ii)=Ti*(1.+(ii-icl2/4.)/icl2)
         enddo
c         if(ledge)write(*,*)'Tr-contours',(zclv2(ii),ii=1,icl2)
      endif

      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      zclv(1)=Trmin
      zclv(2)=Tmax
      icl=-2
      ntype=2+16
      if(.not.lconline)ntype=ntype+32
      call contourl(Tr(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      ntype=2
      zclv(1)=10
      icl=0
      if(lconline)
     $     call contourl(Tr(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv1,icl,zrho,xrho,ntype)
      write(*,*)'Trmin,Trmax',Trmin,Trmax
      call ticrev()
      call gradlegend(Trmin,Tmax,
     $     1.1,0.1,1.1,0.9,0.03,.false.)
      call ticrev()
      call legendline(1.03,0.5,258,'T!dr!d'//char(0))
      call ticnumset(8)
      call axis()
      call axlabels('z','r sin!Aq!@')
      if(larrows)then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              char(ichar('_')+128)//char(0),0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0)) 
         call charsize(0.,0.)
         call color(15)
         call legendline(0.8,0.95,258,'  v=1'//char(0))
      endif
c      call pltend()
c Ttp plot
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      zclv(1)=Trmin
      zclv(2)=Tmax
      icl=-2
      ntype=2+16
      if(.not.lconline)ntype=ntype+32
      call contourl(Ttp(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      ntype=2
      zclv2(1)=10
      icl=0
      if(lconline)
     $     call contourl(Ttp(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv1,icl,zrho,xrho,ntype)
      write(*,*)'Tmin,Tmax',Tmin,Tmax
      call ticrev()
      call gradlegend(Trmin,Tmax,
     $     1.1,0.1,1.1,0.9,0.03,.false.)
      call ticrev()
      call legendline(1.03,0.5,258,'T!d!A`!@!d'//char(0))
      call ticnumset(8)
      call axis()
      call axlabels('z','r sin!Aq!@')
      if(larrows)then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              char(ichar('_')+128)//char(0),0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0))
         call charsize(0.,0.)
         call color(15)
         call legendline(0.8,0.95,258,'  v=1'//char(0))
      endif
      call multiframe(0,0,0)
      call pltend()

      call lautoplot(rcc,Ttpave,nrhere,.false.,.true.)
      call axlabels('r','T')
      call legendline(.8,.95,0,'T!d!A`!@!d')
      call dashset(1)
      call legendline(.8,.9,0,'T!dr!d')
      call polyline(rcc,Trave,nrhere)
      call pltend()
      end

c***************************************************************************
c***************************************************************************
c Contouring of the potential, on distorted mesh.
      subroutine condisphi(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline,lpcic,ledge)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline,lpcic,ledge
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
      character*20 cstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(2*ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      real x2rho(NRFULL+1,0:NTHFULL+1),z2rho(NRFULL+1,0:NTHFULL+1)
c      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
            z2rho(i,j)=zrho(i,j)
            x2rho(i,j)=xrho(i,j)            
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.
         if(lpcic)then
            z2rho(i,1)=zrho(i,0)
            x2rho(i,1)=xrho(i,0)
            z2rho(i,nthhere)=zrho(i,nthhere+1)
            x2rho(i,nthhere)=xrho(i,nthhere+1)
         endif
      enddo
c Adjust the outside angle centers if necessary.
c Should never be necessary.
c      if(tcc(1).eq.1)then
c         tcc(1)=0.25*(3.+tcc(2))
c         tcc(0)=1.
c         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
c         tcc(nthhere+1)=-1.
c      endif

c ledge fixing
      if (ledge) then
         rhomax=rhomax/2.
      endif

      zscale=log(1000.)/(ncont-1.)
      do j=1,ncont
c Logarithmic
         zclv(j)=rhomax + (rhomin-rhomax)*
     $        exp(zscale*float(ncont-j))/exp(zscale*float(ncont-1))
c Linear
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+ rhomin
c         write(*,*) "zcle : ",zclv(j)
      enddo
c      write(*,*) "min,max",rhomin,rhomax

c      write(*,*)'Contours=',zclv
      icl=-ncont

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax
c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,35000,130000,65000,130000)
      ntype=2+16+32
      call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)

c Call a second time for contours, without the highest.
         if(lconline)then
c If this is not very bipolar
            if(abs(rhomax)-abs(rhomin).gt.0.2*abs(rhomax-rhomin))then
c Rational logarithmic contours
               do indx=1,ncont-2,3
                  base=10.**(-2+(indx-1)/3)
                  zclv(ncont+indx)=-base
                  zclv(ncont+indx+1)=-base*2.
                  zclv(ncont+indx+2)=-base*5.
c     New positive contours.
                  zclv(ncont-indx)=base
                  zclv(ncont-(indx+1))=base*2.
                  zclv(ncont-(indx+2))=base*5.
               enddo
               zclv(ncont)=0.
            endif
            write(*,'(a,30f6.2)')'Contours=',zclv
            ntype=2
            icl=(ncont-1)
c            write(*,'(2f8.4)')
c     $           (zrho(nrhere,k),z2rho(nrhere,k),k=0,nthhere+1)
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(ncont),icl,zrho,xrho,ntype)
c Fixed contour plot avoiding the hacked coordinates:
            call color(igreen())
            call contourl(phi(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(ncont),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(iskyblue())
c            call contourl(phi(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
c     $           zclv(1),icl,zrho,xrho,ntype)
            call contourl(phi(1,1),cworka,NRFULL+1,nrhere,nthhere,
     $           zclv(1),icl,z2rho(1,1),x2rho(1,1),ntype)
            call color(15)
         endif
c      endif
      call legendline(0.47,1.07,258,'!Af!@'//char(0))
      call axis()
      call axlabels('z','r sin!Aq!@')
      if(larrows) then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere,it
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              '!A_!@',0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0))
         call charsize(0.,0.)
         write(cstring,'(''  v='',f4.1)') v1
         call color(15)
         call legendline(0.8,0.95,258,cstring(1:8)//char(0))
      endif
c     call pltend()
      end

c***************************************************************************

c***************************************************************************
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
      rave=rmom1/(i2-i1+1)
      rlen=r(i2)-r(i1)
      slope=12.*(rmom1)/(rlen*rlen)/(i2-i1+1)
c      write(*,*)rmom0,rmom1,r(i1),r(i2),rmid,rlen
      end
c***************************************************************************
c Contouring of the charge density, rho, on distorted mesh. And temperature.
c This version temperature plots are separate.
      subroutine condisplay2(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
      character*20 cstring
      character*30 tstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif
c Correct the outside angle centers if necessary.
      if(tcc(1).eq.1)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
      enddo

      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1)) + rhomin
      enddo
      icl=-ncont

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax
c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,25000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rho(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)
c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
            call fitrange(rhomin,rhomax,ncont,ipow,
     $           fac10,delta,first,xlast)
c               write(*,*)'rhomax,fac10,first,delta',
c     $              rhomax,fac10,first,delta
            do j=1,ncont
               zclv(j)=(first+j*delta)
            enddo
            ntype=2
            icl=(ncont-1)
            call contourl(rho(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $           zclv,icl,zrho,xrho,ntype)
            write(*,*)'Density Contours=',zclv
            call fwrite(delta,iwd,1,cstring)
            tstring=' contour spacing: '//cstring
c            call legendline(-.1,-.22,258,tstring)
         endif
c      endif
c Fit closer than boxtitle
      call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
c      call boxtitle('n/n!A!d;!d!@')
      call axis()
      call axlabels('z','r sin!Aq!@')
      call color(12)
      if(ir.le.0.or.ir.ge.100) ir=10
      do j=1,nthhere,it
         do i=1,nrhere,max(nrhere/ir,1)
            vri=vrsum(i,j)/(psum(i,j)+1.e-5)
            vti=vtsum(i,j)/(psum(i,j)+1.e-5)
            size=basesize/v1*sqrt(vri**2+vti**2)
            angle=atan2(vti,vri)+acos(tcc(j))
            call charsize(size,0.3*size)
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $           '!A_!@',0.)
         enddo
      enddo
      call charangl(0.)
      size=basesize
      call charsize(size,0.3*size)
      call legendline(0.8,0.95,258,'!A_!@'//char(0))
      call charsize(0.,0.)
      write(cstring,'(''  v='',f4.1)') v1
      call color(15)
      call legendline(0.8,0.95,258,cstring(1:8)//char(0))
      call pltend()

c Temperature plots.
c      call multiframe(2,1,3)
      Tmax=0.
      Tmin=1000.
      Trmax=0.
      Trmin=1000.
      do i=1,nrhere
         do j=1,nthhere
            Tr(i,j)=(vr2sum(i,j)*psum(i,j)-vrsum(i,j)**2)
     $           /(psum(i,j)**2+1.e-5)
            Ttp(i,j)=0.5*(vtp2sum(i,j)*psum(i,j)-
     $           (vpsum(i,j)**2+vtsum(i,j)**2) )
     $           /(psum(i,j)**2+1.e-5)
            if(Tr(i,j).lt.Trmin)Trmin=Tr(i,j)
            if(Tr(i,j).gt.Trmax)Trmax=Tr(i,j)
            if(Ttp(i,j).lt.Tmin)Tmin=Ttp(i,j)
            if(Ttp(i,j).gt.Tmax)Tmax=Ttp(i,j)
         enddo
         Tr(i,0)=Tr(i,1)
         Tr(i,nthhere+1)=Tr(i,nthhere)
         Ttp(i,0)=Ttp(i,1)
         Ttp(i,nthhere+1)=Ttp(i,nthhere)
      enddo
      Tmax=0.8*Tmax
      
c Tr plot
c      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
c         zclv(j)=Trmin+(Tmax-Trmin)*(1.*(j-1)/float(ncont-1))
c      enddo
c      icl=-ncont
c
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      zclv(1)=Trmin
      zclv(2)=Tmax
      icl=-2
      ntype=2+16
      if(.not.lconline)ntype=ntype+32
      call contourl(Tr(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      ntype=2
      zclv(1)=10
      icl=0
      if(lconline)
     $     call contourl(Tr(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      write(*,*)'Trmin,Trmax',Trmin,Trmax
      call gradlegend(Trmin,Tmax,
     $     .1,1.25,.9,1.25,0.02,.true.)
      call legendline(0.47,1.07,258,'T!dr!d')
c      call ticrev()
c      call gradlegend(Trmin,Tmax,
c     $     1.1,0.1,1.1,0.9,0.03,.false.)
c      call ticrev()
c      call legendline(1.03,0.5,258,'T!dr!d'//char(0))
      call ticnumset(8)
      call axis()
      call axlabels('z','r sin!Aq!@')
      if(larrows)then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              char(ichar('_')+128)//char(0),0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0)) 
         call charsize(0.,0.)
         call color(15)
         call legendline(0.8,0.95,258,'  v=1'//char(0))
      endif
      call pltend()
c Ttp plot
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      zclv(1)=Trmin
      zclv(2)=Tmax
      icl=-2
      ntype=2+16
      if(.not.lconline)ntype=ntype+32
      call contourl(Ttp(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      ntype=2
      zclv(1)=10
      icl=0
      if(lconline)
     $     call contourl(Ttp(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      write(*,*)'Tmin,Tmax',Tmin,Tmax
      call gradlegend(Trmin,Tmax,
     $     .1,1.25,.9,1.25,0.02,.true.)
      call legendline(0.47,1.07,258,'T!d!A`!@!d')
c      call ticrev()
c      call gradlegend(Trmin,Tmax,
c     $     1.1,0.1,1.1,0.9,0.03,.false.)
c      call ticrev()
c      call legendline(1.03,0.5,258,'T!d!A`!@!d'//char(0))
      call ticnumset(8)
      call axis()
      call axlabels('z','r sin!Aq!@')
      if(larrows)then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              char(ichar('_')+128)//char(0),0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0))
         call charsize(0.,0.)
         call color(15)
         call legendline(0.8,0.95,258,'  v=1'//char(0))
      endif
      call multiframe(0,0,0)
c     call pltend()
      end

c***************************************************************************
c***************************************************************************
c Contouring of the charge density, rho, on distorted mesh. Unlabelled.
c Top and bottom.
      subroutine condisunlabel(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline
c Common data:
      include 'piccompost.f'
c      include 'cic/piccompost.f'
c      save
      character*20 cstring
      character*30 tstring
      character cworka(nr*(NTHFULL+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(ncont)
      real xrho(NRFULL+1,0:NTHFULL+1),zrho(NRFULL+1,0:NTHFULL+1)
c      real Tr(NRFULL+1,0:NTHFULL+1),Ttp(NRFULL+1,0:NTHFULL+1)
      real xrhob(NRFULL+1,0:NTHFULL+1)
      save xrho,zrho,xrhob

      external ACCIRCLE
      real basesize
      parameter (basesize=.02)

      if(nthhere.gt.NTHFULL)then
         write(*,*)' Condisplay error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',NRFULL,NTHFULL
         stop
      endif
c Correct the outside angle centers if necessary.
      if(tcc(1).eq.1)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2-1.e-3)
            xrhob(i,j)=-xrho(i,j)
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
      enddo

      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1)) + rhomin
      enddo
      icl=-ncont

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax
c      call multiframe(2,2,3)

      call ticnumset(8)
c top half
      call pltinaspect(-rpmax,rpmax,-rpmax,rpmax)
      call accisgradinit(-25000,00000,25000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rho(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call contourl(rho(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrhob,ntype)
c      call gradlegend(zclv(1),zclv(abs(icl)),
c     $     .1,1.25,.9,1.25,0.02,.true.)
c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
            call fitrange(rhomin,rhomax,ncont,ipow,
     $           fac10,delta,first,xlast)
c               write(*,*)'rhomax,fac10,first,delta',
c     $              rhomax,fac10,first,delta
            do j=1,ncont
               zclv(j)=(first+j*delta)
            enddo
            ntype=2
            icl=(ncont-1)
            call contourl(rho(1,0),cworka,NRFULL+1,nrhere,nthhere+2,
     $           zclv,icl,zrho,xrho,ntype)
            write(*,*)'Density Contours=',zclv
            call fwrite(delta,iwd,1,cstring)
            tstring=' contour spacing: '//cstring
c            call legendline(-.1,-.22,258,tstring)
         endif
c      endif
c Fit closer than boxtitle
c      call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
c      call boxtitle('n/n!A!d;!d!@')
c         call axis()
c         ilength=1
         call charsize(wx2nx(2.)-wx2nx(0.),wx2nx(2.)-wx2nx(0.))
         call accircle(wx2nx(0.),wy2ny(0.))
         call pathfill()
c      call axlabels('z','r sin!Aq!@')
      call color(ibrickred())
      if(ir.le.0.or.ir.ge.100) ir=10
      do j=1,nthhere,it
         do i=1,nrhere,max(nrhere/ir,1)
            vri=vrsum(i,j)/(psum(i,j)+1.e-5)
            vti=vtsum(i,j)/(psum(i,j)+1.e-5)
            size=basesize/v1*sqrt(vri**2+vti**2)
            angle=atan2(vti,vri)+acos(tcc(j))
            call charsize(size,0.3*size)
c top
            call charangl(180.*angle/3.141593)
            call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $           '!A_!@',0.)
c bottom
            call charangl(-180.*angle/3.141593)
            call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrhob(i,j)),
     $           '!A_!@',0.)
         enddo
      enddo
      call charangl(0.)
      size=basesize
c      call charsize(size,0.3*size)
c      call legendline(0.8,0.95,258,'!A_!@'//char(0))
c      call charsize(0.,0.)
c      write(cstring,'(''  v='',f4.1)') v1
c      call color(15)
c      call legendline(0.8,0.95,258,cstring(1:8)//char(0))
      call multiframe(0,0,0)
c     call pltend()
      end

c***************************************************************************
c Data reading subroutine
      subroutine readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,nphere,
     $     phipic,rhopic,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,debyelen,vprobe,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)
      logical lreaddiag,lpcic,ledge
      character*100 string,filename
      real phipic(1000),rhopic(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      include 'piccompost.f'
      real rholocal(0:NRFULL,0:NTHFULL)
      character*256 charin
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      ierr=0
c Read the data file:
c__________________________________________________________________

      open(10,file=filename,status='old',err=101)
c Line for nothing.
      read(10,*)charin
      read(10,'(a)')charin
c      write(*,*)charin
      read(charin,*,err=201,end=201)
     $     dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe,damplen,Bz
     $     ,icolntype,colnwt
 201  continue
      write(*,'(a,a)')'  dt    vd     Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp damplen  Bz...'
      write(*,'(2f7.4,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3,f6.2,f7.3,$)')
     $     dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe,damplen,Bz
      if(icolntype.gt.0)then
         write(*,'(i2,f7.3)',err=212)icolntype,colnwt
      else
         write(*,*)
      endif
 212  read(10,*,err=202)nrTi
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
      read(10,*,err=200)
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
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)
     $     icolntype,colwt,Eneutral,vneutral,Tneutral
 402  close(10)
      write(*,*)'nrhere,nthhere,icolntype,colwt'
      write(*,*)nrhere,nthhere,icolntype,colwt
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
c For now, use rholocal to fix --ds problem, maybe wrong ?
            if(rholocal(k1,k2).gt.rhomax)rhomax=rholocal(k1,k2)
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
         rholocal(k1,0)=rholocal(k1,1)
         rholocal(k1,nthhere+1)=rholocal(k1,nthhere)
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

