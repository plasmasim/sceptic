c Version of July 2005 changes in tempmul.
c Update July 2006. Included in cvs.
c fcol Aug 2006 for collisional use.
      program forcegen
      integer nptmax,nsmax
      parameter (nptmax=100,nsmax=100,nfilemax=20)
      character*100 string,filename(nfilemax),charin,string2,compfile
      character*100 abrstring
      character debchar
      real charge1(nptmax),ffield1(nptmax),felec1(nptmax), fion1(nptmax)
     $     ,ftot1(nptmax),ffl1(nptmax)
      real charge2(nptmax),ffield2(nptmax),felec2(nptmax), fion2(nptmax)
     $     ,ftot2(nptmax),ffl2(nptmax)
      real fcombined(nptmax),fcomb2(nptmax)
c      real fcolk(nptmax),fcol(nptmax),fcoul(nptmax),
      real v(nptmax),vprbarr(nptmax)
      real phit1(nptmax),phit2(nptmax),vtheory(nptmax),qt1(nptmax)
      real colnwtarr(nptmax),fluxarr(nptmax),changlaflux(nptmax)
      real collampe(nptmax),fluxbyoml(nptmax),debyearr(nptmax)
      real chenline(nptmax),chencol(nptmax)
      real xs(nptmax),ys(nptmax)
      integer icoltarr(nptmax)
c      logical ldecompose,ltheory,
      logical louter,lchargeplot,lpotplot
c      logical lasymp,lliboff,lmason,lmason2,lstandard,
      logical lfcolk,lfcomp
      logical ldterm1
c Zobnin variables.
      integer nkapts,nz2,nz3,nzmax
      parameter (nkapts=14,nz2=11,nz3=17,nzmax=40)
      integer xyarr(2,nkapts),xya2(2,nz2),xya3(2,nz3)
      integer nlen(3)
      real roverl(nzmax,3),rlog(nzmax,3),fltptl(nzmax,3)
      logical linit
c Data from Lampe et al POP03 for Ti=.01, a/lambda_s=.015
      parameter (nflampe=16)
      real flampe(nflampe),colflampe(nflampe),flm(nflampe),clm(nflampe)
      data flampe/0.,-25,-127.,-265.,-495.,-715.,-905.,-1195.,-1585.,
     $     -1925.,-2235.,-2535.,-2865.,-3215.,-3585.,-3955./
      data colflampe/0.,10.,50.,107.,207.,327.,437.,647.,967.,
     $     1327.,1697.,2177.,2787.,3567.,4537.,5637./
c upper point scaled positions and values:
      data cl1/5997./fl1/-3585./clv1/0.5/flv1/2./
c
c     Floating potential for spherical collisional collection of Zobnin
c     JETP 91, 483 (2000). For neon, Ti/Te=.01, versus \lambda_s/mfp.
c Traced data:
c a/R_D=0.06
      data xyarr/
     $   1967,2955,2497,3285,2977,3585,3677,4015,4257,4345,4757,4625,
     $   5227,4865,5617,5015,5937,5095,6277,5125,6557,5065,6857,4925,
     $   7067,4775,7357,4515/          
c 0.012
      data xya2/
     $	 1987,3135,2907,3695,4007,4375,4767,4835,5407,5185,5917,5425,
     $	 6287,5565,6597,5645,6957,5715,7327,5745,8107,5735/
c 0.24
      data xya3/
     $   1997,2595,2457,2935,2947,3225,3357,3435,3777,3615,4307,3805,
     $	 4757,3925,5087,3975,5327,3995,5617,3955,5837,3875,6137,3715,
     $ 	 6367,3555,6597,3355,6807,3135,7067,2835,7377,2445/
      data nlen/nkapts,nz2,nz3/
      data linit/.false./
      save
c Default curve to use
      if(it.eq.0)it=1
c Initialize if needed.
      if(.not.linit)then
c axis locator points of tracing
         xm1=3630
         xp1=8340
         y0=7350
         y2=2565
         do i=1,nkapts
            rlog(i,1)=(-1.+2.*(xyarr(1,i)-xm1)/(xp1-xm1))
            roverl(i,1)=10**rlog(i,1)
            fltptl(i,1)=-2.*(xyarr(2,i)-y0)/(y2-y0)
         enddo
         do i=1,nz2
            rlog(i,2)=(-1.+2.*(xya2(1,i)-xm1)/(xp1-xm1))
            roverl(i,2)=10**rlog(i,2)
            fltptl(i,2)=-2.*(xya2(2,i)-y0)/(y2-y0)
         enddo
         do i=1,nz3
            rlog(i,3)=(-1.+2.*(xya3(1,i)-xm1)/(xp1-xm1))
            roverl(i,3)=10**rlog(i,3)
            fltptl(i,3)=-2.*(xya3(2,i)-y0)/(y2-y0)
         enddo
         linit=.true.
      endif

c     Defaults.
      louter=.true.
      lchargeplot=.false.
      lpotplot=.false.
      lfcolk=.false.
      lfcomp=.false.
      ldterm1=.true.
      Ti=1.
      phip=-2.5
      rmtoz=1.
      ifile=0
      finmax=0.
      tempmul=1.
      ipfnum=3
      colnwtmin=1.e20
      ic0=1
      
c     Deal with arguments
      do 1 i=1,iargc()
         call getarg(i,string)
c     write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
            if(string(1:2) .eq. '-o') louter=.false.
            if(string(1:2) .eq. '-c') lchargeplot=.true.
            if(string(1:2) .eq. '-p') lpotplot=.true.
            if(string(1:2) .eq. '-k') lfcolk=.true.
            if(string(1:2) .eq. '-a') ldterm1=.false.
            if(string(1:2) .eq. '-n') ipfnum=-ipfnum
            if(string(1:2) .eq. '-r') read(string(3:),*,err=71,end=71
     $           )rmtoz
            if(string(1:2) .eq. '-y') read(string(3:),*,err=71,end=71
     $           )finmax
            if(string(1:2) .eq. '-v') read(string(3:),*,err=71,end=71
     $           )tempmul
 71         if(string(1:2) .eq. '-?') goto 99
            if(string(1:2) .eq. '-f') then
               lfcomp=.true.
               read(string(3:),'(a)') compfile
            endif
         else
            if(ifile.ge.nfilemax) stop 'Too many input files'
            ifile=ifile+1
            filename(ifile)=string
         endif
 1    continue
 3    continue
      if(i.eq.1)goto 99
      write(*,*)'Read ',ifile,' files', (filename(jf),jf=1,ifile)
      do 100 jf=1,ifile
         open(10,file=filename(jf),status='old',err=200)
         do i=1,nptmax
 11         read(10,fmt='(a)',err=201,end=201)charin
c     write(*,*)charin
c     Section to parse file names to determine parameters.        
            if(charin(1:1).eq.'=')then
               do j=1,nsmax-4
c     write(*,*)charin(j:j)
                  if(charin(j:j).eq.'T' .and. charin(j+4:j+4).eq.'v'
     $                 )then
                     read(charin(j+1:j+1),'(i1)')itemp
                     read(charin(j+3:j+3),'(i1)')itpow
                     if(charin(j+2:j+2).eq.'e')then
                        Ti=float(itemp)*10.**itpow
                     else
                        Ti=float(itemp)/10.**itpow
                     endif
                     if(i.eq.1)write(*,*)'Ti=',Ti,'  tempmul=',tempmul
                     read(charin(j+5:j+7),'(i3)')ivel
                     v(i)=ivel*.01
                     if(i.eq.1) then
                        read(charin(j+15:j+15),'(i1)')ideb
                        read(charin(j+16:j+16),'(a1)')debchar
                        read(charin(j+17:j+17),'(i1)')idebp
                        if(debchar.eq.'e') then
                           debye=ideb*10**idebp
                        else
                           debye=ideb*10.**(-idebp)
                        endif
                     endif
                     goto 12
                  endif
               enddo
c     section for reading either general files with top data, or old
c     summary.
 12            read(10,'(a)')charin
               icolntype=0
               if(istrstr(charin,'dt').ne.0)then
                  read(10,*, err=120,end=120)
     $                 dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen
     $                 ,vprobe,damplen,Bz,icolntype,colnwt
 120              read(10,'(a)')charin
                  write(*,'(a,a)')'  dt    vd     Ti   steps  ',
     $                 'rhoinf phiinf  fave  debyelen Vp'
                  write(*,
     $          '(2f7.4,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3,f6.2,f7.3)')
     $              dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
                  if(icolntype.ne.0)then 
                     write(*,'(i4,f8.3)')icolntype,colnwt
c     Find minimum:
                     if(colnwt.lt.colnwtmin)then
                        if(colnwt .ne. 0.)then
                           colnwtmin=colnwt
                           icmin=i
                        else
                           ic0=i
                        endif
                     endif
                     debyearr(i)=debyelen
                     colnwtarr(i)=colnwt
                     icoltarr(i)=icolntype
                     collampe(i)=colnwt/sqrt(2.*Ti)
     $                 *(debyelen/sqrt(1+1./Ti))
                     ZTei=1./Ti
                     Zmepi=sqrt(4.*3.14159/(1837.*rmtoz))
                     omlf=-omlfloat(0.,ZTei,Zmepi)
                     foml=sqrt(Ti)*(1.-omlf/Ti)/sqrt(2.*3.14159)
                     fluxbyoml(i)=fave/foml
                  else
                     colnwtarr(i)=0.
                  endif
                  vprbarr(i)=vprobe
                  fluxarr(i)=fave
                  read(10,'(a)')charin
               else
                  vprbarr(i)=0.
                  colnwtarr(i)=0.
                  fluxarr(i)=fave
               endif
c Chang and Laframboise collisional ion flux:
               changlaflux(i)=changlaf(Ti,Vprobe,colnwt+1.e-4)
               write(*,*)'colnwt,kappai,expion,changlaflux(i)',
     $              colnwt,kappai,expion,changlaflux(i)
c 
               read(charin,*,end=121)charge1(i),ffield1(i), felec1(i)
     $              ,fion1(i),ftot1(i)
 121           read(10,*)charge2(i),ffield2(i),felec2(i),fion2(i)
     $              ,ftot2(i)
               ffl1(i)=ffield1(i)*debyelen**2
               ffl2(i)=ffield2(i)*debyelen**2               
            else
               goto 11
            endif
c 
            fourpi=4.*3.14159
            if(vprbarr(i).ne.0.)then
               write(*,*)'Using vprobe information ',vprbarr(i),
     $              ' not charge-based estimate',phip
               write(*,*)'OML=',-omlfloat(v(i)/sqrt(2.*Ti),1./Ti, sqrt(4
     $              .*3.14159/(1837.*rmtoz)))
               
               phip=vprbarr(i)
            endif
         enddo
 201     write(*,*)' End of file',jf,'   Contains',i-1,' cases.'
         close(10)
         i=i-1
c--------------------------------------------------
c Start of plotting etc.
         if(colnwtmin.ne.1.e20)then
            colmin=10.**(nint(alog10(colnwtmin)-2.49999))
            do kk=1,i
               if(colnwtarr(kk).lt. colnwtmin) then
                  colnwtarr(kk)=colmin
               endif
            enddo
            colnwtmin=colmin
         endif

         
         if(lchargeplot)then
c........................
c     Plotting of Charge
            if(jf.eq.1)then
               call pfset(ipfnum)
               call pltinit(0.,5.,-80.,-20.)
c     call pltinit(0.,5.,-250.,-0.)
               call charsize(.02,.02)
               call ticset(.015,.015,-.03,-.03,0,0,0,0)
               call axis()
               call axlabels('v!df!D  /(ZT!de!d/m!di!d)!u1/2!u',
     $              'Charge  /(!Ae!@!d0!dr!dp!dT!de!d/e)')
c     call legendline(-.55,1.,257,'  !Al!@!dD!d   T!di!d')
               call legendline(-.6,1.,257,'  !Al!@!dDe!d')
            endif
            call color(jf)
            icr=index(filename(jf),'i')
c     write(*,*)filename(jf)(1:20),icr
            if(icr.ne.0)then
c     write(string,'('' '',f4.0,'' '',f3.1,'' i'')')debye,Ti
               write(string,'('' '',f4.0,'' i'')')debye
            else
c     OML charge estimate only for floating
               ZTei=1./Ti
               Zmepi=sqrt(4.*3.14159/(1837.*rmtoz))
               do k=1,nptmax
                  vtheory(k)=5.*k/nptmax
                  phit1(k)=-omlfloat(vtheory(k)/sqrt(2./ZTei),ZTei,Zmepi
     $                 )
                  vt2t=vtheory(k)**2+3.*Ti
c     xlst=debye*sqrt(vt2t/(3.+vt2t))
                  xlst=debye
                  qt1(k)=phit1(k)*4*3.14159*(1.+1./xlst)
               enddo
               call dashset(1)
               call polyline(vtheory,qt1,nptmax)
               call dashset(0)
c     write(string,'('' '',f4.0,'' '',f3.1,'' f'')')debye,Ti
               write(string,'('' '',f4.0,'' f'')')debye
            endif
            call legendline(-.6,1-.05*jf,-jf,string)
            call winset(.true.)
            call polyline(v,charge1,i)
            call polymark(v,charge1,i,jf)
            call winset(.false.)
            if(jf.eq.ifile)call pltend()
         elseif(lpotplot)then
c........................
c     Plotting of Potential
            xleg=.6
            yleg=.95
            if(jf.eq.1)then
c     analytic form
               ZTei=1./Ti
               Zmepi=sqrt(4.*3.14159/(1837.*rmtoz))
               do k=1,nptmax
                  vtheory(k)=5.*k/nptmax
                  phit1(k)=-omlfloat(vtheory(k)/sqrt(2.),1.,Zmepi)
                  phit2(k)=-omlfloat(vtheory(k)/sqrt(2./ZTei),ZTei,Zmepi
     $                 )
               enddo
               call pfset(ipfnum)
               call minmax(phit2,nptmax,phimin,phimax)
c               write(*,*)phimin,phimax,phimax-phimin
               if(phimax-phimin .lt. .5)phimin=phimax-.5
               call fitinit(0.,5.,phimin,phimax)
c     call pltinit(0.,5.,-2.8,-1.7)
               call charsize(.02,.02)
               call axis()
               call axlabels('   v!df!d   /(ZT!de!d/m!di!d)!u1/2!u',
     $              '      !Af!@!df!d  /(T!de!d/e)')
               call legendline(xleg,yleg+.01,257
     $              ,'  !Al!@!dDe!d   T!di!d')
               call dashset(2)
               call polyline(vtheory,phit1,nptmax)
               call polyline(vtheory,phit2,nptmax)
               call dashset(0)
            endif
            call color(jf)
            if(Ti.lt.0.1) then
               write(string,'('' '',f4.0,'' '',f4.2)')debyelen,Ti
            else
               write(string,'('' '',f4.0,'' '',f3.1)')debyelen,Ti
            endif
c     call legendline(xleg,yleg-.05*jf,-jf,string)
            call legendline(xleg,yleg-.05*jf,jf,string)
c     call polyline(v,vprbarr,i)
            call polymark(v,vprbarr,i,jf)
            if(jf.eq.ifile)call pltend()
         else
c........................
c Default Plotting
            write(*,*)'Debye=',debye
            write(*,*
     $           )'Velocity  Sceptic-Inner  -Outer  Analytic  Anal/Sc-I'
     $           ,' Khrapak   K/Sc-I'
            write(*,'(7f10.4)')(v(j),ftot1(j),ftot2(j),fcombined(j),
     $           fcombined(j)/ftot1(j), fcomb2(j),fcomb2(j)/ftot1(j),
     $           j=1,i)
            if(icolntype.ne.0)then
            write(*,*
     $           )'Velocity    Flux    Sceptic-Inner  -Outer  Colnwt'
               write(*,'(5f10.4)')(v(j),fluxarr(j),
     $           ftot1(j),ftot2(j),colnwtarr(j)
     $              ,j=1,i)
               ZTei=1./Ti
               Zmepi=sqrt(4.*3.14159/(1837.*rmtoz))
               omlf=-omlfloat(0.,ZTei,Zmepi)
c Zero drift form.
               foml=sqrt(Ti)*(1.-omlf/Ti)/sqrt(2.*3.14159)
               write(string2,'('' OML M/Z='',i2)')nint(rmtoz)
               potlabr=-potlka(1./debyelen)
c Electron flux density at potential potlabr in units of sqrt(T_e/m_i).
               fabr=exp(potlabr)/Zmepi/sqrt(2.)
               abrstring=' ABR !Al!@!dDe!d='
               call fwrite(debyelen,iwidth,1,
     $              abrstring(lentrim(abrstring)+1:))
               write(*,*)'rmtoz,omlfloat,potlabr,fabr',
     $              rmtoz,omlf,potlabr,fabr
               do kl=1,nflampe
c Scale into Lampe's units:
                  flampe(kl)=1.+flampe(kl)*flv1/fl1
                  colflampe(kl)=colflampe(kl)*clv1/cl1
                  write(*,*)colflampe(kl),flampe(kl)
c Scale into my units
                  flm(kl)=flampe(kl)*foml
c Unfortunately, Lampe uses lambda_s instead of lambda_De
                  clm(kl)=colflampe(kl)*sqrt(2.*Ti)
     $                 *sqrt(1+1./Ti)/debyelen
                  write(*,*)clm(kl),flm(kl)
               enddo
            endif
            call pfset(ipfnum)
            call minmax(colnwtarr,i,cmin,cmax)
            if(cmin.ne.cmax)then
c We have a number of different collisionalities.
c Force plot:
               call lautomark(colnwtarr,ftot1,i,.true.,.false.,1)
               call axlabels(
     $'Collision Frequency !An!@!dc!d/[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]'
     $              ,'Drag Force /[r!dp!d!u2!un!de!dT!de!d]')
               call pltend()
c Flux plot
               call minmax(fluxarr,i,fmin,fmax)
               fmin=min(fmin,foml*.9)
               fmax=max(fmax,fabr*1.05)
               call pltinit(cmin,cmax,fmin,fmax)
               call fitrange(fmin,fmax,6,nxfac,xfac,xdelta,fymin,fymax)
               call scalewn(cmin,cmax,fymin,fymax,.true.,.false.)
               call axis()
               call axlabels(
     $'Collision Frequency !An!@!dc!d/[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Flux /[n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u4!Ap!@r!dp!d!u2!u]')
c Alternative axis.
               sfac=(debyelen/sqrt(1+1./Ti))/sqrt(2.*Ti)
               call axptset(1.,1.)
               call ticrev()
               call altxaxis(sfac,sfac)
               call jdrwstr(.6,.67,
     $    'Collision Frequency !An!@!dc!d/[v!dti!d/!Al!@!ds!d]',0.)
               call axptset(0.,0.)
               call ticrev()

               do kk=1,i
                  imark=4
                  if(icoltarr(kk).eq.1)imark=5
                  call polymark(colnwtarr(kk),fluxarr(kk),1,imark)
               enddo
               call sortincrease(colnwtarr,fluxarr,icoltarr,i,2,
     $              xs,ys,ns)
               call polyline(xs(2),ys(2),ns-1)
               call sortincrease(colnwtarr,fluxarr,icoltarr,i,1,
     $              xs,ys,ns)
               call polyline(xs,ys,ns)
c               call polymark(colnwtarr,fluxarr,i,1)
               call legendline(.5,.05,5,' SCEPTIC 1')
               call legendline(.5,.1,4,' SCEPTIC 2')
c Compare with Chang and Laframboise.
               call winset(.true.)
               call color(5)
               do kk=1,i
                  imark=1
                  if(icoltarr(kk).eq.1)imark=2
                  call polymark(colnwtarr(kk),changlaflux(kk),1,imark)
               enddo
               call dashset(4)
               call sortincrease(colnwtarr,changlaflux,icoltarr,i,2,
     $              xs,ys,ns)
               call polyline(xs,ys,ns)
               call dashset(3)
               call sortincrease(colnwtarr,changlaflux,icoltarr,i,1,
     $              xs,ys,ns)
               call polyline(xs,ys,ns)
               call dashset(0)
               call legendline(.5,.15,2,' Continuum 1')
               call legendline(.5,.2,1,' Continuum 2')
c               call polymark(colnwtarr,changlaflux,i,imark)
               colmax=10.
               colmin=.1
               do kk=1,nptmax
                  chencol(kk)=colmin*
     $                exp((log(colmax)-log(colmin))*(kk-1.)/(nptmax-1.))
                  chenline(kk)=changlaf(Ti,Vprobe,chencol(kk))
               enddo
               call winset(.true.)
               if(vprbarr(1).eq.vprbarr(i))then
                  call polyline(chencol,chenline,nptmax)
               endif
               call color(15)
c Compare with heuristic collisional effects.
c Zero drift form.
c               foml=sqrt(Ti)*(1.-Vprobe/Ti)/sqrt(2.*3.14159)
c               write(*,*)'foml(0)=',foml
c Finite drift form.
c               foml=sqrt(2.*Ti)*omlux(vd/sqrt(2.*Ti),-vprbarr(ic0)/Ti)/4.
               write(*,*)vprbarr(ic0),Ti,'  OML flux=',foml
c               write(*,*)'OML Flux(0,Ti)',
c     $              sqrt(2.*Ti)*omlux(0.*vd/sqrt(2.*Ti),-Vprobe/Ti)/4.
               colleft=colnwtmin
               call vecw(colleft,foml,0)
               call vecw(colleft*5.,foml,1)
               call drcstr(' OML')
c mark abr
               call vecw(colleft,fabr,0)
               call vecw(colleft*5.,fabr,1)
               call drcstr(abrstring)
c Overload chencol/line
               colmax=.1
               colmin=colnwtmin
               slambda=sqrt(1+debyelen**2/(1+1/(Ti+vd**2)))
               rt=slambda*log(1.-2.*Vprobe/(3.*Ti*(1+slambda)))
               write(*,*)'debyelen,slambda,rt=',debyelen,slambda,rt
c               rt=debyelen
               do kk=1,nptmax
                  chencol(kk)=colmin*
     $                exp((log(colmax)-log(colmin))*(kk-1.)/(nptmax-1.))
                  chenline(kk)=foml+sqrt(1./(4.*3.14159))*
     $                 rt*(rt+1)**2*chencol(kk)
               enddo
c Heuristic collision correction line:
c               call polyline(chencol,chenline,nptmax)
c Add Lampe line:
               call dashset(2)
               call color(ibrickred())
               call polyline(clm(2),flm(2),nflampe-1)
               call legendline(.5,.25,0,' Lampe')
               call color(15)
               call dashset(0)
c work around accis bug:
               call vecw(1.,1.,0)
               call pltend()
c Linear plot intended for Lampe comparisons
               call minmax(fluxarr,i,fmin,fmax)
               call minmax(clm,nflampe,cfmin,cfmax)
               call fitinit(0.,cfmax,min(fmin,flm(1)),1.3*fmax)
               call axis()
               call polymark(colnwtarr,fluxarr,i,1)
               call axlabels(
     $'Collision frequency !An!@!dc!d/[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Flux /[n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u4!Ap!@r!dp!d!u2!u]')
c mark oml
               call vecw(0.,foml,0)
               call vecw(.1*cfmax,foml,1)
               call drcstr(string2)
c mark abr
               call vecw(0.,fabr,0)
               call vecw(.1*cfmax,fabr,1)
               call drcstr(abrstring)
               call color(2)
               call polyline(clm,flm,nflampe)
               call color(15)
               call pltend()
c Plot in Lampe's normalization:
c               call pfpsset(1)
               call fitinit(0.,.5,0.,4.)
               call axis
c               call polymark(collampe,fluxbyoml,i,1)
c Alternative coded plot:
               do kl=1,i
                  mark=4
                  if(debyearr(kl)-66.7.lt.1.) mark=2
                  if(icoltarr(kl).eq.2)
     $                 call polymark(collampe(kl),fluxbyoml(kl),1,mark)
               enddo
               call axlabels(
     $     'Collision Frequency !An!@!dc!d/[v!dti!d/!Al!@!ds!d]',
     $              'Flux /OML-value',)
               call polyline(colflampe,flampe,nflampe)
c               call pfpsset(0)
               call legendline(.5,.25,4,' !Al!@!dDe!d=666.7 r!dp!d')
               call legendline(.5,.2,2,' !Al!@!dDe!d=66.7 r!dp!d')
               call legendline(.5,.15,0,' Lampe et al')
               call pltend()
c Floating potential plot
               call minmax(vprbarr,i,vpmin,vpmax)
c               vpmin=min(vpmin,1.1*omlf)
               call fitinit(0.,.025,vpmin,vpmax)
               call fitrange(-vpmax,-vpmin,5,nxfac,xfac,xdelta,
     $              vymin,vymax)
               write(*,*)vpmin,vpmax,vprbarr(ic0),vymin,vymax
               call scalewn(cmin,cmax,-vymin,-vymax,.true.,.false.)
               call axis()
               call polymark(colnwtarr,vprbarr,i,1)
               call axlabels(
     $'Collision Frequency !An!@!dc!d/[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Sphere Potential /[T!de!d/e]')
               call vecw(colleft,omlf,0)
               call vecw(colleft*5.,omlf,1)
               call drcstr(string2)
               sfac=(debyelen/sqrt(1+1./Ti))/sqrt(2.*Ti)
               call axptset(1.,1.)
               call ticrev()
               call altxaxis(sfac,sfac)
               call jdrwstr(.6,.67,
c     $    'Collision Frequency !An!@!dc!d/[v!dti!d/!Al!@!ds!d]',0.)
     $    '!An!@!dc!d/[v!dti!d/!Al!@!ds!d]',1.2)
               call axptset(0.,0.)
               call ticrev()
               if(rmtoz.eq.20) then
c Zobnin comparison
                  icurve=1
                  if(abs(debyelen-833.) .lt. 10.) icurve=2
                  if(abs(debyelen-41.6) .lt. 2.) icurve=3
                  sfac=sfac*sqrt(3.1415926)/2.
                  call scalewn(sfac*cmin,sfac*cmax,
     $                 -vymin,-vymax,.true.,.false.)
                  call color(ibrickred())
                  call polyline(roverl(1,icurve),
     $                 fltptl(1,icurve),nlen(icurve))
                  call legendline(.1,.1,0,' Zobnin et al')
                  call color(15)
               endif
               call pltend()
            endif
         endif
 100  continue

      call exit(0)
 220     stop 'lfcomp file open error'

 200  write(*,*)' Error opening file',filename(jf)
 99   write(*,*)'Usage: ./fcol [switches] file1, file2, ...\n',
     $     ' produce file by headsum, or head -n2 ;tail -n 3 T*.dat'
      write(*,*)'Switches:'
      write(*,*)'-c chargeplot, -d no decompose, -t no theory,',
     $     ' -p potential plot',
     $     '-o no outer  -r?? ratio of mass to charge.',
     $     '-a don''t add 1 to debyelen (default add to compensate).',
     $     '-ffilename comparison force file.  -yf.f max force (y). ',
     $     '-vf.f set tempmul (1). -n no screen output.'

      end
c**************************************************************************
      function colkforce(vd,phip,Ti,xl2)
c v<<v_ti approximation to collision force, per Khrapak
c according to which F= 8 sqrt(2pi)/3 a^2 n_i M v_Ti vd [ 1+rho(v_Ti)/2a
c + rho^2/4a^2 ln\Lambda]
c where v_Ti=sqrt(T_i/m_i), 
c rho=b_90=q1q2/(4pi\epsilon_0 M v^2)
c The first two terms are just the collection force, 
c the last the orbit force (which I originally omitted).
c However, since, Khrapak defines the collection impact parameter as
c a sqrt(1+2rho_0/a), he really means to define rho in terms of which
c the potential at a is rho_0 v^2 / a. Therefore a consistent definition
c of rho is 
c    rho(vt) = |phip /v_Ti^2|  
c and in fact only thus do we get full agreement
c and a is probe size=1. Normalization has M=1, Te=1 etc.
      colkforce= 8.*sqrt(2.*3.14159)/3.*sqrt(Ti)*vd *(
     $     1.-phip/(2.*Ti)
     $     +(phip/(2.*Ti))**2 *0.5*xl2
     $     )
      end
c ********************************************************************
c Collection force
      function colforce(u,chi)
      if(u.ne.0.)then
      colforce=(u*(2.*u**2+1+2.*chi)*exp(-u**2) +
     $     (4.*u**4+4.*u**2-1-2.*(1-2.*u**2)*chi)*sqrt(3.14159)*0.5
     $     *(1.-erfcc(u)) )/u**2
      else
         colforce=0.
      endif
      end


      function chandrafunc(x)
      if(x.ne.0.)then
         chandrafunc=((1.-erfcc(x))-2.*x*(exp(-x**2)/sqrt(3.14159)))
     $     /(2.*x**2)
      else
         chandrafunc=0.
      endif
      end


      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
c*********************************************************************
      function omlux(Ui,x)
c The functional dependence of OML collection vs drift U, and potl x.
      U=Ui
c More accurate to make the cut-off lower 1.e-2 vs 1.e-5
      if(U.lt.1.e-2) U=1.e-2
      omlux=( (1.+(1.+2.*x)/(2.*U**2))*(1.-erfcc(U))
     $     + exp(-U**2)/(sqrt(3.14159)*U) )*U
      end
c*****************************************************************
      function omlfloat(U,ZTei,Zmepi)
c Floating potential according to OML approximation.
c U is ratio of drift to ion thermal velocity
c ZTei is ZT_e/T_i
c Zmepi is sqrt(Z m_e 4\pi/m_i)
c Value returned is in units of ZT_e
      x1=2.*ZTei
 1    x=x1
      x1=-ZTei*alog( 0.25*(Zmepi/sqrt(ZTei))*omlux(U,x) )
c      write(*,*)x,omlux(u,x)
      if(abs(x1-x)/x.gt.1.e-5) goto 1
      omlfloat=x/ZTei
      end
c******************************************************************
c Return the position (not pointer) of match in string
c Trailing blanks are ignored.
      function istrstr(string,match)
      integer istrstr
      character*(*) string,match
      ls=lentrim(string)
      lm=lentrim(match)
      istrstr=0
      do i=1,ls
         do j=1,lm
            if(string(i+j-1:i+j-1).ne.match(j:j)) goto 101
         enddo
c Here when matched.
         istrstr=i
         goto 102
 101  enddo
 102  return
      end
c******************************************************************
c Obtain the length of a string omitting trailing blanks.
      function lentrim(string)
      character*(*) string
      do i=len(string),1,-1
         if(string(i:i).ne.' ') goto 101
      enddo
      i=0
 101  lentrim=i
      end
c********************************************************************
c Chen and Laframboise collisional ion flux:
c Needed to fix vti ambiguity there's a factor of order unity.
c First guess was low by factor 1.329:
c      kappai=(4/3.)*sqrt(2.*Ti)/(colnwt)
      function changlaf(Ti,Vprobe,colnwt)
      kappai=sqrt(3.141592*2.*Ti)/colnwt
      expion=exp(Vprobe/Ti)
      fluxii=(1.+kappai)*(-Vprobe)/
     $     (Ti*(1.-(1.+kappai*Vprobe/Ti)*expion))
      changlaf=fluxii*Ti/(1+kappai)/(colnwt)
      end
c********************************************************************
c Return ordered arrays for plotting etc. Very inefficient!
      subroutine sortincrease(x,y,m,ni,mc,xs,ys,ns)
c Input arrays, to sort in order of increasing x
      real x(ni),y(ni)
      integer m(ni)
c Use only those points having m(i)=mc
      integer mc
c Output arrays of sorted couples
      real xs(ni),ys(ni)
c With this many selected points:
      integer ns

      ns=0
      xv=1.e-30
      do i=1,ni
         xmin=1.e30
         k=0
         do j=1,ni
            if(x(j).lt.xmin .and. x(j).gt.xv .and. m(j).eq.mc)then
               k=j
               xmin=x(j)
            endif
         enddo
         if(k.eq.0) goto 1
         ns=ns+1
         xs(ns)=x(k)
         ys(ns)=y(k)
         xv=xmin
      enddo
 1    continue
      end
