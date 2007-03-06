c Version of July 2005 changes in tempmul.
c Update July 2006. Included in cvs.
      program forcegen
      integer nptmax,nsmax
      parameter (nptmax=100,nsmax=100,nfilemax=20)
      character*100 string,filename(nfilemax),charin,string2,compfile
      character debchar
      real charge1(nptmax),ffield1(nptmax),felec1(nptmax), fion1(nptmax)
     $     ,ftot1(nptmax),ffl1(nptmax)
      real charge2(nptmax),ffield2(nptmax),felec2(nptmax), fion2(nptmax)
     $     ,ftot2(nptmax),ffl2(nptmax)
      real fcol(nptmax),fcoul(nptmax),fcombined(nptmax),fcomb2(nptmax)
      real fcolk(nptmax)
      real v(nptmax),vprbarr(nptmax)
      real phit1(nptmax),phit2(nptmax),vtheory(nptmax),qt1(nptmax)
      real colnwtarr(nptmax),fluxarr(nptmax),changlaflux(nptmax)
      real chenline(nptmax),chencol(nptmax)
      logical ldecompose,ltheory,louter,lchargeplot,lpotplot
      logical lasymp,lliboff,lmason,lmason2,lstandard,lfcolk,lfcomp
      logical ldterm1
c Data from Lampe et al POP03 for Ti=.01, a/lambda=.015
      parameter (nflampe=14)
      real flampe(nflampe),colflampe(nflampe)
      data flampe/0.,-265.,-495.,-715.,-905.,-1195.,-1585.,
     $     -1925.,-2235.,-2535.,-2865.,-3215.,-3585.,-3955./
      data colflampe/0.,107.,207.,327.,437.,647.,967.,
     $     1327.,1697.,2177.,2787.,3567.,4537.,5637./
c upper point scaled positions and values:
      data cl1/5997./fl1/-3585./clv1/0.5/flv1/2./

c Data from Khrapak PoP05 on lamda_eff/lambda_de vs v/sqrt(Ti/mi)
c for Te/Ti = 100 (case 1) and 10 (case 2).
c This lambda_eff is what to put in the coulomb log expression
c \ln\Lambda = \ln (\lambda_eff/R) - 1/2, where R=Ze^2/T_i(1+u^2) is b_90.
c It is for a point particle.
c      parameter (nlk=10)
c      real ulk(nlk)
c      real xlk1(nlk),xlk2(nlk)
c      data ulk/.1,.3,.6,1.,1.5,2.,3.,6.,10.,30./
c      data xlk1/0.1,0.101,0.108,0.123,0.158,0.210,0.380,0.805,0.995,1./
c      data xlk2/0.298,0.304,0.324,0.367,0.460,0.592,0.901,1.,1.,1./
c     Defaults.
      lstandard=.true.
      lmason=.false.
      lmason2=.false.
      lliboff=.false.
      lasymp=.false.
      ldecompose=.true.
      ltheory=.true.
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
            if(string(1:2) .eq. '-j') lstandard=.false.
            if(string(1:2) .eq. '-t') ltheory=.false.
            if(string(1:2) .eq. '-d') ldecompose=.false.
            if(string(1:2) .eq. '-o') louter=.false.
            if(string(1:2) .eq. '-m') lmason=.true.
            if(string(1:2) .eq. '-s') lasymp=.true.
            if(string(1:2) .eq. '-l') lliboff=.true.
            if(string(1:2) .eq. '-2') lmason2=.true.
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
                     colnwtarr(i)=colnwt
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
               ffl1(i)=ffield1(i)*debye**2
               ffl2(i)=ffield2(i)*debye**2               
            else
               goto 11
            endif
            fourpi=4.*3.14159
            
c     This fits ok for lambda=5
c     vt2=v(i)**2+(8./3.14159)*Ti
            if(lstandard) then 
               vt2=tempmul*Ti +v(i)**2
            else
               vt2=tempmul*Ti +v(i)**2
c Enhancement of transition to zero ion shielding.
     $              +v(i)**2*(v(i)/(.5+0.05*alog(rmtoz)
c Version of Orleans:
c     $                  +.1*(Ti/.01)**.5))**3
c it must be identical to 
c                        +Ti**.5))**3
c The following is a good fit from lambda=5 to 20 but there might be better
c     $                  +.1*(Ti/.01)**(.3+.04*debye)
c This is arguably as good but weights the upramp more than the peak.
     $                  +.1 +(debye/5.)**1*(sqrt(Ti)-.1)
     $              ))**3

               
            endif
            dterm=debye**2*(vt2/(tempmul+vt2))
            if(ldterm1)dterm=dterm+1
            xls=sqrt(dterm)
            
c     The following model is used for the Yukawa form:
c     phi_y = q_eff/4pi exp(-r/xls)/r.
c     
c     Then the first method of matching to the data is that
c     q_eff is chosen to make the model charge inside r=1 equal to the
c     measured charge (charge1). That leads to 
c     qeff=charge1(i)*exp(1./xls)*xls/(1.+xls)
c     charge1(i)=-4 pi r^2 d phi_y/dr
c     =q_eff*exp(-r./xls)(1./r+1/xls)
c     and that then gives the value of phi_y at r=1:
c     phi_y(1) = q_eff/4pi exp(-1/xls) = phip:
            phip=(charge1(i)/fourpi)*xls/(1.+xls)
c     However, the above notwithstanding. if the information for the
c     actual
c     probe potential is available then use that potential itself for
c     phip
c     and determine q_eff by requiring the model potential to match at r
c     =1,
c     namely q_eff = phip exp(1/xls).
c     For large lambda_d these agree quite closely.
c     They should agree to first order in r_p/lambda for a Yukawa form.
            if(vprbarr(i).ne.0.)then
               write(*,*)'Using vprobe information ',vprbarr(i),
     $              ' not charge-based estimate',phip
               write(*,*)'OML=',-omlfloat(v(i)/sqrt(2.*Ti),1./Ti, sqrt(4
     $              .*3.14159/(1837.*rmtoz)))
               
               phip=vprbarr(i)
            endif
            
            qeff=fourpi*phip*exp(1./xls)
c     phi2 is the potential at r=1 from an _unshielded_ charge of
c     magnitude qeff.
            phi2=qeff/fourpi
            
            b90=-phip/vt2
c            b90=-phip/(v(i)**2+2.*Ti+(v(i)/.5)**4*.008/(.03+Ti))
c This is the form that most accords with Khrapak's work:
            b90k=-phip/(v(i)**2+2.*Ti)
            bc2=1 - 2.*phip/vt2
            b902=b90*b90
c This is an obsolete form.
c     b90k=-phi2/(v(i)**2+2.*Ti)
c This is for testing purposes when we want the b90 to fall quickly.
c            b90k=b90
c            bc3=1 - 2.*phip/(v(i)**2+3.*Ti)
c            if(sqrt(bc3).gt.b90k .and. lstandard)then
c Abrupt jump transition for illustration purposes.
c               write(*,*)'sqrt(bc3),b90k',sqrt(bc3),b90k
c               bc2=1.
c               b902=0.
c               dterm=debye**2
c            endif
            
c Effective potential for collection (superceded).
            phicol=-phip/Ti
c Use the calculated OML value of potential at low lambda.
c It happens to give a better result.
            phicol=omlfloat(v(i)/sqrt(2.*Ti),1./Ti,
     $           sqrt(4.*3.14159/(1837.*rmtoz)))/Ti
c            phicol=0.
            fcol(i)=sqrt(3.14159)*Ti*colforce(v(i)/sqrt(2.*Ti),phicol)
c     $        -charge1(i)/(fourpi*Ti))
            
c     Asymptotics
c     Standard form 
            xlnlambda2=max(log((b902+dterm)/(b902+bc2)),0.)
c This modified form is the same as Khrapak but with modified b90/bc2.
            if(.not.lstandard)
c     $           xlnlambda2=2.*max(log((b90+xls)/sqrt(b90**2+bc2)),0.)
     $           xlnlambda2=2.*max(log((xls+b90)/(1.+b90)),0.)
c This modified form shows that it is the bc change that is most important.
c            xlnlambda2=2.*max(log((b90+xls)/sqrt(b90**2+bc3)),0.)
C     Forms moved into functions.
c My asymptotic form
            if(lasymp)xlnlambda2=2.*
     $           xlnmasymp(Ti,v(i),xls,sqrt(bc2),phip)
c     Liboff form
            if(lliboff)xlnlambda2=2.*xlnmliboff(Ti,v(i),xls,sqrt(bc2)
     $           ,phip)
c     Mason form
            if(lmason)xlnlambda2=2.*xlnmmason(Ti,v(i),xls,sqrt(bc2),phip
     $           )
c     Alternate Mason form:
            if(lmason2)xlnlambda2=2.*xlnmmason2(Ti,v(i),xls,sqrt(bc2)
     $           ,phip)
            
c     fcoul(i)=phi2**2*(1./(2.*Ti))*fourpi*xlnlambda2*
            fcoul(i)=phip**2*(1./(2.*Ti))*fourpi*xlnlambda2
     $           * chandrafunc(v(i)/sqrt(2*Ti))
            fcombined(i)=fcoul(i)+fcol(i)
            
c     xlnlambda2=max(log(1.+((dterm)/(b902+bc2))),0.)
c     Kilgore approximate integrated form gives results about 10x too
c     big:
c     xlnlambda2=.898*(log(1.+7.8*xls/b90))**2
c     Kilgore unintegrated approx is not so crazy, about factor 2 too
c     big.
c     xlnlambda2=.9369*(2./3.14159)*log(1+(61.32/4.)*dterm/b902)
c     Khrapak integrated form:
c            dtermk=debye**2/(1.+tempmul/(tempmul*Ti+v(i)**2))
            dtermk=debye**2/(1.+1./(Ti+v(i)**2))
c            dtermk=debye**2
            if(ldterm1)dtermk=dtermk+1
            xlnlambda3=2.*max(log((sqrt(dtermk)+b90k)/(1.+b90k)),0.)
c     Here there is ambiguity between phip and phi2 in Khrapak's papers.
c     phi2 gives larger drag values which are better at lower lambda but
c     worse for higher.
c     fcoul(i)=phi2**2*(1./(2.*Ti))*fourpi*xlnlambda2*
            fcoul(i)=phip**2*(1./(2.*Ti))*fourpi*xlnlambda3
     $           * chandrafunc(v(i)/sqrt(2*Ti))
            fcomb2(i)=fcoul(i)+fcol(i)
            fcolk(i)=colkforce(vd,phip,Ti,xlnlambda3)
            write(*,'(''b90='',f7.4,''  bc='',f7.4,''  charge='',f7.2,'/
     $           /'''  shield='',f7.4,''  xLm3='',f7.4)')b90,sqrt(bc2)
     $           ,charge1(i),sqrt(dterm),xlnlambda3
            if(lfcolk)write(*,*)'fcolk=',fcolk(i)
         enddo
 201     write(*,*)' End of file',jf,'   Contains',i-1,' cases.'
         close(10)
         i=i-1

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
               write(string,'('' '',f4.0,'' '',f4.2)')debye,Ti
            else
               write(string,'('' '',f4.0,'' '',f3.1)')debye,Ti
            endif
c     call legendline(xleg,yleg-.05*jf,-jf,string)
            call legendline(xleg,yleg-.05*jf,jf,string)
c     call polyline(v,vprbarr,i)
            call polymark(v,vprbarr,i,jf)
            if(jf.eq.ifile)call pltend()
         else
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
            endif
c     Plotting
            call pfset(ipfnum)
            call minmax(colnwtarr,i,cmin,cmax)
            if(cmin.ne.cmax)then
c We have a number of different collisionalities.
               call lautomark(colnwtarr,ftot1,i,.true.,.false.,1)
               call axlabels(
     $    'Collision frequency /[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $           'Drag Force /[r!dp!d!u2!un!de!dT!de!d]')
               call pltend()
               call lautomark(colnwtarr,fluxarr,i,.true.,.false.,1)
               call axlabels(
     $    'Collision frequency /[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Flux /[n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u4!Ap!@r!dp!d!u2!u]')
c Compare with Chang and Laframboise.
               call color(5)
               call polymark(colnwtarr,changlaflux,i,5)
               colmax=10.
               colmin=.1
               do kk=1,nptmax
                  chencol(kk)=colmin*
     $                exp((log(colmax)-log(colmin))*(kk-1.)/(nptmax-1.))
                  chenline(kk)=changlaf(Ti,Vprobe,chencol(kk))
               enddo
               call winset(.true.)
               call polyline(chencol,chenline,nptmax)
               call color(15)
c Compare with heuristic collisional effects.
c Zero drift form.
c               foml=sqrt(Ti)*(1.-Vprobe/Ti)/sqrt(2.*3.14159)
c               write(*,*)'foml(0)=',foml
c Finite drift form.
               foml=sqrt(2.*Ti)*omlux(vd/sqrt(2.*Ti),-vprbarr(ic0)/Ti)/4.
               write(*,*)vprbarr(ic0),Ti,'  OML flux=',foml
c               write(*,*)'OML Flux(0,Ti)',
c     $              sqrt(2.*Ti)*omlux(0.*vd/sqrt(2.*Ti),-Vprobe/Ti)/4.
               colleft=colnwtmin
               call vecw(colleft,foml,0)
               call vecw(colleft*5.,foml,1)
               call drcstr(' OML')
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
               call polyline(chencol,chenline,nptmax)
c work around accis bug:
               call vecw(1.,1.,0)
               call pltend()
c Linear plot intended for Lampe comparisons
               ZTei=1./Ti
               Zmepi=sqrt(4.*3.14159/(1837.*rmtoz))
               omlf=-omlfloat(0.,ZTei,Zmepi)
c Zero drift form.
               foml=sqrt(Ti)*(1.-omlf/Ti)/sqrt(2.*3.14159)
               write(*,*)'omlfloat,rmtoz',omlf,rmtoz
               write(string2,'('' OML M/Z='',i2)')nint(rmtoz)
               do kl=1,nflampe
c Scale into Lampe's units:
                  flampe(kl)=1.+flampe(kl)*flv1/fl1
                  colflampe(kl)=colflampe(kl)*clv1/cl1
                  write(*,*)colflampe(kl),flampe(kl)
c Scale into my units
                  flampe(kl)=flampe(kl)*foml
c Unfortunately, Lampe uses lambda_s instead of lambda_De
                  colflampe(kl)=colflampe(kl)*sqrt(2.*Ti)
     $                 *sqrt(1+1./Ti)/debyelen
                  write(*,*)colflampe(kl),flampe(kl)
               enddo
               call minmax(fluxarr,i,fmin,fmax)
               call minmax(colflampe,nflampe,cfmin,cfmax)
               call fitinit(0.,cfmax,min(fmin,flampe(1)),1.3*fmax)
               call axis()
               call polymark(colnwtarr,fluxarr,i,1)
               call axlabels(
     $    'Collision frequency /[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Flux /[n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u4!Ap!@r!dp!d!u2!u]')
               call vecw(0.,foml,0)
               call vecw(.1*cfmax,foml,1)
               call drcstr(string2)
               call color(2)
               call polyline(colflampe,flampe,nflampe)
               call color(15)
               call pltend()
c Floating potential plot
               call minmax(vprbarr,i,vpmin,vpmax)
               call fitinit(0.,.025,vpmax,vprbarr(ic0))
               call axis()
               call polymark(colnwtarr,vprbarr,i,1)
               call axlabels(
     $    'Collision frequency /[(ZT!de!d/m!di!d)!u1/2!u/r!dp!d]',
     $    'Sphere Potential /[T!de!d/e]')
               call vecw(0.,omlf,0)
               call vecw(.005,omlf,1)
               call drcstr(string2)
               call pltend()
            endif
c            call minmax(ftot2,i,fmin,fmax)
            call minmax(ftot1,i,fmin,fmax)
            if(ldecompose)then
c               fmin=-ftot2(i)*.3
               fmin=-fmax*.3
               fmax=fmax*1.2
            else
               fmin=0.
               fmax=fmax*1.1
            endif
            if(finmax.ne.0)fmax=finmax
            call ticnumset(10)
            call fitinit(0.,v(i),fmin,fmax)
            call charsize(.02,.02)
c     call axis()
c Instead of that, do things by hand
            call ticset(.015,.015,-.03,-.03,0,0,0,0)
            call yaxis(1.,0.)
            call xaxis(0.,1.)
            call ticset(.01,.01,-.03,-.03,0,0,0,0)
            call ticlabtog()
            call xaxis(0.,.5)
            call ticlabtog()
            call axis2()
            call axlabels('Drift velocity /[ZT!de!d/m!di!d]!u1/2!u',
     $           'Drag Force /[r!dp!d!u2!un!de!dT!de!d]')
c     charin='!Al!@!dD!d='
c     call fwrite(debye,iwidth,1,charin(12:))
            
            call fwrite(debye,iwidth1,1,string)
            call fwrite(Ti,iwidth2,2,string2)
            charin=' '
            charin(1:12)='!Al!@!dD!d='
            charin(12:11+iwidth1)=string(1:iwidth1)
            charin(12+iwidth1:21+iwidth1)=', T!di!d='
            charin(21+iwidth1:)=string2
            call legendline(.41,.95,257,charin(1:21+iwidth1+iwidth2))
            call winset(.true.)
            call polymark(v,ftot1,i,1)
            call polyline(v,ftot1,i)
            yleg=.95
            call legendline(.0,yleg,1,' Total SCEPTIC')
            call vecw(0.,0.,0)
            call vecw(v(i),0.,1)
            if(ldecompose)then
               call polymark(v,ffl1,i,2)
               call polyline(v,ffl1,i)
               call polymark(v,felec1,i,3)
               call polyline(v,felec1,i)
               call polyline(v,fion1,i)
               call polymark(v,fion1,i,4)
               call legendline(.0,yleg-.05,2,' E-field')
               call legendline(.0,yleg-.1,3,' Electrons')
               call legendline(.0,yleg-.15,4,' Ions')
               call legendline(.4,yleg-.12,0,' Inner Sphere')
               if(louter)then
                  call color(4)
                  call dashset(4)
                  call legendline(.4,yleg-.07,0,' Outer Boundary')
                  call polymark(v,ftot2,i,1)
                  call polyline(v,ftot2,i)
                  call polymark(v,ffl2,i,2)
                  call polymark(v,felec2,i,3)
                  call polyline(v,felec2,i)
                  call polymark(v,fion2,i,4)
                  call polyline(v,fion2,i)
               endif
            endif
            if(ltheory)then
               if(lfcolk)then
                  call color(iblue())
                  call dashset(2)
                  call polyline(v,fcolk,i)
                  call polymark(v,fcolk,i,7)
               endif
               call color(ired())
               call dashset(3)
               call polyline(v,fcol,i)
               yoffs=0.
               if(ldecompose) then
                  call polyline(v,fcoul,i)
               else
                  yoffs=.80
               endif
c     call dashset(1)
               call dashset(0)
               call polymark(v,fcombined,i,11)
               call polyline(v,fcombined,i)
               if(.not.(lasymp.and.lmason.and.lliboff.and.lmason2))then
                  if(lstandard) then
                     call legendline(.4,.1+yoffs,-11,' Standard')
                  else
                     call legendline(.4,.1+yoffs,-11,' Adjusted')
                  endif
               else
                  call legendline(.4,.1+yoffs,-11,' Asymptotic')
               endif
c     call dashset(5)
               call color(iblue())
               call polymark(v,fcomb2,i,10)
               call polyline(v,fcomb2,i)
               call legendline(.4,.05+yoffs,-10,' Khrapak')
            endif
            if(lfcomp)then
               write(*,*)'Doing force comparison; file ',compfile(1:30)
               open(10,file=compfile,status='old',err=220)
               do j=1,nptmax
                  read(10,*,end=221,err=221)debyelen,vdi,vpi,force,error
                  ftot2(j)=force
                  v(j)=vdi
               enddo
 221           close(10)
               call color(6)
               call polymark(v,ftot2,j-1,7)
               call polyline(v,ftot2,j-1)
            endif
            call pltend()
         endif
 100  continue

      call exit(0)
 220     stop 'lfcomp file open error'

 200  write(*,*)' Error opening file',filename(jf)
 99   write(*,*)'Usage: ./forcesum [switches] file1, file2, ...\n',
     $     ' produce file by headsum, or head -n2 ;tail -n 3 T*.dat'
      write(*,*)'Switches:'
      write(*,*)'-c chargeplot, -d no decompose, -t no theory,',
     $     ' -p potential plot',
     $     '-s use asymptotic form,',
     $     ' not standard:log((b902+dterm)/(b902+bc2)',
     $     '-m Mason form -l Liboff form -2 Mason2 form.',
     $     '-o no outer  -r?? ratio of mass to charge.',
     $     '-k plot Khrapak collection force (v<<v_ti approx)',
     $     '-a don''t add 1 to debyelen (default add to compensate).',
     $     '-j use adjusted transition form fit.',
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
c*******************************************************************
      function xlnmmason(Ti,v,xls,bc,phip)
c Log lambda for shielded coulomb potential,
c  per fit to Mason et al Phys Fl 10,1827(1967)
c  extended by incorporating a minimum impact parameter.
c T is ion temperature. v drift velocity. 
c xls=\lambda_s is shielding length.
c bc is the collection impact parameter
c phip is the probe potential, at radius 1.
c
c phi0 (\lambda_s/r) exp(-r/\lambda_s) = potential.
      phi0=phip*exp(1./xls)/xls
c b90L is the equivalent minimum impact parameter.
      b90L=-phi0/(Ti+1.12*v**2)
c Tstar is Mason's parameter: xls*T/phi0, 
c    but with flow- and collection-correction.
      Tstar=xls/sqrt(b90L**2+bc**2)
c a1 and a2 are fitting parameters.
      a1=1.75
      a2=7.0
      xlnmmason=alog(1.+a1*tstar) -(0.268+alog(a1))*tstar/(a2+tstar)
      end
c*******************************************************************
      function xlnmmason2(Ti,v,xls,bc,phip)
c Log lambda for shielded coulomb potential,
c  per fit to Mason et al Phys Fl 10,1827(1967). Alternate cutoff.
c  extended by incorporating a minimum impact parameter.
c T is ion temperature. v drift velocity. 
c xls=\lambda_s is shielding length.
c bc is the collection impact parameter
c phip is the probe potential, at radius 1.
c
c phi0 (\lambda_s/r) exp(-r/\lambda_s) = potential.
      phi0=phip*exp(1./xls)/xls
c b90L is the equivalent minimum impact parameter.
      b90L=-phi0/(Ti+1.12*v**2)
c Tstar is Mason's parameter: xls*T/phi0, 
c    with flow- but without collection-correction.
      Tstar=xls/abs(b90L)
c a1 and a2 are fitting parameters.
      a1=1.75
      a2=7.0
      xlnmmason2=alog(1.+a1*tstar) -(0.268+alog(a1))*tstar/(a2+tstar)
c Collection correction applied after the fact.
      xlnmmason2=xlnmmason2-0.5*alog(1+ bc**2/b90L**2)
      end
c******************************************************************
      function xlnmliboff(Ti,v,xls,bc,phip)
c Log lambda for shielded coulomb potential,
c  per Liboff. Simply made asymptotically correct.
c T is ion temperature. v drift velocity. 
c xls=\lambda_s is shielding length.
c bc is the collection impact parameter
c phip is the probe potential, at radius 1.
c
c phi0 (\lambda_s/r) exp(-r/\lambda_s) = potential.
      phi0=phip*exp(1./xls)/xls
c b90L is the equivalent minimum impact parameter.
      b90L=-phi0/(Ti+1.12*v**2)
c Tstar is Mason's parameter: xls*T/phi0, 
c    but with flow- and collection-correction.
      Tstar=xls/sqrt(b90L**2+bc**2)
c a1 and a2 are fitting parameters.
      xlnmliboff=max(alog(tstar)-0.268,0.)
      end
c*****************************************************************
      function xlnmasymp(Ti,v,xls,bc,phip)
c Log lambda for shielded coulomb potential,
c  per asymptotic approximations.
c T is ion temperature. v drift velocity. 
c xls=\lambda_s is shielding length.
c bc is the collection impact parameter
c phip is the probe potential, at radius 1.
      vt2=v**2+ 3.*Ti
      b902=(phip/vt2)**2
      bc2=bc**2
      a=sqrt(bc2/b902)
      sl=sqrt(xls**2/b902)
c     New asymptotics
      xlnmasymp=((sl+a)/sqrt(1.+(sl+a)**2))
     $        *alog(1.+sl/sqrt(1.+a**2))
     $     *exp(-a/sl)

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
c Need to fix vti ambiguity there's a factor of order unity.
c First guess was low by factor 1.329:
c      kappai=(4/3.)*sqrt(2.*Ti)/(colnwt)
      function changlaf(Ti,Vprobe,colnwt)
      kappai=sqrt(3.141592*2.*Ti)/colnwt
      expion=exp(Vprobe/Ti)
      fluxii=(1.+kappai)*(-Vprobe)/
     $     (Ti*(1.-(1.+kappai*Vprobe/Ti)*expion))
      changlaf=fluxii*Ti/(1+kappai)/(colnwt)
      end
