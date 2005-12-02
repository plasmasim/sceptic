c Plot from summary
C Contour plot from summaries. Based on sumplot Dec 2002
      integer imax,ifilemax,itmax
c Each file contains a column of data at particular debyelength, up to
c imax runs, maximum theta size itmax, up to ifilemax files.
      parameter (imax=100,ifilemax=30,itmax=100)
      real debyearr(imax,ifilemax),vparr(imax,ifilemax)
     $     ,debyel(imax,ifilemax)
      real cost(imax),flux(imax),fit(imax),dummy(imax)
      real gup(imax,ifilemax),gdown(imax,ifilemax)
     $     ,updown(imax,ifilemax),vdj(imax,ifilemax)
      real cosar(itmax,ifilemax,imax),fluxar(itmax,ifilemax,imax)
      character cworka(imax,ifilemax)
      integer numv(ifilemax)
      real Tiv(ifilemax),deb1(ifilemax),debl1(ifilemax)
      real coulf(imax),vpcoul(imax)
      character*100 filename
      character*100 charin
      character*100 carg
      logical lprint,lfitplot,lcos,llam
      integer nfit
      integer iskip
      real fluxmax
c Fit the curves by exponentials, on the up and downstream sides.
c gamma_up|down= cx*exp(vd*(up|down)x)
c And angular dependence is then
c gamma=cx*exp(vd*0.5*((1.-cost)upx + (1.+cost)downx)
      real upx,downx,cx
      data upx,downx,cx/0.64,-.7,.62/
      data fluxmax/0./
c Default
      lfitplot=.false.
      lprint=.false.
      lcos=.false.
      llam=.false.
      vdmax=0.
      iskip=1
      filename='summary'
      gmax=0.
      nfit=0

c Read command line arguments. Switches first, then filenames.
c Some switches require being before the filenames.
      is=0
      i=0
 11   i=i+1
      call getarg(i+is,carg)
      if(carg(1:2) .eq. '-s')then
         read(carg(3:),*)iskip
         i=i-1
         is=is+1
         goto 11
      elseif(carg(1:2) .eq. '-p')then
         lprint=.true.
         i=i-1
         is=is+1
         goto 11
      elseif(carg(1:2) .eq. '-f')then
         lfitplot=.true.
         i=i-1
         is=is+1
         goto 11         
      elseif(carg(1:2) .eq. '-c')then
         lcos=.true.
         i=i-1
         is=is+1
         goto 11         
      elseif(carg(1:2) .eq. '-l')then
         llam=.true.
         i=i-1
         is=is+1
         goto 11
      elseif(carg(1:2) .eq. '-n')then
         read(carg(3:),*)nfit
         i=i-1
         is=is+1
         goto 11         
      elseif(carg(1:2) .eq. '-m')then
         read(carg(3:),*)fluxmax
         i=i-1
         is=is+1
         goto 11         
      endif
      if(carg(1:1) .ne. ' ') then
         filename=carg
      else
         if(i.gt.1)then
            i=i-1
            goto 12
         else
            write(*,*)'You probably need to specify multiple files'
            write(*,*)'Usage ./fluxVp [-s -p -f] sumfile1 ...'
            write(*,*)' -sN skip N, -p print, -f fitplot,',
     $           ' -n? fit ? end points to cos-plot,'
            write(*,*)' -c cos(theta) plots, -l plot vs lambda.'
            write(*,*)' -m??.? specify max flux for -c plots'
            call exit
         endif         
      endif

      
      open(10,file=filename,status='old',err=105)
      call pfset(3)

      numv(i)=0
      do j=1,imax
         read(10,*,end=101,err=101)charin
         read(10,*)dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
         read(10,*)charin
         read(10,*)nrhere,nthhere
c         write(*,*)nrhere,nthhere
         write(*,901)dt,vd,Ti,isteps,fave,debyelen,vprobe
 901     format('dt=',f5.3,' vd=',f5.2,' Ti=',f4.1,' Steps=',i4,
     $        ' Flux=',f7.4,' Lambda=',e7.2,' Vp=',f5.1,$)
         write(*,*)j,i
         deb1(i)=debyelen
         debyearr(j,i)=debyelen
         debl1(i)=log10(debyelen)
         debyel(j,i)=log10(debyelen)
         vparr(j,i)=vprobe
         if(abs(vpmax).lt. abs(vprobe)) vpmax=vprobe
         read(10,*)charin
         if(j.eq.1 .and. lcos)then
            if(fluxmax.eq.0.)then
               if(Ti.ge.5. .or. debyelen.gt.0.2)then
                  call pltinit(-1.,1.,0.,25.)
               else
                  call pltinit(-1.,1.,0.,5.)
               endif
            else
               call pltinit(-1.,1.,0.,fluxmax)
            endif
c     call scalewn(-1.,1.,0.02,2.5,.false.,.true.)
            call charsize(0.02,0.02)
            call axis()
            call axlabels('cos!Aq!@',
     $     'Flux Density  /n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u')
            call winset(.true.)
            call charsize(0.0,0.0)
         endif
         Tiv(i)=Ti
         do it=nthhere,1,-1
            read(10,*)cost(it),flux(it)
            if(lprint)write(*,*)cost(it),flux(it)
            cosar(it,i,j)=cost(it)
            fluxar(it,i,j)=flux(it)
c Variable overall coefficient.
c            cx=sqrt((1.7+Ti)/(2.*3.14159))
            fit(it)=cx*
     $           exp(vd*0.5*((1.-cost(it))*upx+ (1.+cost(it))*downx))
         enddo
         if(mod(j,iskip).eq.0 .and. lcos)then
            call color(mod(j-1,15)+1)
            call fwrite(-vprobe,iwidth,1,charin)
            call labeline(cost,flux,nthhere,charin,iwidth)
            
c     call charsize(.0,.0)
c     call polyline(cost,flux,nthhere)
c     call vecw(-0.999,cx*exp(upx*vd),0)
c     call vecw(.999,cx*exp(downx*vd),1)
            call dashset(5)
            if(lfitplot)call polyline(cost,fit,nthhere)
            call dashset(0)
            call color(15)
         endif
         if(nfit.ge.2)then
            its=nthhere-nfit+1
            call fitlin(cost(its),flux(its),nfit,dummy(its),0,
     $           A,B,siga,sigb,chi2,Q)
            gdown(j,i)=A+B
            if(mod(j,iskip).eq.0 .and. lcos)then
               call color(13)
               call vecw(cost(its),A+B*cost(its),0)
               call vecw(cost(its+nfit-1),A+B*cost(its+nfit-1),1)
               call color(15)
            endif
            its=1
            call fitlin(cost(its),flux(its),nfit,dummy(its),0,
     $           A,B,siga,sigb,chi2,Q)
            gup(j,i)=A-B
            if(mod(j,iskip).eq.0 .and. lcos)then
               call color(13)
               call vecw(cost(its),A+B*cost(its),0)
               call vecw(cost(its+nfit-1),A+B*cost(its+nfit-1),1)
               call color(15)
            endif
         else
            gdown(j,i)=flux(nthhere)
     $           + (flux(nthhere)-flux(nthhere-1))*
     $           (1.-cost(nthhere))/(cost(nthhere)-cost(nthhere-1))
            gup(j,i)=flux(1)
     $           + (flux(1)-flux(2))*
     $           (-1.-cost(1))/(cost(1)-cost(2))
         endif
         write(*,*)'gdown/up',gdown(j,i),gup(j,i)
         updown(j,i)=gup(j,i)/gdown(j,i)
c Convert to calibration factor.
         updown(j,i)=log(updown(j,i))/vd
         vdj(j,i)=vd
         if(gup(j,i).gt.gmax)gmax=gup(j,i)
         if(gdown(j,i).gt.gmax)gmax=gdown(j,i)
         if(vd.gt.vdmax)vdmax=vd
         numv(i)=numv(i)+1
      enddo
 101  continue
      if(lcos)then
         call charsize(0.02,0.02)
         call fwrite(Ti,iwidth,1,carg)
         call jdrwstr(wx2nx(0.),.68,'T!di!d='//carg(1:13),1.)
         call fwrite(Vprobe,iwidth,2,carg)
         call jdrwstr(wx2nx(0.),.64,'!Af!@!dp!d='//carg,1.)
         call fwrite(debyelen,iwidth,2,carg)
         call jdrwstr(wx2nx(0.),.60,'!Al!@!dDe!d='//carg,1.)
         close(10)
         call charsize(0.0,0.0)

c      write(*,*)'debyearr,debyel,vparr:'
c      write(*,*)((debyearr(kj,ki),ki=1,i),kj=1,numv(i))
c      write(*,*)((debyel(kj,ki),ki=1,i),kj=1,numv(i))
c      write(*,*)((vparr(kj,ki),ki=1,i),kj=1,numv(i))
c         write(carg,'(''T!di!d='',f5.1)')Ti
c         call jdrwstr(wx2nx(0.),.6,carg(1:13),1.)
c         call fwrite(debyelen,iwidth,2,charin)
c         call jdrwstr(wx2nx(0.),.57,'!Al!@!dDe!d='//charin(1:iwidth),1.)
         call pltend()
      endif
c Write a table of the results.
      itemp=10*Ti
      write(filename,'(''tab'',i3.3,''.tex'')')itemp
      open(13,file=filename)
      write(13,501)
 501  format('\\begin{tabular}{',$)
      write(13,*)'|c|'
      do kj=1,numv(i)
         write(13,511)'c'
      enddo
      write(13,*)'|'
 511  format(a,$)
      write(13,*)'}\\hline'
 502  format(f6.3,'&',$)
 503  format(f6.3,' \\\\')
 504  format('\\quad{$v_f$}:&',$)
      write(13,504)
      do kj=1,numv(i)-1
         write(13,502)vdj(kj,i)
      enddo
      write(13,503)vdj(numv(i),i)
      write(13,*)'$\\cos\\theta$ &','\\multicolumn{',numv(i)-1,'}{c}',
     $     '{Ion Flux density, $\\Gamma$',
     $     '\\quad($/n_{i\\infty}(ZT_e/m)^{1/2})$}&\\\\'
      
      write(13,*)'\\hline'
      do kk=1,nthhere
         write(13,502)cosar(kk,i,1)
         do kj=1,numv(i)-1
            write(13,502)fluxar(kk,i,kj)
         enddo
         write(13,503)fluxar(kk,i,numv(i))
      enddo
      write(13,*)'\\hline\\end{tabular}'
      close(13)

c While there were non-null command line arguments repeat.
      if(carg(1:1) .eq. ' ') then
         goto 12
      endif
      goto 11
 12   continue

c Calculate the coulflux values
      do k=1,imax
         vpcoul(k)=vparr(1,1)+(vpmax-vparr(1,1))
     $        *(k-1.)/(imax-1.)
         coulf(k)=log(coulflux(3.14159,-vpcoul(k)/Ti,
     $        vd/sqrt(Ti))
     $        /coulflux(0.,-vpcoul(k)/Ti,vd/sqrt(Ti)))/vd
c         coulf(k)=log(coulflux(3.14159,-vpcoul(k),vd)
c     $        /coulflux(0.,-vpcoul(k),vd)/vd
      enddo
c      write(*,*)vpmax,vparr(1,1),vd
c      write(*,*)vpcoul,coulf
c Plot of the upstream and downstream values versus vprobe
c      write(*,*)'i,',i,'  vdj',(vdj(k,i),k=1,numv(i)),(gup(k,i),k=1,numv(i))
c      call pltinit(0.,vpmax,-1.,1.5)
      call pltinit(0.,vpmax,-3.,2.)
c      call scalewn(0.,vpmax,1.,10.,.false.,.true.)
      call axis()
      call axis2() 
c      call axlabels('v!df!d (nondimensional)',
      call axlabels('!Af!@!dprobe!d  /(T!de!d/e)',
     $     'ln(R)/v!df!d')
      call fwrite(Ti,iwidth,1,charin)
      call jdrwstr(0.04,.6,'T!di!d='//charin,1.)
      call color(13)
      call vecw(-5.,1.34,0)
      call vecw(vpmax,1.34,1)
      call jdrwstr(wx2nx(vpmax),wy2ny(1.34)+.02,'!Al!@!dDe!d= 0',-2.)
      call color(12)
      call polyline(vpcoul,coulf,imax)
      call jdrwstr(wx2nx(vpcoul(imax)),wy2ny(coulf(imax))+.03,
     $     '!Al!@!dDe!d= !A;!@',-2.)
c Drop the zero debye length case from plotting points.
      imin=0
      if(debyearr(1,1).lt.1.e-3)imin=1
      do kk=1+imin,i
         k=kk-imin
         call color(k)
         call polymark(vparr(1,kk),updown(1,kk),numv(kk),k)
         call polyline(vparr(1,kk),updown(1,kk),numv(kk))
         write(charin,'('' !Al!@!dDe!d='',f6.2)')debyearr(1,kk)
         call legendline(-0.6,.05*(i-k),k,charin(1:20))
c         if(Tiv(k).gt.8.)call freeflight(vdmax,Tiv(k),0)
      enddo
      call color(15)
      call pltend()

      call multiframe(1,2,3)
      call pltinit(0.,vpmax,0.,gmax)
c      call scalewn(0.,vpmax,.5,10.,.false.,.true.)
      call axis()
      do kk=1+imin,i
c         cxi=sqrt((1.7+Tiv(k))/(2.*3.14159))
         k=kk-imin
         call color(k)
c         call vecw(0.,cxi,0)
c         call vecw(vdmax,gun*cxi/cx,1)
         call polymark(vparr(1,kk),gup(1,kk),numv(kk),k)
         call polyline(vparr(1,kk),gup(1,kk),numv(kk))
         write(charin,'('' !Al!@!dDe!d='',f6.2)')debyearr(1,kk)
         call legendline(-0.1,.95+0.05*(k-i),k,charin(1:20))
c         if(Tiv(k).gt.8.)call freeflight(vdmax,Tiv(k),1)
      enddo
      call color(15)
c      call axlabels('v!df!d (nondimensional)',
c     $     'Flux density (nondimensional)')
      call axlabels('!Af!@!dprobe!d  /(T!de!d/e)',
     $     'Flux Density  /n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u')
      call boxtitle('Upstream')
      call axis2()

      call pltinit(0.,vpmax,0.,gmax)
c      call scalewn(0.,vpmax,.5,10.,.false.,.true.)
      call axis()
      do kk=1+imin,i
         k=kk-imin
         call color(k)
         call polymark(vparr(1,kk),gdown(1,kk),numv(kk),k)
         call polyline(vparr(1,kk),gdown(1,kk),numv(kk))
         write(charin,'('' !Al!@!dDe!d='',f6.2)')debyearr(1,kk)
         call legendline(-0.1,.95+0.05*(k-i),k,charin(1:20))
         if(Tiv(k).gt.8.)call freeflight(vpmax,Tiv(k),-1)
      enddo
      call color(15)
c      call lautomark(vdj(1,i),gdown(1,i),numv(i),.false.,.true.,2)
      call axlabels('!Af!@!dprobe!d  /(T!de!d/e)',
     $     'Flux Density  /n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u')
c      call axlabels('v!df!d (nondimensional)',
c     $     'Flux density (nondimensional)')
      call boxtitle('Downstream')
      call axis2()
      call pltend()
c      write(*,*)(cost(i),flux(i),i=1,nthhere)
c
c Contour plot needs the array to be rectangular.
      call multiframe(0,0,0)
c      write(*,*)i,numv(1)
c     $     ,vparr(1,1),vparr(numv(1),1),debyel(1,1),debyel(1,i)
      call pltinit(vparr(1,1),vparr(numv(1),1),debyel(1,1),debyel(1,i))
      call scalewn(vparr(1,1),vparr(numv(1),1),
     $     debyearr(1,1),debyearr(1,i),
     $     .false.,.true.)
      call axis()
      call scalewn(vparr(1,1),vparr(numv(1),1),debyel(1,1),debyel(1,i),
     $     .false.,.false.)
      call axlabels('Probe Potential','Debye Length')
      call boxtitle('Contours: ln(R)/v!df!d')
      zclv=0.
      call contourl(updown,cworka,imax,numv(1),i,zclv,0,
     $     vparr,debyel,2)
      call pltend()

      if(llam)then
c     Plot flux versus lambda. requires consistent numv
         call multiframe(1,2,3)
         call pltinit(0.,1.,0.,1.)
         call scalewn(deb1(1),deb1(i),0.,gmax,.true.,.false.)
         call axis()
         call axlabels('!Al!@!dDe!d','Flux density')
         call boxtitle('upstream')
         call crossplot(imax,numv(1),i,vparr(1,1),deb1,gup,
     $        ' !Af!@!dp!d=',.95-.05*numv(1),.05)
c     call pltend()
         call pltinit(0.,1.,0.,1.)
         call scalewn(deb1(1),deb1(i),0.,gmax,.true.,.false.)
         call axis()
         call axlabels('!Al!@!dDe!d','Flux density')
         call boxtitle('downstream')
         call crossplot(imax,numv(1),i,vparr(1,1),deb1,gdown,
     $        ' !Af!@!dp!d=',.95-.05*numv(1),.05)
         call pltend()
         call multiframe(0,0,0)

         call pltinit(0.,1.,0.,1.)
         call scalewn(deb1(1),deb1(i),-3.,2.,.true.,.false.)
         call axis()
         call axlabels('!Al!@!dDe!d',
     $        'ln(R)/v!df!d')
         call crossplot(imax,numv(1),i,vparr(1,1),deb1,updown,
     $        ' !Af!@!dp!d=',.08+numv(1)*0.05,-.05)
         call pltend()
      endif
c
c contour each probe potential flux vs costheta and lamdadebye.
c Requires the arrays to have constant dimensions.
c i is the number of files (lambdas)
      do ip=1,numv(1)
c         call pltinit(1.,0.,1.,0.)
c         call scalewn(cosar(1,1,1),cosar(nthhere,1,1),
c     $        deb1(1),deb1(i),.false.,.true.)
c         call axis()
c         call scalewn(cosar(1,1,1),cosar(nthhere,1,1),
c     $        debl1(1),debl1(i),.false.,.true.)
c         call axlabels('cos!Aq!@','!Al!@')
c         write(charin,'(''V!dprobe!d='',f4.0)')vparr(ip,1)
c         call boxtitle(charin)
c         call scalewn(1.,float(nthhere),1.,float(i),.false.,.false.)
c         call contourl(fluxar(1,1,ip),cworka,itmax,nthhere,i,10.,0,
c     $        cosar(1,1,1),debl1,1)
c         call pltend()
 21      call pltinit(0.,1.,0.,1.)
         isw=0
         call hidweb(cosar(1,1,1),debl1,fluxar(1,1,ip),
     $        itmax,nthhere,i,isw)
         call hdprset(-3,scbn(3))
c         call hdproject(0,1,3,1.,1)
         call scalewn(cosar(1,1,1),cosar(nthhere,1,1),
     $        deb1(1),deb1(i),.false.,.true.)
         call color(13)
         call axis()
         call axlabels('cos!Aq!@','log(!Al!@)')
c         call scalewn(cosar(1,1,1),cosar(nthhere,1,1),
c     $        debl1(1),debl1(i),.false.,.false.)
         call contourl(fluxar(1,1,ip),cworka,itmax,nthhere,i,10.,0,
     $        cosar(1,1,1),deb1,1)
         write(charin,'(''V!dprobe!d='',f4.0)')vparr(ip,1)
         call boxtitle(charin(1:20))
         call color(15)
         call hdprset(0,1.)
         call eye3d(isw)
         if(isw.eq.1)goto 21
      enddo
      return
 105  write(*,*)'Error opening file:',filename
      return
      end
c*************************************************************
      subroutine freeflight(vmax,Ti,isw)
c isw +1 upstream, -1 downstream, 0 ratio.
      integer nflight
      parameter(nflight=100)
      real vflight(nflight),rflight(nflight),y(nflight)

      pi=3.14159
      spi=sqrt(pi)
      do i=1,nflight
         vflight(i)=abs(vmax)*(i-1)/float(nflight)
         xi=vflight(i)/sqrt(2.*Ti)
         if(isw.gt.0.)xi=-xi
         y(i)=sqrt(2.*Ti)*(-xi*0.5*erfcc(xi) +exp(-xi**2)/(2.*spi))
         rflight(i)=(xi*0.5*erfcc(-xi) +exp(-xi**2)/(2.*spi))/
     $        (-xi*0.5*erfcc(xi) +exp(-xi**2)/(2.*spi))
      enddo
      call dashset(4)
      if(isw.eq.0)then
         call polyline(vflight,rflight,nflight)
      else
         call polyline(vflight,y,nflight)
      endif
      call dashset(0)
      end

c*******************************************************************
      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END

c*******************************************************************
c Plot versus second index.
      subroutine crossplot(nxL,nx,ny,xa,ya,fxy,string,y0,dy)
      integer nxL,nx,ny
      real xa(nx),ya(ny)
      real fxy(nxL,ny)
      character*(*) string

      character*100 char
      integer maxij
      parameter (maxij=1000)
      real fplot(maxij)

      if(nx.gt.maxij .or. ny.gt.maxij) return

      do i=1,nx
         do j=1,ny
            fplot(j)=fxy(i,j)
         enddo
         call color(i)
         call polyline(ya,fplot,ny)
         call polymark(ya,fplot,ny,i)
         write(char,100)string, xa(i)
 100     format(a,f5.1)
         call termchar(char)
         call legendline(-0.07,y0+dy*i,i,char)
c         write(*,*)nx,ny
c         write(*,*)(fplot(kk),kk=1,ny)
c         write(*,*)(ya(kk),kk=1,ny)
      enddo
      call color(15)
      end
c********************************************************************
      SUBROUTINE FITlin(X,Y,NDATA,SIG,MWT,A,B,SIGA,SIGB,CHI2,Q)
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA)
      SX=0.
      SY=0.
      ST2=0.
      B=0.
      IF(MWT.NE.0) THEN
        SS=0.
        DO 11 I=1,NDATA
          WT=1./(SIG(I)**2)
          SS=SS+WT
          SX=SX+X(I)*WT
          SY=SY+Y(I)*WT
11      CONTINUE
      ELSE
        DO 12 I=1,NDATA
          SX=SX+X(I)
          SY=SY+Y(I)
12      CONTINUE
        SS=FLOAT(NDATA)
      ENDIF
      SXOSS=SX/SS
      IF(MWT.NE.0) THEN
        DO 13 I=1,NDATA
          T=(X(I)-SXOSS)/SIG(I)
          ST2=ST2+T*T
          B=B+T*Y(I)/SIG(I)
13      CONTINUE
      ELSE
        DO 14 I=1,NDATA
          T=X(I)-SXOSS
          ST2=ST2+T*T
          B=B+T*Y(I)
14      CONTINUE
      ENDIF
      B=B/ST2
      A=(SY-SX*B)/SS
      SIGA=SQRT((1.+SX*SX/(SS*ST2))/SS)
      SIGB=SQRT(1./ST2)
      CHI2=0.
      IF(MWT.EQ.0) THEN
        DO 15 I=1,NDATA
          CHI2=CHI2+(Y(I)-A-B*X(I))**2
15      CONTINUE
        Q=1.
        SIGDAT=SQRT(CHI2/(NDATA-2))
        SIGA=SIGA*SIGDAT
        SIGB=SIGB*SIGDAT
      ELSE
        DO 16 I=1,NDATA
          CHI2=CHI2+((Y(I)-A-B*X(I))/SIG(I))**2
16      CONTINUE
c        Q=GAMMQ(0.5*(NDATA-2),0.5*CHI2)
        q=0
      ENDIF
      RETURN
      END
