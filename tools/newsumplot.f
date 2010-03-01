c Plot from summary
      integer imax,ifilemax,itmax
      parameter (imax=100,ifilemax=30,itmax=50)
      real cost(imax),flux(imax),fit(imax)
      real gup(imax,ifilemax),gdown(imax,ifilemax)
     $     ,updown(imax,ifilemax),vdj(imax,ifilemax)
      real cosar(itmax,imax,ifilemax),fluxar(itmax,imax,ifilemax)
      integer numv(ifilemax)
      real Tiv(ifilemax)
      character*100 filename
      character*100 charin
      character*100 carg
      logical lprint,lfitplot
      integer iskip
      real fluxmax
c Fit the curves by exponentials, on the up and downstream sides.
c gamma_up|down= cx*exp(vd*(up|down)x)
c And angular dependence is then
c gamma=cx*exp(vd*0.5*((1.-cost)upx + (1.+cost)downx)
      real upx,downx,cx,gun,gdn
      data upx,downx,cx/0.64,-.7,.62/
      data fluxmax/0./
c Default
      lfitplot=.false.
      lprint=.false.
      vdmax=0.
      iskip=1
      filename='summary'

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
         endif
      endif

      open(10,file=filename,status='old',err=105)

      numv(i)=0
      do j=1,50
         read(10,*,end=101,err=101)charin
         read(10,*)dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe
         read(10,*)charin
         read(10,*)nrhere,nthhere
         write(*,*)nrhere,nthhere
         write(*,901)dt,vd,Ti,isteps,rhoinf,phiinf,fave
 901     format('dt=',f6.3,' vd=',f6.3,' Ti=',f5.2,' Steps=',i4,
     $        ' rho=',f7.1,' phi=',f6.2,' Flux=',f7.4)
         read(10,*)charin
         if(j.eq.1)then
            call pfset(3)
            if(fluxmax.eq.0.)then
               if(Ti.ge.5.)then
                  call pltinit(-1.,1.,0.,6.5)
               else
                  call pltinit(-1.,1.,0.,7.5)
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
c         call winset(.true.)
         write(*,*)'Reading data'
         do it=nthhere,1,-1
            read(10,*)cost(it),flux(it)
            if(lprint)write(*,*)cost(it),flux(it)
            cosar(it,j,i)=cost(it)
            fluxar(it,j,i)=flux(it)
c Variable overall coefficient.
c            cx=sqrt((1.7+Ti)/(2.*3.14159))
            fit(it)=cx*
     $           exp(vd*0.5*((1.-cost(it))*upx+ (1.+cost(it))*downx))
         enddo
         gdown(j,i)=flux(nthhere)
     $        + (flux(nthhere)-flux(nthhere-1))*
     $        (1.-cost(nthhere))/(cost(nthhere)-cost(nthhere-1))
         gup(j,i)=flux(1)
     $        + (flux(1)-flux(2))*
     $        (-1.-cost(1))/(cost(1)-cost(2))
c         write(*,'(''Angle Corrections'',2f8.4,''%;'',2f8.4,''%.'')')
c     $        gup(j,i)-flux(1),100.*(gup(j,i)-flux(1))/flux(1),
c     $        gdown(j,i)-flux(nthhere),
c     $        100.*(gdown(j,i)-flux(nthhere))/flux(nthhere)
         updown(j,i)=gup(j,i)/gdown(j,i)
         vdj(j,i)=vd
         if(vd.gt.vdmax)vdmax=vd
         numv(i)=numv(i)+1
         if(mod(j,iskip).eq.0)then
c         call charsize(.01,.01)
         call color(mod(j-1,15)+1)
         call fwrite(vd,iwidth,1,charin)
         call labeline(cost,flux,nthhere,charin,iwidth)
c         call charsize(.0,.0)
c         call polyline(cost,flux,nthhere)
c         call vecw(-0.999,cx*exp(upx*vd),0)
c         call vecw(.999,cx*exp(downx*vd),1)
         call dashset(5)
         if(lfitplot)call polyline(cost,fit,nthhere)
         call dashset(0)
         call color(15)
         endif
      enddo
 101  continue
      call charsize(0.02,0.02)
      call fwrite(Ti,iwidth,1,carg)
      call jdrwstr(wx2nx(0.),.68,'T!di!d='//carg(1:13),1.)
      call fwrite(Vprobe,iwidth,2,carg)
      call jdrwstr(wx2nx(0.),.64,'!Af!@!dp!d='//carg,1.)
      call fwrite(debyelen,iwidth,2,carg)
      call jdrwstr(wx2nx(0.),.60,'!Al!@!dDe!d='//carg,1.)
      close(10)
      write(*,*)'File ',i,'  ',carg(1:30),'  Ti=',Ti
      call charsize(0.0,0.0)
      call pltend()
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
         write(13,502)cosar(kk,1,i)
         do kj=1,numv(i)-1
            write(13,502)fluxar(kk,kj,i)
         enddo
         write(13,503)fluxar(kk,numv(i),i)
      enddo
      write(13,*)'\\hline\\end{tabular}'
      close(13)

c While there were non-null command line arguments repeat.
      if(carg(1:1) .eq. ' ') then
         goto 12
      endif
      goto 11
 12   continue

c Only use this if we are varying the coefficients.
      cx=.61
      gun=cx*exp(upx*vdmax)
      gdn=cx*exp(downx*vdmax)
c      write(*,*)'i,',i,'  vdj',(vdj(k,i),k=1,numv(i)),(gup(k,i),k=1,numv(i))
      call pltinit(0.,1.,0.,1.)
      call scalewn(0.,vdmax,1.,150.,.false.,.true.)
      call axis()
      call axis2() 
c      call axlabels('v!df!d (nondimensional)',
      call axlabels('v!df!d  /(ZT!de!d/m!di!d)!u1/2!u',
     $     'Upstream/Downstream Ratio')
      call winset(.true.)
      call vecw(0.,1.,0)
      call vecw(vdmax,gun/gdn,1)
      do k=1,i
         call color(k)
         call polymark(vdj(1,k),updown(1,k),numv(k),k)
         write(charin,'('' T!di!d='',f5.2)')Tiv(k)
         call legendline(0.0,0.95-.05*k,k,charin(1:20))
         if(Tiv(k).gt.8.)call freeflight(vdmax,Tiv(k),0)
      enddo
      call color(15)
      call pltend()
      call multiframe(1,2,3)
      call pltinit(0.,1.,0.,1.)
      call scalewn(0.,vdmax,.1,10.,.false.,.true.)
      call axis()
      do k=1,i
c         cxi=sqrt((1.7+Tiv(k))/(2.*3.14159))
         call color(k)
c         call vecw(0.,cxi,0)
c         call vecw(vdmax,gun*cxi/cx,1)
         call polymark(vdj(1,k),gup(1,k),numv(k),k)
         write(charin,'('' T!di!d='',f5.2)')Tiv(k)
         call legendline(0.4,.01+0.05*k,k,charin(1:20))
         if(Tiv(k).gt.8.)call freeflight(vdmax,Tiv(k),1)
      enddo
      call color(15)
c      call axlabels('v!df!d (nondimensional)',
c     $     'Flux density (nondimensional)')
      call axlabels('v!df!d  /(ZT!de!d/m!di!d)!u1/2!u',
     $     'Flux Density  /n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u')
      call boxtitle('Upstream')
      call axis2()
      call winset(.true.)
      call vecw(0.,cx,0)
      call vecw(vdmax,gun,1)
      call winset(.false.)

      call pltinit(0.,1.,0.,1.)
      call scalewn(0.,vdmax,.005,2.,.false.,.true.)
      call axis()
      do k=1,i
         call color(k)
         call polymark(vdj(1,k),gdown(1,k),numv(k),k)
         write(charin,'('' T!di!d='',f5.2)')Tiv(k)
         call legendline(-0.1,.01+.05*k,k,charin(1:20))
         if(Tiv(k).gt.8.)call freeflight(vdmax,Tiv(k),-1)
      enddo
      call color(15)
c      call lautomark(vdj(1,i),gdown(1,i),numv(i),.false.,.true.,2)
      call axlabels('v!df!d  /(ZT!de!d/m!di!d)!u1/2!u',
     $     'Flux Density  /n!di!A;!@!d(ZT!de!d/m!di!d)!u1/2!u')
c      call axlabels('v!df!d (nondimensional)',
c     $     'Flux density (nondimensional)')
      call boxtitle('Downstream')
      call axis2()
      call vecw(0.,cx,0)
      call vecw(vdmax,gdn,1)

      call pltend()
c      write(*,*)(cost(i),flux(i),i=1,nthhere)
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

