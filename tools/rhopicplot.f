c Plot a series of rhopic curves getting the data from files on the
c command line. Those files are written by postproc when it plots rho.
      character*79 charin, ch2
      character*100 filename
      character*100 argstr
      logical lvlabel,lclabel,llx
      integer ir
      parameter (ir=401)
      real rpic(ir),rhopic(ir)
      data filename/' '/

      jj=0
      lvlabel=.false.
      lclabel=.false.
      llx=.false.
      do 1 i=1,iargc()
         call getarg(i,argstr)
         write(*,*)'Argument:',argstr(1:40)

         if(argstr(1:1) .eq. ' ') then
            goto 3
         endif
         if(argstr(1:1) .eq. '-') then
            if(argstr(1:2) .eq. '-x') llx=.true.
            if(argstr(1:2) .eq. '-v') lvlabel=.true.
            if(argstr(1:2) .eq. '-c') lclabel=.true.
c         if(argstr(1:2) .eq. '-i') linfrec=.true.
c         if(argstr(1:2) .eq. '-j') read(argstr(3:),*)jstepth
            if(argstr(1:2) .eq. '-?') goto 51
         else
            filename=argstr
            jj=jj+1
            open(13,file=argstr,status='old')
            read(13,'(a)')charin
            write(*,*)charin
c     the error label is necessary since one can encounter short lines.
            read(13,'(2f8.5,f8.4,f8.3,f8.3,f12.5,f10.5,i4,e13.4)'
     $           ,err=201)
     $           dt,vd,Ti,rmax,fave,debyelen,vprobe,icoln,colnwt
 201        continue
            write(*,*)dt,vd,Ti,rmax,fave,debyelen,vprobe,icoln,colnwt
            read(13,*)nrhere
            read(13,'(2f12.5)')(rpic(j),rhopic(j),j=1,nrhere)
            close(13)
            

            write(*,*)'nrhere=',nrhere
            if(jj.eq.1) then
               call pfset(3)
               if(.not.llx)then 
                  call autoplot(rpic,rhopic,nrhere)
               else
                  call pltinit(0.,1.,0.,1.)
                  call scalewn(1.,200.,1.,200.,.true.,.true.)
                  call axis()
                  call polyline(rpic,rhopic,nrhere)
               endif
               call axlabels('r [/r!dp!d]','Ion density [/n!di!A;!@!d]')
            else
               call color(jj)
               call dashset(jj)
               call polyline(rpic,rhopic,nrhere)
            endif
            imid=nrhere/3
c     call fwrite(debyelen,iwd,2,charin)
c     call jdrwstr(wx2nx(rpic(imid)),wy2ny(rhopic(imid))-.015,
c     $        '!Al!@!dDe!d='//charin,1.)
            if(lvlabel)then
               call fwrite(vd,iwd,2,charin)
               call termchar(charin)
               call legendline(0.4,0.8-.05*i,0,'v!dd!d='//charin)
            endif
            if(lclabel .and. icoln.ne.0)then
               if(colnwt.eq.0) then
                  nindex=0.
                  cx=0.
               else
                  nindex=nint(log10(colnwt)-.4999999)
                  cx=colnwt/10.**nindex
               endif
               call fwrite(cx,iwd1,1,charin)
               call iwrite(nindex,iwd2,ch2)
c               write(*,*)nindex,cx,iwd1,iwd2,charin,ch2
               call legendline(0.6,0.8-.05*i,0,
     $              ' !An!@!dc!d='//charin(1:iwd1)//'x10!u'//
     $              ch2(1:iwd2)//'!u'//char(0))
            endif
c         call jdrwstr(wx2nx(rpic(imid)),wy2ny(rhopic(imid))-.015,
c     $        'v!dd!d='//charin,1.)
         endif
 1    continue
 3    continue
      if(filename(1:1) .eq. ' ') goto 51
      slambda=debyelen/sqrt(1.+1./Ti)
      call winset(.true.)
      call lampecomp(slambda)
      call pltend()
      if(i.ne.1)call exit()
 51   write(*,*)'Usage: rhopicplot [-options] file1 [file2 ...]'
      write(*,*)'-x :log-axes -v :v-legend -c :coln-legend'
      end
c*************************************************************
      subroutine lampecomp(slambda)
c Overplot a line corresponding to the data of figure 2 of Lampe et al
c PoP 2003.

      integer npts
      parameter (npts=100)
      real x(npts),y(npts)
      integer nt,nu
      parameter (nt=14,nu=12)
      real xt(nt),yt(nt),dt(2,nt)
      real xu(nu),yu(nu),du(2,nu)
      data dt/
     $     2595,1890,2655,1935,2970,2355,3195,2640,3435,2880,3660,3075,
     $     3900,3270,4215,3465,4755,3750,5505,4080,6600,4500,7890,4935,
     $     9405,5400,10875,5850/
      data du/
     $     2610,3045,2790,3225,3015,3375,3195,3465,3570,3615,4035,3765,
     $     4575,3915,5295,4095,6615,4380,8355,4740,9855,5025,10875,5220/

c Shielding length. Lampe's x=(r-r_p)/slambda
      if(slambda.eq.0) slambda=1
c x-axis 0-3 at 1.e-4
      xl=2490
      yl=7845
      xr= 10905
      yr= 7785
c y-axis 1.e-4-1.e3 at 0.
      xb= 2490
      yb= 7845 
      xtop=2475
      ytop= 900

      do i=1,nt
         xt(i)=(dt(1,i)-(xb*(dt(2,i)-ytop)+xtop*(yb-dt(2,i)))/(yb-ytop))
     $        *3./(xr-xl)
         xt(i)=xt(i)*slambda+1.
         yt(i)=(dt(2,i) - (yl*(dt(1,i)-xr)+yr*(xl-dt(1,i)))/(xl-xr))
     $        *7./(ytop-yb) - 4.
         yt(i)=10.**yt(i)
      enddo

      do i=1,nu
         xu(i)=(du(1,i)-(xb*(du(2,i)-ytop)+xtop*(yb-du(2,i)))/(yb-ytop))
     $        *3./(xr-xl)
         xu(i)=xu(i)*slambda+1.
         yu(i)=(du(2,i) - (yl*(du(1,i)-xr)+yr*(xl-du(1,i)))/(xl-xr))
     $        *7./(ytop-yb) - 4.
         yu(i)=10.**yu(i) +1.
      enddo

      x0=max(xt(1),xu(1))
      x1=min(xt(nt),xu(nu))
      do i=1,npts
         x(i)=x0+(i-1)*(x1-x0)/(npts-1)
         zt=finterp(xt,yt,nt,x(i),p)
         zu=finterp(xu,yu,nu,x(i),p)
         y(i)=zt+zu
      enddo


      write(*,'(2f12.4)')(xt(i),yt(i),i=1,nt)
      write(*,'(2f12.4)')(xu(i),yu(i),i=1,nu)

      call legendline(.04,.25,0,' SCEPTIC')
      call color(3)
      call dashset(1)
      call polyline(xt,yt,nt)
      call legendline(.04,.2,0,' Trapped')
      call color(4)
      call dashset(3)
      call legendline(.04,.15,0,' Untrapped')
      call polyline(xu,yu,nu)
      call color(2)
      call dashset(2)
      call legendline(.04,.1,0,' Total Lampe')
      call polyline(x,y,npts)
      call dashset(0)
      call pltend()
      call lautoplot(xt,yt,nt,.false.,.true.)
      call axlabels('r/r_p','Density')
      call boxtitle('Tracing of Lampe trapped and untrapped densities')
      call polyline(xu,yu,nu)

      end
c********************************************************************
c Given two monotonic functions Q(x), R(x) 
c on a 1-D grid x=1..nq, solve Q(x)=y for x and interpolate R(x)
c That is, invert Q to give x=Q^-1(y) and give R(x), or simply
c interpolate the value of R when Q=y.
c If x is returned as 1, then the y is outside Q's range.
      function finterp(Q,R,nq,y,x)
      real Q(nq),R(nq)
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql

      finterp=0.
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         xpart=(y-Ql)/(Qr-Ql)
         x=xpart+iql
         finterp=R(iql)+xpart*(R(iqr)-R(iql))
      else
         x=iql
         write(*,*)'****** Error!: finterp coincident points'
      endif
      end
c**********************************************************************
