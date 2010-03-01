c Plot a series of phi radial curves getting the data from files on the
c command line. Those files are written by postproc when it plots phi.
      character*79 charin,ch2
      character*100 filename
      character*100 argstr
      logical lvlabel,lclabel,llx
      integer ir
      parameter (ir=401)
      real rpic(ir),phipic(ir),eta(ir),Er(ir)
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
            read(13,'(2f12.5)')(rpic(j),phipic(j),j=1,nrhere)
            close(13)
            
            jmax=nrhere-1
            jemax=nrhere-1
            do j=1,nrhere-1
               eta(j)=-phipic(j)
               if(j.gt.1)then
                  Er(j)=(phipic(j+1)-phipic(j))/(rpic(j+1)-rpic(j))
     $                 *.5*.25*(rpic(j+1)+rpic(j))**2 +
     $                 (phipic(j)-phipic(j-1))/(rpic(j)-rpic(j-1))
     $                 *.5*.25*(rpic(j)+rpic(j-1))**2 +1.e-6
               else
                  Er(j)=2.*(phipic(j+1)-phipic(j))/(rpic(j+1)-rpic(j))
     $                 *((rpic(j)+rpic(j+1))/2.)**2 -
     $                 (phipic(j+2)-phipic(j+1))/(rpic(j+2)-rpic(j+1))
     $                 *((rpic(j+1)+rpic(j+2))/2.)**2 + 1.e-6
               endif
               write(*,*)j,Er(j)
               if(eta(j).le.1.e-4)jmax=min(jmax,j-1)
               if(Er(j).le.1.e-4)jemax=min(jemax,j-1)
            enddo

            write(*,*)'nrhere=',nrhere
            if(jj.eq.1) then
               call pfset(3)
               if(.not.llx)then
                  call autoplot(rpic,phipic,nrhere)
               else
                  call pltinit(0.,1.,0.,1.)
                  call scalewn(1.,200.,1.e-4,10.,.true.,.true.)
                  call axis()
                  call polyline(rpic,eta,jmax)
c                  call lautoplot(rpic,eta,jmax,.true.,.true.)
               endif
               call axlabels('r [/r!dp!d]',
c     $              'Potential !Af!@!dp!d [/T!de!d/e], Charge')
     $              'Potential: -!Af!@!dp!d [/T!de!d/e]')
               call winset(.true.)
c               call polyline(rpic(1),Er(1),jemax-2)
            else
               call color(jj)
               call dashset(jj)
               if(.not.llx) call polyline(rpic,phipic,nrhere)
               if(llx) call polyline(rpic,eta,jmax)
c               call polyline(rpic(1),Er(1),jemax-2)
            endif
            imid=nrhere/3
c     call fwrite(debyelen,iwd,2,charin)
c     call jdrwstr(wx2nx(rpic(imid)),wy2ny(phipic(imid))-.015,
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
c     write(*,*)nindex,cx,iwd1,iwd2,charin,ch2
               call legendline(0.03,0.4-.05*i,0,
     $              ' !An!@!dc!d='//charin(1:iwd1)//'x10!u'//
     $              ch2(1:iwd2)//'!u'//char(0))
            endif
c         call jdrwstr(wx2nx(rpic(imid)),wy2ny(phipic(imid))-.015,
c     $        'v!dd!d='//charin,1.)
         endif
 1    continue
 3    continue
      if(filename(1:1) .eq. ' ') goto 51
      call pltend()
      if(i.ne.1)call exit()
 51   write(*,*)'Usage: phioutplot file1 [file2 ...]'
      write(*,*)'-x :log-axes -v :v-legend -c :coln-legend'
      end
