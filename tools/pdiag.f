c Diagnose the final particle distributions.
c 8 Oct 05.
c Results with this program using particle data saved from runs with
c -x20 -d.2 -s500 -z1.e-5 show very good agreement between 
c the analytic form of f(v) and f(v) from the particle data.
c That confirms that the particle injection in being done correctly, since
c -z1.e-5 is essentially a zero field case, and so one ought to arrive
c at steady state where the distribution is the injection distribution.
      program pdiag
      include 'piccom.f'
      logical success
      integer nv
      parameter (nv=50)
      real fv1(-nv:nv),fv2(-nv:nv),fv3(-nv:nv),v(-nv:nv)
     $     ,fmaxwell(-nv:nv),fz(-nv:nv)
      character*100 string

c Defaults
      Ti=1.
      rmax=5.
      vd=0.
      myid=0

c Deal with arguments.
      if(iargc().eq.0) goto 51
      do 1 i=1,iargc()
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:2) .eq. '-s') diags=.true.
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
         if(string(1:2) .eq. '-x') read(string(3:),*)rmax
         if(string(1:2) .eq. '-v') read(string(3:),*)vd
         if(string(1:3) .eq. '-pf') then
            lfloat=.true.
         elseif(string(1:3) .eq. '-pi') then
            linsulate=.true.
         elseif(string(1:2) .eq. '-p') then
            read(string(3:),*)vprobe
         endif
         if(string(1:2) .eq. '-m') read(string(3:),*)myid
         if(string(1:2) .eq. '-?') goto 51
 1    continue
 3    continue

      vrange=max(4.,4.+5.*vd/sqrt(2.*Ti))
      do j=-nv,nv
         v(j)=j*vrange/nv
         fv1(j)=0.
         fv2(j)=0.
         fv3(j)=0.
         fmaxwell(j)=exp(-v(j)**2)
      enddo

      call partrd(success)
      if(.not.success)then
         write(*,*)'Failed to read the particles in'
         stop
      endif
      
      vmean=0.
      do i=1,npart
         iv1=min(nv,max(-nv,nint(xp(4,i)*nv/vrange)))
         iv2=min(nv,max(-nv,nint(xp(5,i)*nv/vrange)))
         iv3=min(nv,max(-nv,nint(xp(6,i)*nv/vrange)))
         fv1(iv1)=fv1(iv1)+1
         fv2(iv2)=fv2(iv2)+1
         fv3(iv3)=fv3(iv3)+1
         vmean=vmean+xp(6,i)
      enddo
      vmean=vmean/npart
      write(*,*)'Mean vz=',vmean,' cf vd=',vd

      call pfset(3)
      call minmax(fv1(-nv),2*nv+1,vmin,vmax)
      do j=-nv,nv
         fmaxwell(j)=exp(-v(j)**2/(2.*Ti))*vmax
      enddo
      call autoplot(v(-nv),fv1(-nv),2*nv+1)
      call axlabels('v','f(v)')
      call color(2)
      call polyline(v(-nv),fv2(-nv),2*nv+1)
      call color(3)
      call polyline(v(-nv),fv3(-nv),2*nv+1)
      call dashset(1)
      call color(4)
      call polyline(v(-nv),fmaxwell(-nv),2*nv+1)
      do j=-nv,nv
         fz(j)=1.77*vmax*fvcx(v(j)/sqrt(2.*Ti),vd/sqrt(2.*Ti))
      enddo
      call color(5)
      call polyline(v(-nv),fz(-nv),2*nv+1)
      call pltend()

      call exit(0)
 51   write(*,*)'Usage: pdiag -t -x -v -m'
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
      if(ipartmax.eq.npartmax)then
         npart=ipart
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
 101  write(*,*) 'No particle file to read:', filename
      end
c**********************************************************************
c****************************************************************
c FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vd,fvcx
c Return the normalized distribution function v_n f(v) for constant 
c cx collision frequency at a value of normalized velocity u=v/v_n,
c when the normalized drift velocity is ud= (a/\nu_c) /v_n,
c with v_n = sqrt(2T_n/m). a is acceleration, nu_c collision freq.
      if(ud.lt.0.) then
         v=-u
         vd=-ud
      else
         v=u
         vd=ud
      endif
      if(vd.eq.0.)then
         earg=100
      else
         earg=(0.5/vd)**2-v/vd
      endif
      if(earg.lt.50) then
         fvcx=exp(earg)*erfcc(0.5/vd-v)*0.5/vd
      else
c asymptotic form exp(-v^2)/sqrt(\pi):
         fvcx=exp(-v**2)/1.77245385
      endif
      if(.not.fvcx.ge.0) then
         write(*,*)'fvcx error. u=',u,' ud=',ud,' f=',fvcx,earg
         fvcx=0.
      endif
      end
c****************************************************************
c erfc from NR: (is in randf.f)
      FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
c*****************************************************************
