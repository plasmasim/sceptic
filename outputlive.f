      subroutine outputLive(dt,step,icolntype)

      include 'piccom.f'
      include 'colncom.f'
      logical first
      data first /.true./
      character*35 filename
      integer step
      real dt
      filename=' '
      call nameappendexp(filename,'T',Ti,1)
      call nameappendint(filename,'v',nint(100*vd),3)
      if(r(nr).ge.100)then
         call nameappendint(filename,'R',ifix(r(nr)/10.),2)
      else
         call nameappendint(filename,'r',ifix(r(nr)),2)
      endif
      call nameappendint(filename,'P',ifix(abs(Vprobe)),2)
      if (infdbl) then
         call nameappendexp(filename,'L',1e5,1)
      else
         call nameappendexp(filename,'L',debyelen,1)
      endif
      if(Bz.ne.0.) call nameappendexp(filename,'B',Bz,2)
      if(icolntype.eq.1) call nameappendexp(filename,'c',colnwt,1)
      if(icolntype.eq.2) call nameappendexp(filename,'C',colnwt,1)
      idf=nbcat(filename,'.liv')

      if (first) then
         open(10,file=filename)
         first=.false.
      else
         open(10,file=filename,status='old',access='append')
      endif

      write(10,'(a,a)')' nstep',' dt'
      write(10,*) step,dt
      write(10,'(a)') 'Potential'
      do j=1,NRUSED
         write(10,'(10f8.3)')(phi(j,k),k=1,NTHUSED)
      enddo
      close(10)
      end
