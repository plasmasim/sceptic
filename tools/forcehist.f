      program forcehist
      parameter (nstepmax=10000)
      real zmn(nstepmax,5,2)
      character*200 charin,filename,string
      
c     Deal with arguments
      do 1 i=1,iargc()
         call getarg(i,string)
c     write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
c            if(string(1:2) .eq. '-e') collsum=.true.
c            if(string(1:2) .eq. '-j') lstandard=.false.
c            if(string(1:2) .eq. '-r') read(string(3:),*,err=71,end=71
c     $           )rmtoz
         else
            filename=string
         endif
 1    continue
 3    continue
      if(iargc().le.0)goto 71

      open(9,file=filename)
      read(9,*)istepmax
      read(9,'(a)')charin
      read(9,'(i5,5f12.5)')((i,(zmn(i,j,k),j=1,5)
     $     ,k=1,2),i=1,istepmax)
c      write(*,'(i5,5f12.5)')((i,(zmn(i,j,k),j=1,5)
c     $     ,k=1,2),i=1,istepmax)
      close(9)
      write(*,*)'Finished reading',istepmax

      call multiframe(2,1,1)
      call yautoplot(zmn(1,5,1),istepmax)
      call yautoplot(zmn(1,5,2),istepmax)
      call pltend()


      call exit(0)
 71   write(*,*)'Argument error'

      end
   
