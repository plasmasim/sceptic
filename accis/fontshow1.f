c
      program fontshow
c Plot aligned sets of the three installed fonts.
      character*256 entered
      character*129 str1
      real js,px,py
      integer i,j
c
      write(*,*)' Plot to file? (0:no,1:hp,2:ps,3:eps)'
      read(*,*)i
      call pfset(i)
      call pltinit(0.,1.,0.,1.)
      do 8 k=1,3
	 if(k.eq.3)then
	    noff=1
	 else
	    noff=0
	 endif
	 do 3 j=1,3
	    ii=32*j
	    if(k.le.2)then
	       do 1 i=ii,ii+31
		  entered(i-31:i-31)=char(i+(k-1)*128)
    1	       continue
	    else
	       do 11 i=ii,ii+31
c Control-B obsolete form.
		  entered((i-31)*2-1:(i-31)*2-1)=char(2)
		  entered((i-31)*2:(i-31)*2)=char(i)
   11	       continue
	    endif
    3	 continue
	 do 9 j=1,3
	    js=0.
	    do 10 i=1,32
	       py=0.6-(4*(k-1)+j)*0.043
	       px=0.02+.03*i
	       str1=entered((i+32*(j-1))*(noff+1)-noff:
     $		 (i+32*(j-1))*(noff+1))//char(0)
c	       write(*,*)char(27)//'[H',str1(1:70)
	       call drwstr(px,py,str1)
   10	    continue
    9	 continue
    8 continue
      call pltend()
      stop
      end

