C********************************************************************
c*********************************************************************
c  character handling routines.
c**********************************************************************
      integer function nrdline(chrsfont,charb)
c Convert a line of char*1, and return total length including cr,lf.
      character*1 chrsfont(1)
      character*(*) charb
      character*1 ch
      integer i
      do 2 i=1,300
	 ch=chrsfont(i)
	    charb(i:i)=ch
	 if(ch.eq.char(10))goto 3
    2 continue
    3 nrdline=i
      write(*,*)charb
      return
      end
c********************************************************************
c Version 2, uses triple fonts, 384 characters,(cgi)
C Currently broken because of the input problems.
      subroutine getfont(fontname)
      character*12 fontname
      include 'plotcom.h'
      character*300 charb
      integer i,j,k,length
      integer nrdline

      i=0
      j=0
c Open and read in the whole file. This is the broken bit.
      open(11,FILE=fontname,STATUS='old',form='formatted',
     $  access='sequential',err=99)
      read(11,'(A1)',END=98) (chrsfont(j),j=1,buffer)
      write(*,*)chrsfont
   98 close(11)
      write(*,*)'j=',j
c Read the data bytes length from the header
      j=nrdline(chrsfont(1),charb)
      read(charb,'(i6)')length
      print '('' Length='',i6)',length
c Read the index: 24 lines of 16 integers.
      do 1 i=1,24
	 j=j+nrdline(chrsfont(j+1),charb)
c	 print '(i5)',j
	 read(charb,'(16i6)')(chrsaddr(16*(i-1)+k),k=1,16)
c	 print '(16i6)',(chrsaddr(16*(i-1)+k),k=1,16)
    1 continue
c Rearrange the addresses:
      do 2 i=1,384
	chrsaddr(i)=chrsaddr(i)+j
    2 continue
c Now chrsfont(chrsaddr(i)) is the start of the ith character, with the
c Hershey char no skipped.
      return
   99 print '('' Error opening font file:'',a12)',fontname
      stop
      end
