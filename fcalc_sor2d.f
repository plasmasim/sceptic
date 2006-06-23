
      subroutine testifconverged(eps,delta,lconverged,ierr)
      
      include 'mpif.h'
      include 'piccomsor.f'
      logical lconverged
      real convgd
      convgd=abs(delta)
      
      call MPI_ALLREDUCE(MPI_IN_PLACE,convgd,1,MPI_REAL,
     $     MPI_MAX,icommcart,ierr)
      
      if(convgd.lt.eps) then
         lconverged=.true.
      else
         lconverged=.false.
      endif

      end

c **********************************************

      subroutine blockdefine(ierr)

      include 'piccom.f'
      include 'piccomsor.f'

c     istep and jstep are the length that each process needs to process
c     (For example, istep=iside(1)-2). We need even values, except for
c     the uppermost
      integer istep,jstep

      istep=(NRUSED-2)/idims(1)
      if (mod(istep,2).eq.1) istep=istep-1
      jstep=(NTHUSED-2)/idims(2)
      if (mod(jstep,2).eq.1) jstep=jstep-1
      
c     Calulate iorigin
      iorig(1,1)=1

      do j=1,idims(2)
         jinc=(j-1)*nrsize*jstep
         iorig(1,j)=1+jinc
         do i=1,idims(1)
            iinc=(i-1)*istep
            iorig(i,j)=iorig(1,j)+iinc
         enddo
      enddo

c     Calculate the side length of the blocks
         iside(1,1)=istep+2
         iside(2,1)=jstep+2
         iside(1,2)=1+NRUSED-iorig(idims(1),1)
         iside(2,2)=(1+NTHUSED*nrsize-iorig(1,idims(2)))/nrsize

      end
c **********************************************
      subroutine topologyset(ierr)

      include 'piccom.f'
      include 'piccomsor.f'
      include 'mpif.h'

      logical lreorder
      nprocs=numprocs
      do n=1,2
         lperiod(n)=.false.
      enddo

c     create topology
      lreorder=.true.
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,idims,lperiod,
     $     lreorder,icommcart,ierr)
c      write(*,*) myid,icommcart
      call MPI_COMM_RANK(icommcart,mycartid,ierr)
      call MPI_CART_COORDS(icommcart, mycartid,2,icoords,ierr)

      do n=1,2
         if (icoords(n).eq.idims(n)-1) then
            ibt(n)=2
         else
            ibt(n)=1
         endif
         myside(n)=iside(n,ibt(n))
      enddo

      do n=1,2
         call MPI_CART_SHIFT(icommcart,2-n,1,isdr(n),iddr(n),ierr)
         call MPI_CART_SHIFT(icommcart,2-n,-1,isdl(n),iddl(n),ierr)
      enddo

      end
c *********************************************

      subroutine facecreate(ierr)

      include 'piccom.f'
      include 'piccomsor.f'

c     Length of the current face (How many points of given parity on the
c     given face)
      integer length

c     Indication of wether we are taking care of even (1) or odds (2) 
      integer p
c     Scratch stacks
      integer blen(100),disp(100)
c     nn is the normal direction
      do nn=1,2
         do p=1,2
            do i=1,2
               length=iside(nn,i)/2
               if (mod(iside(nn,i),2).eq.1.and.p.eq.1) then
                  length=length+1
               endif
               do k=1,length
                  blen(k)=1
                  disp(k)=(2*(k-1)+p-1)*(2-nn+nrsize*(nn-1))
               enddo
               write(*,*) nn,p,i
               call MPI_TYPE_INDEXED(length,blen(1),disp(1),
     $              MPI_REAL,iface(nn,p,i),ierr)
               call MPI_TYPE_COMMIT(iface(nn,p,i),ierr)

            enddo
         enddo
      enddo

      end

c *********************************************

      subroutine blockcreate(ierr)

      include 'piccom.f'
      include 'mpif.h'
      include 'piccomsor.f'

      call MPI_TYPE_SIZE(MPI_REAL,iSIZEOFREAL,ierr)

c     The goal here is to produce an array of types ktype
c     Each entry of the array represents a block. We end up
c     with 2^n useful entries [(bulk or uppermost)^n], starting
c     at rank ith0=2**ndims=4
      ith0=4
      ktype(1)=MPI_REAL
      inew=2
      do n=1,2
         iprior=2**(n-1)
         istride=((n-1)*nrsize+(2-n))*iSIZEOFREAL
         do i=1,2
            istep=iside(n,i)-2
            do iold=iprior,iprior+2**(n-1)-1
               call MPI_TYPE_CREATE_HVECTOR(istep,1,istride,
     $              ktype(iold),ktype(inew),ierr)
c     We just need to commit the upper level
               if (n.eq.2) call MPI_TYPE_COMMIT(ktype(inew),ierr)
               inew=inew+1
            enddo
         enddo
      enddo

      end
c *********************************************
      subroutine alltoallcreate(ierr)

      include 'piccom.f'
      include 'mpif.h'
      include 'piccomsor.f'

      integer icoordhere(2)
      integer ioffset
      integer ithi,ithj
      integer cartid

      ioffset=1+nrsize
      ithi=0

      do n=1,2
         if (icoords(n)+1.eq.idims(n)) ithi=ithi+2**(n-1)
      enddo

      do cartid=0,nprocs-1
         call MPI_CART_COORDS(icommcart,cartid,2,icoordhere,ierr)
         ithj=0
         do n=1,2
            if (icoordhere(n)+1.eq.idims(n)) ithj=ithj+2**(n-1)
         enddo

         isdispls(cartid)=(iorig(icoords(1)+1,icoords(2)+1)-1+ioffset)*
     $        isizeofreal
         istypes(cartid)=ktype(ith0+ithi)
        irdispls(cartid)=(iorig(icoordhere(1),icoordhere(2))-1+ioffset)*
     $        isizeofreal
         irtypes(cartid)=ktype(ith0+ithj)


c     Gather to process 0 (fortran index 1)
         if (cartid.eq.0.and.mycartid.ne.cartid) then
            iscounts(cartid)=1
         else
            iscounts(cartid)=0
         endif

         if(mycartid.eq.0.and.mycartid.ne.cartid) then
            ircounts(cartid)=1
         else
            ircounts(np)=0
         endif

      enddo

      end

c *********************************************


      subroutine sor2d(ierr)

      include 'piccom.f'
      include 'mpif.h'
      include 'piccomsor.f'

      integer status(MPI_STATUS_SIZE)
      logical firstcall,lconverged
      data firstcall/.true./

      save

      if (firstcall) then

c     Calculate the parameters for each block
         call blockdefine(ierr)
c     Set the cartesian communicator
         call topologyset(ierr)
c     Create types for the boundary communications for each iteration
         call facecreate(ierr)
c     Create types for the final block gathering
         call blockcreate(ierr)
c     Create the arrays for Alltoall
         call alltoallcreate(ierr)
      endif
      
c     Odd counter
      ko=mod(isor_k,2)+1
c     Even counter
      ke=mod(isor_k+1,2)+1
c     MPI tag
      itag=100

c     Origin of the left face
      iolp=iorig(icoords(1)+1,icoords(2)+1)

      do n=1,2
         if (n.eq.1) then
            iorp=iorig(icoords(1)+2,icoords(2)+1)
            iolm=iolp+1
            iorm=iorp+1
         else
            iorp=iorig(icoords(1)+1,icoords(2)+2)
            iolm=iolp+nrsize
            iorm=iorp+nrsize
         endif
      enddo

      ierr=0
      omega=1.
      relax=0.5
      
      i0=icoords(1)*istep+1
      j0=icoords(2)*jstep+1
c Main iteration      
      do isor_k=1,isor_mi
         


         sor_del=0.

         relax=(omega*debyelen**2+1.)/(debyelen**2+1.)

c Set outer BC (Here phi=0 on the boundary)
         do j=1,NTHUSED
            phi(NRFULL,j)=0
         enddo
         

         call sorrelax(myside(1),myside(2),phi(i0,j0),rho(i0,j0),
     $        apc(i0),bpc(i0),cpc(i0,j0),dpc(i0,j0),fpc(i0,j0),
     $        relax,delta,nrsize)


         call testifconverged(sor_eps,delta,lconverged,ierr)

         if(lconverged.and.isor_k.ge.2)goto 11
c Chebychev acceleration:
         if(isor_k.eq.1)then
            omega=1./(1.-0.5*sor_jac**2)
         else
            omega=1./(1.-0.25*sor_jac**2*omega)
         endif
      enddo

c We finished the loop, implies we did not converge.
      ierr=-isor_mi
      return
 11   continue
      ierr=isor_k

      end

c**********************************************************************
c Do a single relaxation.


      subroutine sorrelax(ni,nj,u,q,a,b,c,d,f,relax,sdel,nrsize)

      include 'piccomsor.f'
c     u and q are the potential and the charge. indice 1 represents
c     the left boundary value (Not to be unpgraded in sorrelax), and
c     we go to nbi or nbj (right boundary)
      real u(*),q(nrsize+1,nj)
      real a(nrsize+1),b(nrsize+1)
      real c(nrsize+1,nj),d(nrsize+1,nj),f(nrsize+1,nj)

      do n=1,2
         if (iddr(n).ne.-1) call MPI_SEND(u(iorp),
     $        1,iface(n,ko,ibt(n)),iddr(n),itag,icommcart,ierr)
         if (isdr(n).ne.-1) call MPI_RECV(u(iolp),1,
     $           iface(n,ko,ibt(n)),isdr(n),itag,icommcart,status,ierr)
         
         
         if(iddl(n).ne.-1) call MPI_SEND(u(iolm),1,
     $        iface(n,ke,ibt(n)),iddl(n),itag,icommcart,ierr)
         if(isdl(n).ne.-1) call MPI_RECV(u(iorm),1,
     $        iface(n,ke,ibt(n)),isdl(n),itag,icommcart,status,ierr)
         
      enddo
      

c Hard-wired red-black relaxation order.
      km=mod(isor_k,2)
c If isork is even, send even data to the right and top, odd to the bottom
c and left

      do j=2,nj-1

         do i=2+mod(j+km,2),ni-1,2

            expu=exp(u(i+j*nrsize))
            dnum= a(i)*u(i+1+j*nrsize)+b(i)*u(i-1+j*nrsize) + c(i,j)
     $           *u(i+(j+1)*nrsize)+d(i,j)*u(i+(j-1)*nrsize)
     $           -f(i,j)*u(i+j*nrsize)+ q(i,j) - expu
            dden=f(i,j) + expu
            delta=relax*dnum/dden

            if(abs(delta).gt.abs(sdel))sdel=delta
            uij=u(i+j*nrsize)+delta
            u(i+j*nrsize)=uij


         enddo

      enddo
      end


c***********************************************************************
c Here is where the subroutine is called from sceptic
c***********************************************************************
      subroutine fcalc_sor2d()

      include 'piccom.f'
      include 'piccomsor.f'

      
      sor_jac=1.-4/max(10,NRUSED)**2
      isor_mi=2.5*(NRUSED+NRUSED)
      sor_eps=1.e-5
         
      idims(1)=3
      idims(2)=2

      call sor2d(ierr)
      
c     write(*,'('':'',i3,$)') isor_k
      
      end
